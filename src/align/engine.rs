// Adapted from c3align (Copyright 2026 Yibin Wang). See LICENSE.c3align.
#![allow(unused_imports, unused_variables, unused_unsafe)]
#[allow(dead_code)]
use anyhow::{anyhow, Result as AnyResult};
use crossbeam::channel;
use needletail::{parse_fastx_file, parse_fastx_stdin, parser::SequenceRecord};
use rammap::align::extend::AlignmentContext;
use rammap::align::index::Index;
use rammap::align::map::{AlignFlags, MapContext, MapOptions};
use rammap::align::pipeline::{align_and_format_query, OutputConfig, ReadInfo};
use rammap::api::apply_preset_str;
use rayon::prelude::*;
use rust_htslib::bam::{
    self, header::HeaderRecord, record::Aux, record::Cigar, record::CigarString,
    record::CigarStringView, Header, HeaderView, Read, Reader, Record, Writer,
};
use std::collections::HashMap;
use std::ffi::{CStr, CString};
use std::io::{BufRead, BufReader, Read as stdRead, Write};
use std::path::Path;
use std::sync::mpsc;
use std::sync::Arc;
use std::thread;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use crate::core::common_writer;
use crate::methalign::{MethConfig, MethRefiner};

// Opt-in aggregate elapsed time; worker stages overlap and are not CPU time.
struct AlignTimings {
    enabled: bool,
    micros: [AtomicU64; 13],
    calls: [AtomicU64; 13],
    started: Instant,
}

impl AlignTimings {
    fn new() -> Self {
        Self {
            enabled: std::env::var("CPHASING_ALIGN_TIMING").as_deref() == Ok("1"),
            micros: std::array::from_fn(|_| AtomicU64::new(0)),
            calls: std::array::from_fn(|_| AtomicU64::new(0)),
            started: Instant::now(),
        }
    }

    fn measure<T>(&self, stage: usize, f: impl FnOnce() -> T) -> T {
        if !self.enabled { return f(); }
        let start = Instant::now();
        let result = f();
        self.micros[stage].fetch_add(start.elapsed().as_micros() as u64, Ordering::Relaxed);
        self.calls[stage].fetch_add(1, Ordering::Relaxed);
        result
    }

    fn report(&self) {
        if !self.enabled { return; }
        let names = ["reference_load", "methylation_load", "index_build",
            "map_and_records", "methylation_score", "realign", "finalize",
            "restore_primary", "refresh_sa", "output_write_and_close",
            "rammap_align_and_format", "header_create", "sam_parse"];
        for (stage, (name, value)) in names.iter().zip(&self.micros).enumerate() {
            log::info!("ALIGN_TIMING stage={} aggregate_seconds={:.6} calls={}", name,
                value.load(Ordering::Relaxed) as f64 / 1e6,
                self.calls[stage].load(Ordering::Relaxed));
        }
        log::info!("ALIGN_TIMING wall_seconds={:.6}; worker stage times overlap; read stages cover SE HTS record path",
            self.started.elapsed().as_secs_f64());
    }
}

// Owned by a Rayon map_init state, never shared between concurrent tasks.
// Its lifetime is bounded by a batch using one reference/header. The direct
// PAF path never allocates a header; gap-rescue calls reuse the same dictionary.
struct MappingWorkspace {
    mapping: MapContext,
    header: Option<HeaderView>,
}

impl MappingWorkspace {
    fn new() -> Self {
        Self { mapping: MapContext::new(), header: None }
    }
}

pub struct RammapCtx {
    timings: Arc<AlignTimings>,
    pub methylation: Option<Arc<MethRefiner>>,
    pub target_names: Vec<String>,
    pub error: Arc<std::sync::Mutex<Option<String>>>,
    pub index: Index,
    pub opt: MapOptions,
    pub sensitive_index: Option<Index>,
    pub sensitive_opt: Option<MapOptions>,
    pub re_map: Option<Arc<GenomeREMap>>,
    pub mapq_calibrate: bool,
    pub candidate_probability: bool,
    pub graph_assignment: bool,
    pub gap_rescue: bool,
    pub min_gap_len: usize,
    pub out_cfg: OutputConfig,
    pub output_cigar: bool,
    pub do_cs: bool,
    pub cs_long: bool,
    pub eqx: bool,
    pub secondary: bool,
    pub soft_clip: bool,
}


fn copy_all_aux_except(rec_src: &Record, rec_dst: &mut Record, skip: &[&[u8; 2]]) -> AnyResult<()> {
    for item in rec_src.aux_iter() {
        let (tag, aux) = item?;
        if skip.iter().any(|t| t.as_slice() == tag) {
            continue;
        }
        // Alignment-derived tags take precedence when already present.
        if rec_dst.aux(tag).is_err() { rec_dst.push_aux(tag, aux)?; }
    }
    Ok(())
}

#[derive(Debug, PartialEq, Eq)]
pub struct SeqMetaData {
    pub name: String,
    pub length: u32,
    pub is_alt: bool,
}

fn get_score(record: &Record) -> i64 {
    match record.aux(b"AS") {
        Ok(Aux::I8(v)) => v as i64,
        Ok(Aux::U8(v)) => v as i64,
        Ok(Aux::I16(v)) => v as i64,
        Ok(Aux::U16(v)) => v as i64,
        Ok(Aux::I32(v)) => v as i64,
        Ok(Aux::U32(v)) => v as i64,
        Ok(Aux::Float(v)) => v as i64,
        _ => 0,
    }
}


fn normalize_homeolog_group(stem: &str) -> String {
    let lower = stem.to_ascii_lowercase();
    let Some(rest) = lower.strip_prefix("chr") else {
        return lower;
    };

    let digit_len = rest
        .bytes()
        .take_while(|b| b.is_ascii_digit())
        .count();
    if digit_len == 0 {
        return lower;
    }

    let digits = &rest[..digit_len];
    let suffix = &rest[digit_len..];
    if let Ok(n) = digits.parse::<u32>() {
        format!("chr{}{}", n, suffix)
    } else {
        lower
    }
}

fn infer_homeolog_label_from_target_name(name: &[u8]) -> Option<(String, char)> {
    let raw = std::str::from_utf8(name).ok()?.trim();
    let suffix = raw.chars().last()?.to_ascii_uppercase();
    if !matches!(suffix, 'A' | 'B' | 'C' | 'D') {
        return None;
    }

    let stem = &raw[..raw.len().saturating_sub(suffix.len_utf8())];
    if stem.is_empty() || !stem.chars().any(|c| c.is_ascii_digit()) {
        return None;
    }

    Some((normalize_homeolog_group(stem), suffix))
}

fn record_homeolog_label(r: &Record, target_names: &[&[u8]]) -> Option<(String, char)> {
    let tid = r.tid();
    if tid < 0 {
        return None;
    }
    infer_homeolog_label_from_target_name(target_names.get(tid as usize)?)
}

fn homeolog_tie_score_diff(query_len: usize) -> i64 {
    std::cmp::max(2, ((query_len as f64) * 0.0015).round() as i64)
}

fn set_i32_aux(rec: &mut Record, tag: &[u8; 2], value: i32) {
    remove_aux_if_present(rec, tag);
    rec.push_aux(tag, Aux::I32(value)).ok();
}

fn promote_record_from_context(primary: &mut Record, candidate: &Record, mapq: u8) {
    let old_flag = primary.flags();
    let mut new_record = Record::from(candidate.clone());
    new_record.set_mapq(mapq);

    let mut new_flag = new_record.flags();
    new_flag &= !0x100;
    new_flag |= old_flag & 0x40;
    new_flag |= old_flag & 0x80;
    new_flag |= old_flag & 0x800;
    new_record.set_flags(new_flag);
    set_i32_aux(&mut new_record, b"hr", 1);

    *primary = new_record;
}

fn promote_record_from_homeolog2(primary: &mut Record, candidate: &Record, mapq: u8, support: i64) {
    let old_flag = primary.flags();
    let mut new_record = Record::from(candidate.clone());
    new_record.set_mapq(mapq);

    let mut new_flag = new_record.flags();
    new_flag &= !0x100;
    new_flag |= old_flag & 0x40;
    new_flag |= old_flag & 0x80;
    new_flag |= old_flag & 0x800;
    new_record.set_flags(new_flag);
    set_i32_aux(&mut new_record, b"hr", 2);
    set_i32_aux(&mut new_record, b"hs", support.min(i32::MAX as i64) as i32);

    *primary = new_record;
}

fn promote_record_from_zero_homeolog(
    primary: &mut Record,
    candidate: &Record,
    mapq: u8,
    support: i64,
) {
    let old_flag = primary.flags();
    let mut new_record = Record::from(candidate.clone());
    new_record.set_mapq(mapq);

    let mut new_flag = new_record.flags();
    new_flag &= !0x100;
    new_flag |= old_flag & 0x40;
    new_flag |= old_flag & 0x80;
    new_flag |= old_flag & 0x800;
    new_record.set_flags(new_flag);
    set_i32_aux(&mut new_record, b"zr", 1);
    set_i32_aux(&mut new_record, b"zs", support.min(i32::MAX as i64) as i32);

    *primary = new_record;
}

pub struct AlignmentUnit {
    pub primary: Vec<Record>,
    pub secondary: Vec<Vec<Record>>,
}

impl AlignmentUnit {
    pub fn new() -> Self {
        AlignmentUnit {
            primary: Vec::new(),
            secondary: Vec::new(),
        }
    }

    pub fn from_records(records: Vec<Record>) -> Self {
        let mut au = AlignmentUnit::new();
        let mut idx = 0;
        for r in records {
            if !r.is_secondary() {
                au.add_primary(r);
                idx += 1;
            } else {
                if idx > 0 {
                    au.add_secondary(r, (idx - 1) as usize);
                }
            }
        }
        au
    }

    pub fn add_primary(&mut self, record: Record) {
        self.primary.push(record);
        self.secondary.push(Vec::new());
    }

    pub fn add_secondary(&mut self, record: Record, idx: usize) {
        if idx < self.secondary.len() {
            self.secondary[idx].push(record);
        }
    }

    pub fn rescue_robust(&mut self, mapq: u8, target_names: &[&[u8]]) {
        let mut anchors = Vec::new();
        for r in &self.primary {
            let target: i32 = r.tid();
            if r.mapq() >= mapq && !r.is_secondary() && target >= 0 {
                anchors.push((target, r.pos()));
            }
        }

        let limit_dist = 200_000;

        if anchors.is_empty() {
            return;
        }

        for (p, s) in self.primary.iter_mut().zip(self.secondary.iter_mut()) {
            if p.mapq() >= mapq {
                continue;
            }
            let p_score = get_score(p);
            let p_tid = p.tid();
            let p_pos = p.pos();

            let mut best_candidate_idx = None;
            let mut best_candidate_score = i64::MIN;
            let mut best_dist = i64::MAX;

            if !anchors.is_empty() {
                if p_tid >= 0 {
                    for &(a_tid, a_pos) in &anchors {
                        if p_tid == a_tid {
                            let dist = (p_pos - a_pos).abs();
                            if dist < limit_dist {
                                best_candidate_idx = Some(0);
                                best_candidate_score = p_score;
                                best_dist = dist;
                                break;
                            }
                        }
                    }
                }

                for (j, r) in s.iter().enumerate() {
                    let r_tid = r.tid();
                    let r_pos = r.pos();
                    if r_tid >= 0 {
                        for &(a_tid, a_pos) in &anchors {
                            if r_tid == a_tid && (r_pos - a_pos).abs() < limit_dist {
                                let dist = (r_pos - a_pos).abs();
                                let s_score = get_score(r);
                                if s_score as f32 >= (p_score as f32 * 0.95) {
                                    if dist < limit_dist {
                                        if s_score > best_candidate_score
                                            || (s_score == best_candidate_score && dist < best_dist)
                                        {
                                            best_candidate_score = s_score;
                                            best_candidate_idx = Some(j + 1);
                                            best_dist = dist;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if let Some(idx) = best_candidate_idx {
                if idx == 0 {
                    p.set_mapq(std::cmp::min(mapq, 10));
                } else {
                    let r = &s[idx - 1];
                    let flag = p.flags();
                    let mut new_record = Record::from(r.clone());
                    new_record.set_mapq(std::cmp::min(mapq, 10));
                    new_record.set_flags((new_record.flags() & !(0x100 | 0x800 | 0x40 | 0x80))
                        | (flag & (0x800 | 0x40 | 0x80)));
                    let mut displaced = p.clone();
                    displaced.set_flags((displaced.flags() & !0x800) | 0x100);
                    displaced.set_mapq(0);
                    *p = new_record;
                    s[idx - 1] = displaced;
                }
            } else if p.mapq() == 0 {
                let mut candidates = Vec::with_capacity(s.len() + 1);
                candidates.push((0, p_score));
                for (j, r) in s.iter().enumerate() {
                    candidates.push((j + 1, get_score(r)));
                }

                candidates.sort_by_key(|&(_, score)| std::cmp::Reverse(score));

                if !candidates.is_empty() {
                    let is_distinct = if candidates.len() == 1 {
                        true
                    } else {
                        let best_score = candidates[0].1;
                        let runner_up_score = candidates[1].1;
                        (best_score - runner_up_score >= 15)
                            || ((runner_up_score as f32) < (best_score as f32 * 0.85))
                    };

                    if is_distinct {
                        let best_idx = candidates[0].0;
                        if best_idx == 0 {
                            p.set_mapq(1);
                        } else {
                            let r = &s[best_idx - 1];
                            let flag = p.flags();
                            let mut new_record = Record::from(r.clone());
                            new_record.set_mapq(1);
                            new_record.set_flags((new_record.flags() & !(0x100 | 0x800 | 0x40 | 0x80))
                        | (flag & (0x800 | 0x40 | 0x80)));
                            let mut displaced = p.clone();
                            displaced.set_flags((displaced.flags() & !0x800) | 0x100);
                            displaced.set_mapq(0);
                            *p = new_record;
                            s[best_idx - 1] = displaced;
                        }
                    }
                }
            }
        }
    }

    pub fn rescue(&mut self, mapq: u8, _target_names: &[&[u8]]) {
        const NEAR_DIST: i64 = 2_00_000;
        const MIN_ANCHOR_MAPQ: u8 = 1;
        const RESCUED_MAPQ: u8 = 1;
        const MAX_RESCUED_MAPQ: u8 = 10;

        const SAME_TID_BONUS: i64 = 20;
        const NEAR_ANCHOR_BONUS: i64 = 30;
        const DIFF_TID_PENALTY: i64 = 5;

        const MIN_SCORE_GAP: i64 = 20;
        const MIN_SCORE_RATIO: f32 = 0.85;
        const MIN_CONTEXT_SUPPORT: i64 = 20;

        let min_anchor_mapq = std::cmp::max(mapq, MIN_ANCHOR_MAPQ);

        let mut anchors: Vec<(i32, i64)> = Vec::new();

        for r in &self.primary {
            let tid = r.tid();

            if tid < 0 {
                continue;
            }

            if r.is_unmapped() || r.is_secondary() {
                continue;
            }

            if r.mapq() >= min_anchor_mapq {
                anchors.push((tid, r.pos()));
            }
        }

        if anchors.is_empty() {
            return;
        }

        fn calc_context_score(
            r: &Record,
            anchors: &[(i32, i64)],
            near_dist: i64,
            same_tid_bonus: i64,
            near_anchor_bonus: i64,
            diff_tid_penalty: i64,
        ) -> i64 {
            let tid = r.tid();
            let pos = r.pos();

            if tid < 0 || r.is_unmapped() {
                return i64::MIN / 4;
            }

            let mut score: i64 = 0;
            let mut has_same_tid_anchor = false;

            for &(a_tid, a_pos) in anchors {
                if tid == a_tid {
                    has_same_tid_anchor = true;
                    score += same_tid_bonus;

                    let dist = (pos - a_pos).abs();

                    if dist <= near_dist {
                        score += near_anchor_bonus;

                        let penalty = (dist / 100_000).min(10);
                        score -= penalty;
                    }
                }
            }

            if !has_same_tid_anchor {
                score -= diff_tid_penalty;
            }

            score
        }

        fn promote_record_from_secondary(primary: &mut Record, secondary: &mut Record, mapq: u8) {
            let old_flag = primary.flags();

            let mut new_record = Record::from(secondary.clone());

            new_record.set_mapq(mapq);

            let mut new_flag = new_record.flags();

            new_flag &= !(0x100 | 0x800); // adopt the original segment role
            new_flag |= old_flag & 0x800;
            new_flag |= old_flag & 0x40; // preserve read1 flag
            new_flag |= old_flag & 0x80; // preserve read2 flag

            new_record.set_flags(new_flag);

            let mut displaced = primary.clone();
            displaced.set_flags((displaced.flags() & !0x800) | 0x100);
            displaced.set_mapq(0);
            *primary = new_record;
            *secondary = displaced;
        }

        for (p, s) in self.primary.iter_mut().zip(self.secondary.iter_mut()) {
            if p.mapq() >= mapq {
                continue;
            }

            if p.is_unmapped() {
                continue;
            }

            let mut candidates: Vec<(usize, i64, i64, i64)> = Vec::with_capacity(s.len() + 1);

            let p_score = get_score(p);
            let p_context_score = calc_context_score(
                p,
                &anchors,
                NEAR_DIST,
                SAME_TID_BONUS,
                NEAR_ANCHOR_BONUS,
                DIFF_TID_PENALTY,
            );

            let p_total_score = p_score + p_context_score;

            candidates.push((0, p_total_score, p_score, p_context_score));

            for (j, r) in s.iter().enumerate() {
                if r.is_unmapped() {
                    continue;
                }

                let r_score = get_score(r);
                let r_context_score = calc_context_score(
                    r,
                    &anchors,
                    NEAR_DIST,
                    SAME_TID_BONUS,
                    NEAR_ANCHOR_BONUS,
                    DIFF_TID_PENALTY,
                );

                let r_total_score = r_score + r_context_score;

                candidates.push((j + 1, r_total_score, r_score, r_context_score));
            }

            if candidates.is_empty() {
                continue;
            }

            candidates.sort_by_key(|&(_, total_score, _, _)| std::cmp::Reverse(total_score));

            let best = candidates[0];

            let best_idx = best.0;
            let best_total_score = best.1;
            let best_align_score = best.2;
            let best_context_score = best.3;

            if best_context_score < MIN_CONTEXT_SUPPORT {
                continue;
            }

            let is_distinct = if candidates.len() == 1 {
                true
            } else {
                let second = candidates[1];
                let second_total_score = second.1;

                let score_gap = best_total_score - second_total_score;

                score_gap >= MIN_SCORE_GAP
                    || (second_total_score as f32) < (best_total_score as f32 * MIN_SCORE_RATIO)
            };

            if !is_distinct {
                continue;
            }

            let rescued_mapq = if best_align_score >= p_score {
                std::cmp::min(mapq, MAX_RESCUED_MAPQ)
            } else {
                RESCUED_MAPQ
            };

            if best_idx == 0 {
                p.set_mapq(rescued_mapq);
            } else {
                let r = &mut s[best_idx - 1];
                promote_record_from_secondary(p, r, rescued_mapq);
            }
        }
    }

    pub fn rescue_homeolog_tie(&mut self, mapq: u8, target_names: &[&[u8]]) {
        const MIN_ANCHOR_MAPQ: u8 = 10;
        const MIN_ANCHOR_SCORE: i64 = 100;
        const RESCUED_MAPQ: u8 = 10;

        let mut votes: HashMap<(String, char), i64> = HashMap::new();
        for r in &self.primary {
            if r.is_unmapped() || r.mapq() < MIN_ANCHOR_MAPQ || get_score(r) < MIN_ANCHOR_SCORE {
                continue;
            }
            if let Some((group, suffix)) = record_homeolog_label(r, target_names) {
                let weight = r.mapq() as i64 + (get_score(r) / 100).max(1);
                *votes.entry((group, suffix)).or_insert(0) += weight;
            }
        }

        if votes.is_empty() {
            return;
        }

        for (p, s) in self.primary.iter_mut().zip(self.secondary.iter()) {
            if p.is_unmapped() || p.mapq() > mapq {
                continue;
            }

            let Some((p_group, _)) = record_homeolog_label(p, target_names) else {
                continue;
            };

            let p_score = get_score(p);
            let p_q = get_query_coords_on_read(p);
            let p_qlen = p_q.1.saturating_sub(p_q.0);
            let mut candidates: Vec<(usize, i64, char, usize)> = Vec::with_capacity(s.len() + 1);
            candidates.push((
                0,
                p_score,
                record_homeolog_label(p, target_names).unwrap().1,
                p_qlen,
            ));

            for (j, r) in s.iter().enumerate() {
                if r.is_unmapped() {
                    continue;
                }
                let Some((group, suffix)) = record_homeolog_label(r, target_names) else {
                    continue;
                };
                if group != p_group {
                    continue;
                }
                let q = get_query_coords_on_read(r);
                let qlen = q.1.saturating_sub(q.0);
                candidates.push((j + 1, get_score(r), suffix, qlen));
            }

            if candidates.len() < 2 {
                continue;
            }

            let top_score = candidates
                .iter()
                .map(|(_, score, _, _)| *score)
                .max()
                .unwrap_or(p_score);

            let mut tied: Vec<(usize, i64, char, i64)> = candidates
                .into_iter()
                .filter_map(|(idx, score, suffix, qlen)| {
                    let tie_score_diff = homeolog_tie_score_diff(qlen);
                    if top_score.saturating_sub(score) <= tie_score_diff {
                        let vote = votes.get(&(p_group.clone(), suffix)).copied().unwrap_or(0);
                        Some((idx, score, suffix, vote))
                    } else {
                        None
                    }
                })
                .collect();

            if tied.len() < 2 || tied.iter().all(|(_, _, _, vote)| *vote == 0) {
                continue;
            }

            tied.sort_by(|a, b| {
                b.3.cmp(&a.3)
                    .then_with(|| b.1.cmp(&a.1))
                    .then_with(|| a.0.cmp(&b.0))
            });

            let best = tied[0];
            let second_vote = tied.get(1).map(|x| x.3).unwrap_or(0);
            if best.3 == 0 || best.3 == second_vote {
                continue;
            }

            if best.0 == 0 {
                p.set_mapq(p.mapq().max(RESCUED_MAPQ));
                remove_aux_if_present(p, b"hr");
                p.push_aux(b"hr", Aux::I32(1)).ok();
            } else {
                promote_record_from_context(p, &s[best.0 - 1], RESCUED_MAPQ);
            }
        }
    }

    pub fn rescue_homeolog_interval(&mut self, mapq: u8, target_names: &[&[u8]]) {
        const MIN_SUPPORT_MAPQ: u8 = 10;
        const MIN_SUPPORT_SCORE: i64 = 100;
        const MIN_SUPPORT_LEN: usize = 300;
        const RESCUED_MAPQ: u8 = 10;

        let mut supports: Vec<Option<(String, char, i64)>> = Vec::with_capacity(self.primary.len());
        for (p, s) in self.primary.iter().zip(self.secondary.iter()) {
            if p.is_unmapped()
                || p.mapq() < MIN_SUPPORT_MAPQ
                || get_score(p) < MIN_SUPPORT_SCORE
            {
                supports.push(None);
                continue;
            }

            let Some((group, suffix)) = record_homeolog_label(p, target_names) else {
                supports.push(None);
                continue;
            };

            let q = get_query_coords_on_read(p);
            let qlen = q.1.saturating_sub(q.0);
            if qlen < MIN_SUPPORT_LEN {
                supports.push(None);
                continue;
            }

            let p_score = get_score(p);
            let second_same_group = s
                .iter()
                .filter_map(|r| {
                    let (r_group, _) = record_homeolog_label(r, target_names)?;
                    (r_group == group).then_some(get_score(r))
                })
                .max()
                .unwrap_or(i64::MIN);
            let uniqueness_gap = p_score.saturating_sub(second_same_group);
            let mut weight = if p.mapq() >= 20 { 20 } else { 10 };
            if qlen >= 1000 {
                weight += 5;
            }
            if uniqueness_gap >= 10 {
                weight += 10;
            }
            supports.push(Some((group, suffix, weight)));
        }

        for (idx, (p, s)) in self.primary.iter_mut().zip(self.secondary.iter()).enumerate() {
            if p.is_unmapped() || p.mapq() > mapq {
                continue;
            }

            let Some((p_group, _)) = record_homeolog_label(p, target_names) else {
                continue;
            };

            let p_score = get_score(p);
            let p_q = get_query_coords_on_read(p);
            let p_qlen = p_q.1.saturating_sub(p_q.0);
            let mut candidates: Vec<(usize, i64, char, usize)> = Vec::with_capacity(s.len() + 1);
            candidates.push((
                0,
                p_score,
                record_homeolog_label(p, target_names).unwrap().1,
                p_qlen,
            ));

            for (j, r) in s.iter().enumerate() {
                if r.is_unmapped() {
                    continue;
                }
                let Some((group, suffix)) = record_homeolog_label(r, target_names) else {
                    continue;
                };
                if group != p_group {
                    continue;
                }
                let q = get_query_coords_on_read(r);
                let qlen = q.1.saturating_sub(q.0);
                candidates.push((j + 1, get_score(r), suffix, qlen));
            }

            if candidates.len() < 2 {
                continue;
            }

            let top_score = candidates
                .iter()
                .map(|(_, score, _, _)| *score)
                .max()
                .unwrap_or(p_score);

            let mut suffix_support: HashMap<char, i64> = HashMap::new();
            for (support_idx, support) in supports.iter().enumerate() {
                if support_idx == idx {
                    continue;
                }
                let Some((group, suffix, weight)) = support else {
                    continue;
                };
                if group == &p_group {
                    *suffix_support.entry(*suffix).or_insert(0) += *weight;
                }
            }

            if suffix_support.is_empty() {
                continue;
            }

            let mut tied: Vec<(usize, i64, char, i64)> = candidates
                .into_iter()
                .filter_map(|(candidate_idx, score, suffix, qlen)| {
                    let tie_score_diff = homeolog_tie_score_diff(qlen);
                    if top_score.saturating_sub(score) <= tie_score_diff {
                        let support = suffix_support.get(&suffix).copied().unwrap_or(0);
                        Some((candidate_idx, score, suffix, support))
                    } else {
                        None
                    }
                })
                .collect();

            if tied.len() < 2 || tied.iter().all(|(_, _, _, support)| *support == 0) {
                continue;
            }

            tied.sort_by(|a, b| {
                b.3.cmp(&a.3)
                    .then_with(|| b.1.cmp(&a.1))
                    .then_with(|| a.0.cmp(&b.0))
            });

            let best = tied[0];
            let second_support = tied.get(1).map(|x| x.3).unwrap_or(0);
            if best.3 == 0 || best.3 == second_support {
                continue;
            }

            if best.0 == 0 {
                p.set_mapq(p.mapq().max(RESCUED_MAPQ));
                set_i32_aux(p, b"hr", 2);
                set_i32_aux(p, b"hs", best.3.min(i32::MAX as i64) as i32);
            } else {
                promote_record_from_homeolog2(p, &s[best.0 - 1], RESCUED_MAPQ, best.3);
            }
        }
    }

    pub fn rescue_zero_mapq_homeolog(&mut self, target_names: &[&[u8]]) {
        const MIN_INTERVALS: usize = 2;
        const MIN_QUERY_LEN: usize = 100;
        const MAPQ5_SUPPORT: i64 = 15;
        const MAPQ3_SUPPORT: i64 = 8;
        const MAPQ1_SUPPORT: i64 = 4;

        if self.primary.is_empty()
            || self
                .primary
                .iter()
                .any(|r| !r.is_unmapped() && r.mapq() > 0)
        {
            return;
        }

        #[derive(Clone)]
        struct ZeroCandidate {
            idx: usize,
            score: i64,
            suffix: char,
            qlen: usize,
        }

        let mut interval_candidates: Vec<Option<(String, Vec<ZeroCandidate>)>> =
            Vec::with_capacity(self.primary.len());
        let mut usable_intervals = 0usize;

        for (p, s) in self.primary.iter().zip(self.secondary.iter()) {
            let Some((group, suffix)) = record_homeolog_label(p, target_names) else {
                interval_candidates.push(None);
                continue;
            };
            let p_q = get_query_coords_on_read(p);
            let p_qlen = p_q.1.saturating_sub(p_q.0);
            if p_qlen < MIN_QUERY_LEN {
                interval_candidates.push(None);
                continue;
            }

            let mut candidates = vec![ZeroCandidate {
                idx: 0,
                score: get_score(p),
                suffix,
                qlen: p_qlen,
            }];

            for (j, r) in s.iter().enumerate() {
                if r.is_unmapped() {
                    continue;
                }
                let Some((r_group, r_suffix)) = record_homeolog_label(r, target_names) else {
                    continue;
                };
                if r_group != group {
                    continue;
                }
                let q = get_query_coords_on_read(r);
                let qlen = q.1.saturating_sub(q.0);
                if qlen < MIN_QUERY_LEN {
                    continue;
                }
                candidates.push(ZeroCandidate {
                    idx: j + 1,
                    score: get_score(r),
                    suffix: r_suffix,
                    qlen,
                });
            }

            if candidates.len() < 2 {
                interval_candidates.push(None);
                continue;
            }

            let top_score = candidates
                .iter()
                .map(|c| c.score)
                .max()
                .unwrap_or_else(|| get_score(p));
            candidates.retain(|c| {
                top_score.saturating_sub(c.score) <= homeolog_tie_score_diff(c.qlen)
            });

            if candidates.len() < 2 {
                interval_candidates.push(None);
                continue;
            }

            usable_intervals += 1;
            interval_candidates.push(Some((group, candidates)));
        }

        if usable_intervals < MIN_INTERVALS {
            return;
        }

        let mut support_by_group: HashMap<(String, char), i64> = HashMap::new();
        let mut interval_supports: Vec<Option<(String, HashMap<char, i64>)>> =
            Vec::with_capacity(interval_candidates.len());
        for item in &interval_candidates {
            let Some((group, candidates)) = item else {
                interval_supports.push(None);
                continue;
            };
            let mut best_by_suffix: HashMap<char, i64> = HashMap::new();
            for c in candidates {
                best_by_suffix
                    .entry(c.suffix)
                    .and_modify(|v| *v = (*v).max(c.score))
                    .or_insert(c.score);
            }
            if best_by_suffix.len() < 2 {
                interval_supports.push(None);
                continue;
            }
            let mean = best_by_suffix.values().sum::<i64>() as f64 / best_by_suffix.len() as f64;
            let mut current_support: HashMap<char, i64> = HashMap::new();
            for (suffix, score) in best_by_suffix {
                let delta = (score as f64 - mean).round() as i64;
                if delta > 0 {
                    *support_by_group.entry((group.clone(), suffix)).or_insert(0) += delta;
                    current_support.insert(suffix, delta);
                }
            }
            interval_supports.push(Some((group.clone(), current_support)));
        }

        if support_by_group.is_empty() {
            for p in &mut self.primary {
                set_i32_aux(p, b"za", 1);
            }
            return;
        }

        for (idx, (p, s)) in self.primary.iter_mut().zip(self.secondary.iter()).enumerate() {
            let Some((group, candidates)) = interval_candidates.get(idx).and_then(|x| x.as_ref()) else {
                continue;
            };

            let mut ranked: Vec<(usize, i64, char, i64)> = candidates
                .iter()
                .map(|c| {
                    let mut support = support_by_group
                        .get(&(group.clone(), c.suffix))
                        .copied()
                        .unwrap_or(0);
                    if let Some((self_group, self_supports)) =
                        interval_supports.get(idx).and_then(|x| x.as_ref())
                    {
                        if self_group == group {
                            support = support
                                .saturating_sub(self_supports.get(&c.suffix).copied().unwrap_or(0));
                        }
                    }
                    (c.idx, c.score, c.suffix, support)
                })
                .collect();
            ranked.sort_by(|a, b| {
                b.3.cmp(&a.3)
                    .then_with(|| b.1.cmp(&a.1))
                    .then_with(|| a.0.cmp(&b.0))
            });

            let best = ranked[0];
            let second_support = ranked.get(1).map(|x| x.3).unwrap_or(0);
            if best.3 == 0 || best.3 == second_support {
                set_i32_aux(p, b"za", 1);
                continue;
            }

            let new_mapq = if best.3 >= MAPQ5_SUPPORT {
                5
            } else if best.3 >= MAPQ3_SUPPORT {
                3
            } else if best.3 >= MAPQ1_SUPPORT {
                1
            } else {
                set_i32_aux(p, b"za", 1);
                continue;
            };

            if best.0 == 0 {
                p.set_mapq(new_mapq);
                set_i32_aux(p, b"zr", 1);
                set_i32_aux(p, b"zs", best.3.min(i32::MAX as i64) as i32);
            } else {
                promote_record_from_zero_homeolog(p, &s[best.0 - 1], new_mapq, best.3);
            }
        }
    }

    pub fn filter_overlapped_segments(&mut self) {
        if self.primary.len() < 2 {
            return;
        }

        let mut indices: Vec<usize> = (0..self.primary.len()).collect();
        indices.sort_by_key(|&idx| std::cmp::Reverse(get_score(&self.primary[idx])));

        let mut kept = Vec::new();
        for idx in indices {
            let r = &self.primary[idx];
            let (s, e) = get_query_coords(r);
            if e == s {
                kept.push(idx);
                continue;
            }

            let mut has_substantial_overlap = false;
            for &k_idx in &kept {
                let kr = &self.primary[k_idx];
                let (ks, ke) = get_query_coords(kr);
                if ke == ks {
                    continue;
                }

                let overlap = std::cmp::min(e, ke).saturating_sub(std::cmp::max(s, ks));
                let len1 = e - s;
                let len2 = ke - ks;
                let min_len = std::cmp::min(len1, len2);
                if min_len > 0 && (overlap as f32 / min_len as f32) > 0.25 {
                    has_substantial_overlap = true;
                    break;
                }
            }

            if !has_substantial_overlap {
                kept.push(idx);
            }
        }

        kept.sort();

        let mut new_primary = Vec::with_capacity(kept.len());
        let mut new_secondary = Vec::with_capacity(kept.len());
        for idx in kept {
            new_primary.push(std::mem::replace(&mut self.primary[idx], Record::new()));
            new_secondary.push(std::mem::replace(&mut self.secondary[idx], Vec::new()));
        }
        self.primary = new_primary;
        self.secondary = new_secondary;
    }

    pub fn into_records(self, keep_secondary: bool) -> Vec<Record> {
        if !keep_secondary {
            self.primary
        } else {
            let mut res = Vec::new();
            for (p, mut s_list) in self.primary.into_iter().zip(self.secondary.into_iter()) {
                res.push(p);
                res.append(&mut s_list);
            }
            res
        }
    }
}

#[derive(Clone)]
pub struct GenomeREMap {
    pub chrom_cuts: HashMap<i32, Vec<u32>>,
}

fn remove_aux_if_present(rec: &mut Record, tag: &[u8; 2]) {
    let _ = rec.remove_aux(tag);
}

impl GenomeREMap {
    pub fn build(
        reference_fasta: &str,
        patterns: &[Vec<u8>],
        tid_map: &HashMap<String, i32>,
    ) -> AnyResult<Self> {
        log::info!("Building Genome Restriction Site Map for annotate...");
        let mut chrom_cuts = HashMap::new();
        if patterns.is_empty() {
            return Ok(GenomeREMap { chrom_cuts });
        }

        let mut reader = parse_fastx_file(reference_fasta)?;
        while let Some(Ok(rec)) = reader.next() {
            let name = String::from_utf8_lossy(rec.id()).to_string();
            if let Some(&tid) = tid_map.get(&name) {
                let seq = rec.seq();
                let mut cuts = Vec::new();
                cuts.push(0);

                let n = seq.len();
                let mut i = 0;
                while i < n {
                    let mut matched_len = 0;
                    for pattern in patterns {
                        let p_len = pattern.len();
                        if i + p_len <= n {
                            let mut matches = true;
                            for j in 0..p_len {
                                if seq[i + j].to_ascii_uppercase()
                                    != pattern[j].to_ascii_uppercase()
                                {
                                    matches = false;
                                    break;
                                }
                            }
                            if matches {
                                matched_len = p_len;
                                break;
                            }
                        }
                    }
                    if matched_len > 0 {
                        cuts.push(i.saturating_add(1) as u32);
                        i += matched_len;
                    } else {
                        i += 1;
                    }
                }
                cuts.push(n as u32);
                chrom_cuts.insert(tid, cuts);
            }
        }
        Ok(GenomeREMap { chrom_cuts })
    }

    pub fn annotate_record(&self, r: &mut Record) -> AnyResult<()> {
        let tid = r.tid();
        if tid < 0 {
            return Ok(());
        }

        if let Some(cuts) = self.chrom_cuts.get(&tid) {
            let start = r.pos() as u32;
            let end = r.cigar().end_pos() as u32;

            let idx = match cuts.binary_search(&start) {
                Ok(i) => i,
                Err(i) => i.saturating_sub(1),
            };
            let end_for_fragment = end.saturating_sub(1);
            let end_idx = match cuts.binary_search(&end_for_fragment) {
                Ok(i) => i,
                Err(i) => i.saturating_sub(1),
            };

            let fragment_id = idx as i32;
            remove_aux_if_present(r, b"re");
            remove_aux_if_present(r, b"rd");
            remove_aux_if_present(r, b"rf");
            remove_aux_if_present(r, b"rs");
            r.push_aux(b"re", Aux::I32(fragment_id))?;

            let left_cut = cuts[idx];
            let right_cut = if idx + 1 < cuts.len() {
                cuts[idx + 1]
            } else {
                cuts[cuts.len() - 1]
            };

            let dist_left = start.saturating_sub(left_cut);
            let dist_right = right_cut.saturating_sub(end);
            let min_dist = std::cmp::min(dist_left, dist_right) as i32;

            r.push_aux(b"rd", Aux::I32(min_dist))?;
            let fragment_span = end_idx.saturating_sub(idx).saturating_add(1) as i32;
            r.push_aux(b"rf", Aux::I32(fragment_span))?;

            let site_score = if fragment_span > 2 {
                0
            } else if min_dist <= 50 {
                2
            } else if min_dist <= 200 {
                1
            } else {
                0
            };
            r.push_aux(b"rs", Aux::I32(site_score))?;
        }
        Ok(())
    }
}

fn annotate_re_records(records: &mut [Record], ctx: &RammapCtx) {
    if let Some(re_map) = ctx.re_map.as_ref() {
        for rec in records {
            let _ = re_map.annotate_record(rec);
        }
    }
}

fn aux_i32(rec: &Record, tag: &[u8; 2]) -> Option<i32> {
    match rec.aux(tag) {
        Ok(Aux::I8(v)) => Some(v as i32),
        Ok(Aux::U8(v)) => Some(v as i32),
        Ok(Aux::I16(v)) => Some(v as i32),
        Ok(Aux::U16(v)) => Some(v as i32),
        Ok(Aux::I32(v)) => Some(v),
        Ok(Aux::U32(v)) => Some(v as i32),
        _ => None,
    }
}

fn aux_float(rec: &Record, tag: &[u8; 2]) -> Option<f64> {
    match rec.aux(tag) {
        Ok(Aux::Float(v)) => Some(v as f64),
        Ok(Aux::I8(v)) => Some(v as f64),
        Ok(Aux::U8(v)) => Some(v as f64),
        Ok(Aux::I16(v)) => Some(v as f64),
        Ok(Aux::U16(v)) => Some(v as f64),
        Ok(Aux::I32(v)) => Some(v as f64),
        Ok(Aux::U32(v)) => Some(v as f64),
        _ => None,
    }
}

fn set_string_aux(rec: &mut Record, tag: &[u8; 2], value: &str) {
    remove_aux_if_present(rec, tag);
    rec.push_aux(tag, Aux::String(value)).ok();
}

fn query_interval_groups(records: &[Record], min_overlap_frac: f32) -> (Vec<Option<(usize, usize)>>, Vec<Vec<usize>>) {
    let intervals: Vec<Option<(usize, usize)>> = records
        .iter()
        .map(|rec| {
            if rec.is_unmapped() {
                None
            } else {
                let q = get_query_coords_on_read(rec);
                (q.1 > q.0).then_some(q)
            }
        })
        .collect();

    let mut groups: Vec<Vec<usize>> = Vec::new();
    for (idx, q_opt) in intervals.iter().enumerate() {
        let Some(q) = q_opt else {
            continue;
        };
        let q_len = q.1 - q.0;
        let mut assigned = false;
        for group in &mut groups {
            let gq = intervals[group[0]].unwrap();
            let g_len = gq.1 - gq.0;
            let overlap = interval_overlap(*q, gq);
            let min_len = std::cmp::min(q_len, g_len).max(1);
            if (overlap as f32 / min_len as f32) >= min_overlap_frac {
                group.push(idx);
                assigned = true;
                break;
            }
        }
        if !assigned {
            groups.push(vec![idx]);
        }
    }

    (intervals, groups)
}

fn adjusted_candidate_score(rec: &Record) -> f64 {
    let mut score = get_score(rec) as f64;
    match aux_i32(rec, b"rs") {
        Some(2) => score += 5.0,
        Some(1) => score += 2.0,
        _ => {}
    }
    if aux_i32(rec, b"rf").is_some_and(|v| v > 2) {
        score -= 5.0;
    }
    if aux_i32(rec, b"gr").is_some_and(|v| v > 0) {
        score -= 5.0;
    }
    score
}

fn ambiguity_label(records: &[Record], group: &[usize], target_names: &[&[u8]]) -> &'static str {
    if group.len() <= 1 {
        return "unique";
    }

    let mut homeolog_group: Option<String> = None;
    let mut suffixes = std::collections::HashSet::new();
    let mut all_homeolog = true;
    for &idx in group {
        if let Some((group_name, suffix)) = record_homeolog_label(&records[idx], target_names) {
            if let Some(prev) = homeolog_group.as_ref() {
                if prev != &group_name {
                    all_homeolog = false;
                    break;
                }
            } else {
                homeolog_group = Some(group_name);
            }
            suffixes.insert(suffix);
        } else {
            all_homeolog = false;
            break;
        }
    }

    if all_homeolog && suffixes.len() > 1 {
        "homeolog"
    } else {
        "multi"
    }
}

fn annotate_candidate_probabilities(records: &mut [Record], ctx: &RammapCtx, target_names: &[&[u8]]) {
    if !ctx.candidate_probability {
        return;
    }

    const TAU: f64 = 5.0;
    let (_, groups) = query_interval_groups(records, 0.5);
    for group in groups {
        if group.is_empty() {
            continue;
        }

        let ac = group.len() as i32;
        let amb = ambiguity_label(records, &group, target_names);
        let scores: Vec<f64> = group
            .iter()
            .map(|&idx| adjusted_candidate_score(&records[idx]))
            .collect();
        let max_score = scores
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let weights: Vec<f64> = scores
            .iter()
            .map(|score| ((*score - max_score) / TAU).exp())
            .collect();
        let weight_sum: f64 = weights.iter().sum();

        for (&idx, weight) in group.iter().zip(weights.iter()) {
            let cp = if weight_sum > 0.0 {
                (*weight / weight_sum) as f32
            } else {
                0.0
            };
            remove_aux_if_present(&mut records[idx], b"cp");
            records[idx].push_aux(b"cp", Aux::Float(cp)).ok();
            set_i32_aux(&mut records[idx], b"ac", ac);
            set_string_aux(&mut records[idx], b"am", amb);
        }
    }
}

fn target_interval_overlap(a: &Record, b: &Record) -> i64 {
    if a.tid() < 0 || b.tid() < 0 || a.tid() != b.tid() {
        return 0;
    }
    let a_start = a.pos();
    let a_end = a.cigar().end_pos();
    let b_start = b.pos();
    let b_end = b.cigar().end_pos();
    std::cmp::min(a_end, b_end).saturating_sub(std::cmp::max(a_start, b_start))
}

fn graph_transition_score(prev: &Record, curr: &Record, target_names: &[&[u8]]) -> f64 {
    if prev.is_unmapped() || curr.is_unmapped() {
        return -50.0;
    }

    let mut score = 0.0;
    if aux_i32(prev, b"rs").unwrap_or(0) > 0 && aux_i32(curr, b"rs").unwrap_or(0) > 0 {
        score += 2.0;
    }

    if prev.tid() == curr.tid() && prev.tid() >= 0 {
        if aux_i32(prev, b"re").is_some()
            && aux_i32(prev, b"re") == aux_i32(curr, b"re")
        {
            score -= 12.0;
        }
        let overlap = target_interval_overlap(prev, curr);
        let prev_len = (prev.cigar().end_pos() - prev.pos()).max(1);
        let curr_len = (curr.cigar().end_pos() - curr.pos()).max(1);
        let min_len = std::cmp::min(prev_len, curr_len).max(1);
        if overlap > 0 && overlap as f64 / min_len as f64 >= 0.5 {
            score -= 25.0;
        }
    }

    if let (Some((prev_group, prev_suffix)), Some((curr_group, curr_suffix))) = (
        record_homeolog_label(prev, target_names),
        record_homeolog_label(curr, target_names),
    ) {
        if prev_group == curr_group && prev_suffix == curr_suffix {
            score += 1.5;
        }
    }

    score
}

fn graph_homeolog_anchor_support(
    records: &[Record],
    layers: &[Vec<usize>],
    intervals: &[Option<(usize, usize)>],
    target_names: &[&[u8]],
) -> (HashMap<(String, char), f64>, Vec<HashMap<(String, char), f64>>) {
    const MIN_ANCHOR_MAPQ: u8 = 10;
    const MIN_ANCHOR_LEN: usize = 300;

    let mut global_support: HashMap<(String, char), f64> = HashMap::new();
    let mut layer_supports: Vec<HashMap<(String, char), f64>> = Vec::with_capacity(layers.len());

    for layer in layers {
        let mut current_support: HashMap<(String, char), f64> = HashMap::new();
        let mut best_score_by_label: HashMap<(String, char), i64> = HashMap::new();
        let mut best_layer_score = i64::MIN;
        let mut second_layer_score = i64::MIN;

        for &idx in layer {
            let score = get_score(&records[idx]);
            if score > best_layer_score {
                second_layer_score = best_layer_score;
                best_layer_score = score;
            } else if score > second_layer_score {
                second_layer_score = score;
            }
            if let Some(label) = record_homeolog_label(&records[idx], target_names) {
                best_score_by_label
                    .entry(label)
                    .and_modify(|v| *v = (*v).max(score))
                    .or_insert(score);
            }
        }

        for &idx in layer {
            let rec = &records[idx];
            if rec.is_unmapped() || rec.mapq() < MIN_ANCHOR_MAPQ {
                continue;
            }
            let Some(label) = record_homeolog_label(rec, target_names) else {
                continue;
            };
            let Some((start, end)) = intervals[idx] else {
                continue;
            };
            let aligned_len = end.saturating_sub(start);
            if aligned_len < MIN_ANCHOR_LEN {
                continue;
            }

            let score = get_score(rec);
            let score_gap = if second_layer_score == i64::MIN {
                20
            } else {
                score.saturating_sub(second_layer_score).max(0)
            };
            let label_gap = best_score_by_label
                .iter()
                .filter_map(|(other_label, other_score)| {
                    (other_label.0 == label.0 && other_label.1 != label.1).then_some(*other_score)
                })
                .max()
                .map(|other_score| score.saturating_sub(other_score).max(0))
                .unwrap_or(score_gap);

            let support = rec.mapq() as f64
                + (aligned_len as f64 / 200.0).min(10.0)
                + (score_gap as f64 / 3.0).min(10.0)
                + (label_gap as f64 / 2.0).min(10.0);
            if support <= 0.0 {
                continue;
            }

            current_support
                .entry(label.clone())
                .and_modify(|v| *v = (*v).max(support))
                .or_insert(support);
        }

        for (label, support) in &current_support {
            *global_support.entry(label.clone()).or_insert(0.0) += *support;
        }
        layer_supports.push(current_support);
    }

    (global_support, layer_supports)
}

fn graph_mapq_for_record(rec: &Record, is_primary: bool) -> u8 {
    let gp = aux_float(rec, b"gp").unwrap_or(0.0).clamp(0.0, 0.9999);
    let gh = aux_i32(rec, b"gh").unwrap_or(0);
    let q = get_query_coords_on_read(rec);
    let aligned_len = q.1.saturating_sub(q.0);

    let mut mq = if gp >= 0.999 {
        40
    } else if gp >= 0.99 {
        30
    } else if gp >= 0.95 {
        20
    } else if gp >= 0.80 {
        10
    } else if gp >= 0.65 && (aligned_len >= 300 || gh >= 10) {
        7
    } else if gp >= 0.50 && (aligned_len >= 500 || gh >= 20) {
        5
    } else if gp >= 0.35 && (aligned_len >= 800 || gh >= 30) {
        3
    } else if gp >= 0.20 && gh >= 40 {
        2
    } else if gp >= 0.50 && aligned_len >= 200 {
        1
    } else {
        0
    };

    if mq > 0 && gh >= 40 {
        mq += 8;
    } else if mq > 0 && gh >= 20 {
        mq += 5;
    } else if mq > 0 && gh >= 10 {
        mq += 2;
    }

    if aligned_len < 100 {
        mq = mq.min(5);
    } else if aligned_len < 200 {
        mq = mq.min(10);
    } else if aligned_len < 500 {
        mq = mq.min(20);
    }

    if is_primary {
        mq.min(60)
    } else {
        mq.min(40)
    }
}

fn graph_protected_layer_assignment(
    records: &[Record],
    layer: &[usize],
    target_names: &[&[u8]],
) -> Option<(usize, (String, char))> {
    layer
        .iter()
        .filter_map(|&idx| {
            let rec = &records[idx];
            if rec.is_unmapped() || rec.is_secondary() || rec.mapq() <= 1 {
                return None;
            }
            let label = record_homeolog_label(rec, target_names)?;
            Some((idx, label))
        })
        .max_by(|(a_idx, _), (b_idx, _)| {
            records[*a_idx]
                .mapq()
                .cmp(&records[*b_idx].mapq())
                .then_with(|| get_score(&records[*a_idx]).cmp(&get_score(&records[*b_idx])))
        })
}

fn graph_preserved_mapq(old_mapq: u8) -> u8 {
    if old_mapq > 1 {
        old_mapq
    } else {
        0
    }
}

fn annotate_graph_fragment_assignment(
    records: &mut [Record],
    ctx: &RammapCtx,
    target_names: &[&[u8]],
) {
    if !ctx.graph_assignment {
        return;
    }

    let (intervals, groups) = query_interval_groups(records, 0.5);
    if groups.is_empty() {
        return;
    }

    let mut layers: Vec<Vec<usize>> = groups
        .into_iter()
        .filter(|group| !group.is_empty())
        .collect();
    layers.sort_by_key(|group| {
        intervals[group[0]]
            .map(|(start, end)| (start + end) / 2)
            .unwrap_or(usize::MAX)
    });

    for layer in &layers {
        let gc = layer.len().min(i32::MAX as usize) as i32;
        for &idx in layer {
            set_i32_aux(&mut records[idx], b"gc", gc);
            set_i32_aux(&mut records[idx], b"ga", 0);
            remove_aux_if_present(&mut records[idx], b"gp");
            remove_aux_if_present(&mut records[idx], b"gs");
        }
    }

    if layers.is_empty() {
        return;
    }

    let (homeolog_support, layer_homeolog_supports) =
        graph_homeolog_anchor_support(records, &layers, &intervals, target_names);

    let mut dp: Vec<Vec<f64>> = Vec::with_capacity(layers.len());
    let mut back: Vec<Vec<Option<usize>>> = Vec::with_capacity(layers.len());
    let protected_assignments: Vec<Option<(usize, (String, char))>> = layers
        .iter()
        .map(|layer| graph_protected_layer_assignment(records, layer, target_names))
        .collect();

    for (layer_idx, layer) in layers.iter().enumerate() {
        let mut layer_scores = Vec::with_capacity(layer.len());
        let mut layer_back = vec![None; layer.len()];

        for (cand_idx, &rec_idx) in layer.iter().enumerate() {
            let mut node_score = adjusted_candidate_score(&records[rec_idx]);
            let mut homeolog_prior = 0.0;
            if let Some(label) = record_homeolog_label(&records[rec_idx], target_names) {
                let protected_label = protected_assignments
                    .get(layer_idx)
                    .and_then(|x| x.as_ref())
                    .map(|(_, label)| label);
                let allow_homeolog_prior = match protected_label {
                    Some(protected) => protected == &label,
                    None => true,
                };
                if allow_homeolog_prior {
                    homeolog_prior = homeolog_support.get(&label).copied().unwrap_or(0.0);
                    if let Some(self_support) = layer_homeolog_supports
                        .get(layer_idx)
                        .and_then(|support| support.get(&label))
                    {
                        homeolog_prior = (homeolog_prior - *self_support).max(0.0);
                    }
                    node_score += homeolog_prior;
                }
            }
            set_i32_aux(
                &mut records[rec_idx],
                b"gh",
                homeolog_prior.round().min(i32::MAX as f64) as i32,
            );
            if layer_idx == 0 {
                layer_scores.push(node_score);
                continue;
            }

            let prev_layer = &layers[layer_idx - 1];
            let mut best_score = f64::NEG_INFINITY;
            let mut best_prev = None;
            for (prev_idx, &prev_rec_idx) in prev_layer.iter().enumerate() {
                let transition =
                    graph_transition_score(&records[prev_rec_idx], &records[rec_idx], target_names);
                let score = dp[layer_idx - 1][prev_idx] + transition + node_score;
                if score > best_score {
                    best_score = score;
                    best_prev = Some(prev_idx);
                }
            }
            layer_scores.push(best_score);
            layer_back[cand_idx] = best_prev;
        }

        dp.push(layer_scores);
        back.push(layer_back);
    }

    let mut best_path: Vec<usize> = vec![0; layers.len()];
    let mut best_last = dp
        .last()
        .and_then(|scores| {
            scores
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(idx, _)| idx)
        })
        .unwrap_or(0);

    for layer_idx in (0..layers.len()).rev() {
        best_path[layer_idx] = best_last;
        if layer_idx > 0 {
            best_last = back[layer_idx][best_last].unwrap_or(0);
        }
    }

    for (layer_idx, protected) in protected_assignments.iter().enumerate() {
        let Some((protected_idx, _)) = protected else {
            continue;
        };
        if let Some(pos) = layers[layer_idx].iter().position(|idx| idx == protected_idx) {
            best_path[layer_idx] = pos;
        }
    }

    const TAU: f64 = 8.0;
    for (layer_idx, layer) in layers.iter().enumerate() {
        let max_score = dp[layer_idx]
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let weights: Vec<f64> = dp[layer_idx]
            .iter()
            .map(|score| ((*score - max_score) / TAU).exp())
            .collect();
        let weight_sum: f64 = weights.iter().sum();

        for (cand_idx, &rec_idx) in layer.iter().enumerate() {
            let gp = if weight_sum > 0.0 {
                (weights[cand_idx] / weight_sum) as f32
            } else {
                0.0
            };
            let graph_score = dp[layer_idx][cand_idx].round() as i32;
            set_i32_aux(
                &mut records[rec_idx],
                b"ga",
                if cand_idx == best_path[layer_idx] { 1 } else { 0 },
            );
            set_i32_aux(&mut records[rec_idx], b"gs", graph_score);
            records[rec_idx].push_aux(b"gp", Aux::Float(gp)).ok();
        }
    }

    let selected: Vec<usize> = layers
        .iter()
        .enumerate()
        .filter_map(|(layer_idx, layer)| layer.get(best_path[layer_idx]).copied())
        .collect();
    if selected.is_empty() {
        return;
    }

    let primary_idx = selected
        .iter()
        .copied()
        .max_by(|&a, &b| {
            records[a]
                .mapq()
                .cmp(&records[b].mapq())
                .then_with(|| get_score(&records[a]).cmp(&get_score(&records[b])))
                .then_with(|| {
                    let qa = get_query_coords_on_read(&records[a]);
                    let qb = get_query_coords_on_read(&records[b]);
                    qa.1.saturating_sub(qa.0).cmp(&qb.1.saturating_sub(qb.0))
                })
        })
        .unwrap_or(selected[0]);

    for (idx, rec) in records.iter_mut().enumerate() {
        if rec.is_unmapped() {
            continue;
        }

        let mut flags = rec.flags();
        if idx == primary_idx {
            flags &= !0x100;
            flags &= !0x800;
            let old_mapq = rec.mapq();
            let graph_mapq = graph_mapq_for_record(rec, true);
            set_i32_aux(rec, b"gm", old_mapq as i32);
            rec.set_mapq(graph_preserved_mapq(old_mapq).max(graph_mapq));
            set_i32_aux(rec, b"gf", 1);
        } else if selected.contains(&idx) {
            flags &= !0x100;
            flags |= 0x800;
            let old_mapq = rec.mapq();
            let graph_mapq = graph_mapq_for_record(rec, false);
            set_i32_aux(rec, b"gm", old_mapq as i32);
            rec.set_mapq(graph_preserved_mapq(old_mapq).max(graph_mapq));
            set_i32_aux(rec, b"gf", 2);
        } else if aux_i32(rec, b"ga") == Some(0) {
            flags |= 0x100;
            flags &= !0x800;
            let old_mapq = rec.mapq();
            set_i32_aux(rec, b"gm", old_mapq as i32);
            rec.set_mapq(0);
            set_i32_aux(rec, b"gf", 0);
        }
        rec.set_flags(flags);
    }
}

fn recalibrate_mapq_records(records: &mut [Record], ctx: &RammapCtx) {
    if !ctx.mapq_calibrate {
        return;
    }

    let (intervals, groups) = query_interval_groups(records, 0.5);

    let mut new_mapqs: Vec<Option<u8>> = vec![None; records.len()];
    for group in groups {
        let mut ranked: Vec<(usize, i64)> = group
            .iter()
            .map(|&idx| (idx, get_score(&records[idx])))
            .collect();
        ranked.sort_by_key(|&(_, score)| std::cmp::Reverse(score));
        let best_score = ranked.first().map(|(_, score)| *score).unwrap_or(0);
        let second_score = ranked.get(1).map(|(_, score)| *score);

        for &(idx, score) in &ranked {
            let rec = &records[idx];
            let original_mapq = rec.mapq() as i32;
            let mut mq = if let Some(second) = second_score {
                let gap = score.saturating_sub(second);
                let ratio = if second > 0 {
                    score as f32 / second as f32
                } else {
                    99.0
                };
                if gap >= 40 || ratio >= 1.15 {
                    40
                } else if gap >= 25 || ratio >= 1.10 {
                    30
                } else if gap >= 15 || ratio >= 1.05 {
                    20
                } else if gap >= 8 {
                    10
                } else {
                    1
                }
            } else {
                (original_mapq + 10).min(60)
            };

            if score < best_score {
                mq = mq.min(10);
            }

            let q = intervals[idx].unwrap();
            let aligned_len = q.1.saturating_sub(q.0);
            if aligned_len < 80 {
                mq -= 10;
            } else if aligned_len >= 200 {
                mq += 5;
            }

            match aux_i32(rec, b"rs") {
                Some(2) => mq += 5,
                Some(1) => mq += 2,
                Some(0) if original_mapq <= 1 => mq -= 5,
                _ => {}
            }
            if aux_i32(rec, b"rf").is_some_and(|v| v > 2) {
                mq -= 10;
            }
            if aux_i32(rec, b"gr").is_some_and(|v| v > 0) {
                mq -= 10;
                mq = mq.min(10);
            }

            new_mapqs[idx] = Some(mq.clamp(0, 60) as u8);
        }
    }

    for (rec, mq_opt) in records.iter_mut().zip(new_mapqs.into_iter()) {
        if let Some(mq) = mq_opt {
            remove_aux_if_present(rec, b"om");
            rec.push_aux(b"om", Aux::I32(rec.mapq() as i32)).ok();
            rec.set_mapq(mq);
        }
    }
}

pub(crate) fn refresh_sa_tags(records: &mut [Record], target_names: &[&[u8]]) {
    let selected: Vec<(usize, String)> = records.iter().enumerate()
        .filter(|(_, r)| !r.is_unmapped() && !r.is_secondary())
        .filter_map(|(i, r)| {
            let name = std::str::from_utf8(target_names.get(r.tid() as usize)?).ok()?;
            Some((i, format!("{},{},{},{},{},{};", name, r.pos() + 1,
                if r.is_reverse() { '-' } else { '+' }, r.cigar(), r.mapq(),
                aux_i32(r, b"NM").unwrap_or(0))))
        }).collect();
    for (i, rec) in records.iter_mut().enumerate() {
        remove_aux_if_present(rec, b"tp");
        if !rec.is_unmapped() {
            rec.push_aux(b"tp", Aux::Char(if rec.is_secondary() { b'S' } else { b'P' })).ok();
        }
        remove_aux_if_present(rec, b"SA");
        if rec.is_unmapped() || rec.is_secondary() { continue; }
        let sa: String = selected.iter().filter(|(j, _)| *j != i).map(|(_, s)| s.as_str()).collect();
        if !sa.is_empty() { rec.push_aux(b"SA", Aux::String(&sa)).ok(); }
    }
}

fn finalize_records(records: &mut [Record], ctx: &RammapCtx, target_names: &[&[u8]]) {
    annotate_re_records(records, ctx);
    annotate_candidate_probabilities(records, ctx, target_names);
    recalibrate_mapq_records(records, ctx);
    annotate_graph_fragment_assignment(records, ctx, target_names);
    // Refresh SA only after the final primary CIGAR/sequence restoration.
}

fn record_to_paf_fast(r: &Record, target_names: &[&[u8]], target_lens: &[u64]) -> Option<String> {
    let tid = r.tid();
    if tid < 0 {
        return None;
    }
    let qname = std::str::from_utf8(r.qname()).unwrap_or("*");

    let cigar = r.cigar();
    let mut qlen = 0;
    let mut qstart = 0;
    let mut aligned_qlen = 0;
    let mut block_len = 0;

    for (i, op) in cigar.iter().enumerate() {
        match op {
            Cigar::SoftClip(n) | Cigar::HardClip(n) => {
                qlen += *n as usize;
                if i == 0 {
                    qstart = *n as usize;
                }
            }
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) => {
                qlen += *n as usize;
                aligned_qlen += *n as usize;
                block_len += *n as usize;
            }
            Cigar::Del(n) => {
                block_len += *n as usize;
            }
            _ => {}
        }
    }
    let qend = qstart + aligned_qlen;
    let is_rev = r.is_reverse();
    let strand = if is_rev { '-' } else { '+' };
    let (paf_qstart, paf_qend) = if is_rev {
        (qlen - qend, qlen - qstart)
    } else {
        (qstart, qend)
    };

    let tname_bytes = target_names.get(tid as usize)?;
    let tname = std::str::from_utf8(tname_bytes).ok()?;
    let tlen = target_lens.get(tid as usize)?;
    let tstart = r.pos() as usize;
    let tend = cigar.end_pos() as usize;

    let nm = match r.aux(b"NM") {
        Ok(Aux::I8(v)) => v as usize,
        Ok(Aux::U8(v)) => v as usize,
        Ok(Aux::I16(v)) => v as usize,
        Ok(Aux::U16(v)) => v as usize,
        Ok(Aux::I32(v)) => v as usize,
        Ok(Aux::U32(v)) => v as usize,
        _ => 0,
    };
    let matches = block_len.saturating_sub(nm);
    let mapq = r.mapq();

    use std::fmt::Write as _;
    let mut cg_str = String::with_capacity(cigar.len() * 4);
    for op in cigar.iter() {
        match op {
            Cigar::Match(n) => {
                let _ = write!(cg_str, "{}M", n);
            }
            Cigar::Ins(n) => {
                let _ = write!(cg_str, "{}I", n);
            }
            Cigar::Del(n) => {
                let _ = write!(cg_str, "{}D", n);
            }
            Cigar::RefSkip(n) => {
                let _ = write!(cg_str, "{}N", n);
            }
            Cigar::SoftClip(n) => {
                let _ = write!(cg_str, "{}S", n);
            }
            Cigar::HardClip(n) => {
                let _ = write!(cg_str, "{}H", n);
            }
            Cigar::Pad(n) => {
                let _ = write!(cg_str, "{}P", n);
            }
            Cigar::Equal(n) => {
                let _ = write!(cg_str, "{}=", n);
            }
            Cigar::Diff(n) => {
                let _ = write!(cg_str, "{}X", n);
            }
        }
    }
    let mut paf_line = String::with_capacity(256 + cg_str.len());
    let _ = write!(
        paf_line,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        qname,
        qlen,
        paf_qstart,
        paf_qend,
        strand,
        tname,
        tlen,
        tstart,
        tend,
        matches,
        block_len,
        mapq
    );

    let _ = write!(paf_line, "\tNM:i:{}", nm);
    if let Ok(ascore) = r.aux(b"AS") {
        let val = match ascore {
            Aux::I8(v) => v as i32,
            Aux::U8(v) => v as i32,
            Aux::I16(v) => v as i32,
            Aux::U16(v) => v as i32,
            Aux::I32(v) => v,
            Aux::U32(v) => v as i32,
            _ => 0,
        };
        let _ = write!(paf_line, "\tAS:i:{}", val);
    }

    for tag in [b"s0", b"m0", b"MA"] {
        if let Some(value) = aux_i32(r, tag) {
            let _ = write!(paf_line, "\t{}:i:{}", std::str::from_utf8(tag).unwrap(), value);
        }
    }
    match aux_i32(r, b"gf") {
        Some(0) => paf_line.push_str("\ttp:A:S"),
        Some(1) | Some(2) => paf_line.push_str("\ttp:A:P"),
        _ => {
            let is_secondary = (r.flags() & 0x100) != 0;
            if is_secondary {
                paf_line.push_str("\ttp:A:S");
            } else {
                paf_line.push_str("\ttp:A:P");
            }
        }
    }

    if !cg_str.is_empty() {
        let _ = write!(paf_line, "\tcg:Z:{}", cg_str);
    }
    if let Ok(Aux::String(cs)) = r.aux(b"cs") {
        let _ = write!(paf_line, "\tcs:Z:{}", cs);
    }

    if let Ok(Aux::I32(re)) = r.aux(b"re") {
        let _ = write!(paf_line, "\tre:i:{}", re);
    }
    if let Ok(Aux::I32(rd)) = r.aux(b"rd") {
        let _ = write!(paf_line, "\trd:i:{}", rd);
    }
    if let Ok(Aux::I32(rf)) = r.aux(b"rf") {
        let _ = write!(paf_line, "\trf:i:{}", rf);
    }
    if let Ok(Aux::I32(rs)) = r.aux(b"rs") {
        let _ = write!(paf_line, "\trs:i:{}", rs);
    }
    if let Ok(Aux::I32(om)) = r.aux(b"om") {
        let _ = write!(paf_line, "\tom:i:{}", om);
    }
    if let Ok(Aux::I32(gr)) = r.aux(b"gr") {
        let _ = write!(paf_line, "\tgr:i:{}", gr);
    }
    if let Ok(Aux::I32(hr)) = r.aux(b"hr") {
        let _ = write!(paf_line, "\thr:i:{}", hr);
    }
    if let Ok(Aux::I32(hs)) = r.aux(b"hs") {
        let _ = write!(paf_line, "\ths:i:{}", hs);
    }
    if let Ok(Aux::I32(zr)) = r.aux(b"zr") {
        let _ = write!(paf_line, "\tzr:i:{}", zr);
    }
    if let Ok(Aux::I32(zs)) = r.aux(b"zs") {
        let _ = write!(paf_line, "\tzs:i:{}", zs);
    }
    if let Ok(Aux::I32(za)) = r.aux(b"za") {
        let _ = write!(paf_line, "\tza:i:{}", za);
    }
    if let Ok(Aux::Float(cp)) = r.aux(b"cp") {
        let _ = write!(paf_line, "\tcp:f:{:.4}", cp);
    }
    if let Ok(Aux::I32(ac)) = r.aux(b"ac") {
        let _ = write!(paf_line, "\tac:i:{}", ac);
    }
    if let Ok(Aux::String(am)) = r.aux(b"am") {
        let _ = write!(paf_line, "\tam:Z:{}", am);
    }
    if let Ok(Aux::I32(ga)) = r.aux(b"ga") {
        let _ = write!(paf_line, "\tga:i:{}", ga);
    }
    if let Ok(Aux::Float(gp)) = r.aux(b"gp") {
        let _ = write!(paf_line, "\tgp:f:{:.4}", gp);
    }
    if let Ok(Aux::I32(gs)) = r.aux(b"gs") {
        let _ = write!(paf_line, "\tgs:i:{}", gs);
    }
    if let Ok(Aux::I32(gc)) = r.aux(b"gc") {
        let _ = write!(paf_line, "\tgc:i:{}", gc);
    }
    if let Ok(Aux::I32(gh)) = r.aux(b"gh") {
        let _ = write!(paf_line, "\tgh:i:{}", gh);
    }
    if let Ok(Aux::I32(gf)) = r.aux(b"gf") {
        let _ = write!(paf_line, "\tgf:i:{}", gf);
    }
    if let Ok(Aux::I32(gm)) = r.aux(b"gm") {
        let _ = write!(paf_line, "\tgm:i:{}", gm);
    }

    Some(paf_line)
}

fn map_raw_sequence_to_records(
    qname: &[u8],
    seq_buf: &[u8],
    qual_buf: &[u8],
    original_record: Option<&Record>,
    index: &Index,
    opt: &MapOptions,
    out_cfg: &OutputConfig,
    soft_clip: bool,
    header: &Header,
    timings: &AlignTimings,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> AnyResult<Vec<Record>> {
    let mm_opt: Option<String> = original_record.and_then(|rec| {
        rec.aux(b"MM").ok().and_then(|a| match a {
            Aux::String(s) => Some(s.to_string()),
            _ => None,
        })
    });

    let ml_opt: Option<Vec<u8>> = original_record.and_then(|rec| {
        rec.aux(b"ML").ok().and_then(|a| match a {
            Aux::ArrayU8(arr) => Some(arr.iter().collect::<Vec<u8>>()),
            _ => None,
        })
    });

    let qname_str = std::str::from_utf8(qname).unwrap_or("*");

    let qual_ascii: Option<Vec<u8>> = if qual_buf.iter().all(|&q| q == 255) {
        None
    } else if original_record.is_some() {
        Some(qual_buf.iter().map(|&q| q.saturating_add(33)).collect())
    } else {
        Some(qual_buf.to_vec())
    };
    let qual_str = qual_ascii
        .as_ref()
        .and_then(|q| std::str::from_utf8(q).ok());

    let read_info = ReadInfo {
        qname: qname_str,
        qseq: seq_buf,
        qual: qual_str,
        comment: None,
        n_seg: 1,
        seg_idx: 0,
    };

    let (sam_block, _stats) = timings.measure(10, || align_and_format_query(
        opt,
        index,
        &read_info,
        aln_ctx,
        &mut map_ctx.mapping,
        None,
        None,
        out_cfg,
    ));

    let local_header_view = map_ctx.header.get_or_insert_with(||
        timings.measure(11, || HeaderView::from_header(header)));
    let mut output_records = Vec::new();
    for line in sam_block.lines() {
        if line.is_empty() || line.as_bytes()[0] == b'@' {
            continue;
        }

        let mut rec = timings.measure(12, || Record::from_sam(local_header_view, line.as_bytes()))?;

        if let Some(orig_rec) = original_record {
            copy_all_aux_except(
                orig_rec,
                &mut rec,
                &[
                    b"NM", b"AS", b"tp", b"MM", b"ML", b"SA", b"MD", b"cs", b"cg",
                ],
            )?;
        }

        if (!rec.is_secondary() && !rec.is_supplementary()) || soft_clip {
            if let Some(mm) = mm_opt.as_ref() {
                rec.push_aux(b"MM", Aux::String(mm))?;
            }
            if let Some(ml) = ml_opt.as_ref() {
                rec.push_aux(b"ML", Aux::ArrayU8((&ml[..]).into()))?;
            }
        }

        output_records.push(rec);
    }

    Ok(output_records)
}

fn interval_overlap(a: (usize, usize), b: (usize, usize)) -> usize {
    std::cmp::min(a.1, b.1).saturating_sub(std::cmp::max(a.0, b.0))
}

fn uncovered_query_gaps(records: &[Record], read_len: usize, min_gap_len: usize) -> Vec<(usize, usize)> {
    let mut intervals: Vec<(usize, usize)> = records
        .iter()
        .filter(|r| !r.is_unmapped() && !r.is_secondary())
        .map(get_query_coords_on_read)
        .filter(|(s, e)| e > s)
        .collect();
    intervals.sort_unstable();

    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (s, e) in intervals {
        if let Some(last) = merged.last_mut() {
            if s <= last.1 {
                last.1 = last.1.max(e);
                continue;
            }
        }
        merged.push((s, e));
    }

    let mut gaps = Vec::new();
    let mut cursor = 0usize;
    for (s, e) in merged {
        if s > cursor && s - cursor >= min_gap_len {
            gaps.push((cursor, s));
        }
        cursor = cursor.max(e);
    }
    if read_len > cursor && read_len - cursor >= min_gap_len {
        gaps.push((cursor, read_len));
    }
    gaps
}

fn retain_gap_rescue_records(
    rescued: Vec<Record>,
    primary_records: &[Record],
    gaps: &[(usize, usize)],
    min_gap_len: usize,
) -> Vec<Record> {
    let existing: Vec<(usize, usize)> = primary_records
        .iter()
        .filter(|r| !r.is_unmapped() && !r.is_secondary())
        .map(get_query_coords_on_read)
        .filter(|(s, e)| e > s)
        .collect();
    let min_overlap = std::cmp::min(20, min_gap_len.max(1));

    rescued
        .into_iter()
        .filter_map(|mut rec| {
            if rec.is_unmapped() {
                return None;
            }
            let q = get_query_coords_on_read(&rec);
            if q.1 <= q.0 {
                return None;
            }
            let gap_overlap = gaps.iter().map(|&g| interval_overlap(q, g)).max().unwrap_or(0);
            if gap_overlap < min_overlap {
                return None;
            }
            let q_len = q.1 - q.0;
            let existing_overlap = existing
                .iter()
                .map(|&iv| interval_overlap(q, iv))
                .max()
                .unwrap_or(0);
            if existing_overlap * 2 > q_len {
                return None;
            }
            let flags = rec.flags();
            rec.set_flags(if flags & 0x100 != 0 {
                flags & !0x800
            } else {
                flags | 0x800
            });
            remove_aux_if_present(&mut rec, b"gr");
            rec.push_aux(b"gr", Aux::I32(1)).ok();
            Some(rec)
        })
        .collect()
}

fn project_gap_record_to_read(
    rec: &mut Record,
    qname: &[u8],
    gap_start: usize,
    gap_end: usize,
    read_len: usize,
) {
    let is_rev = rec.is_reverse();
    let left_clip = if is_rev {
        read_len.saturating_sub(gap_end)
    } else {
        gap_start
    };
    let right_clip = if is_rev {
        gap_start
    } else {
        read_len.saturating_sub(gap_end)
    };

    let mut ops: Vec<Cigar> = Vec::new();
    if left_clip > 0 {
        ops.push(Cigar::HardClip(left_clip as u32));
    }
    ops.extend(rec.cigar().iter().cloned());
    if right_clip > 0 {
        ops.push(Cigar::HardClip(right_clip as u32));
    }
    let cigar = CigarString(ops);
    let seq = rec.seq().as_bytes();
    let qual = rec.qual().to_vec();
    rec.set(qname, Some(&cigar), &seq, &qual);
}

fn apply_gap_rescue_to_records(
    output_records: &mut Vec<Record>,
    qname: &[u8],
    seq_buf: &[u8],
    qual_buf: &[u8],
    original_record: Option<&Record>,
    ctx: &RammapCtx,
    header: &Header,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> AnyResult<Vec<Record>> {
    if !ctx.gap_rescue {
        return Ok(Vec::new());
    }

    let (Some(sensitive_index), Some(sensitive_opt)) =
        (ctx.sensitive_index.as_ref(), ctx.sensitive_opt.as_ref())
    else {
        return Ok(Vec::new());
    };

    let gaps = uncovered_query_gaps(output_records, seq_buf.len(), ctx.min_gap_len);
    if gaps.is_empty() {
        return Ok(Vec::new());
    }

    let mut rescued_all = Vec::new();
    for &(gap_start, gap_end) in &gaps {
        let gap_seq = &seq_buf[gap_start..gap_end];
        let gap_qual = &qual_buf[gap_start..gap_end];
        let rescued = map_raw_sequence_to_records(
            qname,
            gap_seq,
            gap_qual,
            original_record,
            sensitive_index,
            sensitive_opt,
            &ctx.out_cfg,
            ctx.soft_clip,
            header,
            &ctx.timings,
            aln_ctx,
            map_ctx,
        )?;
        for mut rec in rescued {
            project_gap_record_to_read(&mut rec, qname, gap_start, gap_end, seq_buf.len());
            rescued_all.push(rec);
        }
    }

    Ok(retain_gap_rescue_records(
        rescued_all,
        output_records,
        &gaps,
        ctx.min_gap_len,
    ))
}

fn process_raw_sequence(
    qname: &[u8],
    seq_buf: &[u8],
    qual_buf: &[u8],
    original_record: Option<&Record>,
    ctx: &RammapCtx,
    header: &Header,
    preset_str: &str,
    gap_rescue: bool,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> AnyResult<Vec<Record>> {
    log::debug!(
        "secondary: {}, soft_clip: {}, preset: {}, best_n: {}, gap_rescue: {}",
        ctx.secondary,
        ctx.soft_clip,
        preset_str,
        ctx.opt.filtering.best_n,
        ctx.gap_rescue
    );

    let mut output_records = map_raw_sequence_to_records(
        qname,
        seq_buf,
        qual_buf,
        original_record,
        &ctx.index,
        &ctx.opt,
        &ctx.out_cfg,
        ctx.soft_clip,
        header,
        &ctx.timings,
        aln_ctx,
        map_ctx,
    )?;

    if gap_rescue {
        let retained = apply_gap_rescue_to_records(
            &mut output_records,
            qname,
            seq_buf,
            qual_buf,
            original_record,
            ctx,
            header,
            aln_ctx,
            map_ctx,
        )?;
        output_records.extend(retained);
    }

    Ok(output_records)
}

fn format_raw_sequence_paf(
    qname: &[u8],
    seq_buf: &[u8],
    ctx: &RammapCtx,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> String {
    let qname_str = std::str::from_utf8(qname).unwrap_or("*");
    let read_info = ReadInfo {
        qname: qname_str,
        qseq: seq_buf,
        qual: None,
        comment: None,
        n_seg: 1,
        seg_idx: 0,
    };
    let mut out_cfg = ctx.out_cfg.clone();
    out_cfg.output_sam = false;
    out_cfg.do_cigar = ctx.output_cigar || ctx.do_cs;
    out_cfg.do_cs = ctx.do_cs;
    out_cfg.cs_long = ctx.cs_long;
    out_cfg.eqx = ctx.eqx;
    let (mut paf_block, _stats) = align_and_format_query(
        &ctx.opt,
        &ctx.index,
        &read_info,
        aln_ctx,
        &mut map_ctx.mapping,
        None,
        None,
        &out_cfg,
    );
    if !paf_block.is_empty() && !paf_block.ends_with('\n') {
        paf_block.push('\n');
    }
    paf_block
}

fn process_single_read(
    original_record: &Record,
    ctx: &RammapCtx,
    header: &Header,
    preset_str: &str,
    gap_rescue: bool,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> AnyResult<Vec<Record>> {
    process_raw_sequence(
        original_record.qname(),
        &original_record.seq().as_bytes(),
        &original_record.qual(),
        Some(original_record),
        ctx,
        header,
        preset_str,
        gap_rescue,
        aln_ctx,
        map_ctx,
    )
}

fn realign_records(
    records: Vec<Record>,
    keep_secondary: bool,
    mapq_rescue: Option<u8>,
    target_names_ref: &[&[u8]],
    rescue_mode: Option<&str>,
) -> Vec<Record> {
    let rescue_mapq = mapq_rescue.unwrap_or(10);
    let mut au = AlignmentUnit::from_records(records);
    match rescue_mode {
        Some("zero-homeolog") | Some("zero-homoeolog") | Some("zero") => {
            au.rescue_zero_mapq_homeolog(target_names_ref)
        }
        Some("homeolog2") | Some("homoeolog2") => {
            au.rescue_homeolog_interval(rescue_mapq, target_names_ref)
        }
        Some("homeolog") | Some("homoeolog") => {
            au.rescue_homeolog_tie(rescue_mapq, target_names_ref)
        }
        Some("sensitive") | Some("rescue") => au.rescue(rescue_mapq, target_names_ref),
        _ => au.rescue_robust(rescue_mapq, target_names_ref),
    }
    au.into_records(keep_secondary)
}

fn get_query_coords(r: &Record) -> (usize, usize) {
    let cigar = r.cigar();
    if cigar.is_empty() {
        return (0, 0);
    }
    let mut leading_clip = 0;
    let mut aligned_qlen = 0;
    let mut seen_aligned = false;
    for op in cigar.iter() {
        match op {
            Cigar::SoftClip(n) | Cigar::HardClip(n) => {
                if !seen_aligned {
                    leading_clip += *n as usize;
                }
            }
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Ins(n) => {
                seen_aligned = true;
                aligned_qlen += *n as usize;
            }
            _ => {}
        }
    }
    (leading_clip, leading_clip + aligned_qlen)
}

fn query_len_from_cigar(r: &Record) -> usize {
    r.cigar()
        .iter()
        .map(|op| match op {
            Cigar::SoftClip(n)
            | Cigar::HardClip(n)
            | Cigar::Match(n)
            | Cigar::Equal(n)
            | Cigar::Diff(n)
            | Cigar::Ins(n) => *n as usize,
            _ => 0,
        })
        .sum()
}

fn get_query_coords_on_read(r: &Record) -> (usize, usize) {
    let (start, end) = get_query_coords(r);
    if !r.is_reverse() {
        return (start, end);
    }

    let qlen = query_len_from_cigar(r);
    (qlen.saturating_sub(end), qlen.saturating_sub(start))
}

enum WorkerInput {
    SeBam(Vec<Record>),
    SeFastx(Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>),
    PeFastx(Vec<(Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>)>),
}

enum OutputBatch {
    Paf(String),
    Bam(Vec<Record>),
}

fn log_input_read_progress(count: &mut usize, next_log_count: &mut usize, added: usize) {
    *count += added;
    while *count >= *next_log_count {
        log::info!("Processed {} input reads...", *next_log_count);
        *next_log_count += 10_000;
    }
}

struct QueryBaseBatchLimit {
    current_bases: usize,
    max_bases: usize,
    min_records: usize,
}

impl QueryBaseBatchLimit {
    fn new(max_bases: usize, min_records: usize) -> Self {
        Self {
            current_bases: 0,
            max_bases: max_bases.max(1),
            min_records: min_records.max(1),
        }
    }

    fn push(&mut self, query_bases: usize, record_count: usize) -> bool {
        self.current_bases = self.current_bases.saturating_add(query_bases);
        record_count >= self.min_records && self.current_bases >= self.max_bases
    }

    fn reset(&mut self) {
        self.current_bases = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batch_limit_tracks_query_bases_and_forces_progress() {
        let mut limit = QueryBaseBatchLimit::new(100, 10);

        assert!(!limit.push(30, 1));
        assert!(!limit.push(60, 9));
        assert!(limit.push(10, 10));

        limit.reset();

        assert!(!limit.push(0, 1));
        assert!(limit.push(100, 10));
    }
}

#[derive(Debug, Clone)]
pub struct AlignPlusConfig {
    pub methylation: Option<MethConfig>,
    pub reference: String,
    pub queries: Vec<String>,
    pub output: String,
    pub output_bam_format: bool,
    pub output_cigar: bool,
    pub eqx: bool,
    pub cs: Option<String>,
    pub kmer: Option<i16>,
    pub window: Option<i16>,
    pub is_hpc: bool,
    pub batch_size: u64,
    pub soft_clip: bool,
    pub secondary: bool,
    pub mid_occ_frac: Option<(f32, Option<f32>)>,
    pub bounds_of_occurrence: Option<(i32, Option<i32>)>,
    pub max_gap: Option<i32>,
    pub max_gap_ref: Option<i32>,
    pub max_frag_len: Option<i32>,
    pub mask_level: Option<f32>,
    pub bw: Option<(i32, Option<i32>)>,
    pub min_cnt: Option<i32>,
    pub min_chain_score: Option<i32>,
    pub matching_score: Option<i32>,
    pub mismatch_penalty: Option<i32>,
    pub gap_open: Option<(i32, Option<i32>)>,
    pub gap_extension: Option<(i32, Option<i32>)>,
    pub z_drop: Option<(i32, Option<i32>)>,
    pub min_dp_max: Option<i32>,
    pub best_n: Option<i32>,
    pub pri_ratio: Option<f32>,
    pub max_qlen: Option<i32>,
    pub mini_batch_size: i64,
    pub seed: Option<i32>,
    pub porec_gap_rescue: bool,
    pub rescue_k: Option<i16>,
    pub rescue_w: Option<i16>,
    pub min_gap_len: usize,
    pub preset: String,
    pub threads: usize,
    pub realign: bool,
    pub rescue_mode: Option<String>,
    pub mapq_rescue: Option<u8>,
    pub porec_mapq_calibrate: bool,
    pub candidate_probability: bool,
    pub graph_assignment: bool,
    pub re_site: Option<String>,
}

pub(crate) type PrimaryIdentity = (i32, i64, bool, CigarString);

/// A promoted secondary may have SEQ="*". Restore a valid primary record
/// from the original molecule; MM/ML remain relative to its original orientation.
pub(crate) fn restore_primary_read(records: &mut [Record], original: &Record, initial: Option<&PrimaryIdentity>) {
    // The original mapped primary already carries the full molecule and its
    // modification tags. Only reuse it for unaligned input and an unchanged
    // alignment; promoted/clipped records still need full restoration.
    let needs_restore = |rec: &Record| {
        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() { return false; }
        let unchanged = original.is_unmapped() && !original.is_reverse()
            && rec.seq_len() == original.seq_len()
            && initial.map(|(tid, pos, reverse, cigar)| {
                rec.tid() == *tid && rec.pos() == *pos && rec.is_reverse() == *reverse
                    && rec.cigar().iter().eq(cigar.iter())
                    && !cigar.iter().any(|op| matches!(op, Cigar::HardClip(_)))
            }).unwrap_or(false);
        !unchanged
    };
    if !records.iter().any(&needs_restore) { return; }
    let mut original_seq = original.seq().as_bytes();
    let mut original_qual = original.qual().to_vec();
    if original.is_reverse() {
        original_seq = bio::alphabets::dna::revcomp(&original_seq);
        original_qual.reverse();
    }
    for rec in records {
        if !needs_restore(rec) { continue; }
        let cigar = CigarString(rec.cigar().iter().map(|op| match op {
            Cigar::HardClip(n) => Cigar::SoftClip(*n), other => *other,
        }).collect());
        let (seq, qual) = if rec.is_reverse() {
            let mut q = original_qual.clone(); q.reverse();
            (bio::alphabets::dna::revcomp(&original_seq), q)
        } else { (original_seq.clone(), original_qual.clone()) };
        rec.set(original.qname(), Some(&cigar), &seq, &qual);
        for tag in [b"MM", b"ML", b"MN"] {
            rec.remove_aux(tag).ok();
            if let Ok(value) = original.aux(tag) { rec.push_aux(tag, value).ok(); }
        }
    }
}

fn process_se_bam_record(
    rec: Record,
    ctx: &RammapCtx,
    header: &Header,
    preset_str: &str,
    realign: bool,
    mapq_rescue: Option<u8>,
    target_names_ref: &[&[u8]],
    secondary: bool,
    rescue_mode: Option<&str>,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> Vec<Record> {
    let result = (|| -> AnyResult<Vec<Record>> {
        if ctx.methylation.is_some() {
            anyhow::ensure!(rec.is_unmapped(), "--meth-bed requires original unaligned reads; found mapped read {}", String::from_utf8_lossy(rec.qname()));
            anyhow::ensure!(crate::methalign::is_contain_methylation(&rec),
                "read {} is missing MM/ML tags required by --meth-bed", String::from_utf8_lossy(rec.qname()));
        }
        let mut res = ctx.timings.measure(3, || process_single_read(&rec, ctx, header, preset_str, true, aln_ctx, map_ctx))?;
        let initial_primary = ctx.timings.measure(7, || res.iter()
            .find(|r| !r.is_unmapped() && !r.is_secondary() && !r.is_supplementary())
            .map(|r| (r.tid(), r.pos(), r.is_reverse(), CigarString(r.cigar().iter().copied().collect()))));
        if let Some(refiner) = &ctx.methylation {
            res = ctx.timings.measure(4, || refiner.refine(res, &rec, &ctx.target_names));
        }
        if realign {
            res = ctx.timings.measure(5, || realign_records(res, true, mapq_rescue, target_names_ref, rescue_mode));
        }
        ctx.timings.measure(6, || finalize_records(&mut res, ctx, target_names_ref));
        ctx.timings.measure(7, || restore_primary_read(&mut res, &rec, initial_primary.as_ref()));
        ctx.timings.measure(8, || refresh_sa_tags(&mut res, target_names_ref));
        if !secondary { res.retain(|r| !r.is_secondary()); }
        Ok(res)
    })();
    match result {
        Ok(records) => records,
        Err(error) => {
            let mut slot = ctx.error.lock().unwrap();
            if slot.is_none() { *slot = Some(error.to_string()); }
            Vec::new()
        }
    }
}

fn process_se_fastx_record(
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    ctx: &RammapCtx,
    header: &Header,
    preset_str: &str,
    realign: bool,
    mapq_rescue: Option<u8>,
    target_names_ref: &[&[u8]],
    secondary: bool,
    rescue_mode: Option<&str>,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> Vec<Record> {
    let result = process_raw_sequence(id, seq, qual, None, ctx, header,
        preset_str, true, aln_ctx, map_ctx);
    match result {
        Ok(mut res) => {
            if realign { res = realign_records(res, true, mapq_rescue, target_names_ref, rescue_mode); }
            finalize_records(&mut res, ctx, target_names_ref);
            refresh_sa_tags(&mut res, target_names_ref);
            if !secondary { res.retain(|r| !r.is_secondary()); }
            res
        }
        Err(error) => {
            let mut slot = ctx.error.lock().unwrap();
            if slot.is_none() { *slot = Some(error.to_string()); }
            Vec::new()
        }
    }
}

fn process_pe_fastx_record(
    id1: &[u8],
    seq1: &[u8],
    qual1: &[u8],
    id2: &[u8],
    seq2: &[u8],
    qual2: &[u8],
    ctx: &RammapCtx,
    header: &Header,
    preset_str: &str,
    realign: bool,
    mapq_rescue: Option<u8>,
    target_names_ref: &[&[u8]],
    rescue_mode: Option<&str>,
    aln_ctx: &mut AlignmentContext,
    map_ctx: &mut MappingWorkspace,
) -> Vec<Record> {
    let mut res1 = process_raw_sequence(
        id1,
        seq1,
        qual1,
        None,
        ctx,
        header,
        preset_str,
        !realign,
        aln_ctx,
        map_ctx,
    )
    .unwrap_or_default();
    let mut res2 = process_raw_sequence(
        id2,
        seq2,
        qual2,
        None,
        ctx,
        header,
        preset_str,
        !realign,
        aln_ctx,
        map_ctx,
    )
    .unwrap_or_default();

    for r in &mut res1 {
        let f = r.flags();
        r.set_flags(f | 0x40);
    }
    for r in &mut res2 {
        let f = r.flags();
        r.set_flags(f | 0x80);
    }

    let mut pe_res = res1;
    pe_res.extend(res2);

    if realign {
        pe_res = realign_records(
            pe_res,
            true,
            mapq_rescue,
            target_names_ref,
            rescue_mode,
        );
    }

    if realign && ctx.gap_rescue {
        let mut r1_post: Vec<Record> = pe_res
            .iter()
            .filter(|r| (r.flags() & 0x40) != 0)
            .cloned()
            .collect();
        let mut r2_post: Vec<Record> = pe_res
            .iter()
            .filter(|r| (r.flags() & 0x80) != 0)
            .cloned()
            .collect();

        if let Ok(retained) =
            apply_gap_rescue_to_records(
                &mut r1_post,
                id1,
                seq1,
                qual1,
                None,
                ctx,
                header,
                aln_ctx,
                map_ctx,
            )
        {
            pe_res.extend(retained.into_iter().map(|mut r| {
                let f = r.flags();
                r.set_flags(f | 0x40);
                r
            }));
        }
        if let Ok(retained) =
            apply_gap_rescue_to_records(
                &mut r2_post,
                id2,
                seq2,
                qual2,
                None,
                ctx,
                header,
                aln_ctx,
                map_ctx,
            )
        {
            pe_res.extend(retained.into_iter().map(|mut r| {
                let f = r.flags();
                r.set_flags(f | 0x80);
                r
            }));
        }
    }

    let r1_rescued: Vec<Record> = pe_res
        .iter()
        .filter(|r| (r.flags() & 0x40) != 0)
        .cloned()
        .collect();
    let r2_rescued: Vec<Record> = pe_res
        .into_iter()
        .filter(|r| (r.flags() & 0x80) != 0)
        .collect();

    let r1_best = r1_rescued.into_iter().max_by(|a, b| {
        let score_a = get_score(a);
        let score_b = get_score(b);
        if score_a != score_b {
            score_a.cmp(&score_b)
        } else {
            a.mapq().cmp(&b.mapq())
        }
    });

    let r2_best = r2_rescued.into_iter().max_by(|a, b| {
        let score_a = get_score(a);
        let score_b = get_score(b);
        if score_a != score_b {
            score_a.cmp(&score_b)
        } else {
            a.mapq().cmp(&b.mapq())
        }
    });

    let mut out_records = Vec::new();

    match (r1_best, r2_best) {
        (Some(mut r1), Some(mut r2)) => {
            let r1_tid = r1.tid();
            let r1_pos = r1.pos();
            let r1_is_rev = r1.is_reverse();

            let r2_tid = r2.tid();
            let r2_pos = r2.pos();
            let r2_is_rev = r2.is_reverse();

            r2.set_qname(r1.qname());

            let mut f1 = 0x1 | 0x40; // Paired, R1
            if r1_is_rev {
                f1 |= 0x10;
            }
            if r2_is_rev {
                f1 |= 0x20;
            }
            r1.set_flags(f1);
            r1.set_mtid(r2_tid);
            r1.set_mpos(r2_pos);

            let mut f2 = 0x1 | 0x80; // Paired, R2
            if r2_is_rev {
                f2 |= 0x10;
            }
            if r1_is_rev {
                f2 |= 0x20;
            }
            r2.set_flags(f2);
            r2.set_mtid(r1_tid);
            r2.set_mpos(r1_pos);

            out_records.push(r1);
            out_records.push(r2);
        }
        (Some(mut r1), None) => {
            let r1_is_rev = r1.is_reverse();
            let mut f1 = 0x1 | 0x40 | 0x8;
            if r1_is_rev {
                f1 |= 0x10;
            }
            r1.set_flags(f1);
            out_records.push(r1);
        }
        (None, Some(mut r2)) => {
            let r2_is_rev = r2.is_reverse();
            r2.set_qname(id1);
            let mut f2 = 0x1 | 0x80 | 0x8;
            if r2_is_rev {
                f2 |= 0x10;
            }
            r2.set_flags(f2);
            out_records.push(r2);
        }
        _ => {}
    }
    finalize_records(&mut out_records, ctx, target_names_ref);
    refresh_sa_tags(&mut out_records, target_names_ref);
    out_records
}

pub fn align(
    reference: &String,
    input_bams: &Vec<String>,
    output_bam: &String,
    output_bam_format: bool,
    output_cigar: bool,
    eqx: bool,
    cs: Option<String>,
    kv: Option<i16>,
    wv: Option<i16>,
    is_hpc: bool,
    batch_size: u64,
    soft_clip: bool,
    secondary: bool,
    mid_occ_frac: Option<(f32, Option<f32>)>,
    bounds_of_occurrence: Option<(i32, Option<i32>)>,
    max_gap: Option<i32>,
    max_gap_ref: Option<i32>,
    max_frag_len: Option<i32>,
    mask_level: Option<f32>,
    bw: Option<(i32, Option<i32>)>,
    min_cnt: Option<i32>,
    min_chain_score: Option<i32>,
    matching_score: Option<i32>,
    mismatch_penalty: Option<i32>,
    gap_open: Option<(i32, Option<i32>)>,
    gap_extension: Option<(i32, Option<i32>)>,
    z_drop: Option<(i32, Option<i32>)>,
    min_dp_max: Option<i32>,
    best_n: Option<i32>,
    pri_ratio: Option<f32>,
    max_qlen: Option<i32>,
    mini_batch_size: i64,
    seed: Option<i32>,
    porec_gap_rescue: bool,
    rescue_k: Option<i16>,
    rescue_w: Option<i16>,
    min_gap_len: usize,
    preset: &str,
    threads: usize,
    realign: bool,
    rescue_mode: Option<String>,
    mapq_rescue: Option<u8>,
    re_site: Option<String>,
) -> AnyResult<()> {
    let config = AlignPlusConfig {
        methylation: None,
        reference: reference.clone(),
        queries: input_bams.clone(),
        output: output_bam.clone(),
        output_bam_format,
        output_cigar,
        eqx,
        cs,
        kmer: kv,
        window: wv,
        is_hpc,
        batch_size,
        soft_clip,
        secondary,
        mid_occ_frac,
        bounds_of_occurrence,
        max_gap,
        max_gap_ref,
        max_frag_len,
        mask_level,
        bw,
        min_cnt,
        min_chain_score,
        matching_score,
        mismatch_penalty,
        gap_open,
        gap_extension,
        z_drop,
        min_dp_max,
        best_n,
        pri_ratio,
        max_qlen,
        mini_batch_size,
        seed,
        porec_gap_rescue,
        rescue_k,
        rescue_w,
        min_gap_len,
        preset: preset.to_string(),
        threads,
        realign,
        rescue_mode,
        mapq_rescue,
        porec_mapq_calibrate: false,
        candidate_probability: false,
        graph_assignment: false,
        re_site,
    };

    align_with_config(config)
}

pub fn align_with_config(config: AlignPlusConfig) -> AnyResult<()> {
    let AlignPlusConfig {
        methylation,
        reference,
        queries: input_bams,
        output: output_bam,
        output_bam_format,
        output_cigar,
        eqx,
        cs,
        kmer: kv,
        window: wv,
        is_hpc,
        batch_size,
        soft_clip,
        secondary,
        mid_occ_frac,
        bounds_of_occurrence,
        max_gap,
        max_gap_ref,
        max_frag_len,
        mask_level,
        bw,
        min_cnt,
        min_chain_score,
        matching_score,
        mismatch_penalty,
        gap_open,
        gap_extension,
        z_drop,
        min_dp_max,
        best_n,
        pri_ratio,
        max_qlen,
        mini_batch_size,
        seed,
        porec_gap_rescue,
        rescue_k,
        rescue_w,
        min_gap_len,
        preset,
        threads,
        realign,
        rescue_mode,
        mapq_rescue,
        porec_mapq_calibrate,
        candidate_probability,
        graph_assignment,
        re_site,
    } = config;

    anyhow::ensure!(threads > 0, "--threads must be greater than zero");
    anyhow::ensure!(kv.is_none_or(|v| (1..=28).contains(&v)), "k-mer size must be between 1 and 28");
    anyhow::ensure!(wv.is_none_or(|v| v > 0), "minimizer window must be positive");
    anyhow::ensure!(input_bams.iter().all(|p| p != &output_bam), "output must differ from the input");
    anyhow::ensure!(rescue_k.is_none_or(|v| (1..=28).contains(&v)) && rescue_w.is_none_or(|v| v > 0),
        "rescue k-mer size must be 1..28 and rescue window must be positive");
    anyhow::ensure!(min_gap_len > 0 && batch_size > 0, "minimum gap length and index batch size must be positive");
    anyhow::ensure!(reference != output_bam, "output must differ from the reference");
    if let Some(meth) = &methylation {
        anyhow::ensure!(meth.bed != output_bam, "output must differ from --meth-bed");

        anyhow::ensure!(!graph_assignment && !porec_mapq_calibrate,
            "experimental graph/MAPQ recalibration cannot be combined with --meth-bed yet");
        for path in &input_bams {
            anyhow::ensure!(path == "-" || path == "/dev/stdin" || path.starts_with("/dev/fd/")
                || path.to_ascii_lowercase().ends_with(".bam") || path.to_ascii_lowercase().ends_with(".sam"),
                "--meth-bed requires unaligned BAM/SAM with MM/ML tags, not FASTQ/FASTA: {path}");
        }
    }
    for path in &input_bams {
        if path != "-" && path != "/dev/stdin" && !path.starts_with("/dev/fd/") {
            std::fs::File::open(path).map_err(|e| anyhow!("cannot open input {path}: {e}"))?;
        }
    }
    if let Some(meth) = &methylation {
        std::fs::File::open(&meth.bed).map_err(|e| anyhow!("cannot open methylation reference {}: {e}", meth.bed))?;
    }
    let preset_str = preset.as_str();
    let rammap_preset = match preset_str {
        "porec" => "lr:hq",
        "porec:map-ont" => "map-ont",
        "porec:lr:hq" => "lr:hq",
        "cifi" => "map-hifi",
        "hic" => "sr",
        other => other,
    };

    let mut k: usize = kv.map(|v| v as usize).unwrap_or(15);
    let mut w: usize = wv.map(|v| v as usize).unwrap_or(10);
    let mut is_hpc = is_hpc;
    let mut opt = MapOptions::default();

    if let Err(err) = apply_preset_str(&mut opt, &mut k, &mut w, &mut is_hpc, rammap_preset) {
        log::warn!(
            "Unknown rammap preset: {} (from {}). Falling back to map-ont: {}",
            rammap_preset,
            preset_str,
            err
        );
        apply_preset_str(&mut opt, &mut k, &mut w, &mut is_hpc, "map-ont")
            .map_err(|e| anyhow!("Failed to apply default rammap preset map-ont: {}", e))?;
    }

    if let Some(v) = kv {
        k = v as usize;
    }
    if let Some(v) = wv {
        w = v as usize;
    }

    let mut best_n = best_n;
    let mut z_drop = z_drop;
    let mut min_chain_score = min_chain_score;
    let mut max_gap = max_gap;
    let mut min_cnt = min_cnt;
    let mut pri_ratio = pri_ratio;
    let mut max_frag_len = max_frag_len;
    let mut mid_occ_frac_for_index = 2e-4_f32;

    if preset_str == "porec" || preset_str == "porec:lr:hq" || preset_str == "porec:map-ont" || preset_str == "cifi" {
        if best_n.is_none() {
            best_n = Some(500);
        }
        if z_drop.is_none() {
            z_drop = Some((200, Some(64)));
        }
        // if min_chain_score.is_none() {
        //     min_chain_score = Some(15);
        // }
        // if max_gap.is_none() {
        //     max_gap = Some(200);
        // }

        if !realign && methylation.is_none() && pri_ratio.is_none() {
            pri_ratio = Some(1.0);
        }
    }

    if preset_str == "hic" {
        if best_n.is_none() {
            best_n = Some(20);
        }
        if min_chain_score.is_none() {
            min_chain_score = Some(40);
        }
        if min_cnt.is_none() {
            min_cnt = Some(3);
        }
        if max_gap.is_none() {
            max_gap = Some(5000);
        }
        if pri_ratio.is_none() {
            pri_ratio = Some(0.8);
        }
        if max_frag_len.is_none() {
            max_frag_len = Some(800);
        }
    }


    let final_patterns: Vec<Vec<u8>> = re_site
        .as_ref()
        .map(|s| {
            s.split(',')
                .map(|p| p.trim().as_bytes().to_vec())
                .filter(|p| !p.is_empty())
                .collect()
        })
        .unwrap_or_default();

    let timings = Arc::new(AlignTimings::new());
    let mut seqs = timings.measure(0, || -> AnyResult<_> {
        let mut reader = parse_fastx_file(&reference)
            .map_err(|e| anyhow!("Failed to read reference file {}: {}", reference, e))?;
        let mut seqs = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| anyhow!("invalid reference FASTA: {e}"))?;
            seqs.push((
                String::from_utf8_lossy(record.id()).into_owned(),
                record.seq().to_vec(),
            ));
        }

        Ok(seqs)
    })?;
    anyhow::ensure!(!seqs.is_empty(), "reference FASTA contains no sequences");
    let processing_error = Arc::new(std::sync::Mutex::new(None));

    if let Some((v1, v2)) = mid_occ_frac {
        if v1 >= 1.0 {
            opt.seeding.mid_occ = v1 as usize;
        } else {
            mid_occ_frac_for_index = v1;
        }
        if let Some(v) = v2 {
            if v >= 1.0 {
                opt.seeding.max_mid_occ = v as i32;
            }
        }
    }

    if let Some((lo, hi)) = bounds_of_occurrence {
        if lo > 0 {
            opt.seeding.min_mid_occ = lo;
        }
        if let Some(v) = hi {
            if v > 0 {
                opt.seeding.max_mid_occ = v;
            }
        }
    }

    if let Some(n) = best_n {
        opt.filtering.best_n = n;
    }

    if let Some(p) = pri_ratio {
        opt.filtering.pri_ratio = p;
    }

    if let Some(n) = min_chain_score {
        opt.chaining.min_chain_score = n;
    }

    if let Some(n) = min_cnt {
        opt.chaining.min_cnt = n;
    }

    if let Some(n) = max_gap {
        opt.chaining.max_gap = n;
    }

    if let Some(n) = max_gap_ref {
        opt.chaining.max_gap_ref = n;
    }

    if let Some(n) = max_frag_len {
        opt.pairing.max_frag_len = n;
    }

    if let Some(v) = mask_level {
        opt.filtering.mask_level = v;
    }

    if let Some((bw1, bw2)) = bw {
        opt.chaining.bandwidth = bw1;
        opt.chaining.bandwidth_long = bw2.unwrap_or(bw1);
    }

    if let Some(v) = matching_score {
        opt.scoring.match_score = v;
    }

    if let Some(v) = mismatch_penalty {
        opt.scoring.mismatch_penalty = v;
    }

    if let Some((v1, v2)) = gap_open {
        opt.scoring.gap_open = v1;
        opt.scoring.gap_open2 = v2.unwrap_or(v1);
    }

    if let Some((v1, v2)) = gap_extension {
        opt.scoring.gap_extend = v1;
        opt.scoring.gap_extend2 = v2.unwrap_or(v1);
    }

    if let Some((v1, v2)) = z_drop {
        opt.alignment.zdrop = v1;
        opt.alignment.zdrop_inv = v2.unwrap_or(v1);
    }

    if let Some(v) = min_dp_max {
        opt.alignment.min_dp_max = v;
    }

    if let Some(v) = max_qlen {
        opt.filtering.max_qlen = v;
    }

    if let Some(v) = seed {
        opt.filtering.seed = v;
    }

    if mini_batch_size > 0 {
        opt.mini_batch_size = mini_batch_size;
    }

    opt.chaining.chn_pen_gap = (opt.chaining.chain_gap_scale as f64 * 0.01 * (k as f64)) as f32;
    opt.chaining.chn_pen_skip =
        opt.filtering.chain_skip_scale * opt.scoring.match_score as f32 * 0.01;

    if !secondary && !realign && methylation.is_none() && !candidate_probability && !graph_assignment {
        opt.flags.insert(AlignFlags::NO_PRINT_2ND);
    } else {
        opt.flags.remove(AlignFlags::NO_PRINT_2ND);
    }

    if soft_clip {
        opt.flags.insert(AlignFlags::SOFTCLIP);
    } else {
        opt.flags.remove(AlignFlags::SOFTCLIP);
    }

    let batch_size = usize::try_from(batch_size)
        .map_err(|_| anyhow!("batch_size is too large for this platform: {}", batch_size))?;

    let rayon_pool = Arc::new(rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| anyhow!("Failed to build Rayon thread pool: {}", e))?);
    let worker_pool = rayon_pool.clone();

    log::info!(
        "Indexing reference with {} sequences by rammap Index::build (preset={}, k={}, w={}, hpc={}, batch_size={})...",
        seqs.len(),
        rammap_preset,
        k,
        w,
        is_hpc,
        batch_size
    );

    (|| -> AnyResult<()> {

        let enable_gap_rescue = porec_gap_rescue
            && matches!(preset_str, "porec" | "porec:lr:hq" | "porec:map-ont");
        // rammap consumes its input and keeps a packed reference. Retain the
        // ASCII buffers only while another index or methylation needs them.
        let mut index = timings.measure(2, || {
            let index_seqs = if methylation.is_some() || enable_gap_rescue {
                seqs.clone()
            } else {
                std::mem::take(&mut seqs)
            };
            rayon_pool.install(|| Index::build(index_seqs, w, k, is_hpc, batch_size))
        });
        index.index = 0;
        if opt.seeding.mid_occ == 0 {
            opt.seeding.mid_occ = index.cal_mid_occ(
                mid_occ_frac_for_index,
                opt.seeding.min_mid_occ,
                opt.seeding.max_mid_occ,
            );
        }

        let (sensitive_index, sensitive_opt) = if enable_gap_rescue {
            let rescue_k = rescue_k.map(|v| v as usize).unwrap_or(11);
            let rescue_w = rescue_w.map(|v| v as usize).unwrap_or(5);
            let mut rescue_opt = opt.clone();
            rescue_opt.chaining.min_chain_score = rescue_opt.chaining.min_chain_score.min(10);
            rescue_opt.chaining.min_cnt = rescue_opt.chaining.min_cnt.min(2);
            rescue_opt.filtering.best_n = rescue_opt.filtering.best_n.max(50);
            rescue_opt.filtering.pri_ratio = rescue_opt.filtering.pri_ratio.min(0.6);
            rescue_opt.alignment.zdrop = rescue_opt.alignment.zdrop.min(100);
            rescue_opt.alignment.zdrop_inv = rescue_opt.alignment.zdrop_inv.min(50);
            rescue_opt.chaining.chn_pen_gap =
                (rescue_opt.chaining.chain_gap_scale as f64 * 0.01 * (rescue_k as f64)) as f32;
            rescue_opt.chaining.chn_pen_skip =
                rescue_opt.filtering.chain_skip_scale * rescue_opt.scoring.match_score as f32 * 0.01;

            log::info!(
                "Building Pore-C sensitive gap rescue index (k={}, w={}, min_gap_len={})...",
                rescue_k,
                rescue_w,
                min_gap_len
            );
            let mut rescue_index = timings.measure(2, || {
                let rescue_seqs = if methylation.is_some() { seqs.clone() } else { std::mem::take(&mut seqs) };
                rayon_pool.install(|| Index::build(rescue_seqs, rescue_w, rescue_k, is_hpc, batch_size))
            });
            rescue_index.index = 0;
            if rescue_opt.seeding.mid_occ == 0 {
                rescue_opt.seeding.mid_occ = rescue_index.cal_mid_occ(
                    mid_occ_frac_for_index,
                    rescue_opt.seeding.min_mid_occ,
                    rescue_opt.seeding.max_mid_occ,
                );
            }
            (Some(rescue_index), Some(rescue_opt))
        } else {
            (None, None)
        };

        log::info!(
            "Reference indexed. Starting alignment with best_n={} pri_ratio={} mid_occ={} mini_batch_size={}",
            opt.filtering.best_n,
            opt.filtering.pri_ratio,
            opt.seeding.mid_occ,
            opt.mini_batch_size
        );

        let cs_mode = cs.as_deref();
        let do_cs = matches!(cs_mode, Some("short") | Some("long"));
        let cs_long = matches!(cs_mode, Some("long"));

        let out_cfg = OutputConfig {
            do_cigar: true,
            do_cs,
            cs_long,
            do_md: false,
            do_ds: false,
            eqx,
            output_sam: true,
            rg_id: None,
            split_mode: false,
        };

        let mut header = Header::new();
        for seq in &index.seqs {
            let mut r = HeaderRecord::new(b"SQ");
            r.push_tag(b"SN", &seq.name);
            r.push_tag(b"LN", &seq.len.to_string());
            header.push_record(&r);
        }
        let header_view = HeaderView::from_header(&header);

        let tid_map: HashMap<String, i32> = header_view
            .target_names()
            .iter()
            .enumerate()
            .map(|(i, name)| (String::from_utf8_lossy(name).to_string(), i as i32))
            .collect();

        let re_map = if !final_patterns.is_empty() {
            Some(Arc::new(GenomeREMap::build(
                &reference,
                &final_patterns,
                &tid_map,
            )?))
        } else {
            None
        };
        let enable_mapq_calibrate = porec_mapq_calibrate
            && matches!(
                preset_str,
                "porec" | "porec:lr:hq" | "porec:map-ont" | "cifi"
            );
        let enable_candidate_probability = candidate_probability
            && matches!(
                preset_str,
                "porec" | "porec:lr:hq" | "porec:map-ont" | "cifi"
            );
        let enable_graph_assignment = graph_assignment
            && matches!(
                preset_str,
                "porec" | "porec:lr:hq" | "porec:map-ont" | "cifi"
            );

        let meth_refiner = timings.measure(1, || methylation
            .map(|cfg| MethRefiner::from_owned(cfg, seqs)).transpose())?.map(Arc::new);
        let rammap_ctx = RammapCtx {
            timings: timings.clone(),
            methylation: meth_refiner.clone(),
            target_names: index.seqs.iter().map(|seq| seq.name.clone()).collect(),
            error: processing_error.clone(),
            index,
            opt,
            sensitive_index,
            sensitive_opt,
            re_map,
            mapq_calibrate: enable_mapq_calibrate,
            candidate_probability: enable_candidate_probability,
            graph_assignment: enable_graph_assignment,
            gap_rescue: enable_gap_rescue,
            min_gap_len,
            out_cfg,
            output_cigar,
            do_cs,
            cs_long,
            eqx,
            secondary,
            soft_clip,
        };

        if input_bams.len() > 1 {
            log::info!("Starting alignment of input files...");
        }

        let is_paf = if output_bam_format {
            false
        } else if output_bam.ends_with(".bam") || output_bam.ends_with(".sam") {
            false
        } else {
            true
        };
        log::info!("Preparing output writers (is_paf: {})...", is_paf);
        let mut bam_writer = if !is_paf {
            let format = if output_bam.ends_with(".sam") { bam::Format::Sam } else { bam::Format::Bam };
            let mut w = if output_bam == "-" {
                Writer::from_stdout(&header, format)?
            } else {
                Writer::from_path(&output_bam, &header, format)?
            };
            let _ = w.set_threads(threads.min(8));
            Some(w)
        } else {
            None
        };

        let mut paf_writer = if is_paf {
            Some(common_writer(&output_bam))
        } else {
            None
        };
        let target_names_owned: Vec<Vec<u8>> = header_view
            .target_names()
            .iter()
            .map(|name| name.to_vec())
            .collect();
        let target_names: Vec<&[u8]> = target_names_owned
            .iter()
            .map(|name| name.as_slice())
            .collect();
        let target_names_ref = target_names.as_slice();

        let target_lens: Vec<u64> = target_names_ref
            .iter()
            .map(|n| {
                header_view
                    .tid(*n)
                    .and_then(|tid| header_view.target_len(tid))
                    .unwrap_or(0) as u64
            })
            .collect();
        let target_lens_ref = &target_lens;

        let (tx_in, rx_in) = channel::bounded::<WorkerInput>(1);
        let (tx_out, rx_out) = channel::bounded::<OutputBatch>(2);
        let worker_batch_bases = if mini_batch_size > 0 {
            mini_batch_size as usize
        } else {
            500_000_000
        };
        let worker_min_records = threads.saturating_mul(4).clamp(1, 10_000);

        let is_fastq_pe = input_bams.len() == 2
            && {
                let ext1 = Path::new(&input_bams[0])
                    .extension()
                    .map(|e| e.to_string_lossy().to_lowercase());
                let ext2 = Path::new(&input_bams[1])
                    .extension()
                    .map(|e| e.to_string_lossy().to_lowercase());
                let is_fx = |ext: Option<String>| {
                    ext.map(|e| e == "fastq" || e == "fq" || e == "gz" || e == "fasta" || e == "fa")
                        .unwrap_or(false)
                };
                is_fx(ext1) && is_fx(ext2)
            }
            && (preset_str == "hic" || preset_str == "sr");

        thread::scope(|s| {
            let tx_in_clone = tx_in.clone();
            let input_error = processing_error.clone();
            s.spawn(move || {
                let chunk_capacity = worker_min_records;
                let mut input_read_count = 0usize;
                let mut next_input_log_count = 10_000usize;
                if is_fastq_pe {
                    log::info!(
                        "Detected Paired-End FASTX files for Hi-C. Sync-processing R1 and R2..."
                    );
                    let mut reader1 =
                        parse_fastx_file(&input_bams[0]).expect("Failed to open R1 fastq");
                    let mut reader2 =
                        parse_fastx_file(&input_bams[1]).expect("Failed to open R2 fastq");

                    loop {
                        let mut chunk = Vec::with_capacity(chunk_capacity);
                        let mut batch_limit =
                            QueryBaseBatchLimit::new(worker_batch_bases, worker_min_records);
                        loop {
                            match (reader1.next(), reader2.next()) {
                                (Some(Ok(r1)), Some(Ok(r2))) => {
                                    let id1 = r1.id().to_vec();
                                    let seq1 = r1.seq().to_vec();
                                    let qual1 = r1
                                        .qual()
                                        .map(|q| q.to_vec())
                                        .unwrap_or_else(|| vec![255u8; seq1.len()]);

                                    let id2 = r2.id().to_vec();
                                    let seq2 = r2.seq().to_vec();
                                    let qual2 = r2
                                        .qual()
                                        .map(|q| q.to_vec())
                                        .unwrap_or_else(|| vec![255u8; seq2.len()]);

                                    let query_bases = seq1.len().saturating_add(seq2.len());
                                    chunk.push((id1, seq1, qual1, id2, seq2, qual2));
                                    if batch_limit.push(query_bases, chunk.len()) {
                                        break;
                                    }
                                }
                                _ => break,
                            }
                        }
                        if chunk.is_empty() {
                            break;
                        }
                        let chunk_len = chunk.len();
                        if tx_in_clone.send(WorkerInput::PeFastx(chunk)).is_err() {
                            break;
                        }
                        log_input_read_progress(
                            &mut input_read_count,
                            &mut next_input_log_count,
                            chunk_len,
                        );
                    }
                } else {
                    for path in input_bams {
                        log::info!("Starting alignment for: {}", path);
                        let is_stdin = path == "-" || path == "/dev/stdin";
                        let is_pipe = path.starts_with("/dev/fd/");

                        let is_bam = if is_stdin || is_pipe {
                            true
                        } else {
                            let p = Path::new(&path);
                            match p.extension() {
                                Some(ext) => {
                                    let ext_str = ext.to_string_lossy().to_lowercase();
                                    ext_str == "bam" || ext_str == "cram" || ext_str == "sam"
                                }
                                None => false,
                            }
                        };

                        if is_bam {
                            log::info!("Processing HTS file for: {}", path);
                            let mut reader = if is_stdin {
                                Reader::from_stdin().expect("Failed to read from stdin")
                            } else {
                                Reader::from_path(path).expect("Failed to read BAM file")
                            };
                            let _ = reader.set_threads(threads.min(16));

                            let mut chunk = Vec::with_capacity(chunk_capacity);
                            let mut batch_limit =
                                QueryBaseBatchLimit::new(worker_batch_bases, worker_min_records);
                            for r in reader.records() {
                                let rec = match r {
                                    Ok(rec) => rec,
                                    Err(error) => {
                                        *input_error.lock().unwrap() = Some(format!("failed to read input alignment: {error}"));
                                        return;
                                    }
                                };
                                {
                                    let query_bases = rec.seq_len();
                                    chunk.push(rec);
                                    if batch_limit.push(query_bases, chunk.len()) {
                                        let chunk_len = chunk.len();
                                        let batch = std::mem::replace(
                                            &mut chunk,
                                            Vec::with_capacity(chunk_capacity),
                                        );
                                        if tx_in_clone.send(WorkerInput::SeBam(batch)).is_err() {
                                            break;
                                        }
                                        batch_limit.reset();
                                        log_input_read_progress(
                                            &mut input_read_count,
                                            &mut next_input_log_count,
                                            chunk_len,
                                        );
                                    }
                                }
                            }
                            if !chunk.is_empty() {
                                let chunk_len = chunk.len();
                                if tx_in_clone.send(WorkerInput::SeBam(chunk)).is_ok() {
                                    log_input_read_progress(
                                        &mut input_read_count,
                                        &mut next_input_log_count,
                                        chunk_len,
                                    );
                                }
                            }
                        } else if let Ok(mut n_reader) = parse_fastx_file(&path) {
                            log::info!("Detected Fastx format for: {}", &path);
                            let mut chunk = Vec::with_capacity(chunk_capacity);
                            let mut batch_limit =
                                QueryBaseBatchLimit::new(worker_batch_bases, worker_min_records);
                            while let Some(rec) = n_reader.next() {
                                let rec = match rec {
                                    Ok(rec) => rec,
                                    Err(error) => {
                                        *input_error.lock().unwrap() = Some(format!("invalid FASTX input {path}: {error}"));
                                        return;
                                    }
                                };
                                let id = rec.id().to_vec();
                                let seq = rec.seq().to_vec();
                                let qual = rec
                                    .qual()
                                    .map(|q| q.to_vec())
                                    .unwrap_or_else(|| vec![255u8; seq.len()]);
                                let query_bases = seq.len();
                                chunk.push((id, seq, qual));
                                if batch_limit.push(query_bases, chunk.len()) {
                                    let chunk_len = chunk.len();
                                    let batch = std::mem::replace(
                                        &mut chunk,
                                        Vec::with_capacity(chunk_capacity),
                                    );
                                    if tx_in_clone.send(WorkerInput::SeFastx(batch)).is_err() {
                                        break;
                                    }
                                    batch_limit.reset();
                                    log_input_read_progress(
                                        &mut input_read_count,
                                        &mut next_input_log_count,
                                        chunk_len,
                                    );
                                }
                            }
                            if !chunk.is_empty() {
                                let chunk_len = chunk.len();
                                if tx_in_clone.send(WorkerInput::SeFastx(chunk)).is_ok() {
                                    log_input_read_progress(
                                        &mut input_read_count,
                                        &mut next_input_log_count,
                                        chunk_len,
                                    );
                                }
                            }
                        } else {
                            log::info!("Processing input HTS file for: {}", path);
                            let mut reader = if is_stdin {
                                Reader::from_stdin().expect("Failed to read from stdin")
                            } else {
                                Reader::from_path(path).expect("Failed to read BAM file")
                            };
                            let _ = reader.set_threads(threads.min(16));
                            let mut chunk = Vec::with_capacity(chunk_capacity);
                            let mut batch_limit =
                                QueryBaseBatchLimit::new(worker_batch_bases, worker_min_records);
                            for r in reader.records() {
                                let rec = match r {
                                    Ok(rec) => rec,
                                    Err(error) => {
                                        *input_error.lock().unwrap() = Some(format!("failed to read input alignment: {error}"));
                                        return;
                                    }
                                };
                                {
                                    let query_bases = rec.seq_len();
                                    chunk.push(rec);
                                    if batch_limit.push(query_bases, chunk.len()) {
                                        let chunk_len = chunk.len();
                                        let batch = std::mem::replace(
                                            &mut chunk,
                                            Vec::with_capacity(chunk_capacity),
                                        );
                                        if tx_in_clone.send(WorkerInput::SeBam(batch)).is_err() {
                                            break;
                                        }
                                        batch_limit.reset();
                                        log_input_read_progress(
                                            &mut input_read_count,
                                            &mut next_input_log_count,
                                            chunk_len,
                                        );
                                    }
                                }
                            }
                            if !chunk.is_empty() {
                                let chunk_len = chunk.len();
                                if tx_in_clone.send(WorkerInput::SeBam(chunk)).is_ok() {
                                    log_input_read_progress(
                                        &mut input_read_count,
                                        &mut next_input_log_count,
                                        chunk_len,
                                    );
                                }
                            }
                        }
                    }
                }
                drop(tx_in_clone);
            });

            let tx_out_clone = tx_out.clone();
            let rescue_mode_str = rescue_mode.as_deref();
            s.spawn(move || worker_pool.install(|| {
                while let Ok(input) = rx_in.recv() {
                    let results: OutputBatch = match input {
                        WorkerInput::SeBam(recs) => {
                            if is_paf {
                                let paf_strings: String = recs
                                    .into_par_iter()
                                    .map_init(
                                        || (AlignmentContext::new(), MappingWorkspace::new()),
                                        |(aln_ctx, map_ctx), rec| {
                                            let mut acc = String::new();
                                            if rammap_ctx.methylation.is_some() || realign
                                                || rammap_ctx.gap_rescue
                                                || rammap_ctx.re_map.is_some()
                                                || rammap_ctx.mapq_calibrate
                                                || rammap_ctx.candidate_probability
                                                || rammap_ctx.graph_assignment
                                            {
                                                let aligned = process_se_bam_record(
                                                    rec,
                                                    &rammap_ctx,
                                                    &header,
                                                    preset_str,
                                                    realign,
                                                    mapq_rescue,
                                                    target_names_ref,
                                                    secondary,
                                                    rescue_mode_str,
                                                    aln_ctx,
                                                    map_ctx,
                                                );
                                                for r in aligned {
                                                    if let Some(paf_line) = record_to_paf_fast(
                                                        &r,
                                                        target_names_ref,
                                                        target_lens_ref,
                                                    ) {
                                                        acc.push_str(&paf_line);
                                                        acc.push('\n');
                                                    }
                                                }
                                            } else {
                                                acc.push_str(&format_raw_sequence_paf(
                                                    rec.qname(),
                                                    &rec.seq().as_bytes(),
                                                    &rammap_ctx,
                                                    aln_ctx,
                                                    map_ctx,
                                                ));
                                            }
                                            acc
                                        },
                                    )
                                    .reduce(
                                        || String::new(),
                                        |mut a, b| {
                                            a.push_str(&b);
                                            a
                                        },
                                    );
                                OutputBatch::Paf(paf_strings)
                            } else {
                                let flat_recs: Vec<Record> = recs
                                    .into_par_iter()
                                    .map_init(
                                        || (AlignmentContext::new(), MappingWorkspace::new()),
                                        |(aln_ctx, map_ctx), rec| {
                                        process_se_bam_record(
                                            rec,
                                            &rammap_ctx,
                                            &header,
                                            preset_str,
                                            realign,
                                            mapq_rescue,
                                            target_names_ref,
                                            secondary,
                                            rescue_mode_str,
                                            aln_ctx,
                                            map_ctx,
                                        )
                                    })
                                    .flatten()
                                    .collect();
                                OutputBatch::Bam(flat_recs)
                            }
                        }
                        WorkerInput::SeFastx(chunk) => {
                            if is_paf {
                                let paf_strings: String = chunk
                                    .into_par_iter()
                                    .map_init(
                                        || (AlignmentContext::new(), MappingWorkspace::new()),
                                        |(aln_ctx, map_ctx), (id, seq, qual)| {
                                            let mut acc = String::new();
                                            if rammap_ctx.methylation.is_some() || realign
                                                || rammap_ctx.gap_rescue
                                                || rammap_ctx.re_map.is_some()
                                                || rammap_ctx.mapq_calibrate
                                                || rammap_ctx.candidate_probability
                                                || rammap_ctx.graph_assignment
                                            {
                                                let recs = process_se_fastx_record(
                                                    &id,
                                                    &seq,
                                                    &qual,
                                                    &rammap_ctx,
                                                    &header,
                                                    preset_str,
                                                    realign,
                                                    mapq_rescue,
                                                    target_names_ref,
                                                    secondary,
                                                    rescue_mode_str,
                                                    aln_ctx,
                                                    map_ctx,
                                                );
                                                for r in recs {
                                                    if let Some(paf_line) = record_to_paf_fast(
                                                        &r,
                                                        target_names_ref,
                                                        target_lens_ref,
                                                    ) {
                                                        acc.push_str(&paf_line);
                                                        acc.push('\n');
                                                    }
                                                }
                                            } else {
                                                acc.push_str(&format_raw_sequence_paf(
                                                    &id,
                                                    &seq,
                                                    &rammap_ctx,
                                                    aln_ctx,
                                                    map_ctx,
                                                ));
                                            }
                                            acc
                                        },
                                    )
                                    .reduce(
                                        || String::new(),
                                        |mut a, b| {
                                            a.push_str(&b);
                                            a
                                        },
                                    );
                                OutputBatch::Paf(paf_strings)
                            } else {
                                let recs: Vec<Record> = chunk
                                    .into_par_iter()
                                    .map_init(
                                        || (AlignmentContext::new(), MappingWorkspace::new()),
                                        |(aln_ctx, map_ctx), (id, seq, qual)| {
                                        process_se_fastx_record(
                                            &id,
                                            &seq,
                                            &qual,
                                            &rammap_ctx,
                                            &header,
                                            preset_str,
                                            realign,
                                            mapq_rescue,
                                            target_names_ref,
                                            secondary,
                                            rescue_mode_str,
                                            aln_ctx,
                                            map_ctx,
                                        )
                                    })
                                    .flatten()
                                    .collect();
                                OutputBatch::Bam(recs)
                            }
                        }
                        WorkerInput::PeFastx(chunk) => {
                            if is_paf {
                                let paf_strings: String = chunk
                                    .into_par_iter()
                                    .map_init(
                                        || (AlignmentContext::new(), MappingWorkspace::new()),
                                        |(aln_ctx, map_ctx), (id1, seq1, qual1, id2, seq2, qual2)| {
                                            let mut acc = String::new();
                                            let recs = process_pe_fastx_record(
                                                &id1,
                                                &seq1,
                                                &qual1,
                                                &id2,
                                                &seq2,
                                                &qual2,
                                                &rammap_ctx,
                                                &header,
                                                preset_str,
                                                realign,
                                                mapq_rescue,
                                                target_names_ref,
                                                rescue_mode_str,
                                                aln_ctx,
                                                map_ctx,
                                            );
                                            for r in recs {
                                                if let Some(paf_line) = record_to_paf_fast(
                                                    &r,
                                                    target_names_ref,
                                                    target_lens_ref,
                                                ) {
                                                    acc.push_str(&paf_line);
                                                    acc.push('\n');
                                                }
                                            }
                                            acc
                                        },
                                    )
                                    .reduce(
                                        || String::new(),
                                        |mut a, b| {
                                            a.push_str(&b);
                                            a
                                        },
                                    );
                                OutputBatch::Paf(paf_strings)
                            } else {
                                let recs: Vec<Record> = chunk
                                    .into_par_iter()
                                    .map_init(
                                        || (AlignmentContext::new(), MappingWorkspace::new()),
                                        |(aln_ctx, map_ctx), (id1, seq1, qual1, id2, seq2, qual2)| {
                                        process_pe_fastx_record(
                                            &id1,
                                            &seq1,
                                            &qual1,
                                            &id2,
                                            &seq2,
                                            &qual2,
                                            &rammap_ctx,
                                            &header,
                                            preset_str,
                                            realign,
                                            mapq_rescue,
                                            target_names_ref,
                                            rescue_mode_str,
                                            aln_ctx,
                                            map_ctx,
                                        )
                                    })
                                    .flatten()
                                    .collect();
                                OutputBatch::Bam(recs)
                            }
                        }
                    };
                    if tx_out_clone.send(results).is_err() {
                        break;
                    }
                }
                drop(tx_out_clone);
            }));

            drop(tx_in);
            drop(tx_out);

            let mut count = 0;
            while let Ok(batch) = rx_out.recv() {
                match batch {
                    OutputBatch::Paf(paf_block) => {
                        if !paf_block.is_empty() {
                            if let Some(ref mut pw) = paf_writer {
                                timings.measure(9, || pw.write_all(paf_block.as_bytes()))?;
                                count += paf_block.as_bytes().iter().filter(|&&c| c == b'\n').count();
                            }
                        }
                    }
                    OutputBatch::Bam(records) => {
                        if let Some(ref mut bw) = bam_writer {
                            for rec in records {
                                timings.measure(9, || bw.write(&rec))?;
                                count += 1;
                            }
                        }
                    }
                }
            }
            if let Some(error) = processing_error.lock().unwrap().take() {
                return Err(anyhow!(error));
            }
            log::info!("Finished. Total output mappings: {}", count);

            Ok::<(), anyhow::Error>(())
        })?;

        timings.measure(9, || { drop(bam_writer.take()); drop(paf_writer.take()); });
        Ok(())
    })()?;
    timings.report();
    Ok(())
}
