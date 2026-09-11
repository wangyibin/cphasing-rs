#![allow(dead_code)]
#![allow(unused_imports)]
use anyhow::Result as anyResult;
use crossbeam_channel::{Receiver, Sender, bounded, unbounded};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use memchr::memchr_iter;
use noodles::bed::record;
use polars::prelude::*;
use rayon::prelude::*;
use rust_htslib::bam::{
    self, Header, HeaderView, Read, Reader, Record, Writer, header::HeaderRecord, record::Aux,
    record::Cigar, record::CigarString, record::CigarStringView,
};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufRead;
use std::io::{Cursor, Read as StdRead, Seek, Write};
use std::path::Path;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::rc::Rc;
use std::str;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};
use tempfile::NamedTempFile;

use crate::core::BaseTable;
use crate::core::common_reader;
use crate::fastx::Fastx;
use crate::methy::qual_to_prob;

// Preserve the existing permissive bedGraph parsing and start-coordinate semantics.
fn parse_bedgraph_line(line: &str, cov_cutoff: f64) -> Option<(&str, i64)> {
    if line.as_bytes().first() == Some(&b'#') { return None; }
    let mut fields = line.split_ascii_whitespace();
    let ctg = fields.next()?;
    let start = fields.next()?;
    let _end = fields.next();
    let cov: f64 = fields.next()?.parse().ok()?;
    if cov < cov_cutoff { return None; }
    Some((ctg, start.parse().ok()?))
}

pub fn parse_bedgraph(
    bedgraph: &String,
    cov_cutoff: f64,
) -> anyResult<HashMap<String, HashSet<i64>>> {
    let profile = std::env::var("CPHASING_METH_TIMING").as_deref() == Ok("1");
    let mut text_time = Duration::ZERO;
    let mut table_time = Duration::ZERO;
    let mut lines = 0u64;
    let mut retained = 0u64;
    let mut fh = common_reader(bedgraph);
    let mut map: HashMap<String, HashSet<i64>> = HashMap::new();
    let mut line = String::new();
    let mut current_name = String::new();
    let mut current_sites = HashSet::new();
    let mut contig_runs = 0u64;

    loop {
        let text_start = profile.then(Instant::now);
        line.clear();
        let n = fh.read_line(&mut line)?;
        let parsed = if n == 0 { None } else { parse_bedgraph_line(&line, cov_cutoff) };
        if let Some(start) = text_start { text_time += start.elapsed(); }
        if n == 0 { break; }
        lines += 1;
        if let Some((ctg, s)) = parsed {
            let table_start = profile.then(Instant::now);
            if current_name != ctg {
                // Temporarily own the active set to avoid per-row name allocation
                // and contig hashing. Revisit an earlier contig by removing its
                // existing set, preserving duplicates and arbitrarily ordered input.
                if !current_name.is_empty() {
                    map.insert(std::mem::take(&mut current_name), std::mem::take(&mut current_sites));
                }
                (current_name, current_sites) = map.remove_entry(ctg)
                    .unwrap_or_else(|| (ctg.to_owned(), HashSet::new()));
                contig_runs += 1;
            }
            current_sites.insert(s);
            if let Some(start) = table_start { table_time += start.elapsed(); }
            retained += 1;
        }
    }

    let table_start = profile.then(Instant::now);
    if !current_name.is_empty() { map.insert(current_name, current_sites); }
    if let Some(start) = table_start { table_time += start.elapsed(); }
    if profile {
        log::info!("METH_TIMING contig_runs={}", contig_runs);
        log::info!("METH_TIMING stage=text_read_parse seconds={:.9} lines={}", text_time.as_secs_f64(), lines);
        log::info!("METH_TIMING stage=site_table seconds={:.9} retained_rows={}", table_time.as_secs_f64(), retained);
    }
    log::info!(
        "Load {} methylation sites from bedgraph {} after filter by ref_prob_cutoff {}",
        map.values().map(|s| s.len()).sum::<usize>(),
        bedgraph,
        cov_cutoff
    );
    Ok(map)
}

// split records to primary and supplementary with multiple secondary alignments, respectively
pub fn split_records(records: Vec<Record>) -> Vec<Vec<Record>> {
    let mut result: Vec<Vec<Record>> = Vec::new();

    let mut current_group: Vec<Record> = Vec::new();
    for record in records.into_iter() {
        let is_supplementary = record.is_supplementary();
        let is_secondary = record.is_secondary();
        let is_primary = !is_supplementary && !is_secondary;
        if is_primary || is_supplementary {
            if !current_group.is_empty() {
                result.push(std::mem::take(&mut current_group));
            }
            current_group.push(record);
        } else if is_secondary {
            current_group.push(record);
        }
    }

    if !current_group.is_empty() {
        result.push(current_group);
    }

    result
}

pub fn get_as(record: &Record) -> Option<i32> {
    match record.aux(b"AS") {
        Ok(Aux::I8(val)) => Some(val as i32),
        Ok(Aux::I16(val)) => Some(val as i32),
        Ok(Aux::I32(val)) => Some(val),
        Ok(Aux::U8(val)) => Some(val as i32),
        Ok(Aux::U16(val)) => Some(val as i32),
        Ok(Aux::U32(val)) => Some(val as i32),
        _ => None,
    }
}

pub fn is_contain_methylation(record: &Record) -> bool {
    record.aux(b"MM").is_ok() && record.aux(b"ML").is_ok()
}

fn revcomp(bs: &[u8]) -> Vec<u8> {
    fn comp(b: u8) -> u8 {
        match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        }
    }
    let mut v: Vec<u8> = bs.iter().map(|b| comp(*b)).collect();
    v.reverse();
    v
}

fn query_alignment_termini(rec: &bam::Record) -> (i64, i64) {
    let cig = rec.cigar();
    let mut qstart = 0i64;
    let mut qlen_aln = 0i64;
    let mut left = true;
    for op in cig.iter() {
        match op {
            Cigar::HardClip(len) if left => {
                qstart += *len as i64;
            }
            Cigar::SoftClip(len) if left => {
                qstart += *len as i64;
            }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                qlen_aln += *len as i64;
                left = false;
            }
            Cigar::Ins(len) => {
                qlen_aln += *len as i64;
                left = false;
            }
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                left = false;
            }
            Cigar::SoftClip(_) => { /* trailing soft */ }
            Cigar::HardClip(len) => {
                /* trailing hard will add to end */
                let _ = len;
            }
            _ => {
                left = false;
            }
        }
    }
    // trailing hard clip
    let mut right_h = 0i64;
    if let Some(last) = cig.last() {
        if let Cigar::HardClip(len) = last {
            right_h = *len as i64;
        }
    }
    (qstart, qstart + qlen_aln + right_h)
}

fn condense_cigar(cig: &bam::record::CigarStringView) -> String {
    let mut lsoft = 0i64;
    let mut rsoft = 0i64;
    let mut hard_l = 0i64;
    let mut hard_r = 0i64;
    let mut m: i64 = 0;
    let mut ins: i64 = 0;
    let mut del: i64 = 0;
    let ops = cig.iter().collect::<Vec<_>>();
    for (i, op) in ops.iter().enumerate() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => m += *l as i64,
            Cigar::Ins(l) => ins += *l as i64,
            Cigar::Del(l) => del += *l as i64,
            Cigar::SoftClip(l) => {
                if i == 0 {
                    lsoft += *l as i64;
                }
                if i == ops.len() - 1 {
                    rsoft += *l as i64;
                }
            }
            Cigar::HardClip(l) => {
                if i == 0 {
                    hard_l += *l as i64;
                }
                if i == ops.len() - 1 {
                    hard_r += *l as i64;
                }
            }
            _ => {}
        }
    }
    let mut s = String::new();
    if lsoft + hard_l > 0 {
        s.push_str(&format!("{}S", lsoft + hard_l));
    }
    let match_len = if ins > del { m + del } else { m + ins };
    s.push_str(&format!("{}M", match_len.max(0)));
    let indel = del - ins;
    if indel > 0 {
        s.push_str(&format!("{}D", indel));
    } else if indel < 0 {
        s.push_str(&format!("{}I", -indel));
    }
    if rsoft + hard_r > 0 {
        s.push_str(&format!("{}S", rsoft + hard_r));
    }
    s
}

fn reconstruct_sa(rec: &bam::Record, hdr: &HeaderView) -> anyResult<String> {
    let tid = rec.tid();
    let rname = std::str::from_utf8(hdr.tid2name(tid as u32))?.to_string();
    let pos1 = rec.pos() + 1;
    let strand = if rec.is_reverse() { "-" } else { "+" };
    let cigar = condense_cigar(&rec.cigar());
    let mapq = rec.mapq();
    let nm = match rec.aux(b"NM") {
        Ok(Aux::I8(v)) => v as i64,
        Ok(Aux::I16(v)) => v as i64,
        Ok(Aux::I32(v)) => v as i64,
        Ok(Aux::U8(v)) => v as i64,
        Ok(Aux::U16(v)) => v as i64,
        Ok(Aux::U32(v)) => v as i64,
        _ => 0,
    };
    Ok(format!(
        "{},{},{},{},{},{};",
        rname, pos1, strand, cigar, mapq, nm
    ))
}

fn aligned_pairs(rec: &bam::Record) -> Vec<(i64, i64)> {
    let mut v = Vec::new();
    let mut q = 0i64;

    for op in rec.cigar().iter() {
        match op {
            Cigar::SoftClip(l) | Cigar::HardClip(l) => q += *l as i64,
            _ => break,
        }
    }
    let mut r = rec.pos();
    for op in rec.cigar().iter() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                for _ in 0..(*l as i64) {
                    v.push((q, r));
                    q += 1;
                    r += 1;
                }
            }
            Cigar::Ins(l) => {
                q += *l as i64;
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                r += *l as i64;
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {} // _ => {}
        }
    }
    v
}

fn get_5mc_sites_from_read(
    read_seq: &[u8],
    mm: Option<&str>,
    ml: Option<&[u8]>,
    prob_cutoff: u8,
) -> (Vec<i64>, Vec<i16>) {
    let mm_s = match mm {
        Some(s) => s,
        None => return (Vec::new(), Vec::new()),
    };
    let ml_v = ml.unwrap_or(&[]);

    let c_positions: Vec<i64> = memchr_iter(b'C', read_seq).map(|i| i as i64).collect();

    // -1 means unknown; only observed calls (or implicit low-probability
    // calls under the '.' convention) may contribute to the score.
    let mut flags: Vec<i16> = vec![-1; c_positions.len()];
    let mut ml_idx = 0usize;
    for part in mm_s.trim_end_matches(';').split(';') {
        let mut fields = part.split(',');
        let head = fields.next().unwrap_or("");
        if head.len() < 3 { continue; }
        let codes = head[2..].trim_end_matches(['?', '.']);
        let stride = if codes.bytes().all(|b| b.is_ascii_digit()) { 1 } else { codes.len() };
        let methyl_index = if head.starts_with("C+") {
            codes.bytes().position(|b| b == b'm')
        } else { None };
        if methyl_index.is_some() && !head.ends_with('?') { flags.fill(0); }
        let mut next = 0usize;
        for field in fields {
            let Ok(skip) = field.parse::<usize>() else { continue; };
            let Some(pos) = next.checked_add(skip) else { break; };
            if let Some(k) = methyl_index {
                if let Some(flag) = flags.get_mut(pos) {
                    *flag = ml_v.get(ml_idx.saturating_add(k))
                        .map(|&p| i16::from(p >= prob_cutoff)).unwrap_or(-1);
                }
            }
            next = pos.saturating_add(1);
            ml_idx = ml_idx.saturating_add(stride);
        }
    }

    (c_positions, flags)
}

fn count_5mc_consistency_split(
    cigar: &bam::record::CigarStringView,
    c_positions: &[i64],
    meth_flags: &[i16],
    is_fwd: bool,
    mut r: i64,
    ref_bytes: &[u8],
    ref_len: i64,
    site_set: Option<&HashSet<i64>>,
    cpg: bool,
) -> (i32, i32, i32) {
    let mut q: i64 = 0;
    let mut i: usize = 0;
    let mut matches = 0i32;
    let mut read_missing = 0i32;
    let mut ref_missing = 0i32;

    for op in cigar.iter() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                let len = *l as i64;
                let end_q = q + len;

                while i < c_positions.len() && c_positions[i] < q {
                    i += 1;
                }

                while i < c_positions.len() && c_positions[i] < end_q {
                    let cpos = c_positions[i];
                    let offset = cpos - q;
                    let rpos = r + offset;

                    if meth_flags[i] >= 0 && rpos >= 0 && rpos < ref_len {
                        let rb = ref_bytes[rpos as usize];
                        let is_meth = meth_flags[i] != 0;
                        let hit = if is_fwd {
                            if cpg {
                                let r_after = if rpos + 1 < ref_len {
                                    ref_bytes[(rpos + 1) as usize]
                                } else {
                                    b'N'
                                };
                                r_after == b'G'
                                    && rb == b'C'
                                    && site_set.map(|s| s.contains(&rpos)).unwrap_or(false)
                            } else {
                                rb == b'C' && site_set.map(|s| s.contains(&rpos)).unwrap_or(false)
                            }
                        } else {
                            if cpg {
                                let r_prev = if rpos > 0 {
                                    ref_bytes[(rpos - 1) as usize]
                                } else {
                                    b'N'
                                };
                                r_prev == b'C'
                                    && rb == b'G'
                                    && site_set.map(|s| s.contains(&(rpos - 1))).unwrap_or(false)
                            } else {
                                rb == b'G'
                                    && site_set.map(|s| s.contains(&(rpos - 1))).unwrap_or(false)
                            }
                        };

                        if hit == is_meth {
                            matches += 1;
                        } else {
                            if hit && !is_meth {
                                read_missing += 1;
                            } else if !hit && is_meth {
                                ref_missing += 1;
                            }
                        }
                        // println!("is_meth: {}, hit: {}, cpos: {}, q: {}, r: {}, offset: {}, rpos: {}", is_meth, hit, cpos, q, r, offset, rpos);
                    }
                    i += 1;
                }

                q = end_q;
                r += len;
            }
            Cigar::Ins(l) => {
                q += *l as i64;
                while i < c_positions.len() && c_positions[i] < q {
                    i += 1;
                }
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                r += *l as i64;
            }
            Cigar::SoftClip(l) => {
                q += *l as i64;
                while i < c_positions.len() && c_positions[i] < q {
                    i += 1;
                }
            }
            Cigar::HardClip(l) | Cigar::Pad(l) => {
                q += *l as i64;
                while i < c_positions.len() && c_positions[i] < q {
                    i += 1;
                }
            } // _ => {}
        }
    }

    (matches, read_missing, ref_missing)
}

fn get_mm_ml_from_rec(rec: &bam::Record) -> (Option<String>, Option<Vec<u8>>) {
    let mm_owned = match rec.aux(b"MM") {
        Ok(Aux::String(s)) => Some(String::from_utf8_lossy(s.as_bytes()).to_string()),
        _ => None,
    };
    let ml_owned: Option<Vec<u8>> = match rec.aux(b"ML") {
        Ok(Aux::ArrayU8(v)) => Some(v.iter().collect()),
        Ok(Aux::ArrayI32(v)) => Some(v.iter().map(|x| x.clamp(0, 255) as u8).collect()),
        Ok(Aux::ArrayI16(v)) => Some(v.iter().map(|x| (x as i32).clamp(0, 255) as u8).collect()),
        Ok(Aux::ArrayU16(v)) => Some(v.iter().map(|x| (x as i32).clamp(0, 255) as u8).collect()),
        _ => None,
    };
    (mm_owned, ml_owned)
}

fn recalc_score_and_update(
    rec: &mut bam::Record,
    tid_names: &[String],
    match_score: i32,
    ref_miss_penalty: i32,
    read_miss_penalty: i32,
    prob_cutoff: u8,
    fa: &HashMap<String, String>,
    met_sites: &HashMap<String, HashSet<i64>>,
    cpg: bool,
) -> anyResult<()> {
    // let read_seq_raw = rec.seq().as_bytes();
    // let read_seq = if rec.is_reverse() { revcomp(&read_seq_raw) } else { read_seq_raw.to_vec() };

    let read_seq = rec.seq().as_bytes().to_vec();

    let (mm_owned, ml_owned) = get_mm_ml_from_rec(rec);
    let (pos_read, flags_read) = get_5mc_sites_from_read(
        &read_seq,
        mm_owned.as_deref(),
        ml_owned.as_deref(),
        prob_cutoff,
    );

    let mut pairs: Vec<(i64, i16)> = pos_read.into_iter().zip(flags_read.into_iter()).collect();
    if rec.is_reverse() {
        let l = read_seq.len() as i64;
        for (p, _) in pairs.iter_mut() {
            *p = l - 1 - *p;
        }
    }
    pairs.sort_unstable_by_key(|x| x.0);
    let (c_positions_q, meth_flags_q): (Vec<i64>, Vec<i16>) = pairs.into_iter().unzip();

    let rname = match tid_names.get(rec.tid() as usize) {
        Some(s) => s.as_str(),
        None => return Ok(()),
    };
    let Some(ref_seq) = fa.get(rname) else {
        return Ok(());
    };
    let ref_bytes = ref_seq.as_bytes();
    let ref_len = ref_bytes.len() as i64;
    let is_fwd = !rec.is_reverse();
    let site_set = met_sites.get(rname);

    let (mc, mr, mf) = count_5mc_consistency_split(
        &rec.cigar(),
        &c_positions_q,
        &meth_flags_q,
        is_fwd,
        rec.pos(),
        ref_bytes,
        ref_len,
        site_set,
        cpg,
    );

    let mismatches = mr + mf;

    let orig_as = match rec.aux(b"AS") {
        Ok(Aux::I8(v)) => v as i32,
        Ok(Aux::I16(v)) => v as i32,
        Ok(Aux::I32(v)) => v,
        Ok(Aux::U8(v)) => v as i32,
        Ok(Aux::U16(v)) => v as i32,
        Ok(Aux::U32(v)) => v as i32,
        _ => 0,
    };
    // let new_as = orig_as - penalty * mismatches;
    let new_as = orig_as + match_score * mc - read_miss_penalty * mr - ref_miss_penalty * mf;
    rec.push_aux(b"MA", Aux::I32(mismatches as i32)).ok();
    rec.remove_aux(b"AS").ok();
    rec.push_aux(b"AS", Aux::I32(new_as)).ok();

    Ok(())
}

fn recalc_score_and_update2(
    recs: &mut Vec<bam::Record>,
    read_seq: &[u8],
    c_positions_base: &[i64],
    meth_flags_base: &[i16],
    tid_names: &[String],
    match_score: i32,
    ref_miss_penalty: i32,
    read_miss_penalty: i32,
    prob_cutoff: u8,
    fa: &HashMap<String, String>,
    met_sites: &HashMap<String, HashSet<i64>>,
    cpg: bool,
) -> anyResult<()> {
    for record in recs.iter_mut() {
        // println!("c_positions_base: {:?}", c_positions_base);
        let mut pairs: Vec<(i64, i16)> = c_positions_base
            .iter()
            .copied()
            .zip(meth_flags_base.iter().copied())
            .collect();
        if record.is_reverse() {
            let l = read_seq.len() as i64;
            for (p, _) in pairs.iter_mut() {
                *p = l - 1 - *p;
            }
        }
        pairs.sort_unstable_by_key(|x| x.0);
        let (c_positions_q, meth_flags_q): (Vec<i64>, Vec<i16>) = pairs.into_iter().unzip();

        let is_fwd: bool = !record.is_reverse();
        let rname = match tid_names.get(record.tid() as usize) {
            Some(s) => s.as_str(),
            None => return Ok(()),
        };
        let Some(ref_seq) = fa.get(rname) else {
            return Ok(());
        };
        let ref_bytes = ref_seq.as_bytes();
        let ref_len = ref_bytes.len() as i64;
        let site_set = met_sites.get(rname);
        // println!("rname : {}", rname);
        // println!("{:?}", site_set);
        // println!("site_set_len: {}", site_set.map(|s| s.len()).unwrap_or(0));
        // println!("c_positions_q: {:?}", c_positions_q);
        // println!("read: {}, is_fwd: {}", std::str::from_utf8(record.qname()).unwrap_or("N/A"), is_fwd);
        // println!("pos_len: {}, meth_len: {}", c_positions_q.len(), meth_flags_q.len());
        // let filtered_flags = meth_flags_q.iter().filter(|&&f| f != 0).count();
        // println!("filtered flags (meth): {}", filtered_flags);
        let (mc, mr, mf) = count_5mc_consistency_split(
            &record.cigar(),
            &c_positions_q,
            &meth_flags_q,
            !record.is_reverse(),
            record.pos(),
            ref_bytes,
            ref_len,
            site_set,
            cpg,
        );

        // println!("mc: {}, mr: {}, mf: {}", mc, mr, mf);
        let mismatches = mr + mf;

        let orig_as = match record.aux(b"AS") {
            Ok(Aux::I8(v)) => v as i32,
            Ok(Aux::I16(v)) => v as i32,
            Ok(Aux::I32(v)) => v,
            Ok(Aux::U8(v)) => v as i32,
            Ok(Aux::U16(v)) => v as i32,
            Ok(Aux::U32(v)) => v as i32,
            _ => 0,
        };
        // let new_as = orig_as - penalty * mismatches;
        let delta = match_score * mc - read_miss_penalty * mr - ref_miss_penalty * mf;
        record.remove_aux(b"s0").ok();
        record.push_aux(b"s0", Aux::I32(orig_as)).ok();
        record.remove_aux(b"m0").ok();
        record.push_aux(b"m0", Aux::I32(delta)).ok();
        let new_as = orig_as + delta;
        record.push_aux(b"MA", Aux::I32(mismatches as i32)).ok();
        record.remove_aux(b"AS").ok();
        record.push_aux(b"AS", Aux::I32(new_as)).ok();
    }

    Ok(())
}

/// Refine a complete read candidate set without BAM I/O or an inner thread pool.
/// `original` supplies MM/ML in original read coordinates for online alignment.
fn refine_read_records(
    records: Vec<Record>, original: Option<&Record>, tid_names: &[String],
    fa_hash: &HashMap<String, String>, bg_map: &HashMap<String, HashSet<i64>>,
    match_score: i32, ref_miss_penalty: i32, read_miss_penalty: i32,
    prob_cutoff: u8, designate_mapq: u8, is_set_y: bool, cpg: bool,
) -> Vec<Record> {
    let fa_len = fa_hash.len();
    let bg_len = bg_map.len();
    if records.is_empty() {
        return Vec::<Record>::new();
    }
    let source = original.unwrap_or(&records[0]);
    let (mm_owned, ml_owned, read_seq_raw) =
        if fa_len > 0 && bg_len > 0 && !is_set_y && is_contain_methylation(source) {
            let (mm, ml) = get_mm_ml_from_rec(source);
            let seq = source.seq().as_bytes().to_vec();
            if source.is_reverse() {
                (mm, ml, revcomp(&seq))
            } else {
                (mm, ml, seq)
            }
        } else {
            (None, None, Vec::new())
        };
    let (c_positions_base, meth_flags_base) = get_5mc_sites_from_read(
        &read_seq_raw, mm_owned.as_deref(), ml_owned.as_deref(), prob_cutoff,
    );
    let mut groups = split_records(records);

    if fa_len > 0 && bg_len > 0 {
        for group in groups.iter_mut() {
            if group.len() == 1 || group.iter().any(|r| get_as(r).is_none()) {
                continue;
            }

            let primary_mapq = group[0].mapq();
            if primary_mapq >= designate_mapq {
                continue;
            }

            if is_set_y {
                for rec in group.iter_mut() {
                    let _ = recalc_score_and_update(
                        rec,
                        tid_names,
                        match_score,
                        ref_miss_penalty,
                        read_miss_penalty,
                        prob_cutoff,
                        fa_hash,
                        bg_map,
                        cpg,
                    );
                }
            } else {
                if c_positions_base.is_empty() {
                    continue;
                }
                let _ = recalc_score_and_update2(
                    group,
                    &read_seq_raw,
                    &c_positions_base,
                    &meth_flags_base,
                    tid_names,
                    match_score,
                    ref_miss_penalty,
                    read_miss_penalty,
                    prob_cutoff,
                    fa_hash,
                    bg_map,
                    cpg,
                );
            }

            let as_vec2: Vec<i32> = group.iter().filter_map(|r| get_as(r)).collect();
            let (best_recs_idx, max_as2) = as_vec2
                .iter()
                .copied()
                .enumerate()
                .max_by_key(|&(_, v)| v)
                .unwrap_or((0, i32::MIN));
            let secondary_as2 = as_vec2
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != best_recs_idx)
                .map(|(_, &v)| v)
                .max()
                .unwrap_or(i32::MIN);

            if max_as2 > secondary_as2 {
                group[best_recs_idx].set_mapq(designate_mapq);
                if best_recs_idx != 0 {
                    group[best_recs_idx].push_aux(b"RF", Aux::String("Y")).ok();
                    let _primary_flag = group[0].flags();
                    let _primary_flag0 = _primary_flag & 0x800;
                    let _flag = group[best_recs_idx].flags();

                    group[best_recs_idx].set_flags((_flag & !(0x100 | 0x800)) | _primary_flag0);
                    group[0].set_flags((_primary_flag & !0x800) | 0x100);
                    group[0].set_mapq(0);
                    group[best_recs_idx].push_aux(b"tp", Aux::String("P")).ok();
                    group[0].push_aux(b"tp", Aux::String("S")).ok();
                    group.swap(0, best_recs_idx);
                }
            }
        }
    } else {
        for group in groups.iter_mut() {
            if group.len() == 1 || group.iter().any(|r| get_as(r).is_none()) {
                continue;
            }

            let primary_mapq = group[0].mapq();
            if primary_mapq >= designate_mapq {
                continue;
            }

            let as_vec: Vec<i32> = group.iter().filter_map(|r| get_as(r)).collect();
            if as_vec.len() < 2 {
                continue;
            }

            let primary_as = as_vec[0];
            let (best_idx, max_as) = as_vec
                .iter()
                .copied()
                .enumerate()
                .max_by_key(|&(_, v)| v)
                .unwrap_or((0, i32::MIN));
            let secondary_as = as_vec
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != best_idx)
                .map(|(_, &v)| v)
                .max()
                .unwrap_or(i32::MIN);

            if (primary_as > secondary_as) && (best_idx == 0) {
                group[0].set_mapq(designate_mapq);
                group[0].push_aux(b"RF", Aux::String("Y")).ok();
            }
        }
    }

    groups.into_iter().flatten().collect::<Vec<Record>>()
}

fn process_batch_with_rayon(
    batch: Vec<(usize, Vec<Record>)>,
    tid_names: &[String],
    fa_hash: &HashMap<String, String>,
    bg_map: &HashMap<String, HashSet<i64>>,
    match_score: i32,
    ref_miss_penalty: i32,
    read_miss_penalty: i32,
    prob_cutoff: u8,
    designate_mapq: u8,
    is_set_y: bool,
    cpg: bool,
) -> Vec<Record> {
    let target_names: Vec<&[u8]> = tid_names.iter().map(|name| name.as_bytes()).collect();
    batch.into_par_iter().flat_map(|(_, records)| {
        let original = records.iter().find(|r| !r.is_secondary() && !r.is_supplementary()
            && r.seq_len() > 0 && !r.cigar().iter().any(|op| matches!(op, Cigar::HardClip(_))))
            .cloned();
        let mut refined = refine_read_records(records, None, tid_names, fa_hash, bg_map,
            match_score, ref_miss_penalty, read_miss_penalty, prob_cutoff,
            designate_mapq, is_set_y, cpg);
        if let Some(original) = original.as_ref() {
            crate::align::engine::restore_primary_read(&mut refined, original, None);
        }
        crate::align::engine::refresh_sa_tags(&mut refined, &target_names);
        refined
    }).collect()
}

pub fn parse_bam(
    input_bam: &String,
    fa: &Option<String>,
    bg: &Option<String>,
    match_score: i32,
    ref_miss_penalty: i32,
    read_miss_penalty: i32,
    ref_prob_cutoff: f64,
    prob_cutoff: u8,
    designate_mapq: u8,
    is_set_y: bool,
    cpg: bool,
    output_secondary: bool,
    output_bam: &String,
    threads: usize,
) {
    let chunksize: usize = 50_000;
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };

    let fa_hash: HashMap<String, String> = if let Some(fa) = fa {
        let fasta = Fastx::new(fa);
        fasta.get_chrom_seqs().unwrap()
    } else {
        HashMap::new()
    };
    let bg_map: HashMap<String, HashSet<i64>> = if let Some(bg) = bg {
        parse_bedgraph(bg, ref_prob_cutoff).expect("Failed to parse bedgraph")
    } else {
        HashMap::new()
    };

    let header = Header::from_template(bam.header());
    let headerview = HeaderView::from_header(&header);
    let tid_names: Vec<String> = (0..headerview.target_count())
        .map(|i| String::from_utf8_lossy(headerview.tid2name(i as u32)).to_string())
        .collect();
    let _ = bam.set_threads(threads);

    let (wtx, wrx) = crossbeam_channel::bounded::<Vec<Record>>(10_000);

    let out_path = output_bam.clone();
    let header_for_writer = header;
    let writer_threads = threads as i32;
    let writer_handle = std::thread::spawn(move || {
        let mut writer = Writer::from_path(&out_path, &header_for_writer, bam::Format::Bam)
            .expect("Failed to create BAM writer");
        let _ = writer.set_threads(writer_threads.max(1).try_into().unwrap());
        while let Ok(group) = wrx.recv() {
            for rec in group {
                if !output_secondary && rec.is_secondary() {
                    continue;
                }
                writer.write(&rec).expect("Failed to write record");
            }
        }
    });

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} Loaded {pos} records | speed {per_sec} | elapsed {elapsed}",
        )
        .unwrap(),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100));
    pb.set_draw_target(ProgressDrawTarget::stderr_with_hz(10));
    let mut seen: u64 = 0;

    let mut previous_read: Vec<u8> = Vec::new();
    let mut idx: usize = 0;
    let mut batch: Vec<(usize, Vec<Record>)> = Vec::with_capacity(chunksize);
    let mut records_to_process: Vec<Record> = Vec::new();

    for rec in bam.records() {
        seen += 1;
        if seen % 10_000 == 0 {
            pb.set_position(seen);
        }
        let record = rec.expect("Failed to read record");
        if record.is_unmapped() {
            continue;
        }
        let qname = record.qname();
        if previous_read.is_empty() {
            previous_read.extend_from_slice(qname);
            records_to_process.push(record);
            continue;
        }
        if qname != &previous_read[..] {
            batch.push((idx, std::mem::take(&mut records_to_process)));
            if batch.len() >= chunksize {
                let flat = process_batch_with_rayon(
                    std::mem::take(&mut batch),
                    &tid_names,
                    &fa_hash,
                    &bg_map,
                    match_score,
                    ref_miss_penalty,
                    read_miss_penalty,
                    prob_cutoff,
                    designate_mapq,
                    is_set_y,
                    cpg,
                );
                if !flat.is_empty() {
                    wtx.send(flat).unwrap();
                }
            }
            idx += 1;
            previous_read.clear();
            previous_read.extend_from_slice(qname);
        }
        records_to_process.push(record);
    }

    if !records_to_process.is_empty() {
        batch.push((idx, std::mem::take(&mut records_to_process)));
    }
    if !batch.is_empty() {
        let flat = process_batch_with_rayon(
            std::mem::take(&mut batch),
            &tid_names,
            &fa_hash,
            &bg_map,
            match_score,
            ref_miss_penalty,
            read_miss_penalty,
            prob_cutoff,
            designate_mapq,
            is_set_y,
            cpg,
        );
        if !flat.is_empty() {
            wtx.send(flat).unwrap();
        }
    }

    drop(wtx);
    writer_handle.join().expect("Writer thread panicked");
    pb.set_position(seen);
    pb.finish_with_message("Finished processing BAM file.");

    log::info!(
        "Successfully refined alignments and wrote to {}",
        output_bam
    );
}

/// Options shared by the online aligner and the standalone methylation refiner.
#[derive(Debug, Clone)]
pub struct MethConfig {
    pub bed: String,
    pub match_score: i32,
    pub ref_penalty: i32,
    pub read_penalty: i32,
    pub ref_prob_cutoff: f64,
    pub prob_cutoff: u8,
    pub designate_mapq: u8,
    pub cpg: bool,
}

pub struct MethRefiner {
    config: MethConfig,
    reference: HashMap<String, String>,
    sites: HashMap<String, HashSet<i64>>,
}

impl MethRefiner {
    pub fn new(config: MethConfig, sequences: &[(String, Vec<u8>)]) -> anyResult<Self> {
        Self::from_owned(config, sequences.to_vec())
    }

    /// Reuse the FASTA buffers after indexing instead of retaining another copy.
    pub(crate) fn from_owned(config: MethConfig, sequences: Vec<(String, Vec<u8>)>) -> anyResult<Self> {
        anyhow::ensure!(config.ref_prob_cutoff.is_finite()
            && (0.0..=100.0).contains(&config.ref_prob_cutoff),
            "--meth-ref-prob-cutoff must be between 0 and 100");
        anyhow::ensure!((1..=60).contains(&config.designate_mapq),
            "--meth-designate-mapq must be between 1 and 60");
        let sites = parse_bedgraph(&config.bed, config.ref_prob_cutoff)?;
        let reference_start = (std::env::var("CPHASING_METH_TIMING").as_deref() == Ok("1"))
            .then(Instant::now);
        let reference: HashMap<String, String> = sequences.into_iter().map(|(name, seq)| {
            let mut sequence = String::from_utf8(seq).unwrap_or_else(|error| {
                String::from_utf8_lossy(error.as_bytes()).into_owned()
            });
            sequence.make_ascii_uppercase();
            (name, sequence)
        }).collect();
        if let Some(start) = reference_start {
            log::info!("METH_TIMING stage=reference_prepare seconds={:.9} contigs={}",
                start.elapsed().as_secs_f64(), reference.len());
        }
        anyhow::ensure!(sites.keys().any(|name| reference.contains_key(name)),
            "no reference methylation sites remain on matching contigs; check --meth-bed and --meth-ref-prob-cutoff");
        Ok(Self { config, reference, sites })
    }

    pub fn refine(&self, records: Vec<Record>, original: &Record, names: &[String]) -> Vec<Record> {
        let c = &self.config;
        refine_read_records(records, Some(original), names, &self.reference, &self.sites,
            c.match_score, c.ref_penalty, c.read_penalty, c.prob_cutoff,
            c.designate_mapq, false, c.cpg)
    }
}
