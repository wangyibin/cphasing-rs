#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as anyResult;
use crossbeam_channel::{Receiver, Sender, bounded};
use rayon::prelude::*;
use rust_htslib::bam::{
    self, Header, HeaderView, Read, Reader, Record, Writer, header::HeaderRecord, record::Aux,
    record::Cigar, record::CigarString, record::CigarStringView,
};
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;
use std::thread;

use crate::core::BaseTable;
use crate::core::{common_reader, common_writer};

#[derive(Debug, Clone)]
pub struct PAFLine {
    pub query: String,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub query_strand: char,
    pub target: String,
    pub target_length: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub match_n: u32,
    pub alignment_length: u32,
    pub mapq: u8,
    pub tags: Vec<String>,
}

impl PAFLine {
    fn new(fields: Vec<String>) -> Self {
        PAFLine {
            query: fields[0].clone(),
            query_length: fields[1].parse::<u32>().unwrap(),
            query_start: fields[2].parse::<u32>().unwrap(),
            query_end: fields[3].parse::<u32>().unwrap(),
            query_strand: fields[4].parse::<char>().unwrap(),
            target: fields[5].clone(),
            target_length: fields[6].parse::<u64>().unwrap(),
            target_start: fields[7].parse::<u64>().unwrap(),
            target_end: fields[8].parse::<u64>().unwrap(),
            match_n: fields[9].parse::<u32>().unwrap(),
            alignment_length: fields[10].parse::<u32>().unwrap(),
            mapq: fields[11].parse::<u8>().unwrap(),
            tags: fields[12..].to_vec(),
        }
    }

    fn get_tag(&self, tag: &str) -> String {
        let mut value = String::new();
        for t in &self.tags {
            if t.starts_with(tag) {
                value = t.to_string();
                break;
            }
        }
        value
    }

    fn is_secondary(&self) -> bool {
        let tag = self.get_tag("tp:A:");
        if tag == "tp:A:S" {
            return true;
        } else {
            return false;
        }
    }

    fn is_primary(&self) -> bool {
        let tag = self.get_tag("tp:A:");
        if tag == "tp:A:P" {
            return true;
        } else {
            return false;
        }
    }

    fn from(paf_line: &PAFLine) -> Self {
        PAFLine {
            query: paf_line.query.clone(),
            query_length: paf_line.query_length,
            query_start: paf_line.query_start,
            query_end: paf_line.query_end,
            query_strand: paf_line.query_strand,
            target: paf_line.target.clone(),
            target_length: paf_line.target_length,
            target_start: paf_line.target_start,
            target_end: paf_line.target_end,
            match_n: paf_line.match_n,
            alignment_length: paf_line.alignment_length,
            mapq: paf_line.mapq,
            tags: paf_line.tags.clone(),
        }
    }

    fn to_string(&self) -> String {
        let mut fields: Vec<String> = Vec::new();
        fields.push(self.query.clone());
        fields.push(self.query_length.to_string());
        fields.push(self.query_start.to_string());
        fields.push(self.query_end.to_string());
        fields.push(self.query_strand.to_string());
        fields.push(self.target.clone());
        fields.push(self.target_length.to_string());
        fields.push(self.target_start.to_string());
        fields.push(self.target_end.to_string());
        fields.push(self.match_n.to_string());
        fields.push(self.alignment_length.to_string());
        fields.push(self.mapq.to_string());
        for t in &self.tags {
            fields.push(t.to_string());
        }
        fields.join("\t")
    }

    fn score(&self) -> i32 {
        for t in &self.tags {
            if t.starts_with("AS:i:") {
                if let Ok(v) = t[5..].parse::<i32>() {
                    return v;
                }
            }
        }
        self.match_n as i32
    }
}

fn parse_paf_line_from_str(line: &str) -> Option<PAFLine> {
    let mut it = line.split('\t');
    let query = it.next()?.to_string();
    let query_length = it.next()?.parse::<u32>().ok()?;
    let query_start = it.next()?.parse::<u32>().ok()?;
    let query_end = it.next()?.parse::<u32>().ok()?;
    let query_strand = it.next()?.chars().next().unwrap_or('+');
    let target = it.next()?.to_string();
    let target_length = it.next()?.parse::<u64>().ok()?;
    let target_start = it.next()?.parse::<u64>().ok()?;
    let target_end = it.next()?.parse::<u64>().ok()?;
    let match_n = it.next()?.parse::<u32>().ok()?;
    let alignment_length = it.next()?.parse::<u32>().ok()?;
    let mapq = it.next()?.parse::<u8>().ok()?;

    let tags: Vec<String> = it.map(|s| s.to_string()).collect();

    Some(PAFLine {
        query,
        query_length,
        query_start,
        query_end,
        query_strand,
        target,
        target_length,
        target_start,
        target_end,
        match_n,
        alignment_length,
        mapq,
        tags,
    })
}

#[derive(Debug)]
pub struct PAFTable {
    file: String,
}

impl BaseTable for PAFTable {
    fn new(name: &String) -> PAFTable {
        PAFTable { file: name.clone() }
    }

    fn file_name(&self) -> Cow<'_, str> {
        let path = Path::new(&self.file);
        path.file_name().expect("REASON").to_string_lossy()
    }

    fn prefix(&self) -> String {
        let binding = self.file_name().to_string();
        let file_path = Path::new(&binding);
        let file_prefix = file_path.file_stem().unwrap().to_str().unwrap();

        (*file_prefix).to_string()
    }
}

impl PAFTable {
    pub fn parse(&self) -> anyResult<Box<dyn BufRead + Send + 'static>> {
        let reader = common_reader(&self.file);
        Ok(reader)
    }
}

#[derive(Debug, Clone)]
pub struct PAFReadUnit {
    data: Vec<PAFLine>,
}

impl PAFReadUnit {
    fn new() -> Self {
        PAFReadUnit { data: Vec::new() }
    }

    fn clear(&mut self) {
        self.data.clear();
    }
}

impl Default for PAFReadUnit {
    fn default() -> Self {
        PAFReadUnit::new()
    }
}

pub struct PAFAlignmentUnit {
    Primary: Vec<PAFLine>,
    Secondary: Vec<Vec<PAFLine>>,
}

impl PAFAlignmentUnit {
    fn new() -> Self {
        PAFAlignmentUnit {
            Primary: Vec::new(),
            Secondary: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.Primary.clear();
        self.Secondary.clear();
    }

    fn add_primary(&mut self, record: PAFLine) {
        self.Primary.push(record);
        self.Secondary.push(Vec::new());
    }

    fn add_secondary(&mut self, record: PAFLine, idx: usize) {
        // idx: index of primary alignment
        self.Secondary[idx].push(record);
    }

    fn is_empty(&self) -> bool {
        self.Primary.is_empty()
    }

    fn read_id(&self) -> String {
        self.Primary[0].query.clone()
    }

    fn rescue(&mut self, mapq: u8) {
        let mut anchor_counts: HashMap<&str, u32> = HashMap::with_capacity(8);

        for r in &self.Primary {
            if r.mapq >= mapq {
                *anchor_counts.entry(&r.target).or_insert(0) += 1;
            }
        }

        if anchor_counts.is_empty() {
            return;
        }

        let best_anchor = anchor_counts
            .iter()
            .max_by_key(|&(_, count)| count)
            .map(|(k, _)| k.to_string())
            .unwrap();

        for idx in 0..self.Primary.len() {
            if self.Primary[idx].mapq >= mapq {
                continue;
            }
            let (p_score, p_identity, is_target_best) = {
                let p = &self.Primary[idx];
                let score = p.score();
                let ident = if p.alignment_length > 0 {
                    p.match_n as f32 / p.alignment_length as f32
                } else {
                    0.0
                };
                (score, ident, p.target == best_anchor)
            };

            let mut best_cand_idx: Option<usize> = None;
            let mut best_cand_score = i32::MIN;

            if is_target_best {
                best_cand_idx = Some(0);
                best_cand_score = p_score;
            }

            let s_list = &self.Secondary[idx];
            for (i, s) in s_list.iter().enumerate() {
                if s.target != best_anchor {
                    continue;
                }

                let s_score = s.score();
                let s_identity = if s.alignment_length > 0 {
                    s.match_n as f32 / s.alignment_length as f32
                } else {
                    0.0
                };

                if (s_score as f32) >= (p_score as f32) && s_identity >= (p_identity * 0.90) {
                    if best_cand_idx.is_none() || s_score > best_cand_score {
                        best_cand_score = s_score;
                        best_cand_idx = Some(i + 1);
                    }
                }
            }

            if let Some(cidx) = best_cand_idx {
                if cidx == 0 {
                    self.Primary[idx].mapq = mapq.min(10);
                } else {
                    let mut new_p = self.Secondary[idx].swap_remove(cidx - 1);
                    new_p.mapq = mapq.min(10);

                    for t in new_p.tags.iter_mut() {
                        if t.starts_with("tp:A:") {
                            *t = String::from("tp:A:P");
                            break;
                        }
                    }

                    self.Primary[idx] = new_p;
                }
            }
        }
    }
}

pub fn parse_paf_read_unit(read_unit: &mut PAFReadUnit) -> PAFAlignmentUnit {
    let mut idx: u64 = 0;
    let mut au = PAFAlignmentUnit::new();
    for r in read_unit.data.drain(..) {
        if !r.is_secondary() {
            au.add_primary(r);
            idx += 1;
        } else {
            if idx == 0 {
                log::warn!("Secondary alignment missing primary, skipped.");
                continue;
            }
            let idx2 = idx - 1;
            au.add_secondary(r, idx2 as usize);
        }
    }
    au
}

fn process_paf_batch_with_rayon(batch: Vec<(usize, Vec<PAFLine>)>, mapq: u8) -> Vec<PAFLine> {
    batch
        .into_par_iter()
        .map(|(_, records)| {
            let mut ru = PAFReadUnit { data: records };
            let mut au = parse_paf_read_unit(&mut ru);
            au.rescue(mapq);
            au.Primary
        })
        .flatten()
        .collect()
}

pub fn read_paf(input_paf: &String, mapq: u8, output: &String) {
    use std::io::BufRead;
    use std::sync::{
        Arc,
        atomic::{AtomicU64, AtomicUsize, Ordering},
    };
    use std::time::Duration;

    let rdr = common_reader(&input_paf);
    let (wtx, wrx) = crossbeam_channel::bounded::<Vec<u8>>(100_000);
    let out_path = output.clone();

    let writer_handle = std::thread::spawn(move || {
        let wtr = common_writer(&out_path);
        let mut writer = std::io::BufWriter::new(wtr);
        while let Ok(buf) = wrx.recv() {
            let _ = writer.write_all(&buf);
        }
        let _ = writer.flush();
    });

    let chunksize: usize = 10_000;
    let mut batch: Vec<Vec<String>> = Vec::with_capacity(chunksize);
    let mut records_to_process: Vec<String> = Vec::with_capacity(128);

    let mut previous_read_id = String::new();

    let total_reads = Arc::new(AtomicU64::new(0));
    let total_unmapped = Arc::new(AtomicU64::new(0));
    let pending = Arc::new(AtomicUsize::new(0));
    let send_arc = Arc::new(wtx);

    for line_res in rdr.lines() {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => continue,
        };
        if line.is_empty() {
            continue;
        }

        let tab_pos = line.find('\t');
        let read_id = match tab_pos {
            Some(pos) => &line[..pos],
            None => continue,
        };

        if previous_read_id.is_empty() {
            previous_read_id = read_id.to_string();
        }

        if read_id != previous_read_id {
            batch.push(std::mem::take(&mut records_to_process));
            records_to_process = Vec::with_capacity(128);

            if batch.len() >= chunksize {
                let to_proc = std::mem::take(&mut batch);
                let send = send_arc.clone();
                let pend = pending.clone();
                let t_reads = total_reads.clone();
                let t_unmapped = total_unmapped.clone();

                pend.fetch_add(1, Ordering::SeqCst);

                rayon::spawn_fifo(move || {
                    let mut outbuf: Vec<u8> = Vec::with_capacity(to_proc.len() * 300);

                    for raw_lines in to_proc {
                        let mut parsed_records: Vec<PAFLine> = Vec::with_capacity(raw_lines.len());
                        let mut local_unmapped = 0;

                        for raw_line in raw_lines {
                            if let Some(p) = parse_paf_line_from_str(&raw_line) {
                                if p.target == "*" {
                                    local_unmapped += 1;
                                } else {
                                    parsed_records.push(p);
                                }
                            }
                        }

                        t_reads.fetch_add(1, Ordering::Relaxed);
                        if local_unmapped > 0 && parsed_records.is_empty() {
                            t_unmapped.fetch_add(1, Ordering::Relaxed);
                            continue;
                        }

                        let mut ru = PAFReadUnit {
                            data: parsed_records,
                        };
                        if !ru.data.is_empty() {
                            let mut au = parse_paf_read_unit(&mut ru);
                            au.rescue(mapq);
                            for r in au.Primary {
                                let s = r.to_string();
                                outbuf.extend_from_slice(s.as_bytes());
                                outbuf.push(b'\n');
                            }
                        }
                    }
                    let _ = send.send(outbuf);
                    pend.fetch_sub(1, Ordering::SeqCst);
                });
            }

            previous_read_id = read_id.to_string();
        }

        records_to_process.push(line);
    }

    if !records_to_process.is_empty() {
        batch.push(std::mem::take(&mut records_to_process));
    }
    if !batch.is_empty() {
        let to_proc = std::mem::take(&mut batch);
        let send = send_arc.clone();
        let pend = pending.clone();
        let t_reads = total_reads.clone();
        let t_unmapped = total_unmapped.clone();

        pend.fetch_add(1, Ordering::SeqCst);

        rayon::spawn_fifo(move || {
            let mut outbuf: Vec<u8> = Vec::with_capacity(to_proc.len() * 300);

            for raw_lines in to_proc {
                let mut parsed_records: Vec<PAFLine> = Vec::with_capacity(raw_lines.len());
                let mut local_unmapped = 0;

                for raw_line in raw_lines {
                    if let Some(p) = parse_paf_line_from_str(&raw_line) {
                        if p.target == "*" {
                            local_unmapped += 1;
                        } else {
                            parsed_records.push(p);
                        }
                    }
                }

                t_reads.fetch_add(1, Ordering::Relaxed);
                if local_unmapped > 0 && parsed_records.is_empty() {
                    t_unmapped.fetch_add(1, Ordering::Relaxed);
                    continue;
                }

                let mut ru = PAFReadUnit {
                    data: parsed_records,
                };
                if !ru.data.is_empty() {
                    let mut au = parse_paf_read_unit(&mut ru);
                    au.rescue(mapq);
                    for r in au.Primary {
                        let s = r.to_string();
                        outbuf.extend_from_slice(s.as_bytes());
                        outbuf.push(b'\n');
                    }
                }
            }
            let _ = send.send(outbuf);
            pend.fetch_sub(1, Ordering::SeqCst);
        });
    }

    loop {
        if pending.load(Ordering::SeqCst) == 0 {
            break;
        }
        std::thread::sleep(Duration::from_millis(10));
    }

    drop(send_arc);
    let _ = writer_handle.join();

    log::info!(
        "Processed reads: {}, unmapped: {}",
        total_reads.load(Ordering::Relaxed),
        total_unmapped.load(Ordering::Relaxed)
    );
}

pub fn contact_scores_from_contacts(
    contacts_path: &str,
    hv: &HeaderView,
    topk_per_tid: usize,
) -> HashMap<(i32, i32), f32> {
    let mut reader = common_reader(&contacts_path.to_string());

    let mut counts: HashMap<(i32, i32), u64> = HashMap::new();

    for line in reader.lines() {
        let line = match line {
            Ok(s) => s,
            Err(_) => continue,
        };
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let cols: Vec<&str> = line
            .split(|c| c == '\t' || c == ' ')
            .filter(|s| !s.is_empty())
            .collect();
        if cols.len() < 2 {
            continue;
        }

        let mut tids: Vec<i32> = Vec::with_capacity(2);
        for token in &cols {
            if let Some(t) = hv.tid(token.as_bytes()) {
                tids.push(t as i32);
                if tids.len() == 2 {
                    break;
                }
            }
        }
        if tids.len() < 2 {
            continue;
        }

        let mut cnt: u64 = 1;
        for token in cols.iter().rev() {
            if let Ok(v) = u64::from_str(token) {
                cnt = v;
                break;
            }
            if let Ok(vf) = f64::from_str(token) {
                if vf.is_finite() && vf >= 0.0 {
                    cnt = vf as u64;
                    break;
                }
            }
        }

        let (a, b) = if tids[0] <= tids[1] {
            (tids[0], tids[1])
        } else {
            (tids[1], tids[0])
        };
        *counts.entry((a, b)).or_insert(0) += cnt;
    }

    let mut rows: HashMap<i32, Vec<(i32, f32)>> = HashMap::new();
    for ((a, b), c) in counts {
        let w = (c as f32).ln_1p();
        rows.entry(a).or_default().push((b, w));
        if a != b {
            rows.entry(b).or_default().push((a, w));
        }
    }

    let mut scores: HashMap<(i32, i32), f32> = HashMap::new();
    for (tid, mut lst) in rows {
        lst.sort_by(|x, y| y.1.partial_cmp(&x.1).unwrap());
        if lst.len() > topk_per_tid {
            lst.truncate(topk_per_tid);
        }
        let sum = lst.iter().map(|x| x.1).sum::<f32>().max(1e-6);
        for (other, w) in lst {
            let s = w / sum;
            let key = if tid <= other {
                (tid, other)
            } else {
                (other, tid)
            };
            let e = scores.entry(key).or_insert(0.0);

            *e = if *e == 0.0 { s } else { (*e + s) * 0.5 };
        }
    }
    scores
}

#[derive(Clone, Debug)]
pub struct ReadUnit {
    data: Vec<Record>,
}

impl ReadUnit {
    fn new() -> Self {
        ReadUnit { data: Vec::new() }
    }

    fn clear(&mut self) {
        self.data.clear();
    }
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct AlignmentUnit {
    Primary: Vec<Record>,
    Secondary: Vec<Vec<Record>>,
}

#[allow(unused_variables, unused_assignments)]
impl AlignmentUnit {
    fn new() -> Self {
        AlignmentUnit {
            Primary: Vec::new(),
            Secondary: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.Primary.clear();
        self.Secondary.clear();
    }

    fn add_primary(&mut self, record: Record) {
        self.Primary.push(record);
        self.Secondary.push(Vec::new());
    }

    fn add_secondary(&mut self, record: Record, idx: usize) {
        // idx: index of primary alignment
        self.Secondary[idx].push(record);
    }

    fn is_empty(&self) -> bool {
        self.Primary.is_empty()
    }

    fn read_id(&self) -> String {
        String::from_utf8(self.Primary[0].qname().to_vec())
            .unwrap()
            .to_string()
    }
}

fn parse_read_unit(read_unit: &ReadUnit) -> AlignmentUnit {
    let mut idx: u64 = 0;
    let mut au = AlignmentUnit::new();
    for r in &read_unit.data {
        let read_id = String::from_utf8(r.qname().to_vec()).unwrap().to_string();

        if !r.is_secondary() {
            au.add_primary(r.clone());
            idx += 1;
        } else {
            if idx == 0 {
                log::warn!(
                    "Secondary alignment `{:?}` could not found primary, \
                            skipped.\
                            The input bam should be sorted by read name.",
                    read_id
                );
                continue;
            }
            let idx2 = idx - 1;
            au.add_secondary(r.clone(), idx2 as usize);
        }
    }
    au
}

#[derive(Clone, Default)]
pub struct ContactGraph {
    scores: HashMap<(i32, i32), f32>,
}
impl ContactGraph {
    pub fn new() -> Self {
        Self {
            scores: HashMap::new(),
        }
    }
    pub fn with_scores(scores: HashMap<(i32, i32), f32>) -> Self {
        Self { scores }
    }
    #[inline]
    pub fn score(&self, a: i32, b: i32) -> f32 {
        if a < 0 || b < 0 {
            return 0.0;
        }
        *self
            .scores
            .get(&(a, b))
            .or_else(|| self.scores.get(&(b, a)))
            .unwrap_or(&0.0)
    }
}

#[derive(Clone, Default)]
struct PairAlignmentUnit {
    left_prim: Vec<Record>,
    left_sec: Vec<Record>,
    right_prim: Vec<Record>,
    right_sec: Vec<Record>,
}
impl PairAlignmentUnit {
    fn from_read_unit(ru: &ReadUnit) -> Self {
        let mut p = PairAlignmentUnit::default();
        for r in &ru.data {
            if r.is_secondary() {
                if r.is_first_in_template() {
                    p.left_sec.push(r.clone());
                } else if r.is_last_in_template() {
                    p.right_sec.push(r.clone());
                } else {
                    p.left_sec.push(r.clone());
                }
            } else {
                if r.is_first_in_template() {
                    p.left_prim.push(r.clone());
                } else if r.is_last_in_template() {
                    p.right_prim.push(r.clone());
                } else {
                    p.left_prim.push(r.clone());
                }
            }
        }
        p
    }
    fn has_both_mates(&self) -> bool {
        (!self.left_prim.is_empty() || !self.left_sec.is_empty())
            && (!self.right_prim.is_empty() || !self.right_sec.is_empty())
    }

    fn best_pair_with_contacts(
        &self,
        graph: &ContactGraph,
        alpha: f32,
        beta: f32,
        k: usize,
    ) -> Option<(Record, Record)> {
        let mut left: Vec<Record> = self
            .left_prim
            .iter()
            .cloned()
            .chain(self.left_sec.iter().cloned())
            .collect();

        left.sort_by(|a, b| get_score(b).cmp(&get_score(a)));
        if left.len() > k {
            left.truncate(k);
        }

        let mut right: Vec<Record> = self
            .right_prim
            .iter()
            .cloned()
            .chain(self.right_sec.iter().cloned())
            .collect();
        right.sort_by(|a, b| get_score(b).cmp(&get_score(a)));
        if right.len() > k {
            right.truncate(k);
        }

        if left.is_empty() || right.is_empty() {
            return None;
        }

        let max_score_left = get_score(&left[0]);
        let max_score_right = get_score(&right[0]);

        let mut best = None;
        let mut best_score = f32::NEG_INFINITY;

        for l in &left {
            if (get_score(l) as f32) < (max_score_left as f32 * 0.95) {
                continue;
            }

            let tl = l.tid();
            for r in &right {
                if (get_score(r) as f32) < (max_score_right as f32 * 0.95) {
                    continue;
                }

                let tr = r.tid();

                let s = beta * graph.score(tl, tr) + alpha * (l.mapq() as f32 + r.mapq() as f32);
                if s > best_score {
                    best_score = s;
                    best = Some((l.clone(), r.clone()));
                }
            }
        }
        best
    }
}

impl AlignmentUnit {
    fn rescue(&mut self, mapq: u8) {
        let mut high_mapq: HashMap<u64, u32> = HashMap::new();

        for r in &self.Primary {
            let target: u64 = r.tid().try_into().unwrap();

            if r.mapq() >= mapq && !r.is_secondary() {
                *high_mapq.entry(target).or_insert(0) += 1;
            }
        }

        if high_mapq.is_empty() {
            return;
        }

        let mut best_anchor_target = 0;
        let mut max_count = 0;
        for (target, count) in &high_mapq {
            if *count > max_count {
                max_count = *count;
                best_anchor_target = *target;
            }
        }

        for (p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
            if p.mapq() >= mapq {
                continue;
            }

            let p_score = get_score(p);
            let p_tid: u64 = p.tid().try_into().unwrap();

            let mut best_candidate_idx: Option<usize> = None;
            let mut best_candidate_score = 0;

            if p_tid == best_anchor_target {
                best_candidate_idx = Some(0);
                best_candidate_score = p_score;
            } else {
                for (j, r) in s.iter().enumerate() {
                    let t: u64 = r.tid().try_into().unwrap();
                    if t == best_anchor_target {
                        let s_score = get_score(r);
                        if s_score as f32 >= (p_score as f32 * 0.95) {
                            if best_candidate_idx.is_none() || s_score > best_candidate_score {
                                best_candidate_score = s_score;
                                best_candidate_idx = Some(j + 1);
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
                    new_record.set_flags(flag);
                    *p = new_record;
                }
            }
        }
    }
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

pub fn read_bam(
    input_bam: &String,
    mapq: u8,
    output: &String,
    workers: usize,
    contacts: Option<String>,
) {
    let mut bam = Reader::from_path(input_bam).unwrap();
    let _ = bam.set_threads(std::cmp::max(1, workers));
    let in_header = Header::from_template(bam.header());
    let hv = HeaderView::from_header(&in_header);
    let out_header = Header::from_template(&hv);

    let (unit_tx, unit_rx): (Sender<ReadUnit>, Receiver<ReadUnit>) = bounded(1024);
    let (rec_tx, rec_rx): (Sender<Vec<Record>>, Receiver<Vec<Record>>) = bounded(1024);

    let out_path = output.clone();
    let writer_handle = thread::spawn(move || {
        let mut writer = Writer::from_path(&out_path, &out_header, bam::Format::Bam).unwrap();
        for batch in rec_rx.iter() {
            for rec in batch {
                writer.write(&rec).unwrap();
            }
        }
    });
    let graph = if let Some(path) = contacts {
        log::info!("Loading contacts prior from {}", path);
        let scores = contact_scores_from_contacts(&path, &hv, 64);
        ContactGraph::with_scores(scores)
    } else {
        ContactGraph::new()
    };

    // let graph = contact_scores.map(ContactGraph::with_scores).unwrap_or_else(ContactGraph::new);
    let alpha: f32 = 0.05;
    let beta: f32 = 1.0;
    let topk: usize = 4;

    let mut handles = Vec::new();
    for _ in 0..workers {
        let unit_rx = unit_rx.clone();
        let rec_tx = rec_tx.clone();
        let graph = graph.clone();
        let h = thread::spawn(move || {
            while let Ok(ru) = unit_rx.recv() {
                let pair = PairAlignmentUnit::from_read_unit(&ru);
                let out_records: Vec<Record> = if pair.has_both_mates() {
                    if let Some((mut l, mut r)) =
                        pair.best_pair_with_contacts(&graph, alpha, beta, topk)
                    {
                        l.set_mtid(r.tid());
                        l.set_mpos(r.pos());
                        l.set_insert_size(0);
                        r.set_mtid(l.tid());
                        r.set_mpos(l.pos());
                        r.set_insert_size(0);
                        if l.mapq() < mapq {
                            l.set_mapq(mapq.max(1));
                        }
                        if r.mapq() < mapq {
                            r.set_mapq(mapq.max(1));
                        }
                        vec![l, r]
                    } else {
                        let mut au = parse_read_unit(&ru);
                        au.rescue(mapq);
                        au.Primary
                    }
                } else {
                    let mut au = parse_read_unit(&ru);
                    au.rescue(mapq);
                    au.Primary
                };

                if !out_records.is_empty() {
                    let _ = rec_tx.send(out_records);
                }
            }
        });
        handles.push(h);
    }
    drop(rec_tx);

    let mut old_read_id = String::new();
    let mut ru = ReadUnit::new();

    for r in bam.records() {
        let record = r.unwrap();
        if record.is_unmapped() {
            continue;
        }
        let read_id = String::from_utf8(record.qname().to_vec()).unwrap();

        if old_read_id != read_id {
            if !ru.data.is_empty() {
                unit_tx.send(ru).unwrap();
                ru = ReadUnit::new();
            }
            ru.data.push(record);
            old_read_id = read_id;
        } else {
            ru.data.push(record);
        }
    }
    if !ru.data.is_empty() {
        unit_tx.send(ru).unwrap();
    }
    drop(unit_tx);

    for h in handles {
        let _ = h.join();
    }

    // drop(rec_rx);
    let _ = writer_handle.join();
}
