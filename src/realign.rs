#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as anyResult;
use crossbeam_channel::{bounded, Receiver, Sender};
use std::thread;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::str::FromStr;
use std::path::Path;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;

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
    pub cigar: String,
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
            cigar: fields[12].clone(),
            tags: fields[13..].to_vec(),
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
            cigar: paf_line.cigar.clone(),
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
        fields.push(self.cigar.clone());
        for t in &self.tags {
            fields.push(t.to_string());
        }
        fields.join("\t")
    }
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
        PAFReadUnit {
            data: Vec::new(),
        }
    }

    fn clear (&mut self) {
        self.data.clear();
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

    fn clear (&mut self) {
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
        // count map quality >= 1 {"target1": 2, "target2": 1}
        
        let mut high_mapq: HashMap<String, u32> = HashMap::new();
        let mut high_high_mapq: HashMap<String, u32> = HashMap::new();
        for r in &self.Primary {
            let target = r.target.clone();
            if r.mapq >= mapq {
                let count = high_mapq.entry(target.clone()).or_insert(0);
                *count += 1;
            }

            if r.mapq > 1 {
                let count = high_high_mapq.entry(target).or_insert(0);
                *count += 1;
            }
        }
        

        if high_mapq.len() == 0 {
            return;
        }

        {
        for (p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
            if p.mapq > 1 {
                continue;
            }

            let mut res: HashSet<String> = HashSet::new();
            let mut res_record_idx: HashMap<String, usize> = HashMap::new();
            let target = p.target.clone();
    
            if high_high_mapq.contains_key(&target) {
                res_record_idx.insert(target.clone(), 0);
                res.insert(target);
            }

            for (j, r) in s.iter().enumerate() {
                let target = r.target.clone();
                if high_high_mapq.contains_key(&target) {
                    res_record_idx.insert(target.clone(), j+1);
                    res.insert(target);
                }
            }

            let mut max_target = String::new();
            if res.len() == 0 {
                continue;
            } else if res.len() == 1 {
                let target = res.iter().next().unwrap();
                max_target = target.clone();
            } else {
                // max in res 
                let mut max = 0;
                for target in res {
                    let count = high_mapq.get(&target).unwrap();
                    if *count > max {
                        max = *count;
                        max_target = target;
                    }
                }
        
            }

            let record_idx = res_record_idx.get(&max_target).unwrap();

            high_mapq.entry(max_target).and_modify(|x| *x += 1);
            if record_idx == &0 {
                p.mapq = 1;
            } else {
                let mut new_record = PAFLine::from(&s[*record_idx - 1]);
                // println!("{:?}", p);
                // println!("{:?}", new_record);
                new_record.mapq = 2;
                new_record.tags = new_record.tags.iter().map(|x| {
                    if x.starts_with("tp:A:") {
                        "tp:A:P".to_string()
                    } else {
                        x.to_string()
                    }
                }).collect();
                *p = new_record;
            
            }
        }

        
        for (p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
            if p.mapq >= mapq {
                continue;
            }

            let mut res: HashSet<String> = HashSet::new();
            let mut res_record_idx: HashMap<String, usize> = HashMap::new();
            let target = p.target.clone();
    
            if high_mapq.contains_key(&target) {
                res_record_idx.insert(target.clone(), 0);
                res.insert(target);
            }

            for (j, r) in s.iter().enumerate() {
                let target = r.target.clone();
                if high_mapq.contains_key(&target) {
                    res_record_idx.insert(target.clone(), j+1);
                    res.insert(target);
                }
            }

            let mut max_target = String::new();
            if res.len() == 0 {
                continue;
            } else if res.len() == 1 {
                let target = res.iter().next().unwrap();
                max_target = target.clone();
            } else {
                // max in res 
                let mut max = 0;
                for target in res {
                    let count = high_mapq.get(&target).unwrap();
                    if *count > max {
                        max = *count;
                        max_target = target;
                    }
                }
            }
            
            let record_idx = res_record_idx.get(&max_target).unwrap();
            if record_idx == &0 {
                p.mapq = 1;
            } else {
                let mut new_record = PAFLine::from(&s[*record_idx - 1]);
                new_record.mapq = 1;
                new_record.tags = new_record.tags.iter().map(|x| {
                    if x.starts_with("tp:A:") {
                        "tp:A:P".to_string()
                    } else {
                        x.to_string()
                    }
                }).collect();
                *p = new_record;
            }

    }
    }
    }

}


pub fn parse_paf_read_unit(read_unit: &PAFReadUnit) -> PAFAlignmentUnit {
    let mut idx: u64 = 0;
    let mut au = PAFAlignmentUnit::new();
    for r in &read_unit.data {
        let read_id = r.query.clone();
        

        if !r.is_secondary() {
            au.add_primary(r.clone());
            idx += 1;
        } else {
            if idx == 0 {
                log::warn!("Secondary alignment `{:?}` could not found primary, \
                            skipped.\
                            The input bam should be sorted by read name.", read_id);
                continue;
            }
            let idx2 = idx - 1;
            au.add_secondary(r.clone(), idx2 as usize);
        }
        
    }
    au
}


pub fn read_paf(input_paf: &String, mapq: u8, output: &String) {
    let paf = PAFTable::new(input_paf);
    let parse_result = paf.parse();
    let rdr = match parse_result {
        Ok(v) => v,
        Err(error) => panic!("Error: Could not parse input file: {:?}", paf.file_name()),
    };

    let mut total_reads: u64 = 0;
    let mut total_alignments: u64 = 0;
    let mut total_unmapped: u64 = 0;
    let mut old_read_id = String::from("");

    let wtr = common_writer(output);
    let mut writer = BufWriter::new(wtr);
    
    let mut read_unit = PAFReadUnit::new();

    for line in rdr.lines() {
        let fields: Vec<String> = line.unwrap().split('\t').map(|x| x.to_string()).collect();
        
        assert!(fields.len() > 12, "Error: PAF file should have at least 12 columns");
        let paf_line: PAFLine = PAFLine::new(fields);

        if paf_line.target == "*" {
            total_unmapped += 1;
            total_reads += 1;
            continue;
        }

        let read_id = paf_line.query.clone();

        if old_read_id != paf_line.query {
            if old_read_id != "" {
                let mut au = parse_paf_read_unit(&read_unit);
                au.rescue(mapq);
                for r in au.Primary {
                    
                    writer.write(r.to_string().as_bytes()).unwrap();
                    writer.write(b"\n").unwrap();
                }
            }

            total_reads += 1;
            read_unit.clear();
            read_unit.data.push(paf_line);
        } else {
            read_unit.data.push(paf_line);
        }

        old_read_id = read_id;

        

    }
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

        // 从末尾找一个数字作为计数，找不到则按 1 计
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
            let key = if tid <= other { (tid, other) } else { (other, tid) };
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
        ReadUnit {
            data: Vec::new(),
        }
    }

    fn clear (&mut self) {
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

    fn clear (&mut self) {
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
                                .unwrap().to_string()
    }

//     fn rescue(&mut self, mapq: u8) {
        
//         let mut high_mapq: HashMap<u64, u32> = HashMap::new();
//         let mut high_high_mapq: HashMap<u64, u32> = HashMap::new();
//         for r in &self.Primary {
//             let target: u64 = r.tid().try_into().unwrap();
//             if r.mapq() >= mapq {
//                 let count = high_mapq.entry(target).or_insert(0);
//                 *count += 1;
//             } 

//             if r.mapq() > 1 {
//                 let count = high_high_mapq.entry(target).or_insert(0);
//                 *count += 1;
//             }
//         }
        

//         if high_mapq.len() == 0 {
//             return;
//         }
      
//         for (mut p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
//             if p.mapq() >= mapq {
//                 continue;
//             }

//             let mut res: HashSet<u64> = HashSet::new();
//             let mut res_record_idx: HashMap<u64, usize> = HashMap::new();
//             let target: u64 = p.tid().try_into().unwrap();
//             if high_mapq.contains_key(&target) {
//                 res_record_idx.insert(target, 0);
//                 res.insert(target);


//             }
            
//             for (j, r) in s.iter().enumerate() {
//                 let target: u64 = r.tid().try_into().unwrap();
//                 if high_mapq.contains_key(&target) {
//                     res_record_idx.insert(target, j);
//                     res.insert(target);
//                 }
//             }
            
//             let mut max_target = 0;
//             if res.len() == 0 {
//                 continue;
//             } else if res.len() == 1 {
//                 let target = res.iter().next().unwrap();
//                 max_target = *target;
//             } else {
//                 // max in res 
//                 let mut max = 0;
//                 for target in res {
//                     let count = high_mapq.get(&target).unwrap();
//                     if *count > max {
//                         max = *count;
//                         max_target = target;
//                     }
//                 }
            
//             }
//             let record_idx = res_record_idx.get(&max_target).unwrap();
            
//             if record_idx == &0 {
//                 p.set_tid(max_target.try_into().unwrap());
//                 p.set_mapq(1);
//             } else {
//                 let flag = p.flags();
//                 let r = &s[*record_idx - 1];
//                 let mut new_record = Record::from(r.clone());
//                 p = &mut new_record;
//                 p.set_mapq(1);
//                 p.set_flags(flag);
//             } 
            
//         }
    
    
//     }
 

}

fn parse_read_unit(read_unit: &ReadUnit) -> AlignmentUnit {
    let mut idx: u64 = 0;
    let mut au = AlignmentUnit::new();
    for r in &read_unit.data {
        let read_id = String::from_utf8(r.qname().to_vec())
                                .unwrap().to_string();
        

        if !r.is_secondary() {
            au.add_primary(r.clone());
            idx += 1;
        } else {
            if idx == 0 {
                log::warn!("Secondary alignment `{:?}` could not found primary, \
                            skipped.\
                            The input bam should be sorted by read name.", read_id);
                continue;
            }
            let idx2 = idx - 1;
            au.add_secondary(r.clone(), idx2 as usize);
        }
        
    }
    au
}


// pub fn read_bam(input_bam: &String, mapq: u8, output: &String) {
//     let mut bam = Reader::from_path(input_bam).unwrap();
//     let _ = bam.set_threads(8);
//     let bam_header = Header::from_template(bam.header());
//     let bam_header = HeaderView::from_header(&bam_header);

//     let mut total_reads: u64 = 0;
//     let total_alignments: u64 = 0;
//     let mut total_unmapped: u64 = 0;
//     let mut old_read_id = String::from("");

//     let mut read_unit = ReadUnit::new();

//     let header = Header::from_template(&bam_header);
//     let mut writer = Writer::from_path(output, &header, bam::Format::Bam).unwrap();

//     for r in bam.records() {
//         let record = r.unwrap();
//         let read_id = String::from_utf8(record.qname().to_vec())
//                                 .unwrap().to_string();

//         if record.is_unmapped() {
//             total_unmapped += 1;
//             total_reads += 1;
//             continue;
//         }

//         if old_read_id != read_id {
//             if old_read_id != "" {
//                 let mut au = parse_read_unit(&read_unit);
//                 au.rescue(mapq);
//                 for r in au.Primary {
//                     writer.write(&r).unwrap();
//                 }
//             }


//             total_reads += 1;
//             read_unit.clear();
//             read_unit.data.push(record);
//         } else {
//             read_unit.data.push(record);
//         }
      
//         old_read_id = read_id;
    

//     }

// }
// ...existing code...


// 高阶互作先验（稀疏图）
#[derive(Clone, Default)]
pub struct ContactGraph {
    scores: HashMap<(i32, i32), f32>,
}
impl ContactGraph {
    pub fn new() -> Self {
        Self { scores: HashMap::new() }
    }
    pub fn with_scores(scores: HashMap<(i32, i32), f32>) -> Self {
        Self { scores }
    }
    #[inline]
    pub fn score(&self, a: i32, b: i32) -> f32 {
        if a < 0 || b < 0 { return 0.0; }
        *self.scores.get(&(a, b))
            .or_else(|| self.scores.get(&(b, a)))
            .unwrap_or(&0.0)
    }
}

// 将一对读段拆分为 mate1/mate2 的候选集合
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
                if r.is_first_in_template() { p.left_sec.push(r.clone()); }
                else if r.is_last_in_template() { p.right_sec.push(r.clone()); }
                else { p.left_sec.push(r.clone()); }
            } else {
                if r.is_first_in_template() { p.left_prim.push(r.clone()); }
                else if r.is_last_in_template() { p.right_prim.push(r.clone()); }
                else { p.left_prim.push(r.clone()); }
            }
        }
        p
    }
    fn has_both_mates(&self) -> bool {
        (!self.left_prim.is_empty() || !self.left_sec.is_empty()) &&
        (!self.right_prim.is_empty() || !self.right_sec.is_empty())
    }
    fn best_pair_with_contacts(
        &self, graph: &ContactGraph, alpha: f32, beta: f32, k: usize,
    ) -> Option<(Record, Record)> {
        let mut left: Vec<Record> = self.left_prim.iter().cloned()
            .chain(self.left_sec.iter().cloned()).collect();
        left.sort_by(|a,b| b.mapq().cmp(&a.mapq()));
        if left.len() > k { left.truncate(k); }

        let mut right: Vec<Record> = self.right_prim.iter().cloned()
            .chain(self.right_sec.iter().cloned()).collect();
        right.sort_by(|a,b| b.mapq().cmp(&a.mapq()));
        if right.len() > k { right.truncate(k); }

        if left.is_empty() || right.is_empty() { return None; }

        let mut best = None;
        let mut best_score = f32::NEG_INFINITY;
        for l in &left {
            let tl = l.tid();
            for r in &right {
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

// 修正 rescue：次比对索引用 j+1 与取用时 -1 配对
impl AlignmentUnit {
    fn rescue(&mut self, mapq: u8) {
        let mut high_mapq: HashMap<u64, u32> = HashMap::new();
        let mut high_high_mapq: HashMap<u64, u32> = HashMap::new();
        for r in &self.Primary {
            let target: u64 = r.tid().try_into().unwrap();
            if r.mapq() >= mapq { *high_mapq.entry(target).or_insert(0) += 1; }
            if r.mapq() > 1 { *high_high_mapq.entry(target).or_insert(0) += 1; }
        }
        if high_mapq.is_empty() { return; }

        for (p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
            if p.mapq() >= mapq { continue; }

            let mut res: HashSet<u64> = HashSet::new();
            let mut res_record_idx: HashMap<u64, usize> = HashMap::new();
            let t0: u64 = p.tid().try_into().unwrap();
            if high_mapq.contains_key(&t0) { res_record_idx.insert(t0, 0); res.insert(t0); }
            for (j, r) in s.iter().enumerate() {
                let t: u64 = r.tid().try_into().unwrap();
                if high_mapq.contains_key(&t) { res_record_idx.insert(t, j + 1); res.insert(t); }
            }
            if res.is_empty() { continue; }

            let mut max_target = *res.iter().next().unwrap();
            let mut max_cnt = 0u32;
            for t in res {
                let c = *high_mapq.get(&t).unwrap_or(&0);
                if c > max_cnt { max_cnt = c; max_target = t; }
            }
            let record_idx = *res_record_idx.get(&max_target).unwrap_or(&0);
            if record_idx == 0 {
                p.set_tid(max_target.try_into().unwrap());
                p.set_mapq(1);
            } else {
                let flag = p.flags();
                let r = &s[record_idx - 1];
                let mut new_record = Record::from(r.clone());
                new_record.set_mapq(1);
                new_record.set_flags(flag);
                *p = new_record;
            }
        }
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
        let scores = contact_scores_from_contacts(&path, &hv, 64); // 每行保留 top-64
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
                    if let Some((mut l, mut r)) = pair.best_pair_with_contacts(&graph, alpha, beta, topk) {
                   
                        l.set_mtid(r.tid());
                        l.set_mpos(r.pos());
                        l.set_insert_size(0);
                        r.set_mtid(l.tid());
                        r.set_mpos(l.pos());
                        r.set_insert_size(0);
                        if l.mapq() < mapq { l.set_mapq(mapq.max(1)); }
                        if r.mapq() < mapq { r.set_mapq(mapq.max(1)); }
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
        if record.is_unmapped() { continue; }
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

    for h in handles { let _ = h.join(); }

    // drop(rec_rx);
    let _ = writer_handle.join();

}
