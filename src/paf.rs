#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as anyResult;
use crossbeam_channel::{bounded, unbounded, Sender, Receiver};
use crossbeam::{scope};
use std::thread;
use lazy_static::lazy_static;
use std::any::type_name;
use std::borrow::Cow;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::ffi::OsStr;
use std::path::Path;
use std::error::Error;
use std::sync::{Arc, Mutex, atomic::Ordering};
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::io::{ Read, Write, BufReader, BufRead };
use std::fmt::{ Write as FmtWrite };
use serde::{ Deserialize, Serialize};
use std::time::{Instant, Duration};
use rust_lapper::{Interval, Lapper};
use rayon::prelude::*;
use rand::prelude::*;

use crate::bed::Bed3; 
use crate::core::{ common_reader, common_writer, Prof };
use crate::core::BaseTable;
use crate::porec::PoreCRecordPlus as PoreCRecord;


const GROUP_BATCH_SIZE: usize = 100_000;
const RECORD_BATCH_SIZE: usize = 200_000;

const CHUNK_SIZE: usize = 10_000;
const WRITE_CHUNK: usize = 16 << 20; // 1 MB

const PASS_STR : &str = "pass";
const SINGLETON_STR : &str = "singleton";
const LOW_MQ_STR : &str = "low_mq";
const COMPLEX_STR : &str = "complex";



#[derive(Clone, Debug, Deserialize, Serialize)]
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

#[derive(Debug, Clone)]
pub struct ReadSummary {
    pass: u64,
    singleton: u64, 
    low_mq: u64,
    complex: u64,
    mapping: u64,
}

impl ReadSummary {
    fn new() -> ReadSummary {
        let pass: u64 = 0;
        let singleton: u64 = 0;
        let low_mq: u64 = 0;
        let complex: u64 = 0;
        let mapping: u64 = 0;

        ReadSummary {
            pass: pass,
            singleton: singleton,
            low_mq: low_mq,
            complex: complex,
            mapping: mapping,
        }
    }

    fn to_string(&self) -> String {
        let mut result = String::new();
        result.push_str(&format!("{}\t{}\n", "pass", self.pass));
        result.push_str(&format!("{}\t{}\n", "singleton", self.singleton));
        result.push_str(&format!("{}\t{}\n", "low_mq", self.low_mq));
        result.push_str(&format!("{}\t{}\n", "complex", self.complex));
        result.push_str(&format!("{}\t{}\n", "mapping", self.mapping));

        result 
    }

    fn save(&self, output: &String) {
        let mut wtr = common_writer(output);    

        let result: String = self.to_string();

        wtr.write_all(result.as_bytes()).unwrap();

        log::info!("Successful output summary of read mapping `{}`", output);

    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterReason { 
    Pass, 
    Singleton, 
    LowMQ, 
    Complex 
}

impl FilterReason {
    fn as_str(&self) -> &'static str {
        match self {
            FilterReason::Pass      => PASS_STR,
            FilterReason::Singleton => SINGLETON_STR,
            FilterReason::LowMQ     => LOW_MQ_STR,
            FilterReason::Complex   => COMPLEX_STR,
        }
    }
    fn from_str(s: &str) -> Self {
        match s {
            PASS_STR      => FilterReason::Pass,
            SINGLETON_STR => FilterReason::Singleton,
            LOW_MQ_STR    => FilterReason::LowMQ,
            COMPLEX_STR   => FilterReason::Complex,
            _ => panic!("Unknown filter reason"),
        }
    }
}


// #[derive(Debug)]
// pub struct Concatemer {
//     pub records: Vec<PoreCRecord>,
//     pub filter_reasons: Vec<String>,
// }

// impl Concatemer {

//     pub fn new() -> Concatemer{
//         Concatemer {
//             records: Vec::new(),
//             filter_reasons: Vec::new()
//         }
//     }

//     pub fn with_capacity(capacity: usize) -> Concatemer {
//         Concatemer {
//             records: Vec::with_capacity(capacity),
//             filter_reasons: Vec::with_capacity(capacity),
//         }
//     }

//     pub fn push(&mut self, pcr: PoreCRecord) {
//         self.filter_reasons.push(pcr.filter_reason.clone());
//         self.records.push(pcr);
        
//     }

//     pub fn clear(&mut self) {
//         self.records.clear();
//         self.filter_reasons.clear();
//     }

//     pub fn is_singleton(&self) -> bool {
//         self.filter_reasons.iter().filter(|&x| *x == "pass").count() == 1
//     }

//     pub fn parse_singleton(&mut self) {
//         if let Some(i) = self.records.iter().position(|r| r.filter_reason == "pass") {
//             self.records[i].filter_reason = String::from("singleton");
//             self.filter_reasons[i] = String::from("singleton");
//         }
//     }

//     pub fn is_complex(&self, max_order: &u32 ) -> bool {
//         let pass_count = self.filter_reasons.iter().filter(|&x| *x == "pass").count();
//         if pass_count > *max_order as usize {
//             true
//         } else {
//             false
//         }
//     }

//     pub fn parse_complex(&mut self) {
//         for (i, record) in self.records.iter_mut().enumerate() {
//             if record.filter_reason == "pass" {
//                 record.filter_reason = String::from("complex");
//                 self.filter_reasons[i] = String::from("complex");
//             }
//         }
//     }

//     pub fn filter_digest(&mut self, interval_hash: &HashMap<String, Lapper<usize, u8>> ) {
//         for (i, record) in self.records.iter_mut().enumerate() {
//             let is_in_regions: bool =  record.is_in_regions(interval_hash);
//             if is_in_regions {
//                 if record.filter_reason == "pass" {
//                     record.filter_reason = String::from("pass");
//                     self.filter_reasons[i] = String::from("pass");
//                 } else {
//                     record.filter_reason = String::from("low_mq");
//                     self.filter_reasons[i] = String::from("low_mq");
//                 }
              
//             } else {
//                 if record.filter_reason == "pass" {
//                     record.filter_reason = String::from("low_mq");
//                     self.filter_reasons[i] = String::from("low_mq");
//                 }
//             }
//         }
//     }

//     pub fn filter_edges(&mut self, max_length: u64) {
//         for (i, record) in self.records.iter_mut().enumerate() {
//             if record.target_start < max_length || (record.target_length - record.target_end) < max_length {
//                 record.filter_reason = String::from("low_mq");
//                 self.filter_reasons[i] = String::from("low_mq");
//             }
//         }
//     }

//     pub fn stat(&self) -> HashMap<&str, u32> {
//         let mut pass_count: u32 = 0;
//         let mut singleton_count: u32 = 0;
//         let mut low_mq_count: u32 = 0;
//         let mut complex_count: u32 = 0;
//         let mut stat_info: HashMap<&str, u32> = HashMap::new();


//         for filter_reason in self.filter_reasons.iter() {
//             if filter_reason == &String::from("pass") {
//                 pass_count += 1;
//             } else if filter_reason == &String::from("singleton") {
//                 singleton_count +=1;
//             } else if filter_reason == &String::from("low_mq") {
//                 low_mq_count += 1;
//             } else if filter_reason == &String::from("complex") {
//                 complex_count += 1;
//             } else {
//                 panic!("Unknown filter reason.");
//             }
            
//             stat_info.insert(&PASS_STR, pass_count);
//             stat_info.insert(&SINGLETON_STR, singleton_count);
//             stat_info.insert(&LOW_MQ_STR, low_mq_count);
//             stat_info.insert(&COMPLEX_STR, complex_count);

//         }
//         stat_info
//     }

//     pub fn info(&self) -> &str {
//         let stat_info = self.stat();

//         if *stat_info.get("pass").unwrap_or(&0) >= 2 {
//             return "pass";
//         } else if *stat_info.get("low_mq").unwrap_or(&0) >=1 && self.records.len() > 1  {
//             return "low_mq";
//         } else if *stat_info.get("complex").unwrap_or(&0) >= 1 && self.records.len() > 1 {
//             return "complex";
//         } else {
//            return "singleton";
//         }
//     }
// }

fn fast_parse_int(b: &[u8]) -> Option<u64> {
    let mut n = 0;
    for &c in b {
         if c < b'0' || c > b'9' { return None; }
         n = n * 10 + (c - b'0') as u64;
    }
    Some(n)
}

#[derive(Debug)]
pub struct Concatemer {
    pub records:        Vec<PoreCRecord>,
    pub filter_reasons: Vec<FilterReason>,
}

impl Concatemer {
    pub fn with_capacity(cap: usize) -> Concatemer {
        Concatemer {
            records: Vec::with_capacity(cap),
            filter_reasons: Vec::with_capacity(cap),
        }
    }

    pub fn push(&mut self, mut pcr: PoreCRecord) {
        let fr = FilterReason::from_str(&pcr.filter_reason);
        self.filter_reasons.push(fr);
        pcr.filter_reason = fr.as_str().to_string();
        self.records.push(pcr);
    }

    pub fn clear(&mut self) {
        self.records.clear();
        self.filter_reasons.clear();
    }

    pub fn is_singleton(&self) -> bool {
        self.filter_reasons.iter().filter(|&&r| r == FilterReason::Pass).count() == 1
    }

    pub fn parse_singleton(&mut self) {
        if let Some(i) = self.filter_reasons.iter().position(|&r| r == FilterReason::Pass) {
            self.filter_reasons[i] = FilterReason::Singleton;
            self.records[i].filter_reason = SINGLETON_STR.to_string();
        }
    }

    pub fn is_complex(&self, max_order: u32) -> bool {
        let cnt = self.filter_reasons.iter().filter(|&&r| r == FilterReason::Pass).count();
        cnt > (max_order as usize)
    }

    pub fn parse_complex(&mut self) {
        for i in 0..self.filter_reasons.len() {
            if self.filter_reasons[i] == FilterReason::Pass {
                self.filter_reasons[i] = FilterReason::Complex;
                self.records[i].filter_reason = COMPLEX_STR.to_string();
            }
        }
    }

    pub fn filter_digest(&mut self, interval_hash: &HashMap<String, Lapper<usize, u8>>) {
        for (i, r) in self.records.iter_mut().enumerate() {
            let in_region = r.is_in_regions(interval_hash);
            let new_fr = if in_region && self.filter_reasons[i] == FilterReason::Pass {
                FilterReason::Pass
            } else {
                FilterReason::LowMQ
            };
            self.filter_reasons[i] = new_fr;
            r.filter_reason = new_fr.as_str().to_string();
        }
    }

    pub fn filter_edges(&mut self, max_len: u64) {
        for (i, r) in self.records.iter_mut().enumerate() {
            if r.target_start < max_len || (r.target_length - r.target_end) < max_len {
                self.filter_reasons[i] = FilterReason::LowMQ;
                r.filter_reason = LOW_MQ_STR.to_string();
            }
        }
    }

    pub fn stat(&self) -> HashMap<&str, u32> {
        let mut counts = [0u32; 4];
        for &fr in &self.filter_reasons {
            match fr {
                FilterReason::Pass      => counts[0] += 1,
                FilterReason::Singleton => counts[1] += 1,
                FilterReason::LowMQ     => counts[2] += 1,
                FilterReason::Complex   => counts[3] += 1,
            }
        }
        let mut m = HashMap::new();
        m.insert(PASS_STR,      counts[0]);
        m.insert(SINGLETON_STR, counts[1]);
        m.insert(LOW_MQ_STR,    counts[2]);
        m.insert(COMPLEX_STR,   counts[3]);
        m
    }

    pub fn process(
        &mut self,
        is_filter_digest: bool,
        interval_hash: &HashMap<String, Lapper<usize, u8>>,
        is_filter_edges: bool,
        max_edge_length: u64,
        max_order: u32,
    ) -> (&'static str, Vec<usize>)  {
        let mut pass_indices = Vec::new();
        for (i, record) in self.records.iter_mut().enumerate() {
            if self.filter_reasons[i] != FilterReason::Pass {
                continue;
            }
            if is_filter_digest && !record.is_in_regions(interval_hash) {
                self.filter_reasons[i] = FilterReason::LowMQ;
                continue;
            }
            if is_filter_edges && (record.target_start < max_edge_length || (record.target_length - record.target_end) < max_edge_length) {
                self.filter_reasons[i] = FilterReason::LowMQ;
                continue;
            }

            pass_indices.push(i);
        }

        let pass_count = pass_indices.len();
        let final_status = if pass_count == 0 {
            if self.records.len() > 1 { LOW_MQ_STR } else { SINGLETON_STR }
        } else if pass_count == 1 {
            self.filter_reasons[pass_indices[0]] = FilterReason::Singleton;
            SINGLETON_STR
        } else if pass_count > max_order as usize {
            for &idx in &pass_indices {
                self.filter_reasons[idx] = FilterReason::Complex;
            }
            COMPLEX_STR
        } else {
            PASS_STR
        };


        for i in 0..self.records.len() {
            // self.records[i].filter_reason = self.filter_reasons[i].as_str().to_string();
            let s = self.filter_reasons[i].as_str();
            if self.records[i].filter_reason.as_str() != s {
                self.records[i].filter_reason.clear();
                self.records[i].filter_reason.push_str(s);
            }
        }

        (final_status, pass_indices)
    }


    pub fn info(&self) -> &str {
        let s = self.stat();
        if s[PASS_STR] >= 2 {
            PASS_STR
        } else if s[LOW_MQ_STR] >= 1 && self.records.len() > 1 {
            LOW_MQ_STR
        } else if s[COMPLEX_STR] >= 1 && self.records.len() > 1 {
            COMPLEX_STR
        } else {
            SINGLETON_STR
        }
    }
}

fn process_group(
    grp: &mut Vec<PoreCRecord>, 
    ih: &HashMap<String, Lapper<usize, u8>>, 
    is_digest: bool, is_edges: bool, max_len: u64, max_ord: u32,
    buf: &mut String, p: &mut u64, s: &mut u64, l: &mut u64, c: &mut u64,
) {
    let mut concatemer = Concatemer::with_capacity(grp.len() + 2);
    for r in grp.drain(..) { concatemer.push(r); }
    
    let (final_status, pass_indices) = concatemer.process(
        is_digest, ih, is_edges, max_len, max_ord
    );

    if final_status == PASS_STR {
         for &idx in &pass_indices {
             buf.push_str(&concatemer.records[idx].to_string());
             buf.push('\n');
         }
    }
    match final_status {
        PASS_STR => *p += 1,
        SINGLETON_STR => *s += 1,
        LOW_MQ_STR => *l += 1,
        COMPLEX_STR => *c += 1,
        _ => {}
    }
}

impl PAFTable {
    pub fn parse(&self) -> anyResult<csv::Reader<Box<dyn BufRead + Send>>> {
        let input = common_reader(&self.file);
        let rdr = csv::ReaderBuilder::new()
                            .flexible(true)
                            .has_headers(false)
                            .comment(Some(b'#'))
                            .delimiter(b'\t')
                            .from_reader(input);
        
        Ok(rdr)
    }

    pub fn parse2(&self, secondary: bool) -> Result<impl Iterator<Item = PAFLine>, Box<dyn Error>> {
        let  input = common_reader(&self.file);
        let iter = input.lines()
                    .filter_map(move |l| {
                        let line = l.ok().unwrap();
                        if line.is_empty() || line.starts_with('#') {
                            return None;
                        }

                        let mut cols = line.splitn(13,'\t');
                        let mut it = line.splitn(13, '\t');
                        let q          = it.next()?;
                        let qlen_s     = it.next()?;
                        let qs_s       = it.next()?;
                        let qe_s       = it.next()?;
                        let qstrand_s  = it.next()?;
                        let t          = it.next()?;
                        let tlen_s     = it.next()?;
                        let ts_s       = it.next()?;
                        let te_s       = it.next()?;
                        let match_s    = it.next()?;
                        let alen_s     = it.next()?;
                        let mapq_s     = it.next()?;
                        let optional   = it.next().unwrap_or("");

                        // Some(PAFLine{
                        //     query:          cols.next().unwrap().to_owned(),
                        //     query_length:   cols.next().unwrap().parse().ok().unwrap(),
                        //     query_start:    cols.next().unwrap().parse().ok().unwrap(),
                        //     query_end:      cols.next().unwrap().parse().ok().unwrap(),
                        //     query_strand:   cols.next().unwrap().chars().next().unwrap(),
                        //     target:         cols.next().unwrap().to_owned(),
                        //     target_length:  cols.next().unwrap().parse().ok().unwrap(),
                        //     target_start:   cols.next().unwrap().parse().ok().unwrap(),
                        //     target_end:     cols.next().unwrap().parse().ok().unwrap(),
                        //     match_n:        cols.next().unwrap().parse().ok().unwrap(),
                        //     alignment_length: cols.next().unwrap().parse().ok().unwrap(),
                        //     mapq:           cols.next().unwrap().parse().ok().unwrap(),
                        // })
                        if !secondary {
                            if optional.contains("tp:A:S") {
                                return None;
                            }
                        }
                        
                        let query_length  = qlen_s.parse().ok()?;
                        let query_start   = qs_s.parse().ok()?;
                        let query_end     = qe_s.parse().ok()?;
                        let query_strand  = qstrand_s.chars().next()?;
                        let target_length = tlen_s.parse().ok()?;
                        let target_start  = ts_s.parse().ok()?;
                        let target_end    = te_s.parse().ok()?;
                        let match_n       = match_s.parse().ok()?;
                        let alignment_len = alen_s.parse().ok()?;
                        let mapq          = mapq_s.parse().ok()?;

                        Some(PAFLine {
                            query: q.to_owned(),
                            query_length,
                            query_start,
                            query_end,
                            query_strand,
                            target: t.to_owned(),
                            target_length,
                            target_start,
                            target_end,
                            match_n,
                            alignment_length: alignment_len,
                            mapq,
                        })

                    });

        Ok(iter)
    }

    pub fn paf2table_v0_3_0(&self, bed: &String,
            output: &String, min_quality: &u8, 
            min_identity: &f32,  min_length: &u32,
            max_order: &u32, 
            max_edge_length: &u64,
            secondary: bool,
            num_threads: usize 
        ) -> Result<(), Box<dyn Error>> {
            
        let min_quality = *min_quality;
        let min_identity = *min_identity;
        let min_length = *min_length;
        let max_order = *max_order;
        let max_edge_length = *max_edge_length;

        let is_filter_digest = !bed.is_empty();
        let bed_obj = if is_filter_digest { Bed3::new(bed) } else { Bed3::new(&String::from(".tmp.bed")) };
        let interval_hash = Arc::new(bed_obj.to_interval_hash());
        let is_filter_edges = max_edge_length > 0;

        log::info!("Using {} threads for processing.", num_threads);
        
        let (tx_raw, rx_raw) = bounded::<Vec<u8>>(16); 
        let (tx_writer, rx_writer) = bounded::<String>(1024);

        let output_clone = output.clone();
        let writer_handle = thread::spawn(move || {
            let mut writer = common_writer(&output_clone);
            for chunk in rx_writer {
                writer.write_all(chunk.as_bytes()).unwrap();
            }
        });

        let mut handles = Vec::new();
        use std::sync::atomic::{AtomicU64, Ordering};
        let stats_pass = Arc::new(AtomicU64::new(0));
        let stats_sing = Arc::new(AtomicU64::new(0));
        let stats_low  = Arc::new(AtomicU64::new(0));
        let stats_comp = Arc::new(AtomicU64::new(0));
        
        let global_read_idx = Arc::new(AtomicU64::new(0)); 

        for _ in 0..num_threads {
            let rx = rx_raw.clone();
            let tx_w = tx_writer.clone();
            let ih = interval_hash.clone();
            let s_p = stats_pass.clone();
            let s_s = stats_sing.clone();
            let s_l = stats_low.clone();
            let s_c = stats_comp.clone();
            let g_rid = global_read_idx.clone();

            handles.push(thread::spawn(move || {
                let mut local_p=0; let mut local_s=0; let mut local_l=0; let mut local_c=0;
                let mut output_buffer = String::with_capacity(WRITE_CHUNK);
                
    
                fn process_group_vec(
                    groups: &mut Vec<Vec<PoreCRecord>>,
                    start_id: u64,
                    ih: &HashMap<String, Lapper<usize, u8>>,
                    is_fd: bool, is_fe: bool, max_len: u64, max_ord: u32,
                    buf: &mut String,
                    lp: &mut u64, ls: &mut u64, ll: &mut u64, lc: &mut u64
                ) {
                    for (i, grp) in groups.iter_mut().enumerate() {
                        let rid = start_id + i as u64;
                        for r in grp.iter_mut() { r.read_idx = rid; }
                        
                        process_single_group(grp, ih, is_fd, is_fe, max_len, max_ord, buf, lp, ls, ll, lc);
                    }
                    groups.clear();
                }

                fn process_single_group(
                    grp: &mut Vec<PoreCRecord>, 
                    ih: &HashMap<String, Lapper<usize, u8>>, 
                    is_digest: bool, is_edges: bool, max_len: u64, max_ord: u32,
                    buf: &mut String, p: &mut u64, s: &mut u64, l: &mut u64, c: &mut u64,
                ) {
                    let mut concatemer = Concatemer::with_capacity(grp.len() + 2);
                    for r in grp.drain(..) { concatemer.push(r); }
                    
                    let (final_status, pass_indices) = concatemer.process(
                        is_digest, ih, is_edges, max_len, max_ord
                    );

                    if final_status == PASS_STR {
                            for &idx in &pass_indices {
                                buf.push_str(&concatemer.records[idx].to_string());
                                buf.push('\n');
                            }
                    }
                    match final_status {
                        PASS_STR => *p += 1,
                        SINGLETON_STR => *s += 1,
                        LOW_MQ_STR => *l += 1,
                        COMPLEX_STR => *c += 1,
                        _ => {}
                    }
                }

                let mut block_groups: Vec<Vec<PoreCRecord>> = Vec::with_capacity(2048);

                while let Ok(block) = rx.recv() {
                    let mut start = 0;
                    let len = block.len();
                    
                    let mut current_group: Vec<PoreCRecord> = Vec::with_capacity(32);
                    let mut old_q_start = 0; 
                    let mut old_q_end = 0;
                    let mut has_old = false;

                    while start < len {
                        let end = match block[start..].iter().position(|&b| b == b'\n') {
                            Some(p) => start + p,
                            None => len,
                        };
                        let mut line = &block[start..end];
                        start = end + 1; 
                        
                        if !line.is_empty() && line[line.len()-1] == b'\r' { line = &line[..line.len()-1]; }
                        if line.is_empty() || line[0] == b'#' { continue; }

                        let mut iter = line.split(|&b| b == b'\t');
                        let q_bytes = match iter.next() { Some(v) => v, None => continue };

                        let is_switch = if !has_old {
                            false
                        } else {
                            q_bytes != &block[old_q_start..old_q_end]
                        };

                        if is_switch {
                            if !current_group.is_empty() {
                                block_groups.push(std::mem::take(&mut current_group));
                            
                            }
                        }

                        if !has_old || is_switch {
                            let q_len = q_bytes.len();
                            let q_ptr = q_bytes.as_ptr() as usize;
                            let start_ptr = block.as_ptr() as usize;
                            old_q_start = q_ptr - start_ptr;
                            old_q_end = old_q_start + q_len;
                            has_old = true;
                        }

                            let qlen = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u32, None=>continue };
                            let qs   = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u32, None=>continue };
                            let qe   = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u32, None=>continue };
                            let strand_char = match iter.next() { Some(v)=> if v.is_empty() { '+' } else { v[0] as char }, None=>'+' };
                            let t_bytes = match iter.next() { Some(v)=>v, None=>continue };
                            let tlen = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u64, None=>continue };
                            let ts   = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u64, None=>continue };
                            let te   = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u64, None=>continue };
                            let matches = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u32, None=>continue };
                            let alen    = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u32, None=>continue };
                            let mapq    = match iter.next().and_then(fast_parse_int) { Some(v)=>v as u8, None=>continue };

                            if mapq < min_quality { continue; }
                            if !secondary {
                                let mut is_sec = false;
                                for tag in iter { if tag.starts_with(b"tp:A:S") { is_sec = true; break; } }
                                if is_sec { continue; }
                            }

                            let q_cov = qe - qs;
                            let ident = if alen > 0 { matches as f32 / alen as f32 } else { 0.0 };
                            let fr = if q_cov < min_length || ident < min_identity { LOW_MQ_STR } else { PASS_STR };

                            let t_str = unsafe { String::from_utf8_unchecked(t_bytes.to_vec()) };

                        current_group.push(PoreCRecord {
                                read_idx: 0, // Placeholder
                                query_length: qlen, query_start: qs, query_end: qe, query_strand: strand_char,
                                target: t_str, target_length: tlen, target_start: ts, target_end: te,
                                mapq: mapq, identity: ident, filter_reason: fr.to_string(),
                        });
                    }
                    
                    if !current_group.is_empty() {
                            block_groups.push(current_group);
                    }
            
                    let count = block_groups.len() as u64;
                    if count > 0 {
                        let start_id = g_rid.fetch_add(count, Ordering::Relaxed);
                        
        
                        process_group_vec(&mut block_groups, 
                            start_id, &ih, is_filter_digest, 
                            is_filter_edges, 
                            max_edge_length, 
                            max_order,  
                            &mut output_buffer, 
                            &mut local_p, &mut local_s, &mut local_l, &mut local_c);
                    }
                    
                    if output_buffer.len() >= WRITE_CHUNK {
                        tx_w.send(std::mem::take(&mut output_buffer)).unwrap();
                        output_buffer = String::with_capacity(WRITE_CHUNK); 
                    }
                }
                if !output_buffer.is_empty() {
                    tx_w.send(output_buffer).unwrap();
                }
                
                s_p.fetch_add(local_p, Ordering::Relaxed);
                s_s.fetch_add(local_s, Ordering::Relaxed);
                s_l.fetch_add(local_l, Ordering::Relaxed);
                s_c.fetch_add(local_c, Ordering::Relaxed);
            }));
        }
        drop(tx_writer);

        {
            let mut reader = common_reader(&self.file);
            const BLOCK_SIZE: usize = 64 * 1024 * 1024; 
            let mut buffer = vec![0u8; BLOCK_SIZE]; 
            let mut valid_len = 0;

            loop {
                let bytes_read = reader.read(&mut buffer[valid_len..]).unwrap_or(0);
                if bytes_read == 0 {
                    if valid_len > 0 { tx_raw.send(buffer[..valid_len].to_vec()).unwrap(); }
                    break;
                }
                let end_ptr = valid_len + bytes_read;
                let active_slice = &buffer[..end_ptr];
                
                let last_nl = match active_slice.iter().rposition(|&b| b == b'\n') {
                    Some(p) => p,
                    None => {
                            valid_len = end_ptr;
                            if valid_len == buffer.len() { buffer.resize(buffer.len()*2, 0); }
                            continue;
                    }
                };

                let safe_slice = &active_slice[..=last_nl];
                
                let prev_nl = safe_slice[..last_nl].iter().rposition(|&b| b == b'\n').map(|p| p+1).unwrap_or(0);
                let last_line = &safe_slice[prev_nl..last_nl];
                let tab_pos = last_line.iter().position(|&b| b == b'\t').unwrap_or(last_line.len());
                let last_q_bytes = &last_line[..tab_pos];
                
                let mut cut_pos = safe_slice.len(); 
                let mut found_boundary = false;

                let mut scan_end = prev_nl;
                let scan_limit = if safe_slice.len() > 1024*1024 { safe_slice.len() - 1024*1024 } else { 0 };

                while scan_end > scan_limit {
                    if scan_end == 0 { break; } 
                    let newline_idx = scan_end - 1;
                    let start_idx = safe_slice[..newline_idx].iter().rposition(|&b| b == b'\n').map(|p| p+1).unwrap_or(0);
                    let line = &safe_slice[start_idx..newline_idx];
                    let t = line.iter().position(|&b| b == b'\t').unwrap_or(line.len());
                    let q = &line[..t];
                    if q != last_q_bytes {
                        cut_pos = scan_end;
                        found_boundary = true;
                        break;
                    }
                    if start_idx == 0 { break; }
                    scan_end = start_idx;
                }

                if !found_boundary && end_ptr < buffer.len() {
                        if bytes_read < (buffer.len() - valid_len) {
                            tx_raw.send(safe_slice.to_vec()).unwrap();
                            let rem_len = end_ptr - (last_nl + 1);
                            buffer.copy_within((last_nl+1)..end_ptr, 0);
                            valid_len = rem_len;
                            continue;
                        } 
                        valid_len = end_ptr;
                        buffer.resize(buffer.len() * 2, 0);
                        continue;
                }
                
                let chunk_to_send = &buffer[..cut_pos];
                tx_raw.send(chunk_to_send.to_vec()).unwrap(); 

                let tail_len = end_ptr - cut_pos;
                buffer.copy_within(cut_pos..end_ptr, 0);
                valid_len = tail_len;
            }
        }
        drop(tx_raw);

        for h in handles { h.join().unwrap(); }
        writer_handle.join().unwrap();
        
        let mut summary = ReadSummary::new();
        summary.pass = stats_pass.load(Ordering::Relaxed);
        summary.singleton = stats_sing.load(Ordering::Relaxed);
        summary.low_mq = stats_low.load(Ordering::Relaxed);
        summary.complex = stats_comp.load(Ordering::Relaxed);
        summary.mapping = global_read_idx.load(Ordering::Relaxed); 

        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(output).with_extension("").to_str().unwrap().to_string()
        };
        summary.save(&format!("{}.read.summary", output_prefix));

        log::info!("Successful output Pore-C table `{}` (block-io-atomic-ids)", output);
        Ok(())
    }

    pub fn paf2table(&self, bed: &String,
            output: &String, min_quality: &u8, 
            min_identity: &f32,  min_length: &u32,
            max_order: &u32, 
            max_edge_length: &u64,
            secondary: bool,
            threads: usize,
        ) -> Result<(), Box<dyn Error>> {
            
        let min_quality = *min_quality;
        let min_identity = *min_identity;
        let min_length = *min_length;
        let max_order = *max_order;
        let max_edge_length = *max_edge_length;

        let is_filter_digest = !bed.is_empty();
        let bed_obj = if is_filter_digest { Bed3::new(bed) } else { Bed3::new(&String::from(".tmp.bed")) };
        let interval_hash = Arc::new(bed_obj.to_interval_hash());
        let is_filter_edges = max_edge_length > 0;

        let num_threads = threads;
        log::info!("Using {} threads for processing.", num_threads);
        
        let (tx_raw, rx_raw) = bounded::<Vec<u8>>(200); 
        let (tx_writer, rx_writer) = bounded::<String>(2000);

        let output_clone = output.clone();
        let writer_handle = thread::spawn(move || {
            let mut writer = common_writer(&output_clone);
            for chunk in rx_writer {
                writer.write_all(chunk.as_bytes()).unwrap();
            }
        });

        let mut handles = Vec::new();
        use std::sync::atomic::{AtomicU64, Ordering};
        let stats_pass = Arc::new(AtomicU64::new(0));
        let stats_sing = Arc::new(AtomicU64::new(0));
        let stats_low  = Arc::new(AtomicU64::new(0));
        let stats_comp = Arc::new(AtomicU64::new(0));
        let global_read_idx = Arc::new(AtomicU64::new(1)); 

        for _ in 0..num_threads {
            let rx = rx_raw.clone();
            let tx_w = tx_writer.clone();
            let ih = interval_hash.clone();
            let s_p = stats_pass.clone();
            let s_s = stats_sing.clone();
            let s_l = stats_low.clone();
            let s_c = stats_comp.clone();
            let g_rid = global_read_idx.clone();

            handles.push(thread::spawn(move || {
                let mut local_p=0; let mut local_s=0; let mut local_l=0; let mut local_c=0;
                let mut output_buffer = String::with_capacity(WRITE_CHUNK);
                
                fn parse_line_only(line: &[u8], min_q: u8, min_l: u32, min_id: f32, sec: bool) -> Option<PoreCRecord> {
                    let mut iter = line.split(|&b| b == b'\t');
                    let _q_bytes = iter.next()?;
                    
                    let qlen: u32 = fast_parse_int(iter.next()?)? as u32;
                    let qs: u32   = fast_parse_int(iter.next()?)? as u32;
                    let qe: u32   = fast_parse_int(iter.next()?)? as u32;
                    let strand_c  = iter.next().map(|s| if s.is_empty() { '+' } else { s[0] as char }).unwrap_or('+');
                    let t_bytes   = iter.next()?;
                    let tlen: u64 = fast_parse_int(iter.next()?)?;
                    let ts: u64   = fast_parse_int(iter.next()?)?;
                    let te: u64   = fast_parse_int(iter.next()?)?;
                    let mat: u32  = fast_parse_int(iter.next()?)? as u32;
                    let aln: u32  = fast_parse_int(iter.next()?)? as u32;
                    let mapq: u8  = fast_parse_int(iter.next()?)? as u8;
                    
                    if mapq < min_q { return None; }

                    if !sec {
                        let mut is_sec = false;
                        for tag in iter {
                            if tag.starts_with(b"tp:A:S") { is_sec = true; break; }
                        }
                        if is_sec { return None; }
                    }

                    let q_cov = qe - qs;
                    let ident = if aln > 0 { mat as f32 / aln as f32 } else { 0.0 };
                    let fr = if q_cov < min_l || ident < min_id {
                        LOW_MQ_STR
                    } else {
                        PASS_STR
                    };

                    let t_str = unsafe { String::from_utf8_unchecked(t_bytes.to_vec()) };

                    Some(PoreCRecord {
                        read_idx: 0, 
                        query_length: qlen, 
                        query_start: qs, query_end: qe, query_strand: strand_c,
                        target: t_str, target_length: tlen,
                        target_start: ts, target_end: te,
                        mapq: mapq, identity: ident, filter_reason: fr.to_string(),
                    })
                }

                while let Ok(batch_bytes) = rx.recv() {
                    let mut start = 0;
                    let mut current_group: Vec<PoreCRecord> = Vec::with_capacity(32);
                    // 优化 1：彻底消灭 String，改用零分配的 &[u8] 切片判定 Read 切换
                    let mut old_q_bytes: &[u8] = &[];
                    
                    while start < batch_bytes.len() {
                        let end = match batch_bytes[start..].iter().position(|&b| b == b'\n') {
                            Some(p) => start + p,
                            None => batch_bytes.len(),
                        };
                        let line = &batch_bytes[start..end];
                        start = end + 1; 
                        
                        if line.is_empty() || line[0] == b'#' { continue; }
                        
                        let tab_pos = match line.iter().position(|&b| b == b'\t') {
                            Some(p) => p, None => continue,
                        };
                        let q_bytes = &line[..tab_pos];
                        
                        let is_switch = !old_q_bytes.is_empty() && old_q_bytes != q_bytes;
                        if is_switch {
                            if !current_group.is_empty() {
                                // 优化 2：多核并行无锁 fetch_add 抢占分配全局唯一 ID，完全免除主线程调度负担
                                let rid = g_rid.fetch_add(1, Ordering::Relaxed);
                                for r in &mut current_group { r.read_idx = rid; }
                                process_group(&mut current_group, &ih, is_filter_digest, is_filter_edges, max_edge_length, max_order, &mut output_buffer, &mut local_p, &mut local_s, &mut local_l, &mut local_c);
                                
                                if output_buffer.len() >= WRITE_CHUNK {
                                    tx_w.send(std::mem::take(&mut output_buffer)).unwrap();
                                    output_buffer = String::with_capacity(WRITE_CHUNK);
                                }
                            }
                            old_q_bytes = q_bytes;
                        } else if old_q_bytes.is_empty() {
                            old_q_bytes = q_bytes;
                        }

                        if let Some(pcr) = parse_line_only(line, min_quality, min_length, min_identity, secondary) {
                            current_group.push(pcr);
                        }
                    }
                    
                    if !current_group.is_empty() {
                        let rid = g_rid.fetch_add(1, Ordering::Relaxed);
                        for r in &mut current_group { r.read_idx = rid; }
                        process_group(&mut current_group, &ih, is_filter_digest, is_filter_edges, max_edge_length, max_order, &mut output_buffer, &mut local_p, &mut local_s, &mut local_l, &mut local_c);
                    }
                    
                    if output_buffer.len() >= WRITE_CHUNK {
                        tx_w.send(std::mem::take(&mut output_buffer)).unwrap();
                        output_buffer = String::with_capacity(WRITE_CHUNK);
                    }
                }
                
                if !output_buffer.is_empty() {
                    tx_w.send(output_buffer).unwrap();
                }
                s_p.fetch_add(local_p, Ordering::Relaxed);
                s_s.fetch_add(local_s, Ordering::Relaxed);
                s_l.fetch_add(local_l, Ordering::Relaxed);
                s_c.fetch_add(local_c, Ordering::Relaxed);
            }));
        }
        drop(tx_writer);

        let mut reader = common_reader(&self.file);
        const BATCH_SIZE_LIMIT: usize = 16 * 1024 * 1024; 
        let mut buffer = vec![0u8; BATCH_SIZE_LIMIT * 2];
        let mut bytes_in_buffer = 0;

        #[inline(always)]
        fn get_qname(line: &[u8]) -> &[u8] {
            if let Some(pos) = line.iter().position(|&b| b == b'\t') {
                &line[..pos]
            } else {
                line
            }
        }

        loop {
            let bytes_read = reader.read(&mut buffer[bytes_in_buffer..]).unwrap_or(0);
            if bytes_read == 0 && bytes_in_buffer == 0 {
                break;
            }
            bytes_in_buffer += bytes_read;

            let mut chunk_end = bytes_in_buffer;
            if bytes_in_buffer > BATCH_SIZE_LIMIT {
                if let Some(pos) = buffer[BATCH_SIZE_LIMIT..bytes_in_buffer].iter().position(|&b| b == b'\n') {
                    let nl_idx = BATCH_SIZE_LIMIT + pos;
                    let line_start = match buffer[..nl_idx].iter().rposition(|&b| b == b'\n') {
                        Some(p) => p + 1,
                        None => 0,
                    };
                    let q_last = get_qname(&buffer[line_start..nl_idx]);

                    let mut scan_pos = nl_idx + 1;
                    let mut found = false;
                    while scan_pos < bytes_in_buffer {
                        let next_nl = match buffer[scan_pos..bytes_in_buffer].iter().position(|&b| b == b'\n') {
                            Some(p) => scan_pos + p,
                            None => bytes_in_buffer,
                        };
                        let line = &buffer[scan_pos..next_nl];
                        if !line.is_empty() && line[0] != b'#' {
                            let q_cur = get_qname(line);
                            if q_cur != q_last {
                                chunk_end = scan_pos;
                                found = true;
                                break;
                            }
                        }
                        scan_pos = next_nl + 1;
                    }

                    if !found && bytes_in_buffer < buffer.len() {
                        continue;
                    }
                }
            }

            let chunk = &buffer[..chunk_end];
            
            // 优化 3：主线程除极其轻量级的边界定位外做 0 处理，不切行，不计数，秒级直发给下游通道
            tx_raw.send(chunk.to_vec()).unwrap();

            if chunk_end < bytes_in_buffer {
                buffer.copy_within(chunk_end..bytes_in_buffer, 0);
                bytes_in_buffer -= chunk_end;
            } else {
                bytes_in_buffer = 0;
            }
        }
        drop(tx_raw);

        for h in handles { h.join().unwrap(); }
        writer_handle.join().unwrap();
        
        let final_stats = (
            stats_pass.load(Ordering::Relaxed), 
            stats_sing.load(Ordering::Relaxed),
            stats_low.load(Ordering::Relaxed),
            stats_comp.load(Ordering::Relaxed)
        );
        
        let mut summary = ReadSummary::new();
        summary.pass = final_stats.0;
        summary.singleton = final_stats.1;
        summary.low_mq = final_stats.2;
        summary.complex = final_stats.3;

        summary.mapping = global_read_idx.load(Ordering::Relaxed).saturating_sub(1); 

        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(output).with_extension("").to_str().unwrap().to_string()
        };
        summary.save(&format!("{}.read.summary", output_prefix));

        log::info!("Successful output Pore-C table `{}` (raw-bytes-batched)", output);
        Ok(())
    }
    pub fn to_depth(&self, chromsize: &String, window_size: usize, step_size: usize,
                min_mapq: u8, secondary: bool, output: &String) -> Result<(), Box<dyn Error>> {
        let mut chrom_names = Vec::new();
        let mut chrom_sizes = Vec::new();
        let mut chrom_map = HashMap::new();
        {
            let input = common_reader(chromsize);
            let buf = std::io::BufReader::new(input);
            for line in buf.lines() {
                let line = line?;
                if line.trim().is_empty() || line.starts_with('#') { continue; }
                let mut spl = line.split_whitespace();
                let chr = spl.next().unwrap_or("").to_string();
                let size = spl.next().unwrap_or("0").parse::<usize>().unwrap_or(0);
                chrom_map.insert(chr.clone(), chrom_names.len());
                chrom_names.push(chr);
                chrom_sizes.push(size);
            }
        }
        
        let num_chroms = chrom_names.len();
        let chrom_map = Arc::new(chrom_map);

        let (sender, receiver) = bounded::<Vec<u8>>(100);
        
        let mut handles = Vec::new();
        for _ in 0..8 {
            let receiver = receiver.clone();
            let chrom_map = chrom_map.clone();
            handles.push(thread::spawn(move || {
                let mut local_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];
    
                fn fast_parse_u64(bytes: &[u8]) -> Option<u64> {
                    let mut n = 0;
                    for &b in bytes {
                        if b < b'0' || b > b'9' { return None; }
                        n = n * 10 + (b - b'0') as u64;
                    }
                    Some(n)
                }

                while let Ok(data) = receiver.recv() {
                    let mut start_idx = 0;

                    while start_idx < data.len() {
                        let end_idx = match data[start_idx..].iter().position(|&b| b == b'\n') {
                            Some(pos) => start_idx + pos,
                            None => data.len(), 
                        };
                        let line = &data[start_idx..end_idx];
                        start_idx = end_idx + 1;

                        if line.is_empty() || line[0] == b'#' { continue; }


                        let mut iter = line.split(|&b| b == b'\t');
                        
                        if iter.nth(4).is_none() { continue; } 
                        let target_bytes = match iter.next() { Some(v) => v, None => continue };
                        
                        if iter.next().is_none() { continue; }
                        let t_start_bytes = match iter.next() { Some(v) => v, None => continue };
                        let t_end_bytes = match iter.next() { Some(v) => v, None => continue };
                    
                        if iter.nth(1).is_none() { continue; }
                        
                        let mapq_bytes = match iter.next() { Some(v) => v, None => continue };
                        let mapq = match fast_parse_u64(mapq_bytes) { Some(v) => v as u8, None => continue };

                        if mapq < min_mapq { continue; }
                        if !secondary {
                            let mut is_secondary = false;
                            for tag in iter {
                                if tag.starts_with(b"tp:A:S") {
                                    is_secondary = true;
                                    break;
                                }
                            }
                            if is_secondary { continue; }
                        }

                        let t_start = match fast_parse_u64(t_start_bytes) { Some(v) => v, None => continue };
                        let t_end = match fast_parse_u64(t_end_bytes) { Some(v) => v, None => continue };

                        if let Ok(target_str) = std::str::from_utf8(target_bytes) {
                                if let Some(&idx) = chrom_map.get(target_str) {
                                local_events[idx].push((t_start as u32, 1));
                                local_events[idx].push((t_end as u32, -1));
                                }
                        }
                    }
                }
                local_events
            }));
        }
        
        {
            let mut input = common_reader(&self.file); 
            
            const BUF_SIZE: usize = 512 * 1024; 
            let mut buffer = vec![0u8; BUF_SIZE];
            let mut leftover = Vec::with_capacity(1024);

            loop {
                let bytes_read = input.read(&mut buffer).unwrap_or(0);
                if bytes_read == 0 {
                    break;
                }

                let chunk = &buffer[..bytes_read];

                let last_newline = chunk.iter().rposition(|&b| b == b'\n');

                let to_send = if let Some(pos) = last_newline {
                    
                    let mut send_buf = Vec::with_capacity(leftover.len() + pos + 1);
                    send_buf.extend_from_slice(&leftover);
                    send_buf.extend_from_slice(&chunk[..=pos]);
                    
                    leftover.clear();
                    leftover.extend_from_slice(&chunk[pos+1..]);
                    
                    send_buf
                } else {
                    leftover.extend_from_slice(chunk);
                    continue; 
                };

                if !to_send.is_empty() {
                    sender.send(to_send).unwrap();
                }
            }
            if !leftover.is_empty() {
                sender.send(leftover).unwrap();
            }
        }
        
        drop(sender);

        let mut global_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];

        for handle in handles {
            let local_events = handle.join().unwrap();
            for (i, events) in local_events.into_iter().enumerate() {
                global_events[i].extend(events);
            }
        }

        rayon::scope(|s| {
            let (tx, rx) = bounded::<String>(1024);

            s.spawn(|_| {
                let mut wtr = common_writer(output);
                for s in rx {
                    wtr.write_all(s.as_bytes()).unwrap();
                }
            });

            global_events.into_par_iter().enumerate().for_each_with(tx, |tx, (i, mut events)| {
                if events.is_empty() { return; }
                let chrom_len = chrom_sizes[i];
                let chrom_name = &chrom_names[i];

                events.par_sort_unstable_by_key(|e| e.0);

                let mut prefix_sum = Vec::with_capacity(chrom_len + 1);
                prefix_sum.push(0);
                
                let mut current_cov = 0i32;
                let mut prev_pos = 0usize;
                let mut running_sum = 0i64;

                for (pos, delta) in events {
                    let pos = pos as usize;
                    if pos >= chrom_len { break; }
                    
                    if pos > prev_pos {
                        let cov = current_cov as i64;
                        let count = pos - prev_pos;
                        for _ in 0..count {
                            running_sum += cov;
                            prefix_sum.push(running_sum);
                        }
                    }
                    current_cov += delta;
                    prev_pos = pos;
                }
                
                if prev_pos < chrom_len {
                    let cov = current_cov as i64;
                    let count = chrom_len - prev_pos;
                    for _ in 0..count {
                        running_sum += cov;
                        prefix_sum.push(running_sum);
                    }
                }

                let mut output_buf = String::with_capacity(64 * 1024);
                let mut start = 0usize;

                while start < chrom_len {
                    let end = std::cmp::min(start + window_size, chrom_len);
                    let len = end - start;
                    if len == 0 { break; }
                    
                    let sum = prefix_sum[end] - prefix_sum[start];
                    let mean = sum as f64 / len as f64;
                    
                    write!(output_buf, "{}\t{}\t{}\t{:.3}\n", chrom_name, start, end, mean).unwrap();
                // let mut coverage = vec![0i32; chrom_len];
                // let mut current_cov = 0i32;
                // let mut prev_pos = 0usize;

                // for (pos, delta) in events {
                //     let pos = pos as usize;
                //     if pos >= chrom_len { break; }
                    
                //     if pos > prev_pos {
                //         for k in prev_pos..pos {
                //             coverage[k] = current_cov;
                //         }
                //     }
                //     current_cov += delta;
                //     prev_pos = pos;
                // }
                
                // if prev_pos < chrom_len {
                //     for k in prev_pos..chrom_len {
                //         coverage[k] = current_cov;
                //     }
                // }

                // let mut output_buf = String::with_capacity(64 * 1024);
                // let mut start = 0usize;
                
                // let mut window_buf = Vec::with_capacity(window_size);

                // while start < chrom_len {
                //     let end = std::cmp::min(start + window_size, chrom_len);
                //     let len = end - start;
                //     if len == 0 { break; }
                    
    
                //     window_buf.clear();
                //     window_buf.extend_from_slice(&coverage[start..end]);
                    
                //     let mid = len / 2;
                //     let (_, &mut median_val, _) = window_buf.select_nth_unstable(mid);
                    
                //     let median = if len % 2 == 0 {
                //         let max_lower = *window_buf[..mid].iter().max().unwrap();
                //         (median_val as f64 + max_lower as f64) / 2.0
                //     } else {
                //         median_val as f64
                //     };
                    
                    // write!(output_buf, "{}\t{}\t{}\t{:.3}\n", chrom_name, start, end, median).unwrap();
                    if output_buf.len() > 60 * 1024 {
                        tx.send(std::mem::take(&mut output_buf)).unwrap();
                        output_buf = String::with_capacity(64 * 1024);
                    }

                    start += step_size;
                }
                if !output_buf.is_empty() {
                    tx.send(output_buf).unwrap();
                }
            });
        });

        log::info!(
                "Successful output coverage (window {}bp, step {}bp) to `{}`",
                window_size,
                step_size,
                output
            );
        Ok(())
    }


    pub fn downsample_by_vpc_random(
        &self,
        output: &String,
        seed: u64,
        min_quality: u8,
        min_identity: f32,
        min_length: u32,
        min_order: usize,
        max_order: usize,
        target_reads: Option<u64>,
        target_pairs: Option<u64>,
        target_bases: Option<u64>,
        frac: Option<f64>,
        frac_by_pairs: bool,
        keep_comments: bool,
    ) -> anyResult<()> {

        if [
            target_reads.is_some(),
            target_pairs.is_some(),
            target_bases.is_some(),
            frac.is_some(),
        ]
        .into_iter()
        .filter(|x| *x)
        .count()
            != 1
        {
            anyhow::bail!("Specify exactly one of --reads, --pairs, --bases, or --frac");
        }
        if let Some(f) = frac {
            if !(0.0 < f && f <= 1.0) {
                anyhow::bail!("--frac must be in (0, 1]");
            }
        }

        #[inline]
        fn pairs_of(order: usize) -> u64 {
            let n = order as u64;
            n.saturating_mul(n.saturating_sub(1)) / 2
        }

        #[inline]
        fn hash_qname(q: &str) -> u64 {
            let mut h = DefaultHasher::new();
            q.hash(&mut h);
            h.finish()
        }

        log::info!(
            "PAF downsample pass1 scanning... (min_mapq={}, min_ident={}, min_len={}, order in [{}, {}))",
            min_quality,
            min_identity,
            min_length,
            min_order,
            max_order
        );

        // pass1: scan per query order and pairs
        let mut rdr = common_reader(&self.file);
        let mut line = String::new();

        let mut prev_q: Option<String> = None;
        let mut cur_order: usize = 0;
        let mut cur_qlen: u64 = 0;
        let mut per_read: Vec<(u64, u32, u64, u64)> = Vec::new();

        // let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        let mut not_grouped = false;
        let mut not_grouped_example_prev = String::new();
        let mut not_grouped_example_cur = String::new();

        while rdr.read_line(&mut line)? > 0 {
            let t = line.trim_end();
            if t.is_empty() || t.starts_with('#') {
                line.clear();
                continue;
            }

            let mut it = t.split('\t');
            let qname = it.next().unwrap_or("").to_string();

            if let Some(pq) = prev_q.as_ref() {
                if pq != &qname {
                    not_grouped = true;
                    not_grouped_example_prev = pq.clone();
                    not_grouped_example_cur = qname.clone();
     
                    // break;
                }
            }
            // seen.insert(qname.clone());

            let qlen = it.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            let qstart = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let qend = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);

            let _strand = it.next().unwrap_or("+");
            let _tname = it.next();
            let _tlen = it.next();
            let _tstart = it.next();
            let _tend = it.next();

            let nmatch = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let alen = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let mapq = it.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);

            let qcov = qend.saturating_sub(qstart);
            let identity = if alen > 0 { nmatch as f32 / alen as f32 } else { 0.0 };

            let is_eligible_piece =
                mapq >= min_quality && qcov >= min_length && identity >= min_identity;

            match prev_q.as_ref() {
                None => {
                    prev_q = Some(qname);
                    cur_qlen = qlen;
                    cur_order = if is_eligible_piece { 1 } else { 0 };
                }
                Some(pq) if *pq == qname => {
                    if is_eligible_piece {
                        cur_order += 1;
                    }
                }
                Some(pq) => {
                    if target_bases.is_some() || (cur_order >= min_order && cur_order < max_order) {
                        let p = pairs_of(cur_order);
                        per_read.push((hash_qname(pq), cur_order as u32, p, cur_qlen));
                    }
                    prev_q = Some(qname);
                    cur_qlen = qlen; 
                    cur_order = if is_eligible_piece { 1 } else { 0 };
                }
            }

            line.clear();
        }

        if not_grouped {
            log::warn!(
                "Input PAF is NOT grouped by query name (qname appears in multiple blocks). \
Example transition: prev_q='{}' then qname='{}' reappeared. \
Please sort/group by 1st column before downsampling, e.g.:\n\
  zcat in.paf.gz | sort -k1,1 -S 4G --parallel=8 | gzip -c > in.qsorted.paf.gz",
                not_grouped_example_prev,
                not_grouped_example_cur
            );
        }


        // finalize last
        if let Some(pq) = prev_q.take() {
            if target_bases.is_some() || (cur_order >= min_order && cur_order < max_order) {
                let p = pairs_of(cur_order);
                per_read.push((hash_qname(&pq), cur_order as u32, p, cur_qlen));
            }
        }

        let total_reads = per_read.len() as u64;
        let total_pairs: u64 = per_read.iter().map(|x| x.2).sum();
        let total_bases: u64 = per_read.iter().map(|x| x.3).sum();

        log::info!(
            "Eligible reads: {} ; eligible pairs: {} ; eligible bases: {}",
            total_reads,
            total_pairs,
            total_bases
        );



        let mut rng = StdRng::seed_from_u64(seed);
        per_read.shuffle(&mut rng);

        let mut chosen: HashSet<u64> = HashSet::new();

        if let Some(k) = target_reads {
            for (q, _ord, _p, _blen) in per_read.iter().take(k.min(total_reads) as usize) {
                chosen.insert(q.clone());
            }
            log::info!("Selected reads: {}", chosen.len());
        } else if let Some(tb) = target_bases {
            let mut acc: u64 = 0;
            for (q, _ord, _p, blen) in per_read.iter() {
                if acc >= tb {
                    break;
                }
                chosen.insert(q.clone());
                acc = acc.saturating_add(*blen);
            }
            log::info!(
                "Selected reads: {} (approx bases: {} / target {})",
                chosen.len(),
                acc,
                tb
            );
        } else if let Some(tp) = target_pairs {
            let mut acc: u64 = 0;
            for (q, _ord, p, _blen) in per_read.iter() {
                if acc >= tp {
                    break;
                }
                chosen.insert(q.clone());
                acc = acc.saturating_add(*p);
            }
            log::info!(
                "Selected reads: {} (approx pairs: {} / target {})",
                chosen.len(),
                acc,
                tp
            );
        } else if let Some(f) = frac {
            let want_pairs = frac_by_pairs;
            if want_pairs {
                let target_p = ((total_pairs as f64) * f).round() as u64;
                let mut acc: u64 = 0;
                for (q, _ord, p, _blen) in per_read.iter() {
                    if acc >= target_p {
                        break;
                    }
                    chosen.insert(q.clone());
                    acc = acc.saturating_add(*p);
                }
            } else {
                let target_r = ((total_reads as f64) * f).round() as u64;
                for (q, _ord, _p, _blen) in per_read.iter().take(target_r as usize) {
                    chosen.insert(q.clone());
                }
            }
        }

        // pass2: write
        log::info!("PAF downsample pass2 writing output to `{}`...", output);
        let mut rdr2 = common_reader(&self.file);
        let mut wtr = common_writer(output);
        let mut line2 = String::new();

        while rdr2.read_line(&mut line2)? > 0 {
            let t = line2.trim_end();
            if t.is_empty() {
                line2.clear();
                continue;
            }
            if t.starts_with('#') {
                if keep_comments {
                    wtr.write_all(line2.as_bytes())?;
                }
                line2.clear();
                continue;
            }

            let qname = t.split('\t').next().unwrap_or("");
            if chosen.contains(&hash_qname(qname)) {
                wtr.write_all(line2.as_bytes())?;
            }

            line2.clear();
        }

        log::info!(
            "PAF downsample finished. Output reads: {} (seed={})",
            chosen.len(),
            seed
        );
        Ok(())
    }

    pub fn downsample_by_vpc(
        &self,
        output: &String,
        seed: u64,
        min_quality: u8,
        min_identity: f32,
        min_length: u32,
        min_order: usize,
        max_order: usize,
        target_reads: Option<u64>,
        target_pairs: Option<u64>,
        target_bases: Option<u64>,
        frac: Option<f64>,
        frac_by_pairs: bool,
        keep_comments: bool,
    ) -> anyResult<()> {
        use rand::{Rng, SeedableRng};
        use rand::rngs::StdRng;
        use std::io::{BufRead, Write};
    
        if [
            target_reads.is_some(),
            target_pairs.is_some(),
            target_bases.is_some(),
            frac.is_some(),
        ]
        .into_iter()
        .filter(|x| *x)
        .count()
            != 1
        {
            anyhow::bail!("Specify exactly one of --reads, --pairs, --bases, or --frac");
        }
    
        if let Some(f) = frac {
            if !(0.0 < f && f <= 1.0) {
                anyhow::bail!("--frac must be in (0, 1]");
            }
        }
    
        #[inline]
        fn pairs_of(order: usize) -> u64 {
            let n = order as u64;
            n.saturating_mul(n.saturating_sub(1)) / 2
        }
    
        #[inline]
        fn parse_piece(
            line: &str,
            min_quality: u8,
            min_identity: f32,
            min_length: u32,
        ) -> Option<(String, u64, bool)> {
            let t = line.trim_end();
            if t.is_empty() || t.starts_with('#') {
                return None;
            }
    
            let mut it = t.split('\t');
    
            let qname = it.next().unwrap_or("").to_string();
            let qlen = it.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            let qstart = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let qend = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
    
            let _strand = it.next();
            let _tname = it.next();
            let _tlen = it.next();
            let _tstart = it.next();
            let _tend = it.next();
    
            let nmatch = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let alen = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let mapq = it.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
    
            let qcov = qend.saturating_sub(qstart);
            let identity = if alen > 0 {
                nmatch as f32 / alen as f32
            } else {
                0.0
            };
    
            let is_eligible_piece =
                mapq >= min_quality && qcov >= min_length && identity >= min_identity;
    
            Some((qname, qlen, is_eligible_piece))
        }
    
        log::info!(
            "PAF downsample pass1 scanning totals... (min_mapq={}, min_ident={}, min_len={}, order in [{}, {}))",
            min_quality,
            min_identity,
            min_length,
            min_order,
            max_order
        );
    
        // -------------------------
        // pass1: only compute totals
        // -------------------------
        let mut rdr = common_reader(&self.file);
        let mut line = String::new();
    
        let mut prev_q: Option<String> = None;
        let mut cur_order: usize = 0;
        let mut cur_qlen: u64 = 0;
    
        let mut total_reads: u64 = 0;
        let mut total_pairs: u64 = 0;
        let mut total_bases: u64 = 0;
    
        while rdr.read_line(&mut line)? > 0 {
            if let Some((qname, qlen, is_eligible_piece)) =
                parse_piece(&line, min_quality, min_identity, min_length)
            {
                match prev_q.as_ref() {
                    None => {
                        prev_q = Some(qname);
                        cur_qlen = qlen;
                        cur_order = if is_eligible_piece { 1 } else { 0 };
                    }
                    Some(pq) if *pq == qname => {
                        if is_eligible_piece {
                            cur_order += 1;
                        }
                    }
                    Some(_) => {
                        if target_bases.is_some()
                            || (cur_order >= min_order && cur_order < max_order)
                        {
                            total_reads += 1;
                            total_pairs = total_pairs.saturating_add(pairs_of(cur_order));
                            total_bases = total_bases.saturating_add(cur_qlen);
                        }
    
                        prev_q = Some(qname);
                        cur_qlen = qlen;
                        cur_order = if is_eligible_piece { 1 } else { 0 };
                    }
                }
            }
    
            line.clear();
        }
    
        if prev_q.take().is_some() {
            if target_bases.is_some() || (cur_order >= min_order && cur_order < max_order) {
                total_reads += 1;
                total_pairs = total_pairs.saturating_add(pairs_of(cur_order));
                total_bases = total_bases.saturating_add(cur_qlen);
            }
        }
    
        log::info!(
            "Eligible reads: {} ; eligible pairs: {} ; eligible bases: {}",
            total_reads,
            total_pairs,
            total_bases
        );
    
        if total_reads == 0 {
            log::warn!("No eligible reads found. Writing empty output.");
            let _wtr = common_writer(output);
            return Ok(());
        }
    
        // Selection probability.
        //
        // This version avoids storing per_read and chosen.
        // Therefore --reads / --pairs / --bases are approximate in expectation.
        let keep_prob: f64 = if let Some(k) = target_reads {
            (k as f64 / total_reads as f64).min(1.0)
        } else if let Some(tb) = target_bases {
            if total_bases == 0 {
                0.0
            } else {
                (tb as f64 / total_bases as f64).min(1.0)
            }
        } else if let Some(tp) = target_pairs {
            if total_pairs == 0 {
                0.0
            } else {
                (tp as f64 / total_pairs as f64).min(1.0)
            }
        } else if let Some(f) = frac {
            let _ = frac_by_pairs;
            f
        } else {
            unreachable!()
        };
    
        log::info!("Streaming selection probability: {:.6}", keep_prob);
    
        // -------------------------
        // pass2: decide and write by block
        // -------------------------
        log::info!("PAF downsample pass2 streaming write to `{}`...", output);
    
        let mut rng = StdRng::seed_from_u64(seed);
        let mut rdr2 = common_reader(&self.file);
        let mut wtr = common_writer(output);
        let mut line2 = String::new();
    
        let mut prev_q2: Option<String> = None;
        let mut cur_order2: usize = 0;
        let mut cur_qlen2: u64 = 0;
        let mut block_lines: Vec<String> = Vec::new();
    
        let mut out_reads: u64 = 0;
        let mut out_pairs: u64 = 0;
        let mut out_bases: u64 = 0;
    
        let mut flush_block = |
            block_lines: &mut Vec<String>,
            cur_order: usize,
            cur_qlen: u64,
            rng: &mut StdRng,
            wtr: &mut dyn Write,
            out_reads: &mut u64,
            out_pairs: &mut u64,
            out_bases: &mut u64,
        | -> anyResult<()> {
            if block_lines.is_empty() {
                return Ok(());
            }
    
            let eligible =
                target_bases.is_some() || (cur_order >= min_order && cur_order < max_order);
    
            if eligible {
                let keep = rng.r#gen::<f64>() < keep_prob;
                if keep {
                    for l in block_lines.iter() {
                        wtr.write_all(l.as_bytes())?;
                    }
    
                    *out_reads += 1;
                    *out_pairs = out_pairs.saturating_add(pairs_of(cur_order));
                    *out_bases = out_bases.saturating_add(cur_qlen);
                }
            }
    
            block_lines.clear();
            Ok(())
        };
    
        while rdr2.read_line(&mut line2)? > 0 {
            let t = line2.trim_end();
    
            if t.is_empty() {
                line2.clear();
                continue;
            }
    
            if t.starts_with('#') {
                if keep_comments {
                    wtr.write_all(line2.as_bytes())?;
                }
                line2.clear();
                continue;
            }
    
            if let Some((qname, qlen, is_eligible_piece)) =
                parse_piece(&line2, min_quality, min_identity, min_length)
            {
                match prev_q2.as_ref() {
                    None => {
                        prev_q2 = Some(qname);
                        cur_qlen2 = qlen;
                        cur_order2 = if is_eligible_piece { 1 } else { 0 };
                        block_lines.push(line2.clone());
                    }
                    Some(pq) if *pq == qname => {
                        if is_eligible_piece {
                            cur_order2 += 1;
                        }
                        block_lines.push(line2.clone());
                    }
                    Some(_) => {
                        flush_block(
                            &mut block_lines,
                            cur_order2,
                            cur_qlen2,
                            &mut rng,
                            &mut wtr,
                            &mut out_reads,
                            &mut out_pairs,
                            &mut out_bases,
                        )?;
    
                        prev_q2 = Some(qname);
                        cur_qlen2 = qlen;
                        cur_order2 = if is_eligible_piece { 1 } else { 0 };
                        block_lines.push(line2.clone());
                    }
                }
            }
    
            line2.clear();
        }
    
        flush_block(
            &mut block_lines,
            cur_order2,
            cur_qlen2,
            &mut rng,
            &mut wtr,
            &mut out_reads,
            &mut out_pairs,
            &mut out_bases,
        )?;
    
        log::info!(
            "PAF downsample finished. Output reads: {} ; output pairs: {} ; output bases: {} ; seed={}",
            out_reads,
            out_pairs,
            out_bases,
            seed
        );
    
        Ok(())
    }

    pub fn downsample_by_vpc_multi(
        &self,
        output: &String,
        seed: u64,
        min_quality: u8,
        min_identity: f32,
        min_length: u32,
        min_order: usize,
        max_order: usize,
        target_reads: Option<u64>,
        target_pairs: Option<u64>,
        target_bases: Option<u64>,
        frac: Option<f64>,
        frac_by_pairs: bool,
        keep_comments: bool,
    ) -> anyResult<()> {
        use anyhow::Context;
        use crossbeam_channel::{bounded, Receiver, Sender};
        use rand::{Rng, SeedableRng};
        use rand::rngs::StdRng;
        use std::collections::BTreeMap;
        use std::io::{BufRead, Write};
        use std::sync::Arc;
        use std::thread;
    
        #[derive(Debug)]
        struct Block {
            idx: u64,
            qname: String,
            lines: Vec<String>,
        }
    
        #[derive(Debug)]
        struct ResultBlock {
            idx: u64,
            keep: bool,
            lines: Vec<String>,
            order: usize,
            pairs: u64,
            bases: u64,
        }
    
        #[inline]
        fn pairs_of(order: usize) -> u64 {
            let n = order as u64;
            n.saturating_mul(n.saturating_sub(1)) / 2
        }
    
        #[inline]
        fn parse_line(
            line: &str,
            min_quality: u8,
            min_identity: f32,
            min_length: u32,
        ) -> Option<(String, u64, bool)> {
            let t = line.trim_end();
    
            if t.is_empty() || t.starts_with('#') {
                return None;
            }
    
            let mut it = t.split('\t');
    
            let qname = it.next().unwrap_or("").to_string();
            let qlen = it.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            let qstart = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let qend = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
    
            let _strand = it.next();
            let _tname = it.next();
            let _tlen = it.next();
            let _tstart = it.next();
            let _tend = it.next();
    
            let nmatch = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let alen = it.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let mapq = it.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
    
            let qcov = qend.saturating_sub(qstart);
            let identity = if alen > 0 {
                nmatch as f32 / alen as f32
            } else {
                0.0
            };
    
            let eligible =
                mapq >= min_quality && qcov >= min_length && identity >= min_identity;
    
            Some((qname, qlen, eligible))
        }
    
        fn process_block(
            block: Block,
            min_quality: u8,
            min_identity: f32,
            min_length: u32,
            min_order: usize,
            max_order: usize,
            target_bases: bool,
            keep_prob: f64,
            seed: u64,
        ) -> ResultBlock {
            let mut order: usize = 0;
            let mut bases: u64 = 0;
    
            for line in &block.lines {
                if let Some((_q, qlen, eligible_piece)) =
                    parse_line(line, min_quality, min_identity, min_length)
                {
                    bases = qlen;
                    if eligible_piece {
                        order += 1;
                    }
                }
            }
    
            let pairs = pairs_of(order);
    
            let eligible_read =
                target_bases || (order >= min_order && order < max_order);
    
            // Deterministic per-block RNG.
            // This avoids sharing one RNG across threads and keeps output reproducible.
            let mut rng = StdRng::seed_from_u64(seed ^ block.idx.wrapping_mul(0x9E3779B97F4A7C15));
            let keep = eligible_read && rng.r#gen::<f64>() < keep_prob;
    
            ResultBlock {
                idx: block.idx,
                keep,
                lines: block.lines,
                order,
                pairs,
                bases,
            }
        }
    
        if [
            target_reads.is_some(),
            target_pairs.is_some(),
            target_bases.is_some(),
            frac.is_some(),
        ]
        .into_iter()
        .filter(|x| *x)
        .count()
            != 1
        {
            anyhow::bail!("Specify exactly one of --reads, --pairs, --bases, or --frac");
        }
    
        if let Some(f) = frac {
            if !(0.0 < f && f <= 1.0) {
                anyhow::bail!("--frac must be in (0, 1]");
            }
        }
    
        log::info!(
            "PAF downsample pass1 scanning totals... (min_mapq={}, min_ident={}, min_len={}, order in [{}, {}))",
            min_quality,
            min_identity,
            min_length,
            min_order,
            max_order
        );
    
        // -----------------------------
        // pass1: streaming total count
        // -----------------------------
        let mut rdr = common_reader(&self.file);
        let mut line = String::new();
    
        let mut prev_q: Option<String> = None;
        let mut cur_order: usize = 0;
        let mut cur_qlen: u64 = 0;
    
        let mut total_reads: u64 = 0;
        let mut total_pairs: u64 = 0;
        let mut total_bases: u64 = 0;
    
        while rdr.read_line(&mut line)? > 0 {
            if let Some((qname, qlen, eligible_piece)) =
                parse_line(&line, min_quality, min_identity, min_length)
            {
                match prev_q.as_ref() {
                    None => {
                        prev_q = Some(qname);
                        cur_qlen = qlen;
                        cur_order = if eligible_piece { 1 } else { 0 };
                    }
                    Some(pq) if *pq == qname => {
                        if eligible_piece {
                            cur_order += 1;
                        }
                    }
                    Some(_) => {
                        if target_bases.is_some()
                            || (cur_order >= min_order && cur_order < max_order)
                        {
                            total_reads += 1;
                            total_pairs = total_pairs.saturating_add(pairs_of(cur_order));
                            total_bases = total_bases.saturating_add(cur_qlen);
                        }
    
                        prev_q = Some(qname);
                        cur_qlen = qlen;
                        cur_order = if eligible_piece { 1 } else { 0 };
                    }
                }
            }
    
            line.clear();
        }
    
        if prev_q.take().is_some() {
            if target_bases.is_some() || (cur_order >= min_order && cur_order < max_order) {
                total_reads += 1;
                total_pairs = total_pairs.saturating_add(pairs_of(cur_order));
                total_bases = total_bases.saturating_add(cur_qlen);
            }
        }
    
        log::info!(
            "Eligible reads: {} ; eligible pairs: {} ; eligible bases: {}",
            total_reads,
            total_pairs,
            total_bases
        );
    
        if total_reads == 0 {
            log::warn!("No eligible reads found. Writing empty output.");
            let _wtr = common_writer(output);
            return Ok(());
        }
    
        let keep_prob: f64 = if let Some(k) = target_reads {
            (k as f64 / total_reads as f64).min(1.0)
        } else if let Some(tb) = target_bases {
            if total_bases == 0 {
                0.0
            } else {
                (tb as f64 / total_bases as f64).min(1.0)
            }
        } else if let Some(tp) = target_pairs {
            if total_pairs == 0 {
                0.0
            } else {
                (tp as f64 / total_pairs as f64).min(1.0)
            }
        } else if let Some(f) = frac {
            let _ = frac_by_pairs;
            f
        } else {
            unreachable!()
        };
    
        log::info!("Streaming selection probability: {:.6}", keep_prob);
    
        // -----------------------------
        // pass2: producer-consumer
        // -----------------------------
        log::info!("PAF downsample pass2 producer-consumer writing `{}`...", output);
    
        let n_workers = std::thread::available_parallelism()
            .map(|x| x.get())
            .unwrap_or(4)
            .max(1);
    
        log::info!("Using {} worker threads", n_workers);
    
        let (task_tx, task_rx): (Sender<Block>, Receiver<Block>) = bounded(n_workers * 4);
        let (res_tx, res_rx): (Sender<ResultBlock>, Receiver<ResultBlock>) = bounded(n_workers * 4);
    
        let file = self.file.clone();
        let output = output.clone();
    
        // Producer: read file and group by qname block.
        let producer_handle = thread::spawn(move || -> anyResult<()> {
            let mut rdr = common_reader(&file);
            let mut line = String::new();
    
            let mut idx: u64 = 0;
            let mut cur_q: Option<String> = None;
            let mut cur_lines: Vec<String> = Vec::new();
    
            while rdr.read_line(&mut line)? > 0 {
                let t = line.trim_end();
    
                if t.is_empty() {
                    line.clear();
                    continue;
                }
    
                if t.starts_with('#') {
                    // Comments are handled by writer only if they are emitted.
                    // To keep implementation simple, comments are skipped here.
                    line.clear();
                    continue;
                }
    
                let qname = t.split('\t').next().unwrap_or("").to_string();
    
                match cur_q.as_ref() {
                    None => {
                        cur_q = Some(qname);
                        cur_lines.push(line.clone());
                    }
                    Some(q) if *q == qname => {
                        cur_lines.push(line.clone());
                    }
                    Some(q) => {
                        let block = Block {
                            idx,
                            qname: q.clone(),
                            lines: std::mem::take(&mut cur_lines),
                        };
    
                        task_tx
                            .send(block)
                            .context("failed to send block to workers")?;
    
                        idx += 1;
                        cur_q = Some(qname);
                        cur_lines.push(line.clone());
                    }
                }
    
                line.clear();
            }
    
            if let Some(q) = cur_q.take() {
                if !cur_lines.is_empty() {
                    let block = Block {
                        idx,
                        qname: q,
                        lines: cur_lines,
                    };
    
                    task_tx
                        .send(block)
                        .context("failed to send final block to workers")?;
                }
            }
    
            Ok(())
        });
    
        // Workers.
        let mut worker_handles = Vec::new();
        let target_bases_flag = target_bases.is_some();
    
        for worker_id in 0..n_workers {
            let rx = task_rx.clone();
            let tx = res_tx.clone();
    
            let handle = thread::spawn(move || -> anyResult<()> {
                while let Ok(block) = rx.recv() {
                    let res = process_block(
                        block,
                        min_quality,
                        min_identity,
                        min_length,
                        min_order,
                        max_order,
                        target_bases_flag,
                        keep_prob,
                        seed ^ ((worker_id as u64) << 32),
                    );
    
                    tx.send(res)
                        .context("failed to send result block to writer")?;
                }
    
                Ok(())
            });
    
            worker_handles.push(handle);
        }
    
        drop(task_rx);
        drop(res_tx);
    
        // Consumer / writer: preserve original block order.
        let mut wtr = common_writer(&output);
        let mut pending: BTreeMap<u64, ResultBlock> = BTreeMap::new();
        let mut next_idx: u64 = 0;
    
        let mut out_reads: u64 = 0;
        let mut out_pairs: u64 = 0;
        let mut out_bases: u64 = 0;
    
        for res in res_rx.iter() {
            pending.insert(res.idx, res);
    
            while let Some(res) = pending.remove(&next_idx) {
                if res.keep {
                    for l in res.lines {
                        wtr.write_all(l.as_bytes())?;
                    }
    
                    out_reads += 1;
                    out_pairs = out_pairs.saturating_add(res.pairs);
                    out_bases = out_bases.saturating_add(res.bases);
                }
    
                next_idx += 1;
            }
        }
    
        producer_handle
            .join()
            .map_err(|_| anyhow::anyhow!("producer thread panicked"))??;
    
        for h in worker_handles {
            h.join()
                .map_err(|_| anyhow::anyhow!("worker thread panicked"))??;
        }
    
        log::info!(
            "PAF downsample finished. Output reads: {} ; output pairs: {} ; output bases: {} ; seed={}",
            out_reads,
            out_pairs,
            out_bases,
            seed
        );
    
        Ok(())
    }

    pub fn phase_reads(
        &self,
        contigs_group_file: &String,
        output: &String,
        min_mapq: u8,
    ) -> anyResult<()> {
        use std::collections::HashMap;
        use std::io::{BufRead, Write};

        let mut contig_to_group: HashMap<String, String> = HashMap::new();
        {
            let rdr = common_reader(contigs_group_file);
            let buf = std::io::BufReader::new(rdr);
            for line in buf.lines() {
                let line = line?;
                let parts: Vec<&str> = line.trim().split_whitespace().collect();
                if parts.len() >= 2 {
                    // parts[0] = contig name, parts[1] = group name (e.g. Hap1)
                    contig_to_group.insert(parts[0].to_string(), parts[1].to_string());
                }
            }
        }
        log::info!("Loaded {} contig-to-group mappings.", contig_to_group.len());

        let mut rdr = common_reader(&self.file);
        let mut wtr = common_writer(output);
        let mut line = String::new();

        let mut prev_q: Option<String> = None;
        let mut group_scores: HashMap<String, u32> = HashMap::new();

        let mut emit_read = |qname: &str, scores: &mut HashMap<String, u32>, w: &mut dyn Write| -> anyResult<()> {
            if let Some((best_group, _)) = scores.iter().max_by_key(|entry| entry.1) {
                writeln!(w, "{}\t{}", qname, best_group)?;
            } else {
                writeln!(w, "{}\tUnassigned", qname)?;
            }
            scores.clear();
            Ok(())
        };

        log::info!("Scanning PAF file to phase reads...");

        let mut total_reads = 0;
        let mut phased_reads = 0;

        while rdr.read_line(&mut line)? > 0 {
            let t = line.trim_end();
            if t.is_empty() || t.starts_with('#') {
                line.clear();
                continue;
            }

            let mut it = t.split('\t');
            let qname = it.next().unwrap_or("").to_string();
            
            // 跳过不必要的列
            it.next(); // qlen
            it.next(); // qstart
            it.next(); // qend
            it.next(); // strand

            let tname = it.next().unwrap_or("");
            
            it.next(); // tlen
            it.next(); // tstart
            it.next(); // tend
            it.next(); // nmatch
            it.next(); // alen
            let mapq = it.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);

            let is_switch = match prev_q.as_ref() {
                Some(pq) => pq != &qname,
                None => false,
            };

            if is_switch {
                emit_read(prev_q.as_ref().unwrap(), &mut group_scores, &mut wtr)?;
                total_reads += 1;
                if !group_scores.is_empty() { phased_reads += 1; }
                prev_q = Some(qname.clone());
            } else if prev_q.is_none() {
                prev_q = Some(qname.clone());
            }

            if mapq >= min_mapq {
                if let Some(group) = contig_to_group.get(tname) {
                    *group_scores.entry(group.clone()).or_insert(0) += 1;
                }
            }

            line.clear();
        }

        if let Some(pq) = prev_q {
            emit_read(&pq, &mut group_scores, &mut wtr)?;
            total_reads += 1;
            if !group_scores.is_empty() { phased_reads += 1; }
        }

        log::info!(
            "Phasing finished. Total reads processed: {}, successfully assigned: {}",
            total_reads, phased_reads
        );
        log::info!("Output written to: {}", output);

        Ok(())
    }
}


   