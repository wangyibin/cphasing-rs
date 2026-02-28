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
use std::io::{ Write, BufReader, BufRead };
use std::fmt::{ Write as FmtWrite };
use serde::{ Deserialize, Serialize};
use std::time::{Instant, Duration};
use rust_lapper::{Interval, Lapper};
use rayon::prelude::*;

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

    pub fn paf2table(&self, bed: &String,
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

    pub fn paf2table_fast(&self, bed: &String,
            output: &String, min_quality: &u8, 
            min_identity: &f32,  min_length: &u32,
            max_order: &u32, 
            max_edge_length: &u64,
            secondary: bool    
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

        let num_threads = 20; 
        log::info!("Using {} threads for processing.", num_threads);
        
        let (tx_raw, rx_raw) = bounded::<(u64, Vec<u8>)>(100); 
        let (tx_writer, rx_writer) = bounded::<String>(1000);


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

        for _ in 0..num_threads {
            let rx = rx_raw.clone();
            let tx_w = tx_writer.clone();
            let ih = interval_hash.clone();
            let s_p = stats_pass.clone();
            let s_s = stats_sing.clone();
            let s_l = stats_low.clone();
            let s_c = stats_comp.clone();

            handles.push(thread::spawn(move || {
                let mut local_p=0; let mut local_s=0; let mut local_l=0; let mut local_c=0;
                let mut output_buffer = String::with_capacity(WRITE_CHUNK);
                
                
                fn parse_line(line: &[u8], min_q: u8, min_l: u32, min_id: f32, sec: bool) -> Option<(PoreCRecord, String)> {
                    let mut iter = line.split(|&b| b == b'\t');
                    let q_bytes = iter.next()?;
                    
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
                    let q_str = unsafe { String::from_utf8_unchecked(q_bytes.to_vec()) };

                    Some(PoreCRecord {
                            read_idx: 0, 
                            query_length: qlen, 
                            query_start: qs, query_end: qe, query_strand: strand_c,
                            target: t_str, target_length: tlen,
                            target_start: ts, target_end: te,
                            mapq: mapq, identity: ident, filter_reason: fr.to_string(),
                    }).map(|r| (r, q_str))
                }

                while let Ok((start_idx, batch_bytes)) = rx.recv() {
                    let mut current_idx = start_idx;
                    
                    let mut start = 0;
                    let mut current_group: Vec<PoreCRecord> = Vec::with_capacity(32);
                    let mut old_q_str = String::new();
                    
                    while start < batch_bytes.len() {
                        let end = match batch_bytes[start..].iter().position(|&b| b == b'\n') {
                                Some(p) => start + p,
                                None => batch_bytes.len(),
                        };
                        let line = &batch_bytes[start..end];
                        start = end + 1; 
                        
                        if line.is_empty() { continue; }
                    
                        if let Some((mut pcr, q_str)) = parse_line(line, min_quality, min_length, min_identity, secondary) {
                            
                            let is_switch = !old_q_str.is_empty() && old_q_str != q_str;
                            
                            if is_switch {
                        
                                for r in &mut current_group { r.read_idx = current_idx; }
                                
                                process_group(&mut current_group, &ih, is_filter_digest, is_filter_edges, max_edge_length, max_order, &mut output_buffer, &mut local_p, &mut local_s, &mut local_l, &mut local_c);
                                
                                if output_buffer.len() >= WRITE_CHUNK {
                                    tx_w.send(std::mem::take(&mut output_buffer)).unwrap();
                                    output_buffer = String::with_capacity(WRITE_CHUNK);
                                }
                                
                    
                                current_idx += 1;
                            }
                            
                            if old_q_str.is_empty() || is_switch {
                                old_q_str = q_str;
                            }
                            current_group.push(pcr);
                        }
                    }
                    
                    if !current_group.is_empty() {
                            for r in &mut current_group { r.read_idx = current_idx; }
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
        const IO_BUF_SIZE: usize = 1024 * 1024; // 1MB read buffer
        let mut buffer = vec![0u8; IO_BUF_SIZE];
        let mut leftover = Vec::new();

        const BATCH_SIZE_LIMIT: usize = 16 * 1024 * 1024; 
        let mut current_batch = Vec::with_capacity(BATCH_SIZE_LIMIT + IO_BUF_SIZE);
        
        let mut old_query_bytes = Vec::with_capacity(128);
        let mut first = true;
        
        let mut total_groups: u64 = 0;
        let mut batch_start_idx: u64 = 0;

        loop {
            let n = reader.read(&mut buffer).unwrap_or(0);
            if n == 0 { break; }

            let data = if !leftover.is_empty() {
                let mut v = Vec::with_capacity(leftover.len() + n);
                v.extend_from_slice(&leftover);
                v.extend_from_slice(&buffer[..n]);
                leftover.clear();
                Cow::Owned(v)
            } else {
                Cow::Borrowed(&buffer[..n])
            };
            
            let data_slice = &data;
            let mut start_pos = 0;

            while start_pos < data_slice.len() {
                let end_pos = match data_slice[start_pos..].iter().position(|&b| b == b'\n') {
                    Some(p) => start_pos + p,
                    None => {
                        leftover.extend_from_slice(&data_slice[start_pos..]);
                        break; 
                    }
                };
                let mut line = &data_slice[start_pos..end_pos];
                start_pos = end_pos + 1;

                if !line.is_empty() && line[line.len()-1] == b'\r' { line = &line[..line.len()-1]; }
                if line.is_empty() || line[0] == b'#' { continue; }

                let tab_pos = match line.iter().position(|&b| b == b'\t') {
                    Some(p) => p, None => continue,
                };
                let q_bytes = &line[..tab_pos];

                let is_switch = if first {
                    first = false;
                    old_query_bytes.extend_from_slice(q_bytes);
                    false
                } else {
                    q_bytes != old_query_bytes.as_slice()
                };

                if is_switch {
                    total_groups += 1;
                    
                    old_query_bytes.clear();
                    old_query_bytes.extend_from_slice(q_bytes);

                    if current_batch.len() >= BATCH_SIZE_LIMIT {
                        tx_raw.send((batch_start_idx, std::mem::take(&mut current_batch))).unwrap();
            
                        batch_start_idx = total_groups;
                        
                        current_batch = Vec::with_capacity(BATCH_SIZE_LIMIT + IO_BUF_SIZE);
                    }
                }
                
                current_batch.extend_from_slice(line);
                current_batch.push(b'\n');
            }
        }
        
        if !current_batch.is_empty() {
            tx_raw.send((batch_start_idx, current_batch)).unwrap();
            
            total_groups += 1;
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
    
        summary.mapping = total_groups; 

        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(output).with_extension("").to_str().unwrap().to_string()
        };
        summary.save(&format!("{}.read.summary", output_prefix));

        log::info!("Successful output Pore-C table `{}` (raw-bytes-batched)", output);
        Ok(())
    }
    
    pub fn paf2table1(&self, bed: &String,
            output: &String, min_quality: &u8, 
            min_identity: &f32,  min_length: &u32,
            max_order: &u32, 
            max_edge_length: &u64,
            secondary: bool    
        ) -> Result<(), Box<dyn Error>> {
        
        type IvU8 = Interval<usize, u8>;
        let min_quality = *min_quality;
        let min_identity = *min_identity;
        let min_length = *min_length;
        let max_order = *max_order;
        let max_edge_length = *max_edge_length;

        let is_filter_digest = !bed.is_empty();
        let bed = if is_filter_digest { Bed3::new(bed) } else { Bed3::new(&String::from(".tmp.bed")) };
        let interval_hash = bed.to_interval_hash();
        let is_filter_edges = max_edge_length > 0;

        let num_threads = 10;
        log::info!("Using {} threads for processing.", num_threads);
        let (sender, receiver) = bounded::<Vec<PoreCRecord>>(8192); 
        let (writer_sender, writer_receiver) = bounded::<String>(8192); 

        let mut handles = Vec::with_capacity(num_threads);
        // let prof = std::sync::Arc::new(Prof::new());
        

        for _ in 0..num_threads {
            let receiver: Receiver<Vec<PoreCRecord>> = receiver.clone();
            let writer_sender = writer_sender.clone();
            let interval_hash = interval_hash.clone();
            // let prof_cl = prof.clone();

            handles.push(thread::spawn(move || {
                let init_cap = (max_order as usize) + 2;
                let mut concatemer: Concatemer = Concatemer::with_capacity(init_cap);
                let mut local_stats = (0, 0, 0, 0); // pass, singleton, low_mq, complex
                let mut output_buffer = String::with_capacity(WRITE_CHUNK);
                let mut current_read_idx: Option<u64> = None;
                while let Ok(batch) = receiver.recv() {
                    // for pcr in batch {
                    //     match current_read_idx {
                    //         Some(rid) if rid == pcr.read_idx => {
                    //             concatemer.push(pcr);
                    //         }
                    //         Some(_) => {
                    //             let (final_status, pass_indices) = concatemer.process(
                    //                 is_filter_digest,
                    //                 &interval_hash,
                    //                 is_filter_edges,
                    //                 max_edge_length,
                    //                 max_order,
                    //             );
                    //             if final_status == PASS_STR {
                    //                 for &idx in &pass_indices {
                    //                     output_buffer.push_str(&concatemer.records[idx].to_string());
                    //                     output_buffer.push('\n');
                    //                 }
                    //             }
                    //             match final_status {
                    //                 PASS_STR => local_stats.0 += 1,
                    //                 SINGLETON_STR => local_stats.1 += 1,
                    //                 LOW_MQ_STR => local_stats.2 += 1,
                    //                 COMPLEX_STR => local_stats.3 += 1,
                    //                 _ => {}
                    //             }
                    //             concatemer.clear();
                    //             current_read_idx = Some(pcr.read_idx);
                    //             concatemer.push(pcr);
                    //         }
                    //         None => {
                    //             current_read_idx = Some(pcr.read_idx);
                    //             concatemer.push(pcr);
                    //         }
                    //     }
                    // }
                    let mut groups: Vec<Vec<PoreCRecord>> = Vec::new();
                    let mut cur: Vec<PoreCRecord> = Vec::new();
                    let mut cur_idx: Option<u64> = None;
                    for pcr in batch {
                        match cur_idx {
                            Some(rid) if rid == pcr.read_idx => cur.push(pcr),
                            Some(_) => {
                                groups.push(std::mem::take(&mut cur));
                                cur_idx = Some(pcr.read_idx);
                                cur.push(pcr);
                            }
                            None => {
                                cur_idx = Some(pcr.read_idx);
                                cur.push(pcr);
                            }
                        }
                    }
                    if !cur.is_empty() {
                        groups.push(cur);
                    }

                    let results: Vec<(String, (u64, u64, u64, u64))> = groups
                        .into_iter()
                        .map(|grp| {
                            let mut concatemer = Concatemer::with_capacity(grp.len() + 2);
                            for r in grp {
                                concatemer.push(r);
                            }
                            let (final_status, pass_indices) = concatemer.process(
                                is_filter_digest,
                                &interval_hash,
                                is_filter_edges,
                                max_edge_length,
                                max_order,
                            );
     
                            let mut out = String::new();
                            if final_status == PASS_STR {
                                for &idx in &pass_indices {
                                    out.push_str(&concatemer.records[idx].to_string());
                                    out.push('\n');
                                }
                                // out.reserve(pass_indices.len() * 128);
                                // for &idx in &pass_indices {
                                //     concatemer.records[idx].write_to(&mut out);
                                // }
                            }
                            let stats = match final_status {
                                PASS_STR      => (1, 0, 0, 0),
                                SINGLETON_STR => (0, 1, 0, 0),
                                LOW_MQ_STR    => (0, 0, 1, 0),
                                COMPLEX_STR   => (0, 0, 0, 1),
                                _             => (0, 0, 0, 0),
                            };
                            (out, stats)
                        })
                        .collect();
                    for (chunk, (a, b, c, d)) in results {
                        if !chunk.is_empty() {
                            output_buffer.push_str(&chunk);
                        }
                        local_stats.0 += a;
                        local_stats.1 += b;
                        local_stats.2 += c;
                        local_stats.3 += d;

                        if output_buffer.len() >= WRITE_CHUNK {
                            let chunk = std::mem::take(&mut output_buffer);
                            writer_sender.send(chunk).unwrap();
                            output_buffer = String::with_capacity(WRITE_CHUNK);
                        }
                    }


                }

                // if !concatemer.records.is_empty() {
                //     let (final_status, pass_indices) = concatemer.process(
                //         is_filter_digest,
                //         &interval_hash,
                //         is_filter_edges,
                //         max_edge_length,
                //         max_order,
                //     );
                //     if final_status == PASS_STR {
                //         for &idx in &pass_indices {
                //             output_buffer.push_str(&concatemer.records[idx].to_string());
                //             output_buffer.push('\n');
                //         }
                //     }
                //     match final_status {
                //         PASS_STR => local_stats.0 += 1,
                //         SINGLETON_STR => local_stats.1 += 1,
                //         LOW_MQ_STR => local_stats.2 += 1,
                //         COMPLEX_STR => local_stats.3 += 1,
                //         _ => {}
                //     }
                // }
                if !output_buffer.is_empty() {
                    writer_sender.send(output_buffer).unwrap();
                }
                local_stats
            }));
        }
        drop(writer_sender); // Drop the original sender for the writer

        // Writer thread
        let output_clone = output.clone();
        // let prof_w = prof.clone();
        let writer_handle = thread::spawn(move || {
            let mut writer = common_writer(&output_clone);
            while let Ok(data_string) = writer_receiver.recv() {
                // let t0 = Instant::now();
                writer.write_all(data_string.as_bytes()).unwrap();
                // if prof_w.on {
                //     Prof::add(&prof_w.writer_ns, t0.elapsed().as_nanos());
                //     prof_w.out_bytes.fetch_add(data_string.len() as u64, Ordering::Relaxed);
                // }
            }
        });

        let mut rdr = self.parse2(secondary.clone())?;
        let mut old_query: Option<String> = None;
        let mut current_read_records: Vec<PoreCRecord> = Vec::new();
        let mut read_idx: u64 = 0;
        let mut batch_records: Vec<PoreCRecord> = Vec::with_capacity(RECORD_BATCH_SIZE);
        let mut groups_in_batch: usize = 0;

        for record in rdr {

            let is_switch = match old_query.as_deref() {
                Some(prev) => prev != record.query.as_str(),
                None => false,
            };
            if is_switch && old_query.is_some() {
                batch_records.extend(current_read_records.drain(..));
                groups_in_batch += 1;
                read_idx += 1;

                if groups_in_batch >= GROUP_BATCH_SIZE || batch_records.len() >= RECORD_BATCH_SIZE {
                    sender.send(std::mem::take(&mut batch_records)).unwrap();
                    groups_in_batch = 0;
                    batch_records = Vec::with_capacity(RECORD_BATCH_SIZE);
                }
            }
            let length: u32 = record.query_end - record.query_start;
            let identity = if record.alignment_length > 0 {
                record.match_n as f32 / record.alignment_length as f32
            } else {
                0.0
            };
            let fr = if record.mapq < min_quality || length < min_length || identity < min_identity {
                LOW_MQ_STR
            } else {
                PASS_STR
            };
            if old_query.is_none() || is_switch {
                old_query = Some(record.query.clone());
            }
            let pcr = PoreCRecord::from_paf_record(record, read_idx, identity, fr.to_string());
            current_read_records.push(pcr);
        }

        // Send the last group
        if !current_read_records.is_empty() {
            batch_records.extend(current_read_records.drain(..));
            groups_in_batch += 1;
            read_idx += 1;
        }
        if !batch_records.is_empty() {
            sender.send(batch_records).unwrap();
        }
        drop(sender);
        let mut read_pass_count: u64 = 0;
        let mut read_singleton_count: u64 = 0;
        let mut read_low_mq_count: u64 = 0;
        let mut read_complex_count: u64 = 0;

        for handle in handles {
        let stats = handle.join().unwrap();
            read_pass_count += stats.0;
            read_singleton_count += stats.1;
            read_low_mq_count += stats.2;
            read_complex_count += stats.3;
        }
        writer_handle.join().unwrap();
        // if prof.on {
        //     let rd = prof.read_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let wp = prof.worker_proc_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let wf = prof.worker_fmt_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let ww = prof.writer_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let out_mb = prof.out_bytes.load(Ordering::Relaxed) as f64 / (1024.0*1024.0);
        //     log::info!("PROFILE read+group: {:.3}s, worker process: {:.3}s, worker format: {:.3}s, write(compress): {:.3}s, wrote {:.1} MiB",
        //                rd, wp, wf, ww, out_mb);
        // }
        let mut summary: ReadSummary = ReadSummary::new();
        summary.pass = read_pass_count;
        summary.singleton = read_singleton_count;
        summary.low_mq = read_low_mq_count;
        summary.mapping = read_idx;
        summary.complex = read_complex_count;
    
        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(&output).with_extension("").to_str().unwrap().to_string()
        };

        summary.save(&format!("{}.read.summary", output_prefix));

        log::info!("Successful output Pore-C table `{}`", output);
        Ok(())
    }

    
    pub fn paf2table3(&self, bed: &String,
            output: &String, min_quality: &u8, 
            min_identity: &f32,  min_length: &u32,
            max_order: &u32, 
            max_edge_length: &u64,
            secondary: bool    
        ) -> Result<(), Box<dyn Error>> {
        
        type IvU8 = Interval<usize, u8>;
        let min_quality = *min_quality;
        let min_identity = *min_identity;
        let min_length = *min_length;
        let max_order = *max_order;
        let max_edge_length = *max_edge_length;

        let is_filter_digest = !bed.is_empty();
        let bed = if is_filter_digest { Bed3::new(bed) } else { Bed3::new(&String::from(".tmp.bed")) };
        let interval_hash = bed.to_interval_hash();
        let is_filter_edges = max_edge_length > 0;

        let num_threads = 10;
        log::info!("Using {} threads for processing.", num_threads);
        let (sender, receiver) = bounded::<Vec<PoreCRecord>>(8192); // Channel for read groups
        let (writer_sender, writer_receiver) = bounded::<String>(8192); // Channel for output strings

        let mut handles = Vec::with_capacity(num_threads);
        // let prof = std::sync::Arc::new(Prof::new());

        for _ in 0..num_threads {
            let receiver = receiver.clone();
            let writer_sender = writer_sender.clone();
            let interval_hash = interval_hash.clone(); // Clone for thread
            // let prof_cl = prof.clone();

            handles.push(thread::spawn(move || {
                let init_cap = (max_order as usize) + 2;
                let mut concatemer: Concatemer = Concatemer::with_capacity(init_cap);
                let mut local_stats = (0, 0, 0, 0); // pass, singleton, low_mq, complex
                let mut output_buffer = String::with_capacity(WRITE_CHUNK);
                while let Ok(records) = receiver.recv() {
                    for pcr in records {
                        concatemer.push(pcr);
                    }
                    // let t0 = Instant::now();
                    let (final_status, pass_indices) = concatemer.process(
                        is_filter_digest,
                        &interval_hash,
                        is_filter_edges,
                        max_edge_length,
                        max_order,
                    );
                    // if prof_cl.on { Prof::add(&prof_cl.worker_proc_ns, t0.elapsed().as_nanos()); }
                    // let t1 = Instant::now();
                    // for pcr in concatemer.records.iter() {
                    //     if pcr.filter_reason == PASS_STR {
                    //         output_buffer.push_str(&pcr.to_string());
                    //         output_buffer.push('\n');
                    //     }
                    // }
                    if final_status == PASS_STR {
                        for &idx in &pass_indices {
                            let pcr = &concatemer.records[idx];
                            output_buffer.push_str(&pcr.to_string());
                            output_buffer.push('\n');
                        }
                    }
                    // if prof_cl.on { Prof::add(&prof_cl.worker_fmt_ns, t1.elapsed().as_nanos()); }

                    if output_buffer.len() >= WRITE_CHUNK {
                        let chunk = std::mem::take(&mut output_buffer);
                        writer_sender.send(chunk).unwrap();
                        output_buffer = String::with_capacity(WRITE_CHUNK);
                    }

                    match final_status {
                        PASS_STR => local_stats.0 += 1,
                        SINGLETON_STR => local_stats.1 += 1,
                        LOW_MQ_STR => local_stats.2 += 1,
                        COMPLEX_STR => local_stats.3 += 1,
                        _ => {},
                    }
                    concatemer.clear();
                }

                if !output_buffer.is_empty() {
                    writer_sender.send(output_buffer).unwrap();
                }
                local_stats
            }));
        }
        drop(writer_sender); // Drop the original sender for the writer

        // Writer thread
        let output_clone = output.clone();
        // let prof_w = prof.clone();
        let writer_handle = thread::spawn(move || {
            let mut writer = common_writer(&output_clone);
            while let Ok(data_string) = writer_receiver.recv() {
                // let t0 = Instant::now();
                writer.write_all(data_string.as_bytes()).unwrap();
                // if prof_w.on {
                //     Prof::add(&prof_w.writer_ns, t0.elapsed().as_nanos());
                //     prof_w.out_bytes.fetch_add(data_string.len() as u64, Ordering::Relaxed);
                // }
            }
        });

        let mut rdr = self.parse2(secondary.clone())?;
        let mut old_query: Option<String> = None;
        let mut current_read_records: Vec<PoreCRecord> = Vec::new();
        let mut read_idx: u64 = 0;

        for record in rdr {
            // let query_name = record.query.clone();

            // if old_query.as_deref() != Some(&query_name) && old_query.is_some() {
            //     sender.send(current_read_records).unwrap();
            //     current_read_records = Vec::new();
            //     read_idx += 1;
            // }
            let is_switch = match old_query.as_deref() {
                Some(prev) => prev != record.query.as_str(),
                None => false,
            };
            if is_switch && old_query.is_some() {
                sender.send(current_read_records).unwrap();
                current_read_records = Vec::new();
                read_idx += 1;
            }
            let length: u32 = record.query_end - record.query_start;
            let identity = if record.alignment_length > 0 { record.match_n as f32 / record.alignment_length as f32 } else { 0.0 };

            let fr = if record.mapq < min_quality || length < min_length || identity < min_identity {
                LOW_MQ_STR
            } else {
                PASS_STR
            };

            // old_query = Some(query_name);
            if old_query.is_none() {
                old_query = Some(record.query.clone());
            } else if is_switch {
                old_query = Some(record.query.clone());
            }

            let pcr = PoreCRecord::from_paf_record(record, read_idx, identity, fr.to_string());
                current_read_records.push(pcr);
        }

        // Send the last group
        if !current_read_records.is_empty() {
            sender.send(current_read_records).unwrap();
            read_idx += 1;
        }
        drop(sender); 
        let mut read_pass_count: u64 = 0;
        let mut read_singleton_count: u64 = 0;
        let mut read_low_mq_count: u64 = 0;
        let mut read_complex_count: u64 = 0;

        for handle in handles {
        let stats = handle.join().unwrap();
            read_pass_count += stats.0;
            read_singleton_count += stats.1;
            read_low_mq_count += stats.2;
            read_complex_count += stats.3;
        }
        writer_handle.join().unwrap();
        // if prof.on {
        //     let rd = prof.read_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let wp = prof.worker_proc_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let wf = prof.worker_fmt_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let ww = prof.writer_ns.load(Ordering::Relaxed) as f64 / 1e9;
        //     let out_mb = prof.out_bytes.load(Ordering::Relaxed) as f64 / (1024.0*1024.0);
        //     log::info!("PROFILE read+group: {:.3}s, worker process: {:.3}s, worker format: {:.3}s, write(compress): {:.3}s, wrote {:.1} MiB",
        //                rd, wp, wf, ww, out_mb);
        // }
        let mut summary: ReadSummary = ReadSummary::new();
        summary.pass = read_pass_count;
        summary.singleton = read_singleton_count;
        summary.low_mq = read_low_mq_count;
        summary.mapping = read_idx;
        summary.complex = read_complex_count;
       
        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(&output).with_extension("").to_str().unwrap().to_string()
        };

        summary.save(&format!("{}.read.summary", output_prefix));

        log::info!("Successful output Pore-C table `{}`", output);
        Ok(())
    }


    pub fn paf2table2(&self, bed: &String,
                     output: &String, min_quality: &u8, 
                     min_identity: &f32,  min_length: &u32,
                     max_order: &u32, 
                     max_edge_length: &u64,
                     secondary: bool    
                    ) -> Result<(), Box<dyn Error>> {
        

        type IvU8 = Interval<usize, u8>;
        let min_quality = *min_quality;
        let min_identity = *min_identity;
        let min_length = *min_length;
        let max_order = *max_order;
        let max_edge_length = *max_edge_length;
    
        let is_filter_digest = if bed != "" { true } else { false };
       
        let bed = if is_filter_digest {
            Bed3::new(bed)
        } else {
            Bed3::new(&String::from(".tmp.bed"))
        };
        
        let interval_hash = bed.to_interval_hash();

        let parse_result = self.parse2(secondary.clone());
        
        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };
        
        let mut writer = common_writer(output);
    
        
        let mut read_idx: u64 = 0;
        
        let mut filter_reason: &str = "";
        let mut old_query: Option<String> = None;
        let init_cap = (max_order as usize) + 2;
        let mut concatemer: Concatemer = Concatemer::with_capacity(init_cap);
        
        let mut read_mapping_count: u64 = 0;
        let mut read_pass_count: u64 = 0;
        let mut read_singleton_count: u64 = 0;
        let mut read_low_mq_count: u64 = 0;
        let mut read_complex_count: u64 = 0;
        let is_filter_edges = max_edge_length > 0;
        for (line_num, record) in rdr.enumerate() {
            // parse error to continue
            
            if old_query.as_deref() != Some(&record.query) && old_query.is_some() {
                read_idx += 1;  
                
                if is_filter_digest {
                    concatemer.filter_digest(&interval_hash);
                }

                if is_filter_edges {
                    concatemer.filter_edges(max_edge_length);
                }

                if concatemer.is_singleton() {
                    concatemer.parse_singleton();
                }
                
                if concatemer.is_complex(max_order) {
                    concatemer.parse_complex();
                }

                for pcr in concatemer.records.iter() {
                    if pcr.filter_reason == PASS_STR {
                       writer.write_all(pcr.to_string().as_bytes()).unwrap();
                    }
                }
                
                match concatemer.info() {
                    "pass" => read_pass_count += 1,
                    "singleton" => read_singleton_count += 1,
                    "low_mq" => read_low_mq_count += 1,
                    "complex" => read_complex_count += 1,
                    _ => todo!(),
                }
                concatemer.clear();
                
            }
            let length: u32 = record.query_end - record.query_start;
            let identity = record.match_n as f32 / record.alignment_length as f32;
            
            let fr = if record.mapq < min_quality
                 || length   < min_length
                 || identity < min_identity
            {
                LOW_MQ_STR
            } else {
                PASS_STR
            };

            if old_query.as_deref() != Some(&record.query) {
                old_query = Some(record.query.clone());
            }
            let pcr = PoreCRecord::from_paf_record(record, read_idx, identity, fr.to_string());
            
            concatemer.push(pcr);
        }
        
        // parse last record
        if concatemer.is_singleton() {
            concatemer.parse_singleton();
        }
        for pcr in concatemer.records.iter() {
            if pcr.filter_reason == PASS_STR {
                writer.write_all(pcr.to_string().as_bytes()).unwrap();
            }
        }
        match concatemer.info() {
            "pass" => read_pass_count += 1,
            "singleton" => read_singleton_count += 1,
            "low_mq" => read_low_mq_count += 1,
            "complex" => read_complex_count += 1,
            _ => todo!(),
        }
        read_mapping_count = read_idx + 1;


        let mut summary: ReadSummary = ReadSummary::new();
        summary.pass = read_pass_count;
        summary.singleton = read_singleton_count;
        summary.low_mq = read_low_mq_count;
        summary.mapping = read_mapping_count;
        summary.complex = read_complex_count;

        
        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(&output).with_extension("").to_str().unwrap().to_string()
        };

     
        summary.save(&format!("{}.read.summary", output_prefix));

        log::info!("Successful output Pore-C table `{}`", output);
        Ok(())
    }

    // pub fn to_depth(&self, chromsize: &String, window_size: usize, step_size: usize,
    //             min_mapq: u8, secondary: bool, output: &String) -> Result<(), Box<dyn Error>> {
    //     let mut chromsize_map: HashMap<String, usize> = HashMap::new();
    //     {
    //         let input = common_reader(chromsize);
    //         let buf = std::io::BufReader::new(input);
    //         for line in buf.lines() {
    //             let line = line?;
    //             if line.trim().is_empty() || line.starts_with('#') {
    //                 continue;
    //             }

    //             let mut spl = line.split_whitespace();
    //             let chr = spl.next().unwrap_or("").to_string();
    //             let size_str = spl.next().unwrap_or("0");
    //             let size = size_str.parse::<usize>().unwrap_or(0);
    //             chromsize_map.insert(chr, size);
    //         }
    //     }
        
    //     let mut coverage_diff_map: HashMap<String, BTreeMap<usize, i64>> = HashMap::new();
    //     for (chr, _) in chromsize_map.iter() {
    //         coverage_diff_map.insert(chr.clone(), BTreeMap::new());
    //     }

    //     let mut rdr = self.parse2(secondary)?;
    //     let (sender, receiver) = bounded::<Vec<PAFLine>>(100);
        
    //     let mut handles = Vec::new();
    //     let mut coverage_diff_map = Arc::new(Mutex::new(coverage_diff_map));
    //     for _ in 0..8 {
    //         let receiver = receiver.clone();
    //         let coverage_diff_map_clone = Arc::clone(&coverage_diff_map);
    //         handles.push(thread::spawn(move || {
    //             while let Ok(records) = receiver.recv() {
    //                 let mut local_data = HashMap::new();
    //                 for record in records {
    //                     let mapq = record.mapq;
    //                     if mapq < min_mapq {
    //                         continue;
    //                     }
    //                     let chr = record.target.clone();
    //                     let start = record.target_start as usize;
    //                     let end = record.target_end as usize;

    //                     let local_btm = local_data.entry(chr).or_insert_with(BTreeMap::new);
    //                     *local_btm.entry(start).or_insert(0) += 1;
    //                     *local_btm.entry(end).or_insert(0) -= 1;
    //                 }

    //                 let mut coverage_diff_map_locked = coverage_diff_map_clone.lock().unwrap();
    //                 for (chr, local_btm) in local_data {
    //                     let global_btm = coverage_diff_map_locked
    //                         .entry(chr)
    //                         .or_insert_with(BTreeMap::new);
    //                     for (pos, delta) in local_btm {
    //                         *global_btm.entry(pos).or_insert(0) += delta;
    //                     }
    //                 }

    //             }
    //         }))


    //     }
        
    //     let mut batch = Vec::with_capacity(CHUNK_SIZE);
    //     for record in rdr{
    //         batch.push(record);
    //         if batch.len() >= CHUNK_SIZE {
    //             sender.send(batch).unwrap();
    //             batch = Vec::with_capacity(CHUNK_SIZE);
    //         }
    //     }

    //     if !batch.is_empty() {
    //         sender.send(batch).unwrap();
    //     }

    //     drop(sender);

    //     for handle in handles {
    //         handle.join().unwrap();
    //     }

    //     let mut coverage_diff_map = Arc::try_unwrap(coverage_diff_map).unwrap().into_inner().unwrap();
    //     if step_size == window_size {
    //         let mut wtr = common_writer(output);
    //         for (chr, chrom_length) in &chromsize_map {
    //             if let Some(diff) = coverage_diff_map.get_mut(chr) {
    //                 let mut running = 0i64;
    //                 let mut prev_pos = 0usize;
                    
    //                 let mut bin_start = 0usize;
    //                 let mut bin_sum = 0i64;
    //                 let mut bin_count = 0usize;

    //                 for (&pos, &delta) in diff.iter() {
    //                     while prev_pos < pos && prev_pos < *chrom_length {
    //                         bin_sum += running;
    //                         bin_count += 1;
    //                         prev_pos += 1;

    //                         if prev_pos - bin_start == window_size {
    //                             let mean_coverage = bin_sum as f64 / bin_count as f64;
    //                             let line = format!("{}\t{}\t{}\t{:.3}\n", chr, bin_start, prev_pos, mean_coverage);
    //                             wtr.write_all(line.as_bytes())?;

    //                             bin_sum = 0;
    //                             bin_count = 0;
    //                             bin_start = prev_pos;
    //                         }
    //                     }
        
    //                     running += delta;
    //                 }

    //                 while prev_pos < *chrom_length {
    //                     bin_sum += running;
    //                     bin_count += 1;
    //                     prev_pos += 1;

    //                     if prev_pos - bin_start == window_size {
    //                         let mean_coverage = bin_sum as f64 / bin_count as f64;
    //                         let line = format!("{}\t{}\t{}\t{:.3}\n", chr, bin_start, prev_pos, mean_coverage);
    //                         wtr.write_all(line.as_bytes())?;
    //                         bin_sum = 0;
    //                         bin_count = 0;
    //                         bin_start = prev_pos;
    //                     }
    //                 }

    //                 if bin_count > 0 {
    //                     let mean_coverage = bin_sum as f64 / bin_count as f64;
    //                     let line = format!("{}\t{}\t{}\t{:.3}\n", chr, bin_start, prev_pos, mean_coverage);
    //                     wtr.write_all(line.as_bytes())?;
    //                 }
                    
    //             }
    //         }

    //         log::info!(
    //             "Successful output coverage (window {}bp) to `{}`",
    //             window_size,
    //             output
    //         );

    //     } else {

    //         let mut wtr = common_writer(output);

    //         for (chr, &chrom_len) in &chromsize_map {
    //             let diff = match coverage_diff_map.get_mut(chr) {
    //                 Some(d) => d,
    //                 None => continue,
    //             };

    //             let mut coverage_vec = vec![0i64; chrom_len];
    //             let mut running = 0i64;
    //             let mut prev_pos = 0usize;
    //             for (&pos, &delta) in diff.iter() {
    //                 while prev_pos < pos && prev_pos < chrom_len {
    //                     coverage_vec[prev_pos] = running;
    //                     prev_pos += 1;
    //                 }
    //                 running += delta;
    //             }
    
    //             while prev_pos < chrom_len {
    //                 coverage_vec[prev_pos] = running;
    //                 prev_pos += 1;
    //             }

    //             let mut prefix_sum = vec![0i64; chrom_len + 1];
    //             for i in 0..chrom_len {
    //                 prefix_sum[i+1] = prefix_sum[i] + coverage_vec[i];
    //             }

    //             let mut i = 0usize;
    //             while i < chrom_len {
    //                 let w_end = std::cmp::min(i + window_size, chrom_len);
    //                 let total_cov = prefix_sum[w_end] - prefix_sum[i];
    //                 let window_len = w_end - i;
    //                 if window_len == 0 {
    //                     break;
    //                 }
    //                 let mean_cov = total_cov as f64 / window_len as f64;

    //                 let line = format!("{}\t{}\t{}\t{:.3}\n", chr, i, w_end, mean_cov);
    //                 wtr.write_all(line.as_bytes())?;

    //                 i += step_size; 
    //             }
    //         }

    //         log::info!(
    //             "Successful output coverage (window {}bp, step {}bp) to `{}`",
    //             window_size,
    //             step_size,
    //             output
    //         );
    //     }

    //     Ok(())
    // }

    // pub fn to_depth(&self, chromsize: &String, window_size: usize, step_size: usize,
    //             min_mapq: u8, secondary: bool, output: &String) -> Result<(), Box<dyn Error>> {
    //     let mut chrom_names = Vec::new();
    //     let mut chrom_sizes = Vec::new();
    //     let mut chrom_map = HashMap::new();

    //     {
    //         let input = common_reader(chromsize);
    //         let buf = std::io::BufReader::new(input);
    //         for line in buf.lines() {
    //             let line = line?;
    //             if line.trim().is_empty() || line.starts_with('#') {
    //                 continue;
    //             }

    //             let mut spl = line.split_whitespace();
    //             let chr = spl.next().unwrap_or("").to_string();
    //             let size_str = spl.next().unwrap_or("0");
    //             let size = size_str.parse::<usize>().unwrap_or(0);
    //             chrom_map.insert(chr.clone(), chrom_names.len());
    //             chrom_names.push(chr);
    //             chrom_sizes.push(size);
    //         }
    //     }
        
    //     let num_chroms = chrom_names.len();
    //     let chrom_map = Arc::new(chrom_map);

    //     let (sender, receiver) = bounded::<Vec<PAFLine>>(100);
        
    //     let mut handles = Vec::new();
    //     for _ in 0..8 {
    //         let receiver = receiver.clone();
    //         let chrom_map = chrom_map.clone();
    //         handles.push(thread::spawn(move || {
               
    //             let mut local_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];
    //             while let Ok(records) = receiver.recv() {
    //                 for record in records {
    //                     if record.mapq < min_mapq {
    //                         continue;
    //                     }
    //                     if let Some(&idx) = chrom_map.get(&record.target) {
    //                         local_events[idx].push((record.target_start as u32, 1));
    //                         local_events[idx].push((record.target_end as u32, -1));
    //                     }
    //                 }
    //             }
    //             local_events
    //         }));
    //     }
        
    //     let rdr = self.parse2(secondary)?;
    //     let mut batch = Vec::with_capacity(CHUNK_SIZE);
    //     for record in rdr{
    //         batch.push(record);
    //         if batch.len() >= CHUNK_SIZE {
    //             sender.send(batch).unwrap();
    //             batch = Vec::with_capacity(CHUNK_SIZE);
    //         }
    //     }

    //     if !batch.is_empty() {
    //         sender.send(batch).unwrap();
    //     }

    //     drop(sender);

    //     let mut global_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];
    //     for handle in handles {
    //         let local_events = handle.join().unwrap();
    //         for (i, events) in local_events.into_iter().enumerate() {
    //             global_events[i].extend(events);
    //         }
    //     }

    //     let mut wtr = common_writer(output);

    //     for (i, mut events) in global_events.into_iter().enumerate() {
    //         if events.is_empty() { continue; }
    //         let chrom_len = chrom_sizes[i];
    //         let chrom_name = &chrom_names[i];

    //         events.par_sort_unstable_by_key(|e| e.0);

    //         let mut coverage = vec![0i32; chrom_len];
    //         let mut current_cov = 0i32;
    //         let mut prev_pos = 0;

    //         for (pos, delta) in events {
    //             let pos = pos as usize;
    //             if pos >= chrom_len { break; }
                
    //             if pos > prev_pos {
    //                 for k in prev_pos..pos {
    //                     coverage[k] = current_cov;
    //                 }
    //             }
    //             current_cov += delta;
    //             prev_pos = pos;
    //         }
    //         if prev_pos < chrom_len {
    //             for k in prev_pos..chrom_len {
    //                 coverage[k] = current_cov;
    //             }
    //         }


    //         let mut prefix_sum = vec![0i64; chrom_len + 1];
    //         let mut running = 0i64;
    //         for (k, &cov) in coverage.iter().enumerate() {
    //             running += cov as i64;
    //             prefix_sum[k+1] = running;
    //         }

    //         let mut start = 0usize;
    //         while start < chrom_len {
    //             let end = std::cmp::min(start + window_size, chrom_len);
    //             let len = end - start;
    //             if len == 0 { break; }
                
    //             let sum = prefix_sum[end] - prefix_sum[start];
    //             let mean = sum as f64 / len as f64;
                
    //             let line = format!("{}\t{}\t{}\t{:.3}\n", chrom_name, start, end, mean);
    //             wtr.write_all(line.as_bytes())?;

    //             start += step_size;
    //         }
    //     }

    //     log::info!(
    //         "Successful output coverage (window {}bp, step {}bp) to `{}`",
    //         window_size,
    //         step_size,
    //         output
    //     );

    //     Ok(())
    // }

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

                // let mut prefix_sum = Vec::with_capacity(chrom_len + 1);
                // prefix_sum.push(0);
                
                // let mut current_cov = 0i32;
                // let mut prev_pos = 0usize;
                // let mut running_sum = 0i64;

                // for (pos, delta) in events {
                //     let pos = pos as usize;
                //     if pos >= chrom_len { break; }
                    
                //     if pos > prev_pos {
                //         let cov = current_cov as i64;
                //         let count = pos - prev_pos;
                //         for _ in 0..count {
                //             running_sum += cov;
                //             prefix_sum.push(running_sum);
                //         }
                //     }
                //     current_cov += delta;
                //     prev_pos = pos;
                // }
                
                // if prev_pos < chrom_len {
                //     let cov = current_cov as i64;
                //     let count = chrom_len - prev_pos;
                //     for _ in 0..count {
                //         running_sum += cov;
                //         prefix_sum.push(running_sum);
                //     }
                // }

                // let mut output_buf = String::with_capacity(64 * 1024);
                // let mut start = 0usize;

                // while start < chrom_len {
                //     let end = std::cmp::min(start + window_size, chrom_len);
                //     let len = end - start;
                //     if len == 0 { break; }
                    
                //     let sum = prefix_sum[end] - prefix_sum[start];
                //     let mean = sum as f64 / len as f64;
                    
                //     write!(output_buf, "{}\t{}\t{}\t{:.3}\n", chrom_name, start, end, mean).unwrap();
                let mut coverage = vec![0i32; chrom_len];
                let mut current_cov = 0i32;
                let mut prev_pos = 0usize;

                for (pos, delta) in events {
                    let pos = pos as usize;
                    if pos >= chrom_len { break; }
                    
                    if pos > prev_pos {
                        for k in prev_pos..pos {
                            coverage[k] = current_cov;
                        }
                    }
                    current_cov += delta;
                    prev_pos = pos;
                }
                
                if prev_pos < chrom_len {
                    for k in prev_pos..chrom_len {
                        coverage[k] = current_cov;
                    }
                }

                let mut output_buf = String::with_capacity(64 * 1024);
                let mut start = 0usize;
                
                let mut window_buf = Vec::with_capacity(window_size);

                while start < chrom_len {
                    let end = std::cmp::min(start + window_size, chrom_len);
                    let len = end - start;
                    if len == 0 { break; }
                    
    
                    window_buf.clear();
                    window_buf.extend_from_slice(&coverage[start..end]);
                    
                    let mid = len / 2;
                    let (_, &mut median_val, _) = window_buf.select_nth_unstable(mid);
                    
                    let median = if len % 2 == 0 {
                        let max_lower = *window_buf[..mid].iter().max().unwrap();
                        (median_val as f64 + max_lower as f64) / 2.0
                    } else {
                        median_val as f64
                    };
                    
                    write!(output_buf, "{}\t{}\t{}\t{:.3}\n", chrom_name, start, end, median).unwrap();
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
}


   