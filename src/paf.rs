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

    pub fn paf2table4(&self, bed: &String,
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
        rayon::ThreadPoolBuilder::new().num_threads(num_threads).build_global().unwrap();
        let use_rayon = true;
        let mut rdr = self.parse2(secondary.clone())?;
        let mut old_query: Option<String> = None;
        let mut current_group: Vec<PoreCRecord> = Vec::new();
        let mut read_idx: u64 = 0;

        let mut batch_groups: Vec<Vec<PoreCRecord>> = Vec::with_capacity(GROUP_BATCH_SIZE);
        let mut total_records_in_batch = 0usize;

  
        let mut writer = common_writer(output);
        let mut stats = (0u64,0u64,0u64,0u64); // pass,singleton,low_mq,complex

        let flush_batch = |mut groups: Vec<Vec<PoreCRecord>>,
                           is_filter_digest: bool,
                           interval_hash: &HashMap<String,Lapper<usize,u8>>,
                           is_filter_edges: bool,
                           max_edge_length: u64,
                           max_order: u32,
                           stats: &mut (u64,u64,u64,u64),
                           writer: &mut dyn Write| {
            if groups.is_empty() { return; }

            let results: Vec<(String,(u64,u64,u64,u64))> = groups
                .into_par_iter()
                .map(|grp| {
                    let mut c = Concatemer::with_capacity(grp.len()+2);
                    for r in grp { c.push(r); }
                    let (final_status, pass_indices) = c.process(
                        is_filter_digest, interval_hash, is_filter_edges,
                        max_edge_length, max_order
                    );
                    let mut out = String::new();
                    if final_status == PASS_STR {
                        out.reserve(pass_indices.len()*128);
                        for &idx in &pass_indices {
                            out.push_str(&c.records[idx].to_string());
                            out.push('\n');
                        }
                    }
                    let st = match final_status {
                        PASS_STR => (1,0,0,0),
                        SINGLETON_STR => (0,1,0,0),
                        LOW_MQ_STR => (0,0,1,0),
                        COMPLEX_STR => (0,0,0,1),
                        _ => (0,0,0,0)
                    };
                    (out, st)
                }).collect();
            let mut merged = String::with_capacity(results.iter().map(|(s,_ )| s.len()).sum());
            for (chunk,(a,b,c,d)) in results {
                if !chunk.is_empty() {
                    merged.push_str(&chunk);
                }
                stats.0 += a; stats.1 += b; stats.2 += c; stats.3 += d;
            }
            if !merged.is_empty() {
                writer.write_all(merged.as_bytes()).unwrap();
            }
        };

        for record in rdr {
            let is_switch = match old_query.as_deref() {
                Some(prev) => prev != record.query.as_str(),
                None => false,
            };
            if is_switch && old_query.is_some() {
                batch_groups.push(std::mem::take(&mut current_group));
                read_idx += 1;
                if batch_groups.len() >= GROUP_BATCH_SIZE ||
                   total_records_in_batch >= RECORD_BATCH_SIZE {
                    flush_batch(std::mem::take(&mut batch_groups),
                                is_filter_digest,&interval_hash,
                                is_filter_edges,max_edge_length,max_order,
                                &mut stats,&mut writer);
                    total_records_in_batch = 0;
                }
            }
            let length = record.query_end - record.query_start;
            let identity = if record.alignment_length>0 {
                record.match_n as f32 / record.alignment_length as f32
            } else { 0.0 };
            let fr = if record.mapq < min_quality || length < min_length || identity < min_identity {
                LOW_MQ_STR
            } else { PASS_STR };

            if is_switch || old_query.is_none() {
                old_query = Some(record.query.clone());
            }
            let pcr = PoreCRecord::from_paf_record(record, read_idx, identity, fr.to_string());
            current_group.push(pcr);
            total_records_in_batch += 1;
        }
        if !current_group.is_empty() {
            batch_groups.push(std::mem::take(&mut current_group));
            read_idx += 1;
        }
        if !batch_groups.is_empty() {
            flush_batch(std::mem::take(&mut batch_groups),
                        is_filter_digest,&interval_hash,
                        is_filter_edges,max_edge_length,max_order,
                        &mut stats,&mut writer);
        }

        let mut summary = ReadSummary::new();
        summary.pass = stats.0;
        summary.singleton = stats.1;
        summary.low_mq = stats.2;
        summary.complex = stats.3;
        summary.mapping = read_idx;
        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(output).with_extension("").to_str().unwrap().to_string()
        };
        summary.save(&format!("{}.read.summary", output_prefix));
        log::info!("Successful output Pore-C table `{}` (rayon)", output);
        return Ok(());
    }


    pub fn paf2table(&self, bed: &String,
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

    pub fn to_depth(&self, chromsize: &String, window_size: usize, step_size: usize,
                min_mapq: u8, secondary: bool, output: &String) -> Result<(), Box<dyn Error>> {
        let mut chromsize_map: HashMap<String, usize> = HashMap::new();
        {
            let input = common_reader(chromsize);
            let buf = std::io::BufReader::new(input);
            for line in buf.lines() {
                let line = line?;
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }

                let mut spl = line.split_whitespace();
                let chr = spl.next().unwrap_or("").to_string();
                let size_str = spl.next().unwrap_or("0");
                let size = size_str.parse::<usize>().unwrap_or(0);
                chromsize_map.insert(chr, size);
            }
        }
        
        let mut coverage_diff_map: HashMap<String, BTreeMap<usize, i64>> = HashMap::new();
        for (chr, _) in chromsize_map.iter() {
            coverage_diff_map.insert(chr.clone(), BTreeMap::new());
        }

        let mut rdr = self.parse2(secondary)?;
        let (sender, receiver) = bounded::<Vec<PAFLine>>(100);
        
        let mut handles = Vec::new();
        let mut coverage_diff_map = Arc::new(Mutex::new(coverage_diff_map));
        for _ in 0..8 {
            let receiver = receiver.clone();
            let coverage_diff_map_clone = Arc::clone(&coverage_diff_map);
            handles.push(thread::spawn(move || {
                while let Ok(records) = receiver.recv() {
                    let mut local_data = HashMap::new();
                    for record in records {
                        let mapq = record.mapq;
                        if mapq < min_mapq {
                            continue;
                        }
                        let chr = record.target.clone();
                        let start = record.target_start as usize;
                        let end = record.target_end as usize;

                        let local_btm = local_data.entry(chr).or_insert_with(BTreeMap::new);
                        *local_btm.entry(start).or_insert(0) += 1;
                        *local_btm.entry(end).or_insert(0) -= 1;
                    }

                    let mut coverage_diff_map_locked = coverage_diff_map_clone.lock().unwrap();
                    for (chr, local_btm) in local_data {
                        let global_btm = coverage_diff_map_locked
                            .entry(chr)
                            .or_insert_with(BTreeMap::new);
                        for (pos, delta) in local_btm {
                            *global_btm.entry(pos).or_insert(0) += delta;
                        }
                    }

                }
            }))


        }
        
        let mut batch = Vec::with_capacity(CHUNK_SIZE);
        for record in rdr{
            batch.push(record);
            if batch.len() >= CHUNK_SIZE {
                sender.send(batch).unwrap();
                batch = Vec::with_capacity(CHUNK_SIZE);
            }
        }

        if !batch.is_empty() {
            sender.send(batch).unwrap();
        }

        drop(sender);

        for handle in handles {
            handle.join().unwrap();
        }


        // let mut rdr = self.parse()?;
        // for (line_num, rec) in rdr.deserialize().enumerate() {
        //     let record: PAFLine = match rec {
        //         Ok(r) => r,
        //         Err(e) => {
        //             log::warn!("Error parsing line {} in {}: {:?}", line_num, self.file_name(), e);
        //             continue;
        //         }
        //     };
        //     let mapq = record.mapq;
        //     if mapq < min_mapq {
        //         continue;
        //     }
        //     let chr = record.target.clone();
        //     let start = record.target_start as usize;
        //     let end = record.target_end as usize;

        //     if let Some(diff) = coverage_diff_map.get_mut(&chr) {
        //         *diff.entry(start).or_insert(0) += 1;
        //         *diff.entry(end).or_insert(0) -= 1;
        //     }
        // }

        let mut coverage_diff_map = Arc::try_unwrap(coverage_diff_map).unwrap().into_inner().unwrap();
        if step_size == window_size {
            let mut wtr = common_writer(output);
            for (chr, chrom_length) in &chromsize_map {
                if let Some(diff) = coverage_diff_map.get_mut(chr) {
                    let mut running = 0i64;
                    let mut prev_pos = 0usize;
                    
                    let mut bin_start = 0usize;
                    let mut bin_sum = 0i64;
                    let mut bin_count = 0usize;

                    for (&pos, &delta) in diff.iter() {
                        while prev_pos < pos && prev_pos < *chrom_length {
                            bin_sum += running;
                            bin_count += 1;
                            prev_pos += 1;

                            if prev_pos - bin_start == window_size {
                                let mean_coverage = bin_sum as f64 / bin_count as f64;
                                let line = format!("{}\t{}\t{}\t{:.3}\n", chr, bin_start, prev_pos, mean_coverage);
                                wtr.write_all(line.as_bytes())?;

                                bin_sum = 0;
                                bin_count = 0;
                                bin_start = prev_pos;
                            }
                        }
        
                        running += delta;
                    }

                    while prev_pos < *chrom_length {
                        bin_sum += running;
                        bin_count += 1;
                        prev_pos += 1;

                        if prev_pos - bin_start == window_size {
                            let mean_coverage = bin_sum as f64 / bin_count as f64;
                            let line = format!("{}\t{}\t{}\t{:.3}\n", chr, bin_start, prev_pos, mean_coverage);
                            wtr.write_all(line.as_bytes())?;
                            bin_sum = 0;
                            bin_count = 0;
                            bin_start = prev_pos;
                        }
                    }

                    if bin_count > 0 {
                        let mean_coverage = bin_sum as f64 / bin_count as f64;
                        let line = format!("{}\t{}\t{}\t{:.3}\n", chr, bin_start, prev_pos, mean_coverage);
                        wtr.write_all(line.as_bytes())?;
                    }
                    
                }
            }

            log::info!(
                "Successful output coverage (window {}bp) to `{}`",
                window_size,
                output
            );

        } else {

            let mut wtr = common_writer(output);

            for (chr, &chrom_len) in &chromsize_map {
                let diff = match coverage_diff_map.get_mut(chr) {
                    Some(d) => d,
                    None => continue,
                };

                let mut coverage_vec = vec![0i64; chrom_len];
                let mut running = 0i64;
                let mut prev_pos = 0usize;
                for (&pos, &delta) in diff.iter() {
                    while prev_pos < pos && prev_pos < chrom_len {
                        coverage_vec[prev_pos] = running;
                        prev_pos += 1;
                    }
                    running += delta;
                }
    
                while prev_pos < chrom_len {
                    coverage_vec[prev_pos] = running;
                    prev_pos += 1;
                }

                let mut prefix_sum = vec![0i64; chrom_len + 1];
                for i in 0..chrom_len {
                    prefix_sum[i+1] = prefix_sum[i] + coverage_vec[i];
                }

                let mut i = 0usize;
                while i < chrom_len {
                    let w_end = std::cmp::min(i + window_size, chrom_len);
                    let total_cov = prefix_sum[w_end] - prefix_sum[i];
                    let window_len = w_end - i;
                    if window_len == 0 {
                        break;
                    }
                    let mean_cov = total_cov as f64 / window_len as f64;

                    let line = format!("{}\t{}\t{}\t{:.3}\n", chr, i, w_end, mean_cov);
                    wtr.write_all(line.as_bytes())?;

                    i += step_size; 
                }
            }

            log::info!(
                "Successful output coverage (window {}bp, step {}bp) to `{}`",
                window_size,
                step_size,
                output
            );
        }

        Ok(())
    }

    
}


   