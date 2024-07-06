#[warn(unused_assignments)]
use anyhow::Result as anyResult;
use lazy_static::lazy_static;
use std::any::type_name;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::ffi::OsStr;
use std::path::Path;
use std::error::Error;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rust_lapper::{Interval, Lapper};

use crate::bed::Bed3; 
use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;
use crate::porec::PoreCRecordPlus as PoreCRecord;


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

#[derive(Debug)]
pub struct Concatemer {
    pub records: Vec<PoreCRecord>,
    pub filter_reasons: Vec<String>,
}

impl Concatemer {

    pub fn new() -> Concatemer{
        Concatemer {
            records: Vec::new(),
            filter_reasons: Vec::new()
        }
    }

    pub fn push(&mut self, pcr: PoreCRecord) {
        self.filter_reasons.push(pcr.filter_reason.clone());
        self.records.push(pcr);
        
    }

    pub fn clear(&mut self) {
        self.records.clear();
        self.filter_reasons.clear();
    }

    pub fn is_singleton(&self) -> bool {
        self.filter_reasons.iter().filter(|&x| *x == "pass").count() == 1
    }

    pub fn parse_singleton(&mut self) {
        // for (i, record) in self.records.iter_mut().enumerate() {
        //     if record.filter_reason == "pass" {
        //         record.filter_reason = String::from("singleton");
        //         self.filter_reasons[i] = String::from("singleton");
        //     }
        // }

        if let Some(i) = self.records.iter().position(|r| r.filter_reason == "pass") {
            self.records[i].filter_reason = String::from("singleton");
            self.filter_reasons[i] = String::from("singleton");
        }
    }

    pub fn is_complex(&self, max_order: &u32 ) -> bool {
        let pass_count = self.filter_reasons.iter().filter(|&x| *x == "pass").count();
        if pass_count > *max_order as usize {
            true
        } else {
            false
        }
    }

    pub fn parse_complex(&mut self) {
        for (i, record) in self.records.iter_mut().enumerate() {
            if record.filter_reason == "pass" {
                record.filter_reason = String::from("complex");
                self.filter_reasons[i] = String::from("complex");
            }
        }
    }

    pub fn filter_digest(&mut self, interval_hash: &HashMap<String, Lapper<usize, u8>> ) {
        for (i, record) in self.records.iter_mut().enumerate() {
            let is_in_regions: bool =  record.is_in_regions(interval_hash);
            if is_in_regions {
                if record.filter_reason == "pass" {
                    record.filter_reason = String::from("pass");
                    self.filter_reasons[i] = String::from("pass");
                } else {
                    record.filter_reason = String::from("low_mq");
                    self.filter_reasons[i] = String::from("low_mq");
                }
              
            } else {
                if record.filter_reason == "pass" {
                    record.filter_reason = String::from("low_mq");
                    self.filter_reasons[i] = String::from("low_mq");
                }
            }
        }
    }

    pub fn filter_edges(&mut self, max_length: u64) {
        for (i, record) in self.records.iter_mut().enumerate() {
            if record.target_start < max_length || (record.target_length - record.target_end) < max_length {
                record.filter_reason = String::from("low_mq");
                self.filter_reasons[i] = String::from("low_mq");
            }
        }
    }

    pub fn stat(&self) -> HashMap<&str, u32> {
        let mut pass_count: u32 = 0;
        let mut singleton_count: u32 = 0;
        let mut low_mq_count: u32 = 0;
        let mut complex_count: u32 = 0;
        let mut stat_info: HashMap<&str, u32> = HashMap::new();


        for filter_reason in self.filter_reasons.iter() {
            if filter_reason == &String::from("pass") {
                pass_count += 1;
            } else if filter_reason == &String::from("singleton") {
                singleton_count +=1;
            } else if filter_reason == &String::from("low_mq") {
                low_mq_count += 1;
            } else if filter_reason == &String::from("complex") {
                complex_count += 1;
            } else {
                panic!("Unknown filter reason.");
            }
            
            stat_info.insert(&PASS_STR, pass_count);
            stat_info.insert(&SINGLETON_STR, singleton_count);
            stat_info.insert(&LOW_MQ_STR, low_mq_count);
            stat_info.insert(&COMPLEX_STR, complex_count);

        }
        stat_info
    }

    pub fn info(&self) -> &str {
        let stat_info = self.stat();

        if *stat_info.get("pass").unwrap_or(&0) >= 2 {
            return "pass";
        } else if *stat_info.get("low_mq").unwrap_or(&0) >=1 && self.records.len() > 1  {
            return "low_mq";
        } else if *stat_info.get("complex").unwrap_or(&0) >= 1 && self.records.len() > 1 {
            return "complex";
        } else {
           return "singleton";
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

    pub fn paf2table(&self, bed: &String,
                     output: &String, min_quality: &u8, 
                     min_identity: &f32,  min_length: &u32,
                     max_order: &u32, 
                     max_edge_length: &u64,
                    ) -> Result<(), Box<dyn Error>> {
        

        type IvU8 = Interval<usize, u8>;
    
        let is_filter_digest = if bed != "" { true } else { false };
       
        let bed = if is_filter_digest {
            Bed3::new(bed)
        } else {
            Bed3::new(&String::from(".tmp.bed"))
        };
        
        let interval_hash = bed.to_interval_hash();

        let parse_result = self.parse();
        
        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };
        
        let writer = common_writer(output);
        let mut wtr = csv::WriterBuilder::new()
                            .has_headers(false)
                            .delimiter(b'\t')
                            .from_writer(writer);
        
        let mut read_idx: u64 = 0;
        
        let mut filter_reason: &str = "";
        let mut old_query: String = String::from("");
        let mut concatemer: Concatemer  = Concatemer::new();
        
        let mut read_mapping_count: u64 = 0;
        let mut read_pass_count: u64 = 0;
        let mut read_singleton_count: u64 = 0;
        let mut read_low_mq_count: u64 = 0;
        let mut read_complex_count: u64 = 0;
        let is_filter_edges = *max_edge_length > 0;
        for (line_num, line) in rdr.deserialize().enumerate() {
            // parse error to continue
            let record: PAFLine = match line {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Error: Could not parse input file: {:?} at line {}", self.file_name(), line_num);
                    continue;
                },
            };
            if record.query != old_query && old_query != String::from("") {
                read_idx += 1;  
                
                if is_filter_digest {
                    concatemer.filter_digest(&interval_hash);
                }

                if is_filter_edges {
                    concatemer.filter_edges(*max_edge_length);
                }

                if concatemer.is_singleton() {
                    concatemer.parse_singleton();
                }
                
                if concatemer.is_complex(&max_order) {
                    concatemer.parse_complex();
                }

                for pcr in concatemer.records.iter() {
                    if pcr.filter_reason == "pass" {
                        wtr.serialize(&pcr).unwrap();
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
            
            if (&record.mapq < min_quality) || (&length < min_length) || (&identity < min_identity) {
                filter_reason = "low_mq";
            } else {
                filter_reason = "pass";
            }

            old_query = record.query.clone();

            let pcr = PoreCRecord::from_paf_record(record, read_idx, identity, filter_reason.to_string());
            
            concatemer.push(pcr);
        }
        
        // parse last record
        if concatemer.is_singleton() {
            concatemer.parse_singleton();
        }
        for pcr in concatemer.records.iter() {
            if pcr.filter_reason == "pass" {
                wtr.serialize(&pcr).unwrap();
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

        summary.save(&format!("{}.read.summary", self.prefix()));

        log::info!("Successful output Pore-C table `{}`", output);
        Ok(())
    }
}