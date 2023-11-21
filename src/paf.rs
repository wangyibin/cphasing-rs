#[warn(unused_assignments)]
use anyhow::Result as anyResult;
use std::any::type_name;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::ffi::OsStr;
use std::path::Path;
use std::error::Error;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};


use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;
use crate::porec::PoreCRecord;


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
        let pass_count = self.filter_reasons.iter().filter(|&x| *x == "pass").count();
        if pass_count == 1 {
            true
        } else {
            false
        }
    }

    pub fn parse_singleton(&mut self) {
        for (i, record) in self.records.iter_mut().enumerate() {
            if record.filter_reason == "pass" {
                record.filter_reason = String::from("singleton");
                self.filter_reasons[i] = String::from("singleton");
            }
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



    pub fn stat(&self) -> HashMap<String, u32> {
        let mut pass_count: u32 = 0;
        let mut singleton_count: u32 = 0;
        let mut low_mq_count: u32 = 0;
        let mut complex_count: u32 = 0;
        let mut stat_info: HashMap<String, u32> = HashMap::new();

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
            
            stat_info.insert("pass".to_string(), pass_count);
            stat_info.insert("singleton".to_string(), singleton_count);
            stat_info.insert("low_mq".to_string(), low_mq_count);
            stat_info.insert("complex".to_string(), complex_count);

        }
        stat_info
    }

    pub fn info(&self) -> &str {
        let stat_info = self.stat();

        if stat_info["pass"] >= 2 {
            return "pass";
        } else if stat_info["low_mq"] >=1 && self.records.len() > 1  {
            return "low_mq";
        } else if stat_info["complex"] >= 1 && self.records.len() > 1 {
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

    pub fn paf2table(&self, output: &String, min_quality: &u8, 
                     min_identity: &f32,  min_length: &u32,
                     max_order: &u32,
                    ) -> Result<(), Box<dyn Error>> {
        
        let parse_result = self.parse();

        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };
        
        let output = common_writer(output);
        let mut wtr = csv::WriterBuilder::new()
                            .has_headers(false)
                            .delimiter(b'\t')
                            .from_writer(output);
        
        let mut read_idx: u64 = 0;
        
        let mut filter_reason: &str = "";
        let mut old_query: String = String::from("");
        let mut concatemer: Concatemer  = Concatemer::new();
        
        let mut read_mapping_count: u64 = 0;
        let mut read_pass_count: u64 = 0;
        let mut read_singleton_count: u64 = 0;
        let mut read_low_mq_count: u64 = 0;
        let mut read_complex_count: u64 = 0;

        for line in rdr.deserialize() {
            let record: PAFLine = line?;
            if record.query != old_query && old_query != String::from("") {
                read_idx += 1;  
            
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

        log::info!("Done.");
        Ok(())
    }
}