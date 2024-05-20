#[allow(dead_code)]
use anyhow::Result as anyResult;
use crossbeam_channel::unbounded;
use std::thread;
use itertools::{Itertools, Combinations};
use std::cmp::Ordering;
use std::collections::{ BTreeMap, HashMap };
use std::borrow::Cow;
use std::error::Error;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};

use crate::bed::Bed3;
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, binify, ChromSize, ChromSizeRecord };
use crate::pairs::{ PairRecord, PairHeader };
use crate::paf::PAFLine;


#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PoreCRecordPlus {
    pub read_idx: u64,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub query_strand: char,
    pub target: String, 
    #[serde(skip_serializing)]
    pub target_length: u64,
    pub target_start: u64, 
    pub target_end: u64,
    pub mapq: u8,
    pub identity: f32, 
    pub filter_reason: String,
} 




#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PoreCRecord {
    pub read_idx: u64,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub query_strand: char,
    pub target: String, 
    pub target_start: u64, 
    pub target_end: u64,
    pub mapq: u8,
    pub identity: f32, 
    pub filter_reason: String,
} 

impl PartialOrd for PoreCRecordPlus {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.query_start.partial_cmp(&other.query_start)
    }
}

impl PartialEq for PoreCRecordPlus {
    fn eq(&self, other: &Self) -> bool {
        self.query_start == other.query_start
    }
}



impl PoreCRecordPlus {
    pub fn from_paf_record(record: PAFLine, read_idx: u64,
                            identity: f32, filter_reason: String) -> Self {
        PoreCRecordPlus {
            read_idx: read_idx,
            query_length: record.query_length,
            query_start: record.query_start,
            query_end: record.query_end,
            query_strand: record.query_strand,
            target: record.target,
            target_length: record.target_length,
            target_start: record.target_start,
            target_end: record.target_end,
            mapq: record.mapq,
            identity: identity,
            filter_reason: filter_reason
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.read_idx,
            self.query_length,
            self.query_start,
            self.query_end,
            self.query_strand,
            self.target,
            self.target_start,
            self.target_end,
            self.mapq,
            self.identity,
            self.filter_reason,
        )
    }

    pub fn is_in_regions(&self, interval_hash: &HashMap<String, Lapper<usize, u8>>) -> bool {
        let is_in_regions: bool = if let Some(interval) = interval_hash.get(&self.target) {
            let iv_start = interval.count((self.target_start - 1) as usize, self.target_start as usize);
            let iv_end = interval.count((self.target_end - 1) as usize, self.target_end as usize);
            if iv_start > 0 && iv_end > 0 {
                true
            } else {
                false
            }
        } else {
            false
        };
        is_in_regions
    }


}


impl PoreCRecord {
    pub fn from_paf_record(record: PAFLine, read_idx: u64,
                            identity: f32, filter_reason: String) -> PoreCRecord{
        PoreCRecord {
            read_idx: read_idx,
            query_length: record.query_length,
            query_start: record.query_start,
            query_end: record.query_end,
            query_strand: record.query_strand,
            target: record.target,
            target_start: record.target_start,
            target_end: record.target_end,
            mapq: record.mapq,
            identity: identity,
            filter_reason: filter_reason
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.read_idx,
            self.query_length,
            self.query_start,
            self.query_end,
            self.query_strand,
            self.target,
            self.target_start,
            self.target_end,
            self.mapq,
            self.identity,
            self.filter_reason,
        )
    }

    pub fn is_in_regions(&self, interval_hash: &HashMap<String, Lapper<usize, u8>>) -> bool {
        let is_in_regions: bool = if let Some(interval) = interval_hash.get(&self.target) {
            let iv_start = interval.count((self.target_start - 1) as usize, self.target_start as usize);
            let iv_end = interval.count((self.target_end - 1) as usize, self.target_end as usize);
            if iv_start > 0 && iv_end > 0 {
                true
            } else {
                false
            }
        } else {
            false
        };
        is_in_regions
    }


}

impl PartialOrd for PoreCRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.query_start.partial_cmp(&other.query_start)
    }
}

impl PartialEq for PoreCRecord {
    fn eq(&self, other: &Self) -> bool {
        self.query_start == other.query_start
    }
}

#[derive(Debug)]
pub struct Concatemer {
    pub records: Vec<PoreCRecord>,
}

impl Concatemer {

    pub fn new() -> Concatemer{
        Concatemer {
            records: Vec::new(),
        }
    }

    pub fn push(&mut self, pcr: PoreCRecord) {
        self.records.push(pcr);
    }

    pub fn clear(&mut self) {
        self.records.clear();
    }

    pub fn count(&self) -> usize {
        self.records.len()
    }

    pub fn sort(&mut self) {
        // self.records.sort_by(| a, b | a.partial_cmp(&b).unwrap());
        self.records.sort_unstable_by_key(|x| (x.target.clone(), x.target_start));
    }

    pub fn decompose(&mut self) -> Combinations<std::vec::IntoIter<PoreCRecord>> {
        let r = self.records.clone();
        r.into_iter().combinations(2)
        
    }

}

#[derive(Debug, Clone)]
pub struct ConcatemerSummary {
    pub summary: HashMap<u32, u64>,
}

impl ConcatemerSummary {
    pub fn new() -> ConcatemerSummary {
        ConcatemerSummary { summary: HashMap::<u32, u64>::new() }
    }
    
    pub fn count(&mut self, concatemer: &Concatemer) {

        let concatemer_count: u32 = concatemer.count().try_into().unwrap();
        *self.summary.entry(concatemer_count).or_insert(0) += 1;

    }

    pub fn to_string(&self) -> String {
        let mut vec = self.summary.iter().collect::<Vec<_>>();
        vec.sort_by_key(|(key, _)| *key);

        let string = vec.iter()
                        .map(|(key, value)| format!("{}\t{}", key, value))
                        .collect::<Vec<_>>()
                        .join("\n");

        string
    }

    pub fn save(&self, output: &String) {
        let mut wtr = common_writer(output);    

        let result: String = self.to_string();

        wtr.write_all(result.as_bytes()).unwrap();
        log::info!("Successful output summary of concatemer `{}`", output);
    }
    // pub fn collapse(&self) -> HashMap<&str, u64> {
        
    // }
}

#[derive(Debug)]
pub struct PoreCTable {
    file: String,
}

impl BaseTable for PoreCTable {
    fn new(name: &String) -> PoreCTable {
        PoreCTable { file: name.clone() }
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


impl PoreCTable {
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

    pub fn to_pairs(&self, chromsizes: &String, output: &String, min_quality: u8, 
                        min_order: usize, max_order: usize,
                    ) -> Result<(), Box<dyn Error>> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
        };
        log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
        let chromsizes: ChromSize = ChromSize::new(chromsizes);
        let chromsizes_data: Vec<ChromSizeRecord> = chromsizes.to_vec().unwrap();

        let mut ph: PairHeader = PairHeader::new();
        ph.from_chromsizes(chromsizes_data);


        let mut writer = common_writer(output);
        writer.write_all(ph.to_string().as_bytes()).unwrap();

        let mut wtr = csv::WriterBuilder::new()
                            .has_headers(false)
                            .delimiter(b'\t')
                            .from_writer(writer);

        let mut concatemer: Concatemer = Concatemer::new();  
        let mut concatemer_summary: ConcatemerSummary = ConcatemerSummary::new();

        let mut old_read_idx: u64 = 0; 
        let mut flag: bool = false; 
        let mut read_id: u64 = 0;
        
        for (i, line) in rdr.deserialize().enumerate() {
            let record: PoreCRecord = match line {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}", i + 1);
                    continue
                },
            };
            if record.read_idx != old_read_idx && flag == true {
                let order = concatemer.count();
                if (order < max_order) & (order >= min_order) {
                    concatemer.sort();
                    concatemer_summary.count(&concatemer);
                    for pair in concatemer.decompose() {
                        wtr.serialize(PairRecord::from_pore_c_pair(pair, read_id)).unwrap();
                        read_id += 1;
                    }
                    
                    
                }
                
                concatemer.clear();
            }
            flag = true;
            old_read_idx = record.read_idx;
            if record.mapq < min_quality {
                continue
            }
            
            concatemer.push(record);
           
        }

        // process last concatemer
        if (concatemer.count() < max_order) & (concatemer.count() >= min_order) {
            concatemer.sort();
            concatemer_summary.count(&concatemer);
            for pair in concatemer.decompose() {
                wtr.serialize(PairRecord::from_pore_c_pair(pair, read_id)).unwrap();
                read_id += 1;
            }
        }

        log::info!("Successful output pairs `{}`", output);

        concatemer_summary.save(&format!("{}.concatemer.summary", self.prefix()));
        Ok(())
    }

    pub fn to_depth(&mut self, contigsizes: &String, binsize:u64, min_quality: u8, output: &String) {
        use hashbrown::HashMap as BrownHashMap;
        let mut parse_result = self.parse().unwrap();
        let contigsizes = ChromSize { file: contigsizes.to_string() };
        let contigsizes_data = contigsizes.data().unwrap();

        let mut depth: BrownHashMap<String, BTreeMap<u64, u64>> = BrownHashMap::new();
        let bins_db = binify(&contigsizes_data, binsize).unwrap();
    }

    pub fn intersect(&mut self, hcr_bed: &String, invert: bool, output: &String) {
        type IvU8 = Interval<usize, u8>;
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        let writer = common_writer(output);
        let mut wtr = csv::WriterBuilder::new()
                        .has_headers(false)
                        .delimiter(b'\t')
                        .from_writer(writer);

        for (i, line) in self.parse().unwrap().records().enumerate() {
            let record = match line {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}", i + 1);
                    continue
                },
            };

            let target_start = record[6].parse::<usize>().unwrap();
            let target_end = record[7].parse::<usize>().unwrap();

            let is_in_regions = interval_hash.get(&record[5]).map_or(false, |interval|{
                    interval.count(target_start, target_end) > 0 });

            if is_in_regions ^ invert {
                wtr.write_record(&record);
            }

        }
        
       

        // let output = output.to_string();
        // let input = common_reader(&self.file);
        // let mut rdr = csv::ReaderBuilder::new()
        //                     .flexible(true)
        //                     .has_headers(false)
        //                     .comment(Some(b'#'))
        //                     .delimiter(b'\t')
        //                     .from_reader(input);

        // let (sender, receiver) = unbounded();

        // let producer = thread::spawn(move || {
        //     for (i, line) in rdr.records().enumerate() {
        //         let record = match line {
        //             Ok(v) => v,
        //             Err(error) => {
        //                 log::warn!("Could not parse line {}", i + 1);
        //                 continue;
        //             },
        //     };

        //     let target_start = record[6].parse::<usize>().unwrap();
        //     let target_end = record[7].parse::<usize>().unwrap();

        //     let is_in_regions = interval_hash.get(&record[5]).map_or(false, |interval|{
        //         interval.count(target_start, target_end) > 0 });

        //     if is_in_regions ^ invert {
        //         sender.send(record).unwrap();
        //     }
        // }

    // });


    
    // let writer = common_writer(&output);
    // let mut wtr = csv::WriterBuilder::new()
    //             .has_headers(false)
    //             .delimiter(b'\t')
    //             .from_writer(writer);

    // for record in receiver {
    //     wtr.write_record(&record).unwrap();
    // }

    // wtr.flush().unwrap();
    
    log::info!("Successful output intersection porec table into `{}`", output);

    }
}

pub fn merge_porec_tables(input: Vec<&String>, output: &String) {
    let mut wtr = common_writer(output);
    let mut idx = 0;


    for file in input {
        // let mut rdr = PoreCTable::new(file).parse().unwrap();
        // let mut max_idx: u64 = 0;
        // 'inner: for (i, line) in rdr.deserialize().enumerate() {
        //     let line_number = i + 1;
        //     let mut record: PoreCRecord = match line {
        //         Ok(v) => v,
        //         Err(error) => {
        //             log::warn!("Could not parse line {}", line_number);
        //             continue 'inner;
        //         },
        //     };
        //     if record.read_idx > max_idx {
        //         max_idx = record.read_idx.clone();
        //     }
        //     record.read_idx += idx ;
        
        let mut reader = common_reader(file);
        let mut max_idx: u64 = 0;
        for (i, line) in reader.lines().enumerate() {
            let line = line.unwrap();
            let mut line = line.split("\t");
            let mut read_idx: u64 = line.next().unwrap().parse().unwrap();
            
            if read_idx > max_idx {
                max_idx = read_idx;
            }
            
            read_idx += idx as u64;
            let mut record = line.collect::<Vec<&str>>().join("\t");
            record = format!("{}\t{}", read_idx, record);
            wtr.write_all(record.to_string().as_bytes()).unwrap();
            wtr.write_all(b"\n").unwrap();
        }
        idx += max_idx;
    }
    log::info!("Successful output merge porec tables into `{}`", output);
}

