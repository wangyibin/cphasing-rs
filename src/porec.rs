#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as anyResult;
use crossbeam_channel::{ bounded, Sender, Receiver };
use std::thread;
use itertools::{Itertools, Combinations};
use std::cmp::Ordering;
use std::collections::{ BTreeMap, HashMap };
use std::borrow::Cow;
use std::error::Error;
use std::path::Path;
use std::fs::File;
use std::sync::{Arc, Mutex};
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;
use rand::prelude::*;
use rust_lapper::{Interval, Lapper};

use crate::bed::{ Bed3, Bed4 };
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, binify, ChromSize, ChromSizeRecord };
use crate::pairs::{ PairRecord, PairHeader };
use crate::paf::PAFLine;

enum PosValue {
    U32(u32),
    U64(u64),
}


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
            read_idx,
            query_length: record.query_length,
            query_start: record.query_start,
            query_end: record.query_end,
            query_strand: record.query_strand,
            target: record.target,
            target_length: record.target_length,
            target_start: record.target_start,
            target_end: record.target_end,
            mapq: record.mapq,
            identity,
            filter_reason
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
    
    pub fn write_to(&self, buf: &mut String) {
        use std::fmt::Write as _;

        let _ = write!(buf, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.read_idx, self.query_length, self.query_start, self.query_end, self.query_strand,
            self.target, self.target_start, self.target_end, self.mapq, self.identity,
            self.filter_reason
        );
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
            read_idx,
            query_length: record.query_length,
            query_start: record.query_start,
            query_end: record.query_end,
            query_strand: record.query_strand,
            target: record.target,
            target_start: record.target_start,
            target_end: record.target_end,
            mapq: record.mapq,
            identity,
            filter_reason
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

#[derive(Debug, Default)]
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
        // self.records.sort_unstable_by_key(|x| (x.target, x.target_start));
        self.records.sort_unstable_by(|a, b| {
            match a.target.cmp(&b.target) {
                std::cmp::Ordering::Equal => a.target_start.cmp(&b.target_start),
                other => other,
            }
        });

    }

    // pub fn decompose(&mut self) -> Combinations<std::vec::IntoIter<PoreCRecord>> {
    //     let r = self.records.clone();
    //     r.into_iter().combinations(2)
        
    // }
    pub fn decompose(&self) -> Combinations<std::slice::Iter<'_, PoreCRecord>> {
        self.records.iter().combinations(2)
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

        vec.iter()
            .map(|(key, value)| format!("{}\t{}", key, value))
            .collect::<Vec<_>>()
            .join("\n")
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

    pub fn parse2(&mut self)  -> anyResult<Box<dyn BufRead + Send + 'static>> {
        let input = common_reader(&self.file);
      
        Ok(input)
    }

    pub fn to_pairs_pqs(&self, chromsizes: &String, output: &String, 
        chunksize: usize, min_quality: u8, 
        min_order: usize, max_order: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;
        use crate::pqs::_README as _readme;
        use crate::pqs::_METADATA;

        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
        };
        log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
        
        let _ = std::fs::create_dir_all(output);
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));

        // copy chromsizes to output
        let _ = std::fs::copy(chromsizes, format!("{}/_contigsizes", output));

        let contigsizes = ChromSize::new(chromsizes);
        let contigsizes_data = contigsizes.to_vec().unwrap();
        let max_contig_size = contigsizes_data.iter().map(|x| x.size).max().unwrap();
      
        let pos_type = if max_contig_size < 4294967295 {
            DataType::UInt32
        } else {
            DataType::UInt64
        };

        let pos_type_string = match pos_type {
            DataType::UInt32 => "UInt32",
            DataType::UInt64 => "UInt64",
            _ => "UInt32",
        };

        let mut wtr = common_writer(format!("{}/_readme", output).as_str());
        wtr.write_all(_readme.as_bytes()).unwrap();
        wtr.flush().unwrap();

        let create_date_time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

        let mut wtr = common_writer(format!("{}/_metadata", output).as_str());
        let mut _metadata = _METADATA.to_string();
        _metadata = _metadata.replace("REPLACE", &create_date_time);
        _metadata = _metadata.replace("CHUNKSIZE", &chunksize.to_string());
        _metadata = _metadata.replace("pos_type_lower", pos_type_string.to_lowercase().as_str());
        _metadata = _metadata.replace("pos_type", pos_type_string);
        wtr.write_all(_metadata.as_bytes()).unwrap();
        wtr.flush().unwrap();
        

        let mut concatemer_summary: ConcatemerSummary = ConcatemerSummary::new();

        let (sender, receiver) = bounded::<(usize, Vec<Concatemer>)>(100);
        let mut handles = vec![];
        

        #[inline]
        fn mid_u64(a: u64, b: u64) -> u64 {
            a.saturating_add(b) / 2
        }
        #[inline]
        fn mid_u32(a: u64, b: u64) -> u32 {
            (a.saturating_add(b) / 2) as u32
        }
        #[inline]
        fn enc_strand(c: char) -> i8 {
            match c {
                '+' => 1,
                '-' => -1,
                _ => 0,
            }
        }


        match max_contig_size {
            0..=4294967295 => {
                for _ in 0..8 {
                    let receiver: Receiver<_> = receiver.clone();
                    let output = output.clone();
                    handles.push(thread::spawn(move || {
                        while let Ok((chunk_id, records)) = receiver.recv() {
                            
                            let mut read_idx_vec: Vec<u64> = Vec::new();
                            let mut chrom1_vec: Vec<String> = Vec::new();
                            let mut pos1_vec: Vec<u32> = Vec::new();
                            let mut chrom2_vec: Vec<String> = Vec::new();
                            let mut pos2_vec: Vec<u32> = Vec::new();
                            let mut strand1_vec: Vec<String> = Vec::new();
                            let mut strand2_vec: Vec<String> = Vec::new();
                            let mut mapq_vec: Vec<u8> = Vec::new();
                            let mut read_idx = 0;
                            for mut concatemer in records {
                                concatemer.sort();
                                for pair in concatemer.decompose() {
                                    let (record1, record2) = (pair[0], pair[1]);
                                    read_idx_vec.push(read_idx);
                                    chrom1_vec.push(record1.target.clone());
                                    pos1_vec.push((record1.target_start as u32 + record1.target_end as u32) / 2);
                                    chrom2_vec.push(record2.target.clone());
                                    pos2_vec.push((record2.target_start as u32 + record2.target_end as u32) / 2);
                                    strand1_vec.push(record1.query_strand.to_string());
                                    strand2_vec.push(record2.query_strand.to_string());
                                    mapq_vec.push(std::cmp::min(record1.mapq, record2.mapq));
                                    
                                    read_idx += 1;
                                }
        
                            }
        
                            let df = DataFrame::new(vec![
                                Series::new("read_idx", read_idx_vec),
                                Series::new("chrom1", chrom1_vec),
                                Series::new("pos1", pos1_vec),
                                Series::new("chrom2", chrom2_vec),
                                Series::new("pos2", pos2_vec),
                                Series::new("strand1", strand1_vec),
                                Series::new("strand2", strand2_vec),
                                Series::new("mapq", mapq_vec),
                            ]).unwrap();
        
                            let mut df = df.lazy().with_column(
                                col("read_idx").cast(DataType::String)
                            ).with_column(
                                col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).with_column(
                                col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).with_column(
                                col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).with_column(
                                col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).collect().unwrap();
        
                            let file = format!("{}/q0/{}.parquet", output, chunk_id);
                            let mut file = File::create(file).unwrap();
                            ParquetWriter::new(&mut file).finish(&mut df).unwrap();
                            
                            let mut df = df.lazy().filter(
                                col("mapq").gt_eq(1)
                            ).collect().unwrap();
                
                            let file = format!("{}/q1/{}.parquet", output, chunk_id);
                            let mut file = File::create(file).unwrap();
                            ParquetWriter::new(&mut file)
                                .finish(&mut df)
                                .unwrap();
                        }
                    }))
                }
            },
            _ => {
                for _ in 0..8 {
                    let receiver: Receiver<_> = receiver.clone();
                    let output = output.clone();
                    handles.push(thread::spawn(move || {
                        while let Ok((chunk_id, records)) = receiver.recv() {
                            let mut read_idx_vec: Vec<u64> = Vec::new();
                            let mut chrom1_vec: Vec<String> = Vec::new();
                            let mut pos1_vec: Vec<u64> = Vec::new();
                            let mut chrom2_vec: Vec<String> = Vec::new();
                            let mut pos2_vec: Vec<u64> = Vec::new();
                            let mut strand1_vec: Vec<String> = Vec::new();
                            let mut strand2_vec: Vec<String> = Vec::new();
                            let mut mapq_vec: Vec<u8> = Vec::new();
                            let mut read_idx = 0;
                            for mut concatemer in records {
                                concatemer.sort();
                                for pair in concatemer.decompose() {
                                    let (record1, record2) = (pair[0], pair[1]);
                                    read_idx_vec.push(read_idx);
                                    chrom1_vec.push(record1.target.clone());
                                    pos1_vec.push(record1.target_start + record1.target_end / 2);
                                    chrom2_vec.push(record2.target.clone());
                                    pos2_vec.push(record2.target_start + record2.target_end / 2);
                                    strand1_vec.push(record1.query_strand.to_string());
                                    strand2_vec.push(record2.query_strand.to_string());
                                    mapq_vec.push(std::cmp::min(record1.mapq, record2.mapq));
                                    
                                    read_idx += 1;
                                }
        
                            }
        
                            let df = DataFrame::new(vec![
                                Series::new("read_idx", read_idx_vec),
                                Series::new("chrom1", chrom1_vec),
                                Series::new("pos1", pos1_vec),
                                Series::new("chrom2", chrom2_vec),
                                Series::new("pos2", pos2_vec),
                                Series::new("strand1", strand1_vec),
                                Series::new("strand2", strand2_vec),
                                Series::new("mapq", mapq_vec),
                            ]).unwrap();
        
                            let mut df = df.lazy().with_column(
                                col("read_idx").cast(DataType::String)
                            ).with_column(
                                col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).with_column(
                                col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).with_column(
                                col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).with_column(
                                col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                            ).collect().unwrap();
        
                            let file = format!("{}/q0/{}.parquet", output, chunk_id);
                            let mut file = File::create(file).unwrap();
                            ParquetWriter::new(&mut file).finish(&mut df).unwrap();
                            
                            let mut df = df.lazy().filter(
                                col("mapq").gt_eq(1)
                            ).collect().unwrap();
                
                           
                            let file = format!("{}/q1/{}.parquet", output, chunk_id);
                            let mut file = File::create(file).unwrap();
                            ParquetWriter::new(&mut file)
                                .finish(&mut df)
                                .unwrap();
                        }
                    }))
                }
            }
        }
        
        
        let mut concatemer: Concatemer = Concatemer::new();
        let mut batch = Vec::with_capacity(chunksize);
        let mut first_iteration = true;
        let mut previous_read_idx: u64 = 0;
        let mut record_count = 0;
        let mut chunk_id = 0 as usize;
        let mut line_count = 0;
        let mut total_pair_count = 0;   
        for (i, line) in rdr.deserialize().enumerate() {
            let record: PoreCRecord = match line {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}: {:?}", i + 1, error);
                    continue
                },
            }; 
            
            if !first_iteration && record.read_idx != previous_read_idx {
                let order = concatemer.count();
                if (order < max_order) && (order >= min_order) {
                    concatemer_summary.count(&concatemer);
                    batch.push(std::mem::take(&mut concatemer));
                    
                    record_count += order * (order - 1) / 2;
                } else {
                    concatemer.clear();
                }
                
            }
            first_iteration = false;
            previous_read_idx = record.read_idx;
            
            if record.mapq < min_quality {
                continue
            }
            line_count += 1;
            concatemer.push(record); 
            if record_count >= chunksize{
                total_pair_count += record_count;
                sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                record_count = 0;
                chunk_id += 1;
                line_count = 0;
            }
            
        }
        
        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        
        drop(sender);
        for handle in handles {
            handle.join().unwrap();
        }
    
        log::info!("Processed total {} lines", line_count);
        log::info!("Generated total {} pairs", total_pair_count);
        log::info!("Successful output pairs `{}`", output);
        
        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(&output).with_extension("").to_str().unwrap().to_string()
        };

        concatemer_summary.save(&format!("{}.concatemer.summary", output_prefix));

        Ok(())

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

        let mut read_id: u64 = 0;
        
        let mut first_iteration = true;
        for (i, line) in rdr.deserialize().enumerate() {
            let record: PoreCRecord = match line {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}", i + 1);
                    continue
                },
            };
            if !first_iteration && record.read_idx != old_read_idx {
                let order = concatemer.count();
                if (order < max_order) && (order >= min_order) {
                    concatemer.sort();
                    concatemer_summary.count(&concatemer);
                    for pair in concatemer.decompose() {
                        wtr.serialize(PairRecord::from_pore_c_pair(pair, read_id)).unwrap();
                        read_id += 1;
                    }
                      
                }
                
                concatemer.clear();
            }
            first_iteration = false;
            old_read_idx = record.read_idx;
            if record.mapq < min_quality {
                continue
            }
            
            concatemer.push(record);
           
        }

        // process last concatemer
        if (concatemer.count() < max_order) & (concatemer.count() >= min_order) {
            concatemer.sort();

            let mut batch = Vec::with_capacity(max_order);
            concatemer_summary.count(&concatemer);
            for pair in concatemer.decompose() {
                batch.push(PairRecord::from_pore_c_pair(pair, read_id));
                read_id += 1;
            }

            if !batch.is_empty() {
                for record in batch {
                    wtr.serialize(record).unwrap();
                }
            }
        }

        log::info!("Successful output pairs `{}`", output);
        
        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(&output).with_extension("").to_str().unwrap().to_string()
        };

        concatemer_summary.save(&format!("{}.concatemer.summary", output_prefix));
        Ok(())
    }

    
    // pub fn to_pairs(&self, chromsizes: &String, output: &String, min_quality: u8, 
    //                 min_order: usize, max_order: usize) -> Result<(), Box<dyn Error>> {
    //     let parse_result = self.parse();
    //     let mut rdr = match parse_result {
    //         Ok(v) => v,
    //         Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
    //     };
    //     log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
    //     let chromsizes: ChromSize = ChromSize::new(chromsizes);
    //     let chromsizes_data: Vec<ChromSizeRecord> = chromsizes.to_vec().unwrap();
    
    //     let mut ph: PairHeader = PairHeader::new();
    //     ph.from_chromsizes(chromsizes_data);
    
    //     let mut writer = common_writer(output);
    //     writer.write_all(ph.to_string().as_bytes()).unwrap();
        
    //     let wtr = Arc::new(Mutex::new(writer));

    //     let mut concatemer: Arc<Mutex<Concatemer>> = Arc::new(Mutex::new(Concatemer::new()));  
    //     let mut concatemer_summary: Arc<Mutex<ConcatemerSummary>> = Arc::new(Mutex::new(ConcatemerSummary::new()));
    
    //     let mut old_read_idx: u64 = 0; 
    //     let mut flag: bool = false; 
    //     let mut read_id: u64 = 0;
        
    //     let mut first_iteration = true;
    //     let (sender, receiver) = bounded::<PoreCRecord>(1000);
    //     let mut handles = vec![];
    //     for _ in 0..4 {
    //         let receiver = receiver.clone();
    //         let concatemer = Arc::clone(&concatemer);
    //         let concatemer_summary = Arc::clone(&concatemer_summary);
    //         handles.push(thread::spawn(move || {
    //             while let Ok(record) = receiver.recv() {
    //                 let mut concatemer = concatemer.lock().unwrap();
    //                 let mut concatemer_summary = concatemer_summary.lock().unwrap();
    //                 if !first_iteration && record.read_idx != old_read_idx {
    //                     let order = concatemer.count();
    //                     if (order < max_order) && (order >= min_order) {
    //                         concatemer.sort();
    //                         concatemer_summary.count(&concatemer);
    //                         for pair in concatemer.decompose() {
    //                             let mut wtr = wtr.lock().unwrap();
    //                             writeln!(wtr, "{}", PairRecord::from_pore_c_pair(pair, read_id)).unwrap();
    //                             read_id += 1;
    //                         }
    //                     }
    //                     concatemer.clear();
    //                 }
    //                 first_iteration = false;
    //                 old_read_idx = record.read_idx;
    //                 if record.mapq < min_quality {
    //                     continue
    //                 }
    //                 concatemer.push(record);
    //             }
    //         }));
    //     }
    
    //     for (i, line) in rdr.deserialize().enumerate() {
    //         let record: PoreCRecord = match line {
    //             Ok(v) => v,
    //             Err(error) => {
    //                 log::warn!("Could not parse line {}", i + 1);
    //                 continue
    //             },
    //         };
    //         sender.send(record).unwrap();
    //     }
    
    //     drop(sender);
    //     for handle in handles {
    //         handle.join().unwrap();
    //     }
    
    //     let mut concatemer = Arc::try_unwrap(concatemer).unwrap().into_inner().unwrap();
    //     let mut concatemer_summary = Arc::try_unwrap(concatemer_summary).unwrap().into_inner().unwrap();
    
    //     // process last concatemer
    //     if (concatemer.count() < max_order) & (concatemer.count() >= min_order) {
    //         concatemer.sort();
    
    //         let mut batch = Vec::with_capacity(max_order);
    //         concatemer_summary.count(&concatemer);
    //         for pair in concatemer.decompose() {
    //             batch.push(PairRecord::from_pore_c_pair(pair, read_id));
    //             read_id += 1;
    //         }
    
    //         if !batch.is_empty() {
    //             let mut wtr = wtr.lock().unwrap();
    //             for record in batch {
    //                 writeln!(wtr, "{}", record);
    //             }
    //         }
    //     }
    
    //     log::info!("Successful output pairs `{}`", output);
        
    //     let output_prefix = if output == "-" {
    //         Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
    //     } else {
    //         Path::new(&output).with_extension("").to_str().unwrap().to_string()
    //     };
    
    //     concatemer_summary.save(&format!("{}.concatemer.summary", output_prefix));
    //     Ok(())
    // }

    // pub fn to_depth(&mut self, contigsizes: &String, binsize:u32, min_quality: u8, output: &String) {
    //     use hashbrown::HashMap as BrownHashMap;
    //     let mut parse_result = self.parse().unwrap();
    //     let contigsizes = ChromSize { file: contigsizes.to_string() };
    //     let contigsizes_data = contigsizes.data().unwrap();

    //     let mut depth: BrownHashMap<String, BTreeMap<u32, u32>> = BrownHashMap::new();
    //     let bins_db = binify(&contigsizes_data, binsize).unwrap();
    //     // incomplete
    // }

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
                let _ = wtr.write_record(&record);
            }

        }
        
    
    log::info!("Successful output intersection porec table into `{}`", output);

    }

    pub fn intersect_multi_threads(&mut self, hcr_bed: &String, invert: bool, output: &String) {
        type IvU8 = Interval<usize, u8>;
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        let wtr = common_writer(output);
        

        let (sender, receiver) = bounded::<Vec<String>>(1000);

        let mut handles = vec![];
        let wtr = Arc::new(Mutex::new(wtr));
        
        for _ in 0..8 {
            let interval_hash = interval_hash.clone();
            let wtr = Arc::clone(&wtr);
            let receiver = receiver.clone();
            handles.push(thread::spawn(move || {
                while let Ok(records) = receiver.recv() {
                    // let mut data = vec![];
                    // for record in records {
                    //     let target_start = record[6].parse::<usize>().unwrap();
                    //     let target_end = record[7].parse::<usize>().unwrap();

                    //     let is_in_regions = interval_hash.get(&record[5]).map_or(false, |interval|{
                    //         interval.count(target_start, target_end) > 0 });

                    //     if is_in_regions ^ invert {
                    //         let record = record.iter().join("\t");
                    //         if !record.is_empty() {
                    //             data.push(record);
                    //         }
                            
                    //     }
                    // }
                    let data = records.par_iter().filter_map(|record| {
                        let record = record.split("\t").collect::<Vec<_>>();
                        let target_start = record[6].parse::<usize>().unwrap();
                        let target_end = record[7].parse::<usize>().unwrap();

                        let is_in_regions = interval_hash.get(record[5]).map_or(false, |interval|{
                            interval.count(target_start, target_end) > 0 });

                        if is_in_regions ^ invert {
                            Some(record.iter().join("\t"))
                        } else {
                            None
                        }
                    }).collect::<Vec<_>>();

                    if !data.is_empty() {
                        let mut wtr = wtr.lock().unwrap();
                        let data = data.join("\n") + "\n";
                        wtr.write_all(data.as_bytes()).unwrap();
                    }
                }
            }));
        }

        let batch_size = 10_000;
        let mut batch = Vec::with_capacity(batch_size);
        for (idx, record) in self.parse2().unwrap().lines().enumerate() {
            let record = match record {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}", idx + 1);
                    continue
                },
            };
            batch.push(record);
            if batch.len() == batch_size {
                sender.send(std::mem::take(&mut batch)).unwrap();
            }
        }

        drop(sender);

        for handle in handles {
            handle.join().unwrap();
        }

        log::info!("Successful output intersection porec table into `{}`", output);

    }

    pub fn break_contigs(&mut self, break_bed: &String, output: &String) {
        type IvString = Interval<usize, String>;
        let bed = Bed4::new(break_bed);
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
            
            let is_break_contig = interval_hash.contains_key(&record[5]);
            
            if is_break_contig {
                
                let interval = interval_hash.get(&record[5]).unwrap();
                let res = interval.find(target_start, target_end).collect::<Vec<_>>();
                
                if res.len() > 0 {
                    let break_contig_length = res[0].stop - res[0].start + 1;
                    let new_target_start = target_start - res[0].start + 1;
                    let new_target_end = target_end - res[0].start + 1;
                    if new_target_end <= break_contig_length {
                        let mut new_record = csv::StringRecord::new();
                        new_record.push_field(&record[0]);
                        new_record.push_field(&record[1]);
                        new_record.push_field(&record[2]);
                        new_record.push_field(&record[3]);
                        new_record.push_field(&record[4]);
                        // let target = record[5].to_string();
                        // let new_target = format!("{}:{}-{}", target, res[0].start, res[0].stop);
                        let new_target = &res[0].val;
                        new_record.push_field(&new_target);
                        
                    
                        new_record.push_field(&new_target_start.to_string());
                        new_record.push_field(&new_target_end.to_string());
                        new_record.push_field(&record[8]);
                        new_record.push_field(&record[9]);
                        new_record.push_field(&record[10]);
                    
                    
                        let _ = wtr.write_record(&new_record);
                    }
                   

                } else {
                    let _ = wtr.write_record(&record);
                }
               
               
            } else {
                let _ = wtr.write_record(&record);
            }

        }
        log::info!("Successful output contigs corrected porec table into `{}`", output);
    }

    pub fn dup(&mut self, collapsed_list: &String, seed: usize, output: &String) {
        let reader = common_reader(collapsed_list);
        let mut collapsed_contigs: HashMap<String, Vec<String>> = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let s: Vec<&str> = record.split("\t").collect();
            if s.len() != 2 {
                log::warn!("Invalid record: {}", record);
                continue;
            }
            let contig1 = s[0].to_string();
            let contig2 = s[1].to_string();
            collapsed_contigs.entry(contig1.clone()).or_insert(vec![contig1]).push(contig2);
        }
        
        let seed_bytes = seed.to_ne_bytes();
        let mut seed_array = [0u8; 32];
        for (i, byte) in seed_bytes.iter().enumerate() {
            seed_array[i] = *byte;
        }

        let reader = common_reader(&self.file);

        let mut wtr = common_writer(output);
       
        let output_counts = 0;
        let mut rng = StdRng::from_seed(seed_array);

        for record in reader.lines() {
            let record = record.unwrap(); {
                let mut s: Vec<&str> = record.trim().split("\t").collect();
                let mut contig = s[5].to_string();

                if collapsed_contigs.contains_key(&contig) {
                    let idx = rng.gen_range(0..collapsed_contigs.get(&contig).unwrap().len());
                    contig = collapsed_contigs.get(&contig).unwrap()[idx].clone();
                }

                s[5] = &contig;
                writeln!(wtr, "{}", s.join("\t")).unwrap();
                
                
            }
        }
        
    }

    pub fn split(&mut self, output: &String) {
        let reader = common_reader(&self.file);
        let mut wtr = common_writer(output);
        
        let (sender, receiver) = bounded::<Vec<(u32, String)>>(100);   




        log::info!("Successful output split porec table into `{}`", output);
    }
}

pub fn merge_porec_tables(input: Vec<&String>, output: &String) {
    let mut wtr = common_writer(output);
    let mut idx = 0;

    for file in input {
        let reader = common_reader(file);
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
        idx += max_idx + 1;
    }
    log::info!("Successful output merge porec tables into `{}`", output);
}


// #[derive(Debug)]
// pub struct PoreCParquet {
//     file: String,
// }

// impl PoreCParquet {
//     pub fn new(name: &String) -> PoreCParquet {
//         PoreCParquet { file: name.clone() }
//     }

//     pub fn file_name(&self) -> Cow<'_, str> {
//         let path = Path::new(&self.file);
//         path.file_name().expect("REASON").to_string_lossy()
//     }    

//     pub fn prefix(&self) -> String {
//         let binding = self.file_name().to_string();
//         let file_path = Path::new(&binding);
//         let file_prefix = file_path.file_stem().unwrap().to_str().unwrap();

//         (*file_prefix).to_string()
//     }

//     pub fn to_porec(&self, output: &String) {
//         let mut wtr = common_writer(output);
//         let parquet = parquet::read(&self.file).unwrap();
//         let schema = parquet.schema();
//         let mut record = vec![];
//         for row in parquet.rows() {
//             for (i, column) in row.iter().enumerate() {
//                 let column = column.as_any();
//                 let column = column.downcast_ref::<parquet::column::String>().unwrap();
//                 let value = column.get(0);
//                 record.push(value);
//             }
//             writeln!(wtr, "{}", record.join("\t"));
//             record.clear();
//         }
//         log::info!("Successful output parquet file into `{}`", output);
//     }
// }