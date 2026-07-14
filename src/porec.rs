#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as anyResult;
use coitrees::{COITree, IntervalNode, IntervalTree};
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
use std::sync::atomic::{AtomicU64, Ordering as AtomOrdering };
use std::io::{ Write, BufReader, BufRead };
use std::fmt::Write as FmtWrite;
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

    // pub fn to_pairs_pqs(&self, chromsizes: &String, output: &String, 
    //     chunksize: usize, min_quality: u8, 
    //     min_order: usize, max_order: usize,
    // ) -> anyResult<()> {
    //     use polars::prelude::*;
    //     use crate::pqs::_README as _readme;
    //     use crate::pqs::_METADATA;

    //     let parse_result = self.parse();
    //     let mut rdr = match parse_result {
    //         Ok(v) => v,
    //         Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
    //     };
    //     log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
        
    //     let _ = std::fs::create_dir_all(output);
    //     let _ = std::fs::create_dir_all(format!("{}/q0", output));
    //     let _ = std::fs::create_dir_all(format!("{}/q1", output));

    //     // copy chromsizes to output
    //     let _ = std::fs::copy(chromsizes, format!("{}/_contigsizes", output));

    //     let contigsizes = ChromSize::new(chromsizes);
    //     let contigsizes_data = contigsizes.to_vec().unwrap();
    //     let max_contig_size = contigsizes_data.iter().map(|x| x.size).max().unwrap();
      
    //     let pos_type = if max_contig_size < 4294967295 {
    //         DataType::UInt32
    //     } else {
    //         DataType::UInt64
    //     };

    //     let pos_type_string = match pos_type {
    //         DataType::UInt32 => "UInt32",
    //         DataType::UInt64 => "UInt64",
    //         _ => "UInt32",
    //     };

    //     let mut wtr = common_writer(format!("{}/_readme", output).as_str());
    //     wtr.write_all(_readme.as_bytes()).unwrap();
    //     wtr.flush().unwrap();

    //     let create_date_time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

    //     let mut wtr = common_writer(format!("{}/_metadata", output).as_str());
    //     let mut _metadata = _METADATA.to_string();
    //     _metadata = _metadata.replace("REPLACE", &create_date_time);
    //     _metadata = _metadata.replace("CHUNKSIZE", &chunksize.to_string());
    //     _metadata = _metadata.replace("pos_type_lower", pos_type_string.to_lowercase().as_str());
    //     _metadata = _metadata.replace("pos_type", pos_type_string);
    //     wtr.write_all(_metadata.as_bytes()).unwrap();
    //     wtr.flush().unwrap();
        

    //     let mut concatemer_summary: ConcatemerSummary = ConcatemerSummary::new();

    //     let (sender, receiver) = bounded::<(usize, Vec<Concatemer>)>(100);
    //     let mut handles = vec![];
        

    //     #[inline]
    //     fn mid_u64(a: u64, b: u64) -> u64 {
    //         a.saturating_add(b) / 2
    //     }
    //     #[inline]
    //     fn mid_u32(a: u64, b: u64) -> u32 {
    //         (a.saturating_add(b) / 2) as u32
    //     }
    //     #[inline]
    //     fn enc_strand(c: char) -> i8 {
    //         match c {
    //             '+' => 1,
    //             '-' => -1,
    //             _ => 0,
    //         }
    //     }


    //     match max_contig_size {
    //         0..=4294967295 => {
    //             for _ in 0..8 {
    //                 let receiver: Receiver<_> = receiver.clone();
    //                 let output = output.clone();
    //                 handles.push(thread::spawn(move || {
    //                     while let Ok((chunk_id, records)) = receiver.recv() {
                            
    //                         let mut read_idx_vec: Vec<u64> = Vec::new();
    //                         let mut chrom1_vec: Vec<String> = Vec::new();
    //                         let mut pos1_vec: Vec<u32> = Vec::new();
    //                         let mut chrom2_vec: Vec<String> = Vec::new();
    //                         let mut pos2_vec: Vec<u32> = Vec::new();
    //                         let mut strand1_vec: Vec<String> = Vec::new();
    //                         let mut strand2_vec: Vec<String> = Vec::new();
    //                         let mut mapq_vec: Vec<u8> = Vec::new();
    //                         let mut read_idx = 0;
    //                         for mut concatemer in records {
    //                             concatemer.sort();
    //                             for pair in concatemer.decompose() {
    //                                 let (record1, record2) = (pair[0], pair[1]);
    //                                 read_idx_vec.push(read_idx);
    //                                 chrom1_vec.push(record1.target.clone());
    //                                 pos1_vec.push((record1.target_start as u32 + record1.target_end as u32) / 2);
    //                                 chrom2_vec.push(record2.target.clone());
    //                                 pos2_vec.push((record2.target_start as u32 + record2.target_end as u32) / 2);
    //                                 strand1_vec.push(record1.query_strand.to_string());
    //                                 strand2_vec.push(record2.query_strand.to_string());
    //                                 mapq_vec.push(std::cmp::min(record1.mapq, record2.mapq));
                                    
    //                                 read_idx += 1;
    //                             }
        
    //                         }
        
    //                         let df = DataFrame::new(vec![
    //                             Series::new("read_idx".into(), read_idx_vec).into(),
    //                             Series::new("chrom1".into(), chrom1_vec).into(),
    //                             Series::new("pos1".into(), pos1_vec).into(),
    //                             Series::new("chrom2".into(), chrom2_vec).into(),
    //                             Series::new("pos2".into(), pos2_vec).into(),
    //                             Series::new("strand1".into(), strand1_vec).into(),
    //                             Series::new("strand2".into(), strand2_vec).into(),
    //                             Series::new("mapq".into(), mapq_vec).into(),
    //                         ]).unwrap();
        
    //                         let mut df = df.lazy().with_column(
    //                             col("read_idx").cast(DataType::String)
    //                         ).with_column(
    //                             col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).with_column(
    //                             col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).with_column(
    //                             col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).with_column(
    //                             col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).collect().unwrap();
        
    //                         let file = format!("{}/q0/{}.parquet", output, chunk_id);
    //                         let mut file = File::create(file).unwrap();
    //                         ParquetWriter::new(&mut file).finish(&mut df).unwrap();
                            
    //                         let mut df = df.lazy().filter(
    //                             col("mapq").gt_eq(1)
    //                         ).collect().unwrap();
                
    //                         let file = format!("{}/q1/{}.parquet", output, chunk_id);
    //                         let mut file = File::create(file).unwrap();
    //                         ParquetWriter::new(&mut file)
    //                             .finish(&mut df)
    //                             .unwrap();
    //                     }
    //                 }))
    //             }
    //         },
    //         _ => {
    //             for _ in 0..8 {
    //                 let receiver: Receiver<_> = receiver.clone();
    //                 let output = output.clone();
    //                 handles.push(thread::spawn(move || {
    //                     while let Ok((chunk_id, records)) = receiver.recv() {
    //                         let mut read_idx_vec: Vec<u64> = Vec::new();
    //                         let mut chrom1_vec: Vec<String> = Vec::new();
    //                         let mut pos1_vec: Vec<u64> = Vec::new();
    //                         let mut chrom2_vec: Vec<String> = Vec::new();
    //                         let mut pos2_vec: Vec<u64> = Vec::new();
    //                         let mut strand1_vec: Vec<String> = Vec::new();
    //                         let mut strand2_vec: Vec<String> = Vec::new();
    //                         let mut mapq_vec: Vec<u8> = Vec::new();
    //                         let mut read_idx = 0;
    //                         for mut concatemer in records {
    //                             concatemer.sort();
    //                             for pair in concatemer.decompose() {
    //                                 let (record1, record2) = (pair[0], pair[1]);
    //                                 read_idx_vec.push(read_idx);
    //                                 chrom1_vec.push(record1.target.clone());
    //                                 pos1_vec.push(record1.target_start + record1.target_end / 2);
    //                                 chrom2_vec.push(record2.target.clone());
    //                                 pos2_vec.push(record2.target_start + record2.target_end / 2);
    //                                 strand1_vec.push(record1.query_strand.to_string());
    //                                 strand2_vec.push(record2.query_strand.to_string());
    //                                 mapq_vec.push(std::cmp::min(record1.mapq, record2.mapq));
                                    
    //                                 read_idx += 1;
    //                             }
        
    //                         }
        
    //                         let df = DataFrame::new(vec![
    //                             Series::new("read_idx".into(), read_idx_vec).into(),
    //                             Series::new("chrom1".into(), chrom1_vec).into(),
    //                             Series::new("pos1".into(), pos1_vec).into(),
    //                             Series::new("chrom2".into(), chrom2_vec).into(),
    //                             Series::new("pos2".into(), pos2_vec).into(),
    //                             Series::new("strand1".into(), strand1_vec).into(),
    //                             Series::new("strand2".into(), strand2_vec).into(),
    //                             Series::new("mapq".into(), mapq_vec).into(),
    //                         ]).unwrap();
        
    //                         let mut df = df.lazy().with_column(
    //                             col("read_idx").cast(DataType::String)
    //                         ).with_column(
    //                             col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).with_column(
    //                             col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).with_column(
    //                             col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).with_column(
    //                             col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
    //                         ).collect().unwrap();
        
    //                         let file = format!("{}/q0/{}.parquet", output, chunk_id);
    //                         let mut file = File::create(file).unwrap();
    //                         ParquetWriter::new(&mut file).finish(&mut df).unwrap();
                            
    //                         let mut df = df.lazy().filter(
    //                             col("mapq").gt_eq(1)
    //                         ).collect().unwrap();
                
                           
    //                         let file = format!("{}/q1/{}.parquet", output, chunk_id);
    //                         let mut file = File::create(file).unwrap();
    //                         ParquetWriter::new(&mut file)
    //                             .finish(&mut df)
    //                             .unwrap();
    //                     }
    //                 }))
    //             }
    //         }
    //     }
        
        
    //     let mut concatemer: Concatemer = Concatemer::new();
    //     let mut batch = Vec::with_capacity(chunksize);
    //     let mut first_iteration = true;
    //     let mut previous_read_idx: u64 = 0;
    //     let mut record_count = 0;
    //     let mut chunk_id = 0 as usize;
    //     let mut line_count = 0;
    //     let mut total_pair_count = 0;   
    //     for (i, line) in rdr.deserialize().enumerate() {
    //         let record: PoreCRecord = match line {
    //             Ok(v) => v,
    //             Err(error) => {
    //                 log::warn!("Could not parse line {}: {:?}", i + 1, error);
    //                 continue
    //             },
    //         }; 
            
    //         if !first_iteration && record.read_idx != previous_read_idx {
    //             let order = concatemer.count();
    //             if (order < max_order) && (order >= min_order) {
    //                 concatemer_summary.count(&concatemer);
    //                 batch.push(std::mem::take(&mut concatemer));
                    
    //                 record_count += order * (order - 1) / 2;
    //             } else {
    //                 concatemer.clear();
    //             }
                
    //         }
    //         first_iteration = false;
    //         previous_read_idx = record.read_idx;
            
    //         if record.mapq < min_quality {
    //             continue
    //         }
    //         line_count += 1;
    //         concatemer.push(record); 
    //         if record_count >= chunksize{
    //             total_pair_count += record_count;
    //             sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
    //             record_count = 0;
    //             chunk_id += 1;
    //             line_count = 0;
    //         }
            
    //     }
        
    //     if !batch.is_empty() {
    //         sender.send((chunk_id, batch)).unwrap();
    //     }
        
    //     drop(sender);
    //     for handle in handles {
    //         handle.join().unwrap();
    //     }
    
    //     log::info!("Processed total {} lines", line_count);
    //     log::info!("Generated total {} pairs", total_pair_count);
    //     log::info!("Successful output pairs `{}`", output);
        
    //     let output_prefix = if output == "-" {
    //         Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
    //     } else {
    //         Path::new(&output).with_extension("").to_str().unwrap().to_string()
    //     };

    //     concatemer_summary.save(&format!("{}.concatemer.summary", output_prefix));

    //     Ok(())

    // }

    pub fn to_pairs_pqs(&mut self, chromsizes: &String, output: &String, 
        chunksize: usize, min_quality: u8, 
        min_order: usize, max_order: usize, threads: usize
    ) -> anyResult<()> {
        use polars::prelude::*;
        use crate::pqs::_README as _readme;
        use crate::pqs::_METADATA;

        polars::enable_string_cache();

        let mut rdr = self.parse2().expect("Failed to open input file");
        log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
     
        let _ = std::fs::create_dir_all(output);
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));
        std::fs::copy(chromsizes, format!("{}/_contigsizes", output))?;

        let contigsizes = ChromSize::new(chromsizes);
        let contigsizes_data = contigsizes.to_vec().unwrap();
        let max_contig_size = contigsizes_data.iter().map(|x| x.size).max().unwrap_or(0);
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
        
        // let q0_total = Arc::new(AtomicU64::new(0));
        // let q1_total = Arc::new(AtomicU64::new(0));
        // let mut concatemer_summary: ConcatemerSummary = ConcatemerSummary::new();
        
        // #[inline]
        // fn mid_u64(a: u64, b: u64) -> u64 {
        //     a.saturating_add(b) / 2
        // }
        // #[inline]
        // fn mid_u32(a: u64, b: u64) -> u32 {
        //     (a.saturating_add(b) / 2) as u32
        // }
        // #[inline]
        // fn enc_strand(c: char) -> i8 {
        //     match c {
        //         '+' => 1,
        //         '-' => -1,
        //         _ => 0,
        //     }
        // }


        // let (sender, receiver) = bounded::<(usize, Vec<Concatemer>, u64)>(100);
        // let mut handles = vec![];

        // let num_workers = 8;
        // for _ in 0..num_workers {
        //     let rx = receiver.clone();
        //     let output = output.clone();
        //     let q0_total = Arc::clone(&q0_total);
        //     let q1_total = Arc::clone(&q1_total);

        //     handles.push(thread::spawn(move || {
        //         while let Ok((chunk_id, records, mut current_id)) = rx.recv() {
      
        //             let est_capacity = records.iter().map(|c| c.count() * (c.count() - 1) / 2).sum();
        //             let mut read_idx_vec = Vec::with_capacity(est_capacity);
        //             let mut chrom1_vec = Vec::with_capacity(est_capacity);
        //             let mut pos1_vec = Vec::with_capacity(est_capacity);
        //             let mut chrom2_vec = Vec::with_capacity(est_capacity);
        //             let mut pos2_vec = Vec::with_capacity(est_capacity);
        //             let mut strand1_vec = Vec::with_capacity(est_capacity);
        //             let mut strand2_vec = Vec::with_capacity(est_capacity);
        //             let mut mapq_vec = Vec::with_capacity(est_capacity);

        //             for mut concatemer in records {
        //                 concatemer.sort();
        //                 let recs = &concatemer.records;
        //                 let n = recs.len();
        //                 for i in 0..n {
        //                     for j in i + 1..n {
        //                         let r1 = &recs[i];
        //                         let r2 = &recs[j];
                                
        //                         current_id += 1;
        //                         read_idx_vec.push(current_id);
        //                         chrom1_vec.push(r1.target.clone());
        //                         chrom2_vec.push(r2.target.clone());
                                
                 
        //                         pos1_vec.push(((r1.target_start + r1.target_end) >> 1) as u64);
        //                         pos2_vec.push(((r2.target_start + r2.target_end) >> 1) as u64);
              
        //                         strand1_vec.push(if r1.query_strand == '+' { "+" } else { "-" });
        //                         strand2_vec.push(if r2.query_strand == '+' { "+" } else { "-" });
                                
        //                         mapq_vec.push(std::cmp::min(r1.mapq, r2.mapq));
        //                     }
        //                 }
        //             }

        //             let mut df = DataFrame::new(vec![
        //                 Series::new("read_idx".into(), read_idx_vec).into(),
        //                 Series::new("chrom1".into(), chrom1_vec).into(),
        //                 Series::new("pos1".into(), pos1_vec).into(),
        //                 Series::new("chrom2".into(), chrom2_vec).into(),
        //                 Series::new("pos2".into(), pos2_vec).into(),
        //                 Series::new("strand1".into(), strand1_vec).into(),
        //                 Series::new("strand2".into(), strand2_vec).into(),
        //                 Series::new("mapq".into(), mapq_vec).into(),
        //             ]).unwrap();

        //             let pos_dtype = if max_contig_size < 4294967295 { DataType::UInt32 } else { DataType::UInt64 };
                    
        //             let mut df = df.lazy()
        //                 .with_columns([
        //                     col("read_idx").cast(DataType::String),
        //                     col("pos1").cast(pos_dtype.clone()),
        //                     col("pos2").cast(pos_dtype),
        //                     col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
        //                     col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
        //                     col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
        //                     col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
        //                 ])
        //                 .collect().unwrap();

        //                 q0_total.fetch_add(df.height() as u64, AtomOrdering::Relaxed);
        //             let path0 = format!("{}/q0/{}.parquet", output, chunk_id);
        //             ParquetWriter::new(File::create(path0).unwrap()).finish(&mut df).unwrap();
                    
        //             let mut df_q1 = df.lazy().filter(col("mapq").gt_eq(1)).collect().unwrap();
        //             q1_total.fetch_add(df_q1.height() as u64, AtomOrdering::Relaxed);

        //             if df_q1.height() > 0 {
        //                 let path1 = format!("{}/q1/{}.parquet", output, chunk_id);
        //                 ParquetWriter::new(File::create(path1).unwrap()).finish(&mut df_q1).unwrap();
        //             }
        //         }
                  
        //     }));
        // }

        let q0_total = Arc::new(AtomicU64::new(0));
        let q1_total = Arc::new(AtomicU64::new(0));
        let concatemer_summary = Arc::new(Mutex::new(ConcatemerSummary::new()));

        let (sender, receiver) = bounded::<(usize, Vec<String>, u64)>(100);
        let mut handles = vec![];

        let num_workers = threads;
        for _ in 0..num_workers {
            let rx = receiver.clone();
            let output = output.clone();
            let q0_total = Arc::clone(&q0_total);
            let q1_total = Arc::clone(&q1_total);
            let summary_lock = Arc::clone(&concatemer_summary);

            handles.push(thread::spawn(move || {
                while let Ok((chunk_id, lines, mut current_id)) = rx.recv() {
                    let mut records: Vec<Concatemer> = Vec::with_capacity(lines.len() / 4);
                    let mut concatemer = Concatemer::new();
                    let mut previous_read_idx = u64::MAX;

                    for line in lines {
                        let trimmed = line.trim_end();
                        if trimmed.is_empty() || trimmed.starts_with('#') {
                            continue;
                        }

                        let mut parts = trimmed.split('\t');
                        let read_idx = parts.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);

                        if previous_read_idx != u64::MAX && read_idx != previous_read_idx {
                            let order = concatemer.count();
                            if order >= min_order && order < max_order {
                                records.push(std::mem::take(&mut concatemer));
                            } else {
                                concatemer.clear();
                            }
                        }

                        let q_len      = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
                        let q_start    = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
                        let q_end      = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
                        let q_strand   = parts.next().and_then(|s| s.chars().next()).unwrap_or('+');
                        let target     = parts.next().unwrap_or("").to_string();
                        let t_start    = parts.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
                        let t_end      = parts.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
                        let mapq       = parts.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);

                        concatemer.push(PoreCRecord {
                            read_idx, query_length: q_len, query_start: q_start, query_end: q_end,
                            query_strand: q_strand, target, target_start: t_start, target_end: t_end,
                            mapq, identity: 0.0, filter_reason: "".to_string()
                        });

                        previous_read_idx = read_idx;
                    }

                    let last_order = concatemer.count();
                    if last_order >= min_order && last_order < max_order {
                        records.push(concatemer);
                    }

                    let mut local_sum = HashMap::new();
                    for c in &records {
                        *local_sum.entry(c.count() as u32).or_insert(0) += 1;
                    }
                    {
                        let mut global = summary_lock.lock().unwrap();
                        for (k, v) in local_sum {
                            *global.summary.entry(k).or_insert(0) += v;
                        }
                    }

                    let est_capacity = records.iter().map(|c| c.count() * (c.count() - 1) / 2).sum();
                    if est_capacity == 0 { continue; }

                    let mut read_idx_vec = Vec::with_capacity(est_capacity);
                    let mut chrom1_vec = Vec::with_capacity(est_capacity);
                    let mut chrom2_vec = Vec::with_capacity(est_capacity);
                    
                    let use_u32 = max_contig_size < 4294967295;
                    let mut pos1_u32 = Vec::with_capacity(if use_u32 { est_capacity } else { 0 });
                    let mut pos2_u32 = Vec::with_capacity(if use_u32 { est_capacity } else { 0 });
                    let mut pos1_u64 = Vec::with_capacity(if use_u32 { 0 } else { est_capacity });
                    let mut pos2_u64 = Vec::with_capacity(if use_u32 { 0 } else { est_capacity });

                    let mut strand1_vec = Vec::with_capacity(est_capacity);
                    let mut strand2_vec = Vec::with_capacity(est_capacity);
                    let mut mapq_vec = Vec::with_capacity(est_capacity);

                    for concatemer in &mut records {
                        concatemer.sort();
                        let recs = &concatemer.records;
                        let n = recs.len();
                        for i in 0..n {
                            let r1 = &recs[i];
                            let r1_pos = (r1.target_start + r1.target_end) >> 1;
                            let r1_strand = if r1.query_strand == '+' { "+" } else { "-" };

                            for j in i + 1..n {
                                let r2 = &recs[j];
                                let r2_pos = (r2.target_start + r2.target_end) >> 1;
                                
                                current_id += 1;
                                read_idx_vec.push(current_id); 
                                
                                chrom1_vec.push(r1.target.as_str());
                                chrom2_vec.push(r2.target.as_str());
                                
                                if use_u32 {
                                    pos1_u32.push(r1_pos as u32);
                                    pos2_u32.push(r2_pos as u32);
                                } else {
                                    pos1_u64.push(r1_pos as u64);
                                    pos2_u64.push(r2_pos as u64);
                                }
              
                                strand1_vec.push(r1_strand);
                                strand2_vec.push(if r2.query_strand == '+' { "+" } else { "-" });
                                
                                mapq_vec.push(std::cmp::min(r1.mapq, r2.mapq));
                            }
                        }
                    }

                    let pos1_series = if use_u32 {
                        Series::new("pos1".into(), pos1_u32)
                    } else {
                        Series::new("pos1".into(), pos1_u64)
                    };
                    let pos2_series = if use_u32 {
                        Series::new("pos2".into(), pos2_u32)
                    } else {
                        Series::new("pos2".into(), pos2_u64)
                    };

                    let s_read_idx = Series::new("read_idx".into(), read_idx_vec)
                        .cast(&DataType::String)
                        .unwrap();
                    let s_chrom1 = Series::new("chrom1".into(), chrom1_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_chrom2 = Series::new("chrom2".into(), chrom2_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_strand1 = Series::new("strand1".into(), strand1_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_strand2 = Series::new("strand2".into(), strand2_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_mapq = Series::new("mapq".into(), mapq_vec);

                    let mut df = DataFrame::new(vec![
                        s_read_idx.into(),
                        s_chrom1.into(),
                        pos1_series.into(),
                        s_chrom2.into(),
                        pos2_series.into(),
                        s_strand1.into(),
                        s_strand2.into(),
                        s_mapq.into(),
                    ]).unwrap();

                    q0_total.fetch_add(df.height() as u64, AtomOrdering::Relaxed);
                    let path0 = format!("{}/q0/{}.parquet", output, chunk_id);
                    ParquetWriter::new(File::create(path0).unwrap()).finish(&mut df).unwrap();
                    
                    let mapq_col = df.column("mapq".into()).unwrap().as_materialized_series().u8().unwrap();
                    let mask = mapq_col.gt_eq(1);
                    let mut df_q1 = df.filter(&mask).unwrap();
                    
                    q1_total.fetch_add(df_q1.height() as u64, AtomOrdering::Relaxed);

                    if df_q1.height() > 0 {
                        let path1 = format!("{}/q1/{}.parquet", output, chunk_id);
                        ParquetWriter::new(File::create(path1).unwrap()).finish(&mut df_q1).unwrap();
                    }
                }
            }));
        }

        let mut line_buf = String::new();
        let mut chunk_lines = Vec::with_capacity(chunksize * 2);
        let mut previous_read_idx: u64 = 0;
        let mut current_chunk_lines_len = 0;
        let mut global_pair_offset = 0u64;
        let mut chunk_id = 0;
        let mut first_iteration = true;

        while rdr.read_line(&mut line_buf)? > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let mut read_idx = 0u64;
            if let Some(pos) = trimmed.find('\t') {
                if let Ok(idx) = trimmed[..pos].parse::<u64>() {
                    read_idx = idx;
                }
            }

            if !first_iteration && read_idx != previous_read_idx {
                if current_chunk_lines_len >= chunksize {
                    sender.send((chunk_id, std::mem::take(&mut chunk_lines), global_pair_offset)).unwrap();
                    global_pair_offset += current_chunk_lines_len as u64;
                    chunk_id += 1;
                    current_chunk_lines_len = 0;
                    chunk_lines = Vec::with_capacity(chunksize * 2);
                }
            }

            first_iteration = false;
            previous_read_idx = read_idx;
            current_chunk_lines_len += 1;
            chunk_lines.push(std::mem::take(&mut line_buf)); // 利用 take 免内存拷贝
        }

        if !chunk_lines.is_empty() {
            sender.send((chunk_id, chunk_lines, global_pair_offset)).unwrap();
        }
        
        drop(sender);
        for handle in handles { handle.join().unwrap(); }

        {
            let q0_n = q0_total.load(AtomOrdering::Relaxed);
            let q1_n = q1_total.load(AtomOrdering::Relaxed);

            let mut wtr = common_writer(format!("{}/_metadata_counts", output).as_str());
            writeln!(wtr, "q0_records\t{}", q0_n).unwrap();
            writeln!(wtr, "q1_records\t{}", q1_n).unwrap();
            wtr.flush().unwrap();
        }

        let output_prefix = if output == "-" {
            Path::new(&self.file).with_extension("").to_str().unwrap().to_string()
        } else {
            Path::new(&output).with_extension("").to_str().unwrap().to_string()
        };
        
        let final_summary = concatemer_summary.lock().unwrap();
        final_summary.save(&format!("{}.concatemer.summary", output_prefix));

        Ok(())
    }


    // pub fn to_pairs_pqs(&mut self, chromsizes: &String, output: &String, 
    //     chunksize: usize, min_quality: u8, 
    //     min_order: usize, max_order: usize,
    // ) -> anyResult<()> {
    //     use polars::prelude::*;
    //     use crate::pqs::_README as _readme;
    //     use crate::pqs::_METADATA;

    //     polars::enable_string_cache();

    //     let mut rdr = self.parse2().expect("Failed to open input file");
    //     log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
        
    //     let _ = std::fs::create_dir_all(output);
    //     let _ = std::fs::create_dir_all(format!("{}/q0", output));
    //     let _ = std::fs::create_dir_all(format!("{}/q1", output));
    //     std::fs::copy(chromsizes, format!("{}/_contigsizes", output))?;

    //     let contigsizes = ChromSize::new(chromsizes);
    //     let contigsizes_data = contigsizes.to_vec().unwrap();
    //     let max_contig_size = contigsizes_data.iter().map(|x| x.size).max().unwrap_or(0);
        
    //     let pos_dtype = if max_contig_size < 4294967295 { DataType::UInt32 } else { DataType::UInt64 };
    //     let pos_type_string = if max_contig_size < 4294967295 { "UInt32" } else { "UInt64" };

    //     let mut wtr = common_writer(format!("{}/_readme", output).as_str());
    //     wtr.write_all(_readme.as_bytes()).unwrap();
    //     wtr.flush().unwrap();

    //     let create_date_time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

    //     let mut wtr = common_writer(format!("{}/_metadata", output).as_str());
    //     let mut _metadata = _METADATA.to_string();
    //     _metadata = _metadata.replace("REPLACE", &create_date_time);
    //     _metadata = _metadata.replace("CHUNKSIZE", &chunksize.to_string());
    //     _metadata = _metadata.replace("pos_type_lower", pos_type_string.to_lowercase().as_str());
    //     _metadata = _metadata.replace("pos_type", pos_type_string);
    //     wtr.write_all(_metadata.as_bytes()).unwrap();
    //     wtr.flush().unwrap();
        
    //     let mut concatemer_summary: ConcatemerSummary = ConcatemerSummary::new();

    //     let (sender, receiver) = bounded::<(usize, Vec<Concatemer>, u64)>(64); 
    //     let mut handles = vec![];

    //     let num_workers = 12; 
    //     for _ in 0..num_workers {
    //         let rx = receiver.clone();
    //         let output = output.clone();
    //         let pos_dtype = pos_dtype.clone();

    //         handles.push(thread::spawn(move || {
                
    //                 let est_capacity = 1_000;
                
    //                 let mut read_idx_vec = Vec::with_capacity(est_capacity);
    //                 let mut chrom1_vec = Vec::with_capacity(est_capacity);
    //                 let mut pos1_vec = Vec::with_capacity(est_capacity);
    //                 let mut chrom2_vec = Vec::with_capacity(est_capacity);
    //                 let mut pos2_vec = Vec::with_capacity(est_capacity);
    //                 let mut strand1_vec = Vec::with_capacity(est_capacity);
    //                 let mut strand2_vec = Vec::with_capacity(est_capacity);
    //                 let mut mapq_vec = Vec::with_capacity(est_capacity);
    //             while let Ok((chunk_id, records, mut current_id)) = rx.recv() {
    //                 read_idx_vec.clear();
    //                 chrom1_vec.clear();
    //                 pos1_vec.clear();
    //                 chrom2_vec.clear();
    //                 pos2_vec.clear();
    //                 strand1_vec.clear();
    //                 strand2_vec.clear();
    //                 mapq_vec.clear();

    //                 for mut concatemer in records {
    //                     // In-place sort is cheaper than clone sort
    //                     concatemer.sort();
    //                     let recs = &concatemer.records;
    //                     let n = recs.len();
                        
    //                     for i in 0..n {
    //                         let r1 = &recs[i]; 
            
    //                         let r1_pos = ((r1.target_start + r1.target_end) >> 1) as u64;
    //                         let r1_s = if r1.query_strand == '+' { "+" } else { "-" };
    //                         let r1_mq = r1.mapq;

    //                         for j in i + 1..n {
    //                             let r2 = &recs[j];
                                
    //                             current_id += 1;
                                
    //                             read_idx_vec.push(current_id);
                            
    //                             chrom1_vec.push(r1.target.clone()); 
    //                             chrom2_vec.push(r2.target.clone());


    //                             pos1_vec.push(r1_pos);
    //                             pos2_vec.push(((r2.target_start + r2.target_end) >> 1) as u64);
            
    //                             strand1_vec.push(r1_s);
    //                             strand2_vec.push(if r2.query_strand == '+' { "+" } else { "-" });
                                
    //                             mapq_vec.push(std::cmp::min(r1.mapq, r2.mapq));
    //                         }
    //                     }
    //                 }
    //                 if read_idx_vec.is_empty() { continue; }
                 
    //                 let s_idx = Series::new("read_idx".into(), &read_idx_vec);
    //                 let s_idx = s_idx.cast(&DataType::String).unwrap(); 

    //                 let df = DataFrame::new(vec![
    //                     s_idx.into(), 
    //                     Series::new("chrom1".into(), &chrom1_vec).into(), 
    //                     Series::new("pos1".into(), &pos1_vec).into(), 
    //                     Series::new("chrom2".into(), &chrom2_vec).into(), 
    //                     Series::new("pos2".into(), &pos2_vec).into(), 
    //                     Series::new("strand1".into(), &strand1_vec).into(), 
    //                     Series::new("strand2".into(), &strand2_vec).into(), 
    //                     Series::new("mapq".into(), &mapq_vec).into()
    //                 ]).unwrap();

                   
    //                 let mut df_final = df.lazy()
    //                     .with_columns([
    //                         col("pos1").cast(pos_dtype.clone()),
    //                         col("pos2").cast(pos_dtype.clone()),
    //                         col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
    //                         col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
    //                         col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
    //                         col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)),
    //                     ])
    //                     .collect().unwrap();

                    
    //                 let path0 = format!("{}/q0/{}.parquet", output, chunk_id);
    //                 let file0 = std::fs::File::create(path0).unwrap();
                
    //                 ParquetWriter::new(file0).finish(&mut df_final).unwrap();
                    
    //                 let mut df_q1 = df_final.lazy()
    //                     .filter(col("mapq").gt_eq(1))
    //                     .collect().unwrap();
                    
    //                 if df_q1.height() > 0 {
    //                     let path1 = format!("{}/q1/{}.parquet", output, chunk_id);
    //                     let file1 = std::fs::File::create(path1).unwrap();
    //                     ParquetWriter::new(file1).finish(&mut df_q1).unwrap();
    //                 }
    //             }
    //         }));
    //     }

    //     let mut line_buf = String::new();
    //     let mut concatemer = Concatemer::new();

    //     let mut batch = Vec::with_capacity(chunksize + 100); 
        
    //     let mut first_iteration = true;
    //     let mut previous_read_idx: u64 = 0;
    //     let mut current_chunk_pairs = 0;
    //     let mut global_pair_offset = 0u64;
    //     let mut chunk_id = 0;

    //     while rdr.read_line(&mut line_buf)? > 0 {
    //         let trimmed = line_buf.trim_end();
    //         if trimmed.is_empty() { line_buf.clear(); continue; }
    //         if trimmed.as_bytes()[0] == b'#' { line_buf.clear(); continue; }

    //         let mut parts = trimmed.split('\t');

    //         let read_idx_str = parts.next().unwrap_or("0");
    //         let read_idx = match read_idx_str.parse::<u64>() {
    //             Ok(v) => v,
    //             Err(_) => { line_buf.clear(); continue; }
    //         };

    //         if !first_iteration && read_idx != previous_read_idx {
    //             let order = concatemer.count();
    //             if order >= min_order && order < max_order {
    //                 concatemer_summary.count(&concatemer);
    //                 let pairs = (order * (order - 1) / 2) as usize;
    //                 current_chunk_pairs += pairs;
    //                 batch.push(std::mem::take(&mut concatemer));
                    
    //                 if current_chunk_pairs >= chunksize {
            
    //                     sender.send((chunk_id, std::mem::take(&mut batch), global_pair_offset)).unwrap();
                        
    //                     global_pair_offset += current_chunk_pairs as u64;
    //                     current_chunk_pairs = 0;
    //                     chunk_id += 1;
                     
    //                     batch = Vec::with_capacity(chunksize / 2);
    //                 }
    //             } else {
    //                 concatemer.clear();
    //             }
    //         }


    //         let q_len = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
    //         let q_start = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
    //         let q_end = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
    //         let q_strand = parts.next().map(|s| s.as_bytes().get(0).copied().unwrap_or(b'+') as char).unwrap_or('+');
    //         let target = parts.next().unwrap_or("").to_string(); 
    //         let t_start = parts.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
    //         let t_end = parts.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
    //         let mapq = parts.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);

    //         if mapq >= min_quality {
    //             concatemer.push(PoreCRecord {
    //                 read_idx, query_length: q_len, query_start: q_start,
    //                 query_end: q_end, query_strand: q_strand, target,
    //                 target_start: t_start, target_end: t_end, mapq,
    //                 identity: 0.0, filter_reason: "".to_string(),
    //             });
    //         }

    //         first_iteration = false;
    //         previous_read_idx = read_idx;
    //         line_buf.clear();
    //     }

        
    //     if !batch.is_empty() {
    //             let order = concatemer.count();
    //             if order >= min_order && order < max_order {
    //                 concatemer_summary.count(&concatemer);
    //                 batch.push(concatemer);
    //             }
            
    //         if !batch.is_empty() {
    //             sender.send((chunk_id, batch, global_pair_offset)).unwrap();
    //         }
    //     } else {
               
    //             let order = concatemer.count();
    //             if order >= min_order && order < max_order {
    //                 concatemer_summary.count(&concatemer);
    //                 let b = vec![concatemer];
    //                 sender.send((chunk_id, b, global_pair_offset)).unwrap();
    //             }
    //     }
        
    //     drop(sender);
    //     for handle in handles { handle.join().unwrap(); }
        
    //     let output_prefix = Path::new(output).with_extension(""); 
        
    //     concatemer_summary.save(&format!("{}/concatemer.summary", output));

    //     Ok(())
    // }

    // pub fn to_pairs(&self, chromsizes: &String, output: &String, min_quality: u8, 
    //                     min_order: usize, max_order: usize,
    //                 ) -> Result<(), Box<dyn Error>> {
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
        
    //     let mut wtr = csv::WriterBuilder::new()
    //                         .has_headers(false)
    //                         .delimiter(b'\t')
    //                         .from_writer(writer);

    //     let mut concatemer: Concatemer = Concatemer::new();  
    //     let mut concatemer_summary: ConcatemerSummary = ConcatemerSummary::new();

    //     let mut old_read_idx: u64 = 0; 

    //     let mut read_id: u64 = 0;
        
    //     let mut first_iteration = true;
    //     for (i, line) in rdr.deserialize().enumerate() {
    //         let record: PoreCRecord = match line {
    //             Ok(v) => v,
    //             Err(error) => {
    //                 log::warn!("Could not parse line {}", i + 1);
    //                 continue
    //             },
    //         };
    //         if !first_iteration && record.read_idx != old_read_idx {
    //             let order = concatemer.count();
    //             if (order < max_order) && (order >= min_order) {
    //                 concatemer.sort();
    //                 concatemer_summary.count(&concatemer);
    //                 for pair in concatemer.decompose() {
    //                     wtr.serialize(PairRecord::from_pore_c_pair(pair, read_id)).unwrap();
    //                     read_id += 1;
    //                 }
                      
    //             }
                
    //             concatemer.clear();
    //         }
    //         first_iteration = false;
    //         old_read_idx = record.read_idx;
    //         if record.mapq < min_quality {
    //             continue
    //         }
            
    //         concatemer.push(record);
           
    //     }

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
    //             for record in batch {
    //                 wtr.serialize(record).unwrap();
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

    // pub fn to_pairs(&self, chromsizes: &String, output: &String, min_quality: u8, 
    //                     min_order: usize, max_order: usize,
    //                 ) -> Result<(), Box<dyn Error>> {
    //     let parse_result = self.parse();
    //     let mut rdr = match parse_result {
    //         Ok(v) => v,
    //         Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
    //     };
    //     log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
        
    //     let chromsizes_obj = ChromSize::new(chromsizes);
    //     let chromsizes_data = chromsizes_obj.to_vec().unwrap();
    //     let mut ph = PairHeader::new();
    //     ph.from_chromsizes(chromsizes_data);

    //     let writer = common_writer(output);
    //     let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, writer);
    //     writer.write_all(ph.to_string().as_bytes()).unwrap();
        

    //     let (sender, receiver) = bounded::<(usize, Vec<Concatemer>)>(200);
    //     let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

    //     let num_workers = rayon::current_num_threads();
    //     for _ in 0..num_workers {
    //         let rx = receiver.clone();
    //         let tx = out_sender.clone();
    //         thread::spawn(move || {
    //             let mut local_buf = Vec::with_capacity(1024 * 1024);
    //             while let Ok((chunk_id, batch)) = rx.recv() {
    //                 local_buf.clear();
    //                 for mut concatemer in batch {
    //                     concatemer.sort();
    //                     for pair in concatemer.decompose() {
    //                         let (r1, r2) = (pair[0], pair[1]);
    //                         let line = format!(
    //                             ".\t{}\t{}\t{}\t{}\t{}\t{}\t.\t.\n",
    //                             r1.target, r1.target_start, r2.target, r2.target_start,
    //                             r1.query_strand, r2.query_strand
    //                         );
    //                         local_buf.extend_from_slice(line.as_bytes());
    //                     }
    //                 }
    //                 tx.send((chunk_id, local_buf.clone())).unwrap();
    //             }
    //         });
    //     }
    //     drop(out_sender);

    //     let write_handle = thread::spawn(move || {
    //         let mut pending = BTreeMap::new();
    //         let mut next_chunk = 0;
    //         while let Ok((chunk_id, data)) = out_receiver.recv() {
    //             pending.insert(chunk_id, data);
    //             while let Some(data) = pending.remove(&next_chunk) {
    //                 writer.write_all(&data).unwrap();
    //                 next_chunk += 1;
    //             }
    //         }
    //         writer.flush().unwrap();
    //     });

    //     let mut concatemer = Concatemer::new();
    //     let mut batch = Vec::with_capacity(5000);
    //     let mut old_read_idx: u64 = 0;
    //     let mut first_iteration = true;
    //     let mut chunk_id = 0;
    //     let mut concatemer_summary = ConcatemerSummary::new();

    //     for (i, line) in rdr.deserialize().enumerate() {
    //         let record: PoreCRecord = match line {
    //             Ok(v) => v,
    //             Err(_) => continue,
    //         };

    //         if !first_iteration && record.read_idx != old_read_idx {
    //             let order = concatemer.count();
    //             if order >= min_order && order < max_order {
    //                 concatemer_summary.count(&concatemer);
    //                 batch.push(std::mem::take(&mut concatemer));
    //                 if batch.len() >= 5000 {
    //                     sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
    //                     chunk_id += 1;
    //                 }
    //             } else {
    //                 concatemer.clear();
    //             }
    //         }
    //         first_iteration = false;
    //         old_read_idx = record.read_idx;
    //         if record.mapq >= min_quality {
    //             concatemer.push(record);
    //         }
    //     }

    //     if !batch.is_empty() || concatemer.count() > 0 {
    //         if concatemer.count() >= min_order && concatemer.count() < max_order {
    //             batch.push(concatemer);
    //         }
    //         sender.send((chunk_id, batch)).unwrap();
    //     }

    //     drop(sender);
    //     write_handle.join().unwrap();

    //     let output_prefix = Path::new(output).with_extension("");
    //     concatemer_summary.save(&format!("{}.concatemer.summary", output_prefix.to_str().unwrap()));
        
    //     log::info!("Successful output pairs `{}`", output);
    //     Ok(())
    // }

    pub fn to_pairs(&mut self, chromsizes: &String, output: &String, min_quality: u8, 
                        min_order: usize, max_order: usize,
                    ) -> Result<(), Box<dyn Error>> {
        let mut rdr = self.parse2().expect("Failed to open input file");
        log::info!("Only retain concatemer that order in the range of [{}, {})", min_order, max_order);
        
        let chromsizes_obj = ChromSize::new(chromsizes);
        let chromsizes_data = chromsizes_obj.to_vec().unwrap();
        let mut ph = PairHeader::new();
        ph.from_chromsizes(chromsizes_data);

        let writer = common_writer(output);
        let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, writer);
        writer.write_all(ph.to_string().as_bytes()).unwrap();

        let (sender, receiver) = bounded::<(usize, Vec<Concatemer>, u64)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let num_workers = rayon::current_num_threads();
        for _ in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            thread::spawn(move || {
                let mut local_buf = Vec::with_capacity(1024 * 1024);
                while let Ok((chunk_id, batch, mut current_id)) = rx.recv() {
                    local_buf.clear();
                    for mut concatemer in batch {
                        concatemer.sort();
                        // for pair in concatemer.decompose() {
                        //     read_idx += 1;
                        //     let (r1, r2) = (pair[0], pair[1]);
                        //     let mapq = std::cmp::min(r1.mapq, r2.mapq);
                        //     let pos1 = (r1.target_start + r1.target_end) / 2;
                        //     let pos2 = (r2.target_start + r2.target_end) / 2;
                        //     let line = format!(
                        //         "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        //         read_idx,
                        //         r1.target, pos1, r2.target, pos2,
                        //         r1.query_strand, r2.query_strand, mapq
                        //     );
                        //     local_buf.extend_from_slice(line.as_bytes());
                        // }
                        let records = &concatemer.records;
                        let n = records.len();
                        for i in 0..n {
                            for j in i + 1..n {
                                current_id += 1;
                                let r1 = &records[i];
                                let r2 = &records[j];
                                let mapq = std::cmp::min(r1.mapq, r2.mapq);
                                let pos1 = (r1.target_start + r1.target_end) / 2;
                                let pos2 = (r2.target_start + r2.target_end) / 2;
                                
                                let line = format!(
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                    current_id,
                                    r1.target, pos1, r2.target, pos2,
                                    r1.query_strand, r2.query_strand, mapq
                                );
                                local_buf.extend_from_slice(line.as_bytes());
                            }
                        }
                      
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            });
        }
        drop(out_sender);

        let mut global_pair_counter: u64 = 0;
        let write_handle = thread::spawn(move || {
            let mut pending = BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    writer.write_all(&data).unwrap();
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut line_buf = String::new();
        let mut concatemer = Concatemer::new();
        let mut batch = Vec::with_capacity(5000);
        let mut old_read_idx: u64 = 0;
        let mut first_iteration = true;
        let mut chunk_id = 0;
        let mut concatemer_summary = ConcatemerSummary::new();

        let mut global_pair_counter: u64 = 0;
        while rdr.read_line(&mut line_buf)? > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let mut parts = trimmed.split('\t');
            
            let read_idx = match parts.next() {
                Some(s) => s.parse::<u64>().unwrap_or(0),
                None => { line_buf.clear(); continue; }
            };

            if !first_iteration && read_idx != old_read_idx {
                let order = concatemer.count();
                if order >= min_order && order < max_order {
                    concatemer_summary.count(&concatemer);
                    batch.push(std::mem::take(&mut concatemer));
                    if batch.len() >= 5000 {
                        let pairs_in_batch: u64 = batch.iter()
                            .map(|c| { let n = c.count() as u64; n * (n - 1) / 2 })
                            .sum();
                        sender.send((chunk_id, std::mem::take(&mut batch), global_pair_counter)).unwrap();
                        
                        global_pair_counter += pairs_in_batch;
                        chunk_id += 1;
                    }
                } else {
                    concatemer.clear();
                }
            }

            let q_len = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let q_start = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let q_end = parts.next().and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
            let q_strand = parts.next().and_then(|s| s.chars().next()).unwrap_or('+');
            let target_str = parts.next().unwrap_or("");
            let t_start = parts.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            let t_end = parts.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);
            let mapq = parts.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
            let identity = parts.next().and_then(|s| s.parse::<f32>().ok()).unwrap_or(0.0);
            let filter_reason = parts.next().unwrap_or("");

            first_iteration = false;
            old_read_idx = read_idx;

            if mapq >= min_quality {
                concatemer.push(PoreCRecord {
                    read_idx,
                    query_length: q_len,
                    query_start: q_start,
                    query_end: q_end,
                    query_strand: q_strand,
                    target: target_str.to_string(),
                    target_start: t_start,
                    target_end: t_end,
                    mapq,
                    identity,
                    filter_reason: filter_reason.to_string(),
                });
            }

            line_buf.clear();
        }

        if !batch.is_empty() {
            sender.send((chunk_id, batch, global_pair_counter)).unwrap();
        }

        drop(sender);
        write_handle.join().unwrap();

        let output_prefix = Path::new(output).with_extension("");
        concatemer_summary.save(&format!("{}.concatemer.summary", output_prefix.to_str().unwrap()));
        
        log::info!("Successful output pairs `{}`", output);
        Ok(())
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
        
        for _ in 0..10 {
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

    pub fn intersect_multi_threads_coitree(&mut self, hcr_bed: &String, invert: bool, output: &String) {
        log::info!("Building COI-Trees from BED...");
        let bed = Bed3::new(hcr_bed);
        let lapper_map = bed.to_interval_hash();

        let mut tree_map: HashMap<String, COITree<u8, u32>> = HashMap::new();
        for (chrom, lapper) in lapper_map {
            let nodes: Vec<IntervalNode<u8, u32>> = lapper.intervals.into_iter()
                .map(|iv| {
                    IntervalNode::new(iv.start as i32, iv.stop as i32, iv.val)
                })
                .collect();
            
            tree_map.insert(chrom, COITree::new(&nodes));
        }
        let interval_hash = Arc::new(tree_map);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let num_workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
        log::info!("Intersecting Bed in parallel using {} workers...", num_workers);

        let mut handles = vec![];
        for _ in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let ih = Arc::clone(&interval_hash);

            handles.push(thread::spawn(move || {
                let mut local_buf = Vec::with_capacity(2 * 1024 * 1024);
                let mut tab_indices = [0usize; 8];

                while let Ok((chunk_id, records)) = rx.recv() {
                    local_buf.clear();
                    for record in records {
                        let bytes = record.as_bytes();
                        let len = bytes.len();
                        
                        let mut tab_count = 0;
                        for idx in 0..len {
                            if bytes[idx] == b'\t' {
                                if tab_count < 8 {
                                    tab_indices[tab_count] = idx;
                                    tab_count += 1;
                                } else {
                                    break;
                                }
                            }
                        }

                        if tab_count < 8 { continue; } 

                        let target = &record[tab_indices[4] + 1..tab_indices[5]];
                        let start_str = &record[tab_indices[5] + 1..tab_indices[6]];
                        let end_str = &record[tab_indices[6] + 1..tab_indices[7]];

                        let target_start = match start_str.parse::<i32>() {
                            Ok(val) => val,
                            Err(_) => continue,
                        };
                        let target_end = match end_str.parse::<i32>() {
                            Ok(val) => val,
                            Err(_) => continue,
                        };

                        let is_in_regions = ih.get(target).map_or(false, |tree| {
                            let mut has_overlap = false;
                            tree.query(target_start, target_end, |_| {
                                has_overlap = true;
                            });
                            has_overlap
                        });

                        if is_in_regions ^ invert {
                            local_buf.extend_from_slice(record.as_bytes());
                            local_buf.push(b'\n');
                        }
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender);

        let out_path = output.clone();
        let write_handle = thread::spawn(move || {
            let wtr_file = common_writer(&out_path);
            let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, wtr_file);
            let mut pending = BTreeMap::new();
            let mut next_chunk = 0;

            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    if !data.is_empty() {
                        writer.write_all(&data).unwrap();
                    }
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut rdr = self.parse2().expect("Failed to open porec table for reading");
        let batch_size = 10_000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;
        let mut line_buf = String::new();

        while rdr.read_line(&mut line_buf).unwrap_or(0) > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }
            batch.push(std::mem::take(&mut line_buf));
            if batch.len() >= batch_size {
                sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                batch = Vec::with_capacity(batch_size);
                chunk_id += 1;
            }
        }
        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles { h.join().unwrap(); }
        write_handle.join().unwrap();

        log::info!("Successful output intersection porec table into `{}`", output);
    }


    pub fn break_contigs(&mut self, break_bed: &String, output: &String, threads: usize) {
        type IvString = Interval<usize, String>;
        let bed = Bed4::new(break_bed);
   
        let interval_hash = Arc::new(bed.to_interval_hash());

        let wtr_file = common_writer(output);
        let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, wtr_file);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let num_workers = threads;
        let mut handles = vec![];

        for _ in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let ih = Arc::clone(&interval_hash);

            handles.push(thread::spawn(move || {
                let mut local_buf = Vec::with_capacity(1024 * 1024);
                while let Ok((chunk_id, batch)) = rx.recv() {
                    local_buf.clear();
                    for line in batch {
                        let trimmed = line.trim_end();
                        let fields: Vec<&str> = trimmed.split('\t').collect();
                        
                        if fields.len() < 8 {
                            local_buf.extend_from_slice(line.as_bytes());
                            local_buf.push(b'\n');
                            continue;
                        }

                        let target = fields[5];
                        let mut processed = false;

                        if let Some(intervals) = ih.get(target) {
                            if let (Ok(t_start), Ok(t_end)) = (fields[6].parse::<usize>(), fields[7].parse::<usize>()) {
                                if let Some(iv) = intervals.find(t_start, t_end).next() {
                                    let break_contig_length = iv.stop - iv.start + 1;
                                    
                                    if t_start >= iv.start {
                                        let new_target_start = t_start - iv.start + 1;
                                        let new_target_end = t_end - iv.start + 1;

                                        if new_target_end <= break_contig_length {
                                            for i in 0..5 {
                                                local_buf.extend_from_slice(fields[i].as_bytes());
                                                local_buf.push(b'\t');
                                            }
                                            local_buf.extend_from_slice(iv.val.as_bytes());
                                            local_buf.push(b'\t');
                           
                                            let _ = write!(local_buf, "{}\t{}\t", new_target_start, new_target_end);
                                            for i in 8..=10 {
                                                if i < fields.len() {
                                                    local_buf.extend_from_slice(fields[i].as_bytes());
                                                }
                                                if i < 10 { local_buf.push(b'\t'); }
                                            }
                                            local_buf.push(b'\n');
                                            processed = true;
                                        }
                                    }
                                }
                            }
                        }

                        if !processed {
                            local_buf.extend_from_slice(line.as_bytes());
                            local_buf.push(b'\n');
                        }
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender); 

        let write_handle = thread::spawn(move || {
            let mut pending = std::collections::BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    writer.write_all(&data).unwrap();
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut rdr = self.parse2().expect("Failed to open table for reading");
        let batch_size = 5000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;

        for line_res in rdr.lines() {
            if let Ok(line) = line_res {
                batch.push(line);
                if batch.len() >= batch_size {
                    sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                    batch = Vec::with_capacity(batch_size);
                    chunk_id += 1;
                }
            }
        }
        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles { h.join().unwrap(); }
        write_handle.join().unwrap();

        log::info!("Successful output contigs corrected porec table into `{}`", output);
    }

    pub fn chr_porec_to_contig_porec(
        &mut self,
        contig_bed: &String,
        output: &String,
        threads: usize,
    ) -> anyResult<()> {
        #[derive(Debug, Clone)]
        struct ContigInterval {
            start: u64,
            end: u64,
            contig: String,
        }

        log::info!("Loading contig mapping BED from `{}`...", contig_bed);
        let f = common_reader(contig_bed);
        let rdr = std::io::BufReader::new(f);
        let mut chrom_to_intervals: HashMap<String, Vec<ContigInterval>> = HashMap::new();

        for line in rdr.lines().flatten() {
            let s = line.trim();
            if s.is_empty() || s.starts_with('#') { continue; }
            let fields: Vec<&str> = s.split_whitespace().collect();
            if fields.len() < 4 { continue; }
            let chrom = fields[0].to_string();
            let start: u64 = fields[1].parse().unwrap_or(0);
            let end: u64 = fields[2].parse().unwrap_or(0);
            let contig = fields[3].to_string();

            chrom_to_intervals.entry(chrom)
                .or_default()
                .push(ContigInterval { start, end, contig });
        }

        for ivs in chrom_to_intervals.values_mut() {
            ivs.sort_by_key(|iv| iv.start);
        }

        let chrom_intervals = Arc::new(chrom_to_intervals);
        let wtr_file = common_writer(output);
        let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, wtr_file);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        log::info!("Converting chromosome-level Pore-C alignments to contig-level in parallel...");
        let mut handles = vec![];

        for _ in 0..threads {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let ci = Arc::clone(&chrom_intervals);

            handles.push(thread::spawn(move || {
                let mut local_buf = Vec::with_capacity(1024 * 1024);
                while let Ok((chunk_id, batch)) = rx.recv() {
                    local_buf.clear();
                    for line in batch {
                        let trimmed = line.trim_end();
                        let fields: Vec<&str> = trimmed.split('\t').collect();
                        if fields.len() < 8 {
                            continue;
                        }

                        let target = fields[5];
                        let t_start: u64 = fields[6].parse().unwrap_or(0);
                        let t_end: u64 = fields[7].parse().unwrap_or(0);

                        if t_start == 0 || t_end == 0 || t_start > t_end { continue; }

                        let mid = (t_start + t_end) / 2;
                        let mid_0based = match mid.checked_sub(1) {
                            Some(v) => v,
                            None => continue,
                        };

                        let mut mapped = false;
                        if let Some(ivs) = ci.get(target) {
                            let idx_res = ivs.binary_search_by(|iv| {
                                if mid_0based < iv.start {
                                    std::cmp::Ordering::Greater
                                } else if mid_0based >= iv.end {
                                    std::cmp::Ordering::Less
                                } else {
                                    std::cmp::Ordering::Equal
                                }
                            });

                            if let Ok(idx) = idx_res {
                                let iv = &ivs[idx];
                                if t_start >= iv.start && t_end <= iv.end {
                                    let new_start = t_start - iv.start;
                                    let new_end = t_end - iv.start;

                                    for i in 0..5 {
                                        local_buf.extend_from_slice(fields[i].as_bytes());
                                        local_buf.push(b'\t');
                                    }
                                    local_buf.extend_from_slice(iv.contig.as_bytes());
                                    local_buf.push(b'\t');
                                    let _ = write!(local_buf, "{}\t{}\t", new_start, new_end);

                                    for i in 8..fields.len() {
                                        local_buf.extend_from_slice(fields[i].as_bytes());
                                        if i < fields.len() - 1 { local_buf.push(b'\t'); }
                                    }
                                    local_buf.push(b'\n');
                                    mapped = true;
                                }
                            }
                        }

                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender);

        let write_handle = thread::spawn(move || {
            let mut pending = std::collections::BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    writer.write_all(&data).unwrap();
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut rdr = self.parse2().expect("Failed to open table for reading");
        let batch_size = 5000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;

        for line_res in rdr.lines() {
            if let Ok(line) = line_res {
                batch.push(line);
                if batch.len() >= batch_size {
                    sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                    batch = Vec::with_capacity(batch_size);
                    chunk_id += 1;
                }
            }
        }
        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles { h.join().unwrap(); }
        write_handle.join().unwrap();

        log::info!("Successfully generated contig-level porec table into `{}`", output);
        Ok(())
    }

    pub fn dup_v0_3_0(&mut self, collapsed_list: &String, seed: usize, output: &String) {
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
        
        let mut seed_array = [0u8; 32];
        let seed_bytes = seed.to_ne_bytes();
        seed_array[..8].copy_from_slice(&seed_bytes);
        let mut rng = StdRng::from_seed(seed_array);

        let reader = common_reader(&self.file);
        let mut wtr = common_writer(output);
        
        let mut current_read_idx = String::new();
        let mut read_decisions: HashMap<String, String> = HashMap::new();

        for record in reader.lines() {
            let record = record.unwrap();
            let trimmed = record.trim();
            if trimmed.is_empty() { continue; }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 6 {
                writeln!(wtr, "{}", trimmed).unwrap();
                continue;
            }

            let read_idx = fields[0];
            let contig = fields[5];

            if read_idx != current_read_idx {
                current_read_idx = read_idx.to_string();
                read_decisions.clear();
            }

            let mut final_contig = contig;

            if let Some(candidates) = collapsed_contigs.get(contig) {
                final_contig = read_decisions.entry(contig.to_string()).or_insert_with(|| {
                    let idx = rng.gen_range(0..candidates.len());
                    candidates[idx].clone()
                });
            }

            for (i, field) in fields.iter().enumerate() {
                if i > 0 { write!(wtr, "\t").unwrap(); }
                if i == 5 {
                    write!(wtr, "{}", final_contig).unwrap();
                } else {
                    write!(wtr, "{}", field).unwrap();
                }
            }
            writeln!(wtr).unwrap();
        }
    }

    pub fn dup(&mut self, collapsed_list: &String, seed: usize, output: &String, threads: usize) {
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
        
        let mut seed_array = [0u8; 32];
        let seed_bytes = seed.to_ne_bytes();
        seed_array[..8].copy_from_slice(&seed_bytes);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let collapsed_contigs = Arc::new(collapsed_contigs);
        let num_workers = threads;
        let mut handles = vec![];

        log::info!("Duplicating collapsed contigs in parallel using {} workers...", num_workers);

        for worker_id in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let collapsed = Arc::clone(&collapsed_contigs);
            
            let mut thread_seed = seed_array;
            let thread_offset = (worker_id as u64).to_ne_bytes();
            thread_seed[8..16].copy_from_slice(&thread_offset);
            
            handles.push(thread::spawn(move || {
                let mut rng = StdRng::from_seed(thread_seed);
                let mut current_read_idx = String::new();
                let mut read_decisions: HashMap<String, String> = HashMap::new();
                let mut local_buf = Vec::with_capacity(1024 * 1024);

                while let Ok((chunk_id, batch)) = rx.recv() {
                    local_buf.clear();
                    for line in batch {
                        let trimmed = line.trim_end();
                        if trimmed.is_empty() { continue; }

                        let mut tab_indices = [0usize; 6];
                        let mut tab_count = 0;
                        for (idx, &b) in trimmed.as_bytes().iter().enumerate() {
                            if b == b'\t' {
                                tab_indices[tab_count] = idx;
                                tab_count += 1;
                                if tab_count == 6 {
                                    break;
                                }
                            }
                        }

                        if tab_count < 6 {
                            local_buf.extend_from_slice(trimmed.as_bytes());
                            local_buf.push(b'\n');
                            continue;
                        }

                        let read_idx = &trimmed[0..tab_indices[0]];
                        let contig = &trimmed[tab_indices[4] + 1..tab_indices[5]];

                        if read_idx != current_read_idx {
                            current_read_idx.clear();
                            current_read_idx.push_str(read_idx);
                            read_decisions.clear();
                        }

                        let mut final_contig = contig;
                        if let Some(candidates) = collapsed.get(contig) {
                            final_contig = read_decisions.entry(contig.to_string()).or_insert_with(|| {
                                let idx = rng.gen_range(0..candidates.len());
                                candidates[idx].clone()
                            });
                        }

                        local_buf.extend_from_slice(trimmed[0..tab_indices[4] + 1].as_bytes());
                        local_buf.extend_from_slice(final_contig.as_bytes());
                        local_buf.extend_from_slice(trimmed[tab_indices[5]..].as_bytes());
                        local_buf.push(b'\n');
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender);

        let out_path = output.clone();
        let write_handle = thread::spawn(move || {
            let mut wtr = common_writer(&out_path);
            let mut pending = BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    wtr.write_all(&data).unwrap();
                    next_chunk += 1;
                }
            }
            wtr.flush().unwrap();
        });

        let mut rdr = self.parse2().expect("Failed to open input file for reading");
        let batch_size = 10000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;
        let mut line_buf = String::new();
        let mut previous_read_idx = String::new();
        let mut first_iteration = true;

        while rdr.read_line(&mut line_buf).unwrap_or(0) > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let read_idx = trimmed.split('\t').next().unwrap_or("");
            if !first_iteration && read_idx != previous_read_idx {
                if batch.len() >= batch_size {
                    sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                    batch = Vec::with_capacity(batch_size);
                    chunk_id += 1;
                }
            }

            first_iteration = false;
            previous_read_idx.clear();
            previous_read_idx.push_str(read_idx);
            batch.push(std::mem::take(&mut line_buf));
        }

        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles { h.join().unwrap(); }
        write_handle.join().unwrap();

        log::info!("Successfully executed parallelized dup to `{}`", output);
    }

    pub fn split(&mut self, output: &String) {
        let reader = common_reader(&self.file);
        let mut wtr = common_writer(output);
        
        let (sender, receiver) = bounded::<Vec<(u32, String)>>(100);   




        log::info!("Successful output split porec table into `{}`", output);
    }

    pub fn to_pe_pair_reads_single(
        &mut self,
        fasta_path: &str,
        output_prefix: &str,
        read_len: Option<usize>, 
    ) -> anyResult<()> {
        log::info!("Loading reference genome from `{}`...", fasta_path);
        let mut genome: HashMap<String, Vec<u8>> = HashMap::new();
        let mut reader = common_reader(&fasta_path.to_string());
        let mut line = String::new();
        let mut current_ctg = String::new();
        let mut current_seq = Vec::new();

        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim();
            if trimmed.starts_with('>') {
                if !current_ctg.is_empty() {
                    genome.insert(current_ctg.clone(), current_seq.clone());
                    current_seq.clear();
                }
                current_ctg = trimmed[1..].split_whitespace().next().unwrap().to_string();
            } else {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
            line.clear();
        }
        if !current_ctg.is_empty() {
            genome.insert(current_ctg, current_seq);
        }

        let out_r1 = format!("{}_R1.fa.gz", output_prefix);
        let out_r2 = format!("{}_R2.fa.gz", output_prefix);
        let mut wtr_r1 = common_writer(&out_r1);
        let mut wtr_r2 = common_writer(&out_r2);

        let mut rdr = self.parse2().expect("Failed to open porec table for reading");
        let mut concatemer = Concatemer::new();
        let mut old_read_idx: u64 = u64::MAX;
        let mut line_buf = String::new();
        let mut pair_id = 0;

        let revcomp = |seq: &[u8]| -> Vec<u8> {
            seq.iter()
                .rev()
                .map(|&c| match c {
                    b'A' | b'a' => b'T',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    b'T' | b't' => b'A',
                    b'N' | b'n' => b'N',
                    _ => b'N',
                })
                .collect()
        };

        log::info!("Generating paired-end reads...");

        while rdr.read_line(&mut line_buf)? > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let mut parts = trimmed.split('\t');
            let read_idx = parts.next().unwrap_or("0").parse::<u64>().unwrap_or(0);

            if read_idx != old_read_idx && old_read_idx != u64::MAX {
                concatemer.sort();
                for pair in concatemer.decompose() {
                    let r1 = &pair[0];
                    let r2 = &pair[1];
                    pair_id += 1;

                    let process_read = |r: &PoreCRecord| -> Option<Vec<u8>> {
                        let seq = genome.get(&r.target)?;
                        // 1-based to 0-based
                        let start = r.target_start.saturating_sub(1) as usize;
                        let end = (r.target_end as usize).min(seq.len());
                        if start >= end {
                            return None;
                        }
                        
                        let mut fragment = seq[start..end].to_vec();
                        if r.query_strand == '-' {
                            fragment = revcomp(&fragment);
                        }
                        
                        if let Some(len) = read_len {
                            fragment.truncate(len);
                        }
                        Some(fragment)
                    };

                    if let (Some(seq1), Some(seq2)) = (process_read(&r1), process_read(&r2)) {
                        writeln!(wtr_r1, ">read_{} 1\n{}", pair_id, String::from_utf8_lossy(&seq1)).unwrap();
                        writeln!(wtr_r2, ">read_{} 2\n{}", pair_id, String::from_utf8_lossy(&seq2)).unwrap();
                    }
                }
                concatemer.clear();
            }

            let q_len = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_start = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_end = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_strand = parts.next().and_then(|s| s.chars().next()).unwrap_or('+');
            let target = parts.next().unwrap_or("").to_string();
            let t_start = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let t_end = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let mapq = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);

            concatemer.push(PoreCRecord {
                read_idx, query_length: q_len, query_start: q_start, query_end: q_end,
                query_strand: q_strand, target, target_start: t_start, target_end: t_end,
                mapq, identity: 0.0, filter_reason: "".to_string()
            });

            old_read_idx = read_idx;
            line_buf.clear();
        }


        concatemer.sort();
        for pair in concatemer.decompose() {
            let r1 = &pair[0];
            let r2 = &pair[1];
            pair_id += 1;
            let process_read = |r: &PoreCRecord| -> Option<Vec<u8>> {
                let seq = genome.get(&r.target)?;
                let start = r.target_start.saturating_sub(1) as usize;
                let end = (r.target_end as usize).min(seq.len());
                if start >= end {
                    return None;
                }
                
                let mut fragment = seq[start..end].to_vec();
                if r.query_strand == '-' {
                    fragment = revcomp(&fragment);
                }
                
                if let Some(len) = read_len {
                    fragment.truncate(len);
                }
                Some(fragment)
            };

            if let (Some(seq1), Some(seq2)) = (process_read(&r1), process_read(&r2)) {
                writeln!(wtr_r1, ">read_{} 1\n{}", pair_id, String::from_utf8_lossy(&seq1)).unwrap();
                writeln!(wtr_r2, ">read_{} 2\n{}", pair_id, String::from_utf8_lossy(&seq2)).unwrap();
            }
        }

        log::info!("Successfully generated PE reads to {}_R[12].fa", output_prefix);
        Ok(())
    }


    pub fn to_pe_pair_reads(
        &mut self,
        fasta_path: &str,
        output_prefix: &str,
        read_len: Option<usize>, 
    ) -> anyResult<()> {
        log::info!("Loading reference genome from `{}`...", fasta_path);
        let mut genome_map: HashMap<String, Vec<u8>> = HashMap::new();
        let mut reader = common_reader(&fasta_path.to_string());
        let mut line = String::new();
        let mut current_ctg = String::new();
        let mut current_seq = Vec::new();

        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim();
            if trimmed.starts_with('>') {
                if !current_ctg.is_empty() {
                    genome_map.insert(current_ctg.clone(), current_seq.clone());
                    current_seq.clear();
                }
                current_ctg = trimmed[1..].split_whitespace().next().unwrap().to_string();
            } else {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
            line.clear();
        }
        if !current_ctg.is_empty() {
            genome_map.insert(current_ctg, current_seq);
        }
        
        let genome = Arc::new(genome_map);

        let out_r1 = format!("{}_R1.fa.gz", output_prefix);
        let out_r2 = format!("{}_R2.fa.gz", output_prefix);

        log::info!("Generating paired-end reads in parallel...");

        let (tx_work, rx_work) = bounded::<Vec<Concatemer>>(200);
        let (tx_write, rx_write) = bounded::<(Vec<u8>, Vec<u8>)>(200);
        
        let pair_id_counter = Arc::new(AtomicU64::new(1));
        let num_workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
        
        let mut handles = Vec::with_capacity(num_workers);

        for _ in 0..num_workers {
            let rx = rx_work.clone();
            let tx = tx_write.clone();
            let r#gen = Arc::clone(&genome);
            let counter = Arc::clone(&pair_id_counter);
            
            handles.push(thread::spawn(move || {
                let mut buf_r1 = Vec::with_capacity(2 * 1024 * 1024);
                let mut buf_r2 = Vec::with_capacity(2 * 1024 * 1024);

                let revcomp = |seq: &[u8]| -> Vec<u8> {
                    seq.iter()
                        .rev()
                        .map(|&c| match c {
                            b'A' | b'a' => b'T',
                            b'C' | b'c' => b'G',
                            b'G' | b'g' => b'C',
                            b'T' | b't' => b'A',
                            b'N' | b'n' => b'N',
                            _ => b'N',
                        })
                        .collect()
                };

                while let Ok(batch) = rx.recv() {
                    for mut concatemer in batch {
                        concatemer.sort();
                        for pair in concatemer.decompose() {
                            let r1 = &pair[0];
                            let r2 = &pair[1];

                            let process_read = |r: &PoreCRecord| -> Option<Vec<u8>> {
                                let seq = r#gen.get(&r.target)?;
                                let start = r.target_start.saturating_sub(1) as usize;
                                let end = (r.target_end as usize).min(seq.len());
                                if start >= end {
                                    return None;
                                }
                                
                                let mut fragment = seq[start..end].to_vec();
                                if r.query_strand == '-' {
                                    fragment = revcomp(&fragment);
                                }
                                
                                if let Some(len) = read_len {
                                    fragment.truncate(len);
                                }
                                Some(fragment)
                            };

                            if let (Some(seq1), Some(seq2)) = (process_read(&r1), process_read(&r2)) {
                                let pid = counter.fetch_add(1, AtomOrdering::Relaxed);
                                use std::io::Write;
                                let _ = writeln!(&mut buf_r1, ">read_{}/1\n{}", pid, unsafe { std::str::from_utf8_unchecked(&seq1) });
                                let _ = writeln!(&mut buf_r2, ">read_{}/2\n{}", pid, unsafe { std::str::from_utf8_unchecked(&seq2) });
                            }
                        }
                    }
                    if buf_r1.len() >= 1024 * 1024 {
                        tx.send((buf_r1.clone(), buf_r2.clone())).unwrap();
                        buf_r1.clear();
                        buf_r2.clear();
                    }
                }
                
                if !buf_r1.is_empty() {
                    tx.send((buf_r1, buf_r2)).unwrap();
                }
            }));
        }
        drop(tx_write); // Drop the original sender so writer knows when to stop

        // Writer thread
        let write_handle = thread::spawn(move || {
            let mut wtr_r1 = common_writer(&out_r1);
            let mut wtr_r2 = common_writer(&out_r2);
            while let Ok((chunk_r1, chunk_r2)) = rx_write.recv() {
                wtr_r1.write_all(&chunk_r1).unwrap();
                wtr_r2.write_all(&chunk_r2).unwrap();
            }
            wtr_r1.flush().unwrap();
            wtr_r2.flush().unwrap();
        });

        // Reader (Main) thread
        let mut rdr = self.parse2().expect("Failed to open porec table for reading");
        let mut concatemer = Concatemer::new();
        let mut old_read_idx: u64 = u64::MAX;
        let mut line_buf = String::new();
        
        let batch_size = 2000;
        let mut batch = Vec::with_capacity(batch_size);

        while rdr.read_line(&mut line_buf)? > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let mut parts = trimmed.split('\t');
            let read_idx = parts.next().unwrap_or("0").parse::<u64>().unwrap_or(0);

            if read_idx != old_read_idx && old_read_idx != u64::MAX {
                batch.push(std::mem::take(&mut concatemer));
                if batch.len() >= batch_size {
                    tx_work.send(std::mem::take(&mut batch)).unwrap();
                }
            }

            let q_len      = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_start    = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_end      = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_strand   = parts.next().and_then(|s| s.chars().next()).unwrap_or('+');
            let target     = parts.next().unwrap_or("").to_string();
            let t_start    = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let t_end      = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let mapq       = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);

            concatemer.push(PoreCRecord {
                read_idx, query_length: q_len, query_start: q_start, query_end: q_end,
                query_strand: q_strand, target, target_start: t_start, target_end: t_end,
                mapq, identity: 0.0, filter_reason: "".to_string()
            });

            old_read_idx = read_idx;
            line_buf.clear();
        }

        // Flush remaining
        if concatemer.count() > 0 {
            batch.push(std::mem::take(&mut concatemer));
        }
        if !batch.is_empty() {
            tx_work.send(batch).unwrap();
        }
        drop(tx_work);

        // Wait for workers and writer
        for h in handles {
            h.join().unwrap();
        }
        write_handle.join().unwrap();

        log::info!("Successfully generated PE reads to {}_R[12].fa.gz", output_prefix);
        Ok(())
    }


    pub fn downsample(
        &mut self,
        output: &String,
        seed: u64,
        min_quality: u8,
        min_order: usize,
        max_order: usize,
        target_reads: Option<u64>,
        target_pairs: Option<u64>,
        frac: Option<f64>,
        frac_by_pairs: bool,
    ) -> anyResult<()> {
        if [target_reads.is_some(), target_pairs.is_some(), frac.is_some()]
            .into_iter()
            .filter(|x| *x)
            .count() != 1
        {
            anyhow::bail!("Specify exactly one of --reads, --pairs, or --frac");
        }
        if let Some(f) = frac {
            if !(0.0 < f && f <= 1.0) {
                anyhow::bail!("--frac must be in (0, 1]");
            }
        }

        log::info!(
            "Downsample pass1 scanning... (min_mapq={}, order in [{}, {}))",
            min_quality,
            min_order,
            max_order
        );

        let mut rdr = self.parse2().expect("Failed to open input file");
        let mut line = String::new();

        let mut prev_read: Option<u64> = None;
        let mut cur_order: usize = 0;

        let mut per_read: Vec<(u64, u32, u64)> = Vec::new();

        #[inline]
        fn pairs_of(order: usize) -> u64 {
            let n = order as u64;
            n.saturating_mul(n.saturating_sub(1)) / 2
        }

        while rdr.read_line(&mut line)? > 0 {
            let t = line.trim_end();
            if t.is_empty() || t.starts_with('#') {
                line.clear();
                continue;
            }

            let mut it = t.split('\t');
            let read_idx = it.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);

       
            let mut mapq: u8 = 0;
            for _ in 0..7 {
                if it.next().is_none() {
                    break;
                }
            }
            if let Some(mq_s) = it.next() {
                mapq = mq_s.parse::<u8>().unwrap_or(0);
            }

            let is_eligible_piece = mapq >= min_quality;

            match prev_read {
                None => {
                    prev_read = Some(read_idx);
                    cur_order = if is_eligible_piece { 1 } else { 0 };
                }
                Some(pr) if pr == read_idx => {
                    if is_eligible_piece {
                        cur_order += 1;
                    }
                }
                Some(pr) => {
                    // finalize previous read
                    if cur_order >= min_order && cur_order < max_order {
                        let p = pairs_of(cur_order);
                        per_read.push((pr, cur_order as u32, p));
                    }
                    // start new
                    prev_read = Some(read_idx);
                    cur_order = if is_eligible_piece { 1 } else { 0 };
                }
            }

            line.clear();
        }
        // finalize last
        if let Some(pr) = prev_read {
            if cur_order >= min_order && cur_order < max_order {
                let p = pairs_of(cur_order);
                per_read.push((pr, cur_order as u32, p));
            }
        }

        let total_reads = per_read.len() as u64;
        let total_pairs: u64 = per_read.iter().map(|x| x.2).sum();

        log::info!(
            "Eligible reads: {} ; eligible virtual pairs: {}",
            total_reads,
            total_pairs
        );

        let (want_reads, want_pairs): (Option<u64>, Option<u64>) = if let Some(n) = target_reads {
            (Some(n.min(total_reads)), None)
        } else if let Some(p) = target_pairs {
            (None, Some(p.min(total_pairs)))
        } else {
            // frac
            let f = frac.unwrap();
            if frac_by_pairs {
                (None, Some(((total_pairs as f64) * f).round() as u64))
            } else {
                (Some(((total_reads as f64) * f).round() as u64), None)
            }
        };

        let mut rng = StdRng::seed_from_u64(seed);
        per_read.shuffle(&mut rng);

        let mut chosen: std::collections::HashSet<u64> = std::collections::HashSet::new();

        if let Some(k) = want_reads {
            for (rid, _ord, _p) in per_read.iter().take(k as usize) {
                chosen.insert(*rid);
            }
            log::info!("Selected reads: {}", chosen.len());
        } else if let Some(tp) = want_pairs {
            let mut acc: u64 = 0;
            for (rid, _ord, p) in per_read.iter() {
                if acc >= tp {
                    break;
                }
                chosen.insert(*rid);
                acc = acc.saturating_add(*p);
            }
            log::info!(
                "Selected reads: {} (approx pairs: {} / target {})",
                chosen.len(),
                acc,
                tp
            );
        }

        log::info!("Downsample pass2 writing output to `{}`...", output);
        let mut rdr2 = self.parse2().expect("Failed to open input file (pass2)");
        let mut wtr = common_writer(output);

        let mut line2 = String::new();
        while rdr2.read_line(&mut line2)? > 0 {
            let t = line2.trim_end();
            if t.is_empty() {
                line2.clear();
                continue;
            }
            if t.starts_with('#') {
                wtr.write_all(line2.as_bytes())?;
                line2.clear();
                continue;
            }

            let read_idx = t
                .split('\t')
                .next()
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(u64::MAX);

            if chosen.contains(&read_idx) {
                wtr.write_all(line2.as_bytes())?;
            }

            line2.clear();
        }

        log::info!(
            "Downsample finished. Output reads: {} (seed={})",
            chosen.len(),
            seed
        );
        Ok(())
    }

    pub fn to_depth(
        &mut self, 
        chromsize: &String, 
        window_size: usize, 
        step_size: usize,
        min_mapq: u8, 
        output: &String
    ) -> Result<(), Box<dyn Error>> {
        let mut chrom_names = Vec::new();
        let mut chrom_sizes = Vec::new();
        let mut chrom_map = HashMap::new();

        log::info!("Loading chromsizes from `{}`...", chromsize);
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
                chrom_map.insert(chr.clone(), chrom_names.len());
                chrom_names.push(chr);
                chrom_sizes.push(size);
            }
        }
        
        let num_chroms = chrom_names.len();
        let chrom_map = Arc::new(chrom_map);

        let (sender, receiver) = bounded::<Vec<String>>(200);
        let num_workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(8);
        
        let mut handles = Vec::new();
        for _ in 0..num_workers {
            let receiver = receiver.clone();
            let chrom_map = chrom_map.clone();
            handles.push(thread::spawn(move || {
                let mut local_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];
                while let Ok(lines) = receiver.recv() {
                    for line in lines {
                        if line.is_empty() || line.starts_with('#') {
                            continue;
                        }
                        // pore-c columns: read_idx, q_len, q_start, q_end, q_strand, target, t_start, t_end, mapq ...
                        let mut parts = line.split('\t');
                        let target = parts.nth(5).unwrap_or(""); 
                        let t_start_str = parts.next().unwrap_or("0");
                        let t_end_str = parts.next().unwrap_or("0");
                        let mapq_str = parts.next().unwrap_or("0");

                        let mapq = mapq_str.parse::<u8>().unwrap_or(0);
                        if mapq < min_mapq {
                            continue;
                        }

                        if let Some(&idx) = chrom_map.get(target) {
                            let start = t_start_str.parse::<u32>().unwrap_or(0);
                            let end = t_end_str.parse::<u32>().unwrap_or(0);
                            local_events[idx].push((start, 1));
                            local_events[idx].push((end, -1));
                        }
                    }
                }
                local_events
            }));
        }

        log::info!("Parsing Pore-C table to calculate coverage depth...");
        let mut rdr = self.parse2().expect("Failed to open porec table for reading");
        let mut line_buf = String::new();
        let batch_size = 10_000;
        let mut batch = Vec::with_capacity(batch_size);

        while rdr.read_line(&mut line_buf)? > 0 {
            if !line_buf.trim().is_empty() {
                batch.push(line_buf.clone());
            }
            if batch.len() >= batch_size {
                sender.send(std::mem::take(&mut batch)).unwrap();
            }
            line_buf.clear();
        }

        if !batch.is_empty() {
            sender.send(batch).unwrap();
        }

        drop(sender);

        let mut global_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];
        for handle in handles {
            let local_events = handle.join().unwrap();
            for (i, events) in local_events.into_iter().enumerate() {
                global_events[i].extend(events);
            }
        }

        log::info!("Aggregating coverage across parallel windows...");
        
        let output_clone = output.clone();
        // Using Rayon Scope to handle parallel generation and sequential output file writing
        rayon::scope(|s| {
            let (tx, rx) = bounded::<(usize, String)>(1024);

            s.spawn(move |_| {
                use std::io::Write;
                let writer = common_writer(&output_clone);
                let mut buf_writer = std::io::BufWriter::with_capacity(1024 * 1024, writer);
                let mut pending = std::collections::BTreeMap::new();
                let mut next_idx = 0;
                
                while let Ok((idx, data)) = rx.recv() {
                    pending.insert(idx, data);
                    while let Some(d) = pending.remove(&next_idx) {
                        buf_writer.write_all(d.as_bytes()).unwrap();
                        next_idx += 1;
                    }
                }
                buf_writer.flush().unwrap();
            });

            global_events.into_par_iter().enumerate().for_each_with(tx, |tx, (i, mut events)| {
                if events.is_empty() { 
                    tx.send((i, String::new())).unwrap();
                    return; 
                }

                let chrom_len = chrom_sizes[i];
                let chrom_name = &chrom_names[i];

                events.par_sort_unstable_by_key(|e| e.0);

                let mut coverage = vec![0i32; chrom_len];
                let mut current_cov = 0i32;
                let mut prev_pos = 0;

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

                let mut prefix_sum = vec![0i64; chrom_len + 1];
                let mut running = 0i64;
                for (k, &cov) in coverage.iter().enumerate() {
                    running += cov as i64;
                    prefix_sum[k+1] = running;
                }

                let mut chunk_out = String::with_capacity(4096);
                use std::fmt::Write;

                let mut start = 0usize;
                while start < chrom_len {
                    let end = std::cmp::min(start + window_size, chrom_len);
                    let len = end - start;
                    if len == 0 { break; }
                    
                    let sum = prefix_sum[end] - prefix_sum[start];
                    let mean = sum as f64 / len as f64;
                    
                    writeln!(&mut chunk_out, "{}\t{}\t{}\t{:.3}", chrom_name, start, end, mean).unwrap();

                    start += step_size;
                }
                
                tx.send((i, chunk_out)).unwrap();
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
