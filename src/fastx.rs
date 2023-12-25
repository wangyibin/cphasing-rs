use anyhow::Result as AnyResult;
use bio::io::fastq::{Reader, Record, Writer};
use bio::pattern_matching::horspool::Horspool;
use indexmap::IndexMap;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::io;
use std::path::Path;
use std::ops::AddAssign;

use seq_io::prelude::*;
use seq_io::fastq;
use seq_io::fastx::Reader as FastxReader;
use seq_io::parallel::{ read_process_fastx_records, read_process_recordsets };

use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;

pub struct Fastx {
    pub file: String,
}

impl BaseTable for Fastx {
    fn new(name: &String) -> Fastx {
        Fastx { file: name.clone() }
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

impl Fastx {
    pub fn get_chrom_size(&self) -> AnyResult<HashMap<String, u64>>  {
        let reader = common_reader(&self.file);
        let reader = FastxReader::new(reader);
        let mut chrom_size: HashMap<String, u64> = HashMap::new();

        read_process_fastx_records(reader, 4, 2,
            |record, length| { // runs in worker
                *length = record.seq_lines()
                                .fold(0, |l, seq| l + seq.len());
            },
            |record, length| { // runs in main thread
                chrom_size.insert(record.id().unwrap().to_owned(), *length as u64); 
                None::<()>
            }).unwrap();
    
        Ok(chrom_size)
    }

    pub fn get_chrom_seqs (&self) -> AnyResult<HashMap<String, String>> {
        let reader = common_reader(&self.file);
        let reader = FastxReader::new(reader);
        let mut chrom_seqs: HashMap<String, String> = HashMap::new();

        read_process_fastx_records(reader, 4, 2,
            |record, seq| { // runs in worker
                *seq = record.seq_lines()
                                .fold(String::new(), |s, seq| s + &String::from_utf8(seq.to_vec()).unwrap());
            },
            |record, seq| { // runs in main thread
                chrom_seqs.insert(record.id().unwrap().to_owned(), seq.to_owned()); 
                None::<()>
            }).unwrap();
    
        Ok(chrom_seqs)
    }

    pub fn count_re(&self, patterns: &String) -> AnyResult<IndexMap<String, u64>> {
        
        let mut chrom_count: IndexMap<String, u64> = IndexMap::new();
        let pattern_vec = patterns.split(",").collect::<Vec<&str>>();
    
        // convert to uppercase
        let pattern_vec_uppercase = pattern_vec.iter().map(|x| x.to_uppercase()).collect::<Vec<String>>();
        for pattern in &pattern_vec_uppercase {
            log::info!("Counting restriction site `{}` in `{}`", pattern, self.file);
        }
        // convert to lowercase
        let pattern_vec_lowercase = pattern_vec.iter().map(|x| x.to_lowercase()).collect::<Vec<String>>();
        // merge two vectors
        let pattern_vec = pattern_vec_uppercase.iter().chain(pattern_vec_lowercase.iter()).collect::<Vec<&String>>();
        
        

        for pattern in pattern_vec {
            log::set_max_level(log::LevelFilter::Off);
            let reader = common_reader(&self.file);
            log::set_max_level(log::LevelFilter::Info);
            let reader = FastxReader::new(reader);
            read_process_fastx_records(reader, 4, 2,
                |record, count| { // runs in worker
                    *count = record.seq().windows(pattern.len()).filter(|&x| x == pattern.as_bytes()).count();
                                    
                },
                |record, count| { // runs in main thread
                    

                    chrom_count.entry(record.id().unwrap().to_owned()).or_insert(0).add_assign(*count as u64);
                    None::<()>
                }).unwrap();
        }
        
    
        Ok(chrom_count)
    }

    pub fn digest(&self, patterns: &String, slope: i64) -> AnyResult<HashMap<String, Vec<Vec<i64>>>> {
        
        let pattern_vec = patterns.split(",").collect::<Vec<&str>>();
        
        
        let mut positions: HashMap<String, Vec<Vec<i64>>> = HashMap::new();

        for pattern in pattern_vec {
            let mut pattern_uppercase = pattern.to_owned();
            pattern_uppercase.make_ascii_uppercase();
            log::info!("Identifing restriction site positions `{}` in `{}`", pattern_uppercase, self.file);
            let horspool = Horspool::new(pattern_uppercase.as_bytes());
            let reader = common_reader(&self.file);
            let reader = bio::io::fasta::Reader::new(reader);
            for result in reader.records() {
                let record = result?;
                let name = record.id().to_owned();
                let seq = record.seq();
                // seq uppercase

                let mut seq = seq.to_owned();
                seq.make_ascii_uppercase();
                let seq_length = seq.len();
                let seq_str = std::str::from_utf8(&seq).unwrap();
                let mut occ: Vec<usize> = horspool.find_all(&seq).collect();
                let mut len_vec = vec![seq_length];

                if occ.len() >= 1 {
                    let end_pos = seq_length;
                    for (i, pos) in occ.into_iter().enumerate() {
                        
                        let pos_mid = if pos >= seq_length {
                            pos as i64
                        } else {
                            (pos + pattern.len() / 2) as i64
                        };
                        let mut pos_start = pos_mid - slope;
                        let mut pos_end = pos_mid + slope;
                        if pos_start < 0 {
                            pos_start = 0
                        }
                        if pos_end > end_pos as i64 {
                            pos_end = end_pos as i64
                        }
                        
                        positions.entry(name.clone()).or_insert(Vec::new()).push(vec![pos_start, pos_end]);

                    }
                }
            }
        }
        
        Ok(positions)
    }
}

// split fastq into several files by record number
pub fn split_fastq(input_fastq: &String, output_prefix: &String, 
             record_num: usize) -> Result<(), Box<dyn std::error::Error>> {
    let buf = common_reader(input_fastq);
    let fastq = Reader::new(buf);
    let mut i = 0;
    let mut j = 0;
    log::info!("write {} records to {}", &record_num, format!("{}_{}.fastq.gz", output_prefix, j));
    let writer = common_writer(&format!("{}_{}.fastq.gz", output_prefix, j));
    let mut wtr = Writer::new(writer);
    for r in fastq.records() {
        let record = r?;
        wtr.write_record(&record)?;
        i += 1;
        if i == record_num {
            j += 1;
            let writer = common_writer(&format!("{}_{}.fastq.gz", output_prefix, j));
            wtr = Writer::new(writer);
            log::info!("write {} records to {}", i, format!("{}_{}.fastq.gz", output_prefix, j));
            i = 0;
        }
    }


    
    Ok(())
}

