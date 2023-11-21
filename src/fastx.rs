use anyhow::Result as AnyResult;
use bio::io::fastq::{Reader, Record, Writer};
use std::collections::HashMap;
use std::error::Error;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::io;
use std::path::Path;

use seq_io::prelude::*;
use seq_io::fastx::Reader as FastxReader;
use seq_io::parallel::read_process_fastx_records;

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

    pub fn count_re(&self, motif: &String) -> AnyResult<HashMap<String, u64>> {
        log::info!("Counting motif `{}` in `{}`", motif, self.file);
        let reader = common_reader(&self.file);
        let reader = FastxReader::new(reader);
        let mut chrom_count: HashMap<String, u64> = HashMap::new();

        read_process_fastx_records(reader, 4, 2,
            |record, count| { // runs in worker
                *count = record.seq().windows(motif.len()).filter(|&x| x == motif.as_bytes()).count();
                                
            },
            |record, count| { // runs in main thread
                chrom_count.insert(record.id().unwrap().to_owned(), *count as u64); 
                None::<()>
            }).unwrap();
    
        Ok(chrom_count)
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

