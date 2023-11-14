#![allow(unused_variables, unused_imports)]
use bio::io::fastq;
use bio::pattern_matching::horspool::Horspool;
use std::io;
use std::io::*;
use std::io::prelude::*;
use std::error::Error;
use std::path::Path;
use crate::core::{common_reader, common_writer};

pub fn cut_site(fastq: String, pattern: &[u8], output: String) -> Result<()>{
    let fastq_path = Path::new(&fastq);
    
    let buffered = common_reader(&fastq);
    let reader = fastq::Reader::new(buffered);
    
    let write_buffered = common_writer(&output);
    let mut writer = fastq::Writer::new(write_buffered);

    // let pattern = b"GATCGATC";
    let horspool = Horspool::new(pattern);

    for result in reader.records() {
        let record = result.expect("Error during fastq record parsing");
        let seq = record.seq();
        let seq_length = seq.len();
        let seq_str = std::str::from_utf8(&seq).unwrap();
        let qual = record.qual();
        let mut occ: Vec<usize> = horspool.find_all(seq).collect();
        let mut len_vec = vec![seq_length];
        occ.append(&mut len_vec);

        if occ.len() <= 1 {
            writer.write(record.id(), record.desc(), record.seq(), record.qual()).unwrap();
        } else {
            let mut start_pos = 0;
            let end_pos = seq_length;
            for (i, pos) in occ.into_iter().enumerate() {
                
                let end_pos = if pos >= seq_length {
                    pos
                } else {
                    pos + pattern.len() / 2
                };
 
      
                let cut_seq = &seq_str[start_pos..end_pos];
                let cut_qual = &qual[start_pos..end_pos];
                let new_name = record.id().to_owned() + "/" + &i.to_string();
                writer.write(&new_name, None, cut_seq.as_bytes(), cut_qual).unwrap();

                start_pos += cut_seq.len();
            }
        }
    }
    Ok(())
}