#![allow(unused_variables, unused_imports)]
use anyhow::Result;
use bio::io::fastq;
use bio::pattern_matching::horspool::Horspool;
use needletail::parse_fastx_file;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::error::Error;
use std::io;
use std::io::prelude::*;
use std::io::*;
use std::path::Path;

use crate::core::{common_reader, common_writer};

pub fn cut_site(fastq: String, pattern: &[u8], output: String) -> Result<()> {
    let mut reader = parse_fastx_file(&fastq).expect("Error opening input file");
    let mut writer = common_writer(&output);

    let horspool = Horspool::new(pattern);
    let p_len = pattern.len();
    let half_p_len = p_len / 2;
    ThreadPoolBuilder::new()
        .num_threads(8)
        .build_global()
        .unwrap();

    let chunk_size = 10000;

    loop {
        let mut chunk = Vec::with_capacity(chunk_size);
        for _ in 0..chunk_size {
            if let Some(record) = reader.next() {
                let seq_record = record?;
                chunk.push((
                    seq_record.id().to_vec(),
                    seq_record.seq().to_vec(),
                    seq_record.qual().map(|q| q.to_vec()),
                ));
            } else {
                break;
            }
        }

        if chunk.is_empty() {
            break;
        }

        let results: Vec<Vec<u8>> = chunk
            .into_par_iter()
            .map(|(id, seq, qual)| {
                let mut out_buf = Vec::new();
                let seq_length = seq.len();

                let mut occ: Vec<usize> = horspool.find_all(&seq).collect();
                occ.push(seq_length);

                let mut start_pos = 0;
                for pos in occ {
                    let end_pos = if pos >= seq_length {
                        pos
                    } else {
                        pos + half_p_len
                    };
                    let cut_seq = &seq[start_pos..end_pos];

                    match qual {
                        Some(ref q) => {
                            let cut_qual = &q[start_pos..end_pos];
                            out_buf.push(b'@');
                            out_buf.extend_from_slice(&id);
                            out_buf.push(b'\n');
                            out_buf.extend_from_slice(cut_seq);
                            out_buf.extend_from_slice(b"\n+\n");
                            out_buf.extend_from_slice(cut_qual);
                            out_buf.push(b'\n');
                        }
                        None => {
                            out_buf.push(b'>');
                            out_buf.extend_from_slice(&id);
                            out_buf.push(b'\n');
                            out_buf.extend_from_slice(cut_seq);
                            out_buf.push(b'\n');
                        }
                    }
                    start_pos = end_pos;
                }
                out_buf
            })
            .collect();

        for buf in results {
            writer.write_all(&buf)?;
        }
    }
    writer.flush()?;
    Ok(())
}
