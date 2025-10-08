
use rust_htslib::bam::{ Read, Reader };
use noodles::fasta::{ self as fasta, fai };
use std::env;
use std::collections::{ HashMap, HashSet };
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::aligner::mod_seq_complement;
use crate::bed::{ 
    BedCpG,
    BedCpGRecord,
    BedMethylSimple, 
    BedMethylSimpleRecord,
    BedIterator,
};
use crate::methalign::parse_bedgraph;
use crate::fastx::Fastx;
use crate::core::BaseTable;
use crate::core::{ common_reader, common_writer };


pub struct ModRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub frac: f32,
    pub score: u32,
}

pub fn qual_to_prob(qual: f32) -> f32 {
    let mut qual = qual as f32;
    qual = (qual + 0.5f32) / 256f32; 
    qual
}

pub fn modbam2fastq(input_bam: &String, min_prob: f32,
                    output: &String) -> Result<(), Box<dyn std::error::Error>>  {
    
    let mut bam = Reader::from_path(input_bam).unwrap();
    let mut wtr = common_writer(output);

    for r in bam.records() {
        let record = r.unwrap();
        if record.is_secondary() {
            continue;
        }
        let read_id = record.qname();
        let sequence = record.seq();
        let mut new_sequence: Vec<u8> = sequence.as_bytes();
        if record.is_reverse() {
            // new_sequence.reverse();
            for i in 0..new_sequence.len() {
                new_sequence[i] = match new_sequence[i] {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'G' => b'C',
                    b'C' => b'G',
                    _ => new_sequence[i],
                };
            }
        }
        
        let mut new_quality: Vec<u8> = record.qual().to_vec();

        // let mut remove_g_idx: Vec<usize> = Vec::new();

        if let Ok(mods) = record.basemods_iter() {
            
            for (idx, res) in mods.enumerate() {
                if let Ok( (position, m) ) = res {
                    // let position2 = if m.strand == 0 {
                    //     position + 1
                    // } else {
                    //     position - 1
                    // };

                    // if position2 >= new_sequence.len() as i32 {
                    //     continue;
                    // }

                    let qual = qual_to_prob(m.qual as f32);
                    if qual < min_prob {
                        continue;
                    }
                    
                    new_sequence[position as usize] = m.modified_base as u8;
                    
                   
                   
                    // if new_sequence[position2 as usize] == b'G' 
                    //     || new_sequence[position2 as usize] == b'g' {
                    //     new_sequence[position as usize] = m.modified_base as u8;
                    //     // remove_g_idx.push(position2 as usize);
                    // }
                }
                
            }

            if record.is_reverse() {
                new_sequence.reverse();
                new_quality.reverse();
            }

            // remove_g_idx.sort();
            // for idx in remove_g_idx.iter().rev() {
            //     new_sequence.remove(*idx);
            //     new_quality.remove(*idx);
            // }
            // println!("{}", String::from_utf8(new_sequence.clone()).unwrap());
            new_quality.iter_mut().for_each(|x| *x += 33);
            write!(wtr, "@{}\n{}\n+\n{}\n", String::from_utf8(read_id.to_vec()).unwrap(),
                    String::from_utf8(new_sequence).unwrap(), 
                    String::from_utf8(new_quality).unwrap()).unwrap();

            
        };
    }

    
    Ok(())
}

pub fn modify_fasta(fasta_path: &String, bedgraph: &String, 
        min_frac: f64, 
        output: &String) -> Result<(), Box<dyn std::error::Error>> {
    
    let fasta_path = Path::new(fasta_path);
    let fa_hash: HashMap<String, String> = if let Some(fa) = Some(fasta_path) {
       
        let fa = fa.to_str().unwrap().to_string();
        let fasta = Fastx::new(&fa);
        fasta.get_chrom_seqs().unwrap()
    } else {
        HashMap::new()
    };
    // Create output writer
    let mut wtr = common_writer(output);
    

    // Read bedgraph file 
    let bg_map: HashMap<String, HashSet<i64>> = if let Some(bg) = Some(bedgraph) {
        parse_bedgraph(bedgraph, min_frac).expect("Failed to parse bedgraph")
    } else {
        HashMap::new()
    };


    log::info!("Start to modify fasta file");
    
    for chrom in bg_map.keys() {
        let sites = bg_map.get(chrom).unwrap();
        let record = fa_hash.get(chrom).unwrap();
        let mut seq = record.as_bytes().to_vec();
        log::info!("Processing {}", chrom);
        
        // for (idx, c) in seq.iter_mut().enumerate() {
        //     if sites.contains(&idx.try_into().unwrap()) {
        //         if idx == seq.len() - 1 {
        //             continue;
        //         }
        //         if seq[idx+1] == b'G' || seq[idx+1] == b'g' {
        //             *c = b'm';
        //         }
              
        //     }
        // }
        for idx in 0..seq.len() {
            if sites.contains(&(idx as i64)) {
                if idx == seq.len() - 1 {
                    continue;
                }
                if seq[idx + 1] == b'G' || seq[idx + 1] == b'g' {
                    seq[idx] = b'm'; 
                }
            }
        }

        // seq.retain(|&c| c != b'\0');
        

        let mod_seq = String::from_utf8(seq).unwrap();
        write!(wtr, ">{}\n{}\n", chrom, mod_seq).unwrap();

        
        log::info!("Finish to modify {}", chrom);
 
    }
    
    Ok(())
}