
use rust_htslib::bam::{ Read, Reader };
use noodles::fasta::{ self as fasta, fai };
use std::env;
use std::collections::{ HashMap, HashSet };
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use crate::bed::{ 
    BedCpG,
    BedCpGRecord,
    BedMethylSimple, 
    BedMethylSimpleRecord };
use crate::core::BaseTable;
use crate::core::{ common_reader, common_writer };


pub struct ModRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub frac: f32,
    pub score: u32,
}

fn qual_to_prob(qual: f32) -> f32 {
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
        let read_id = record.qname();
        let sequence = record.seq();
        let mut new_sequence: Vec<u8> = sequence.as_bytes();
        
        let mut new_quality: Vec<u8> = record.qual().to_vec();

        let mut remove_g_idx: Vec<usize> = Vec::new();

        if let Ok(mods) = record.basemods_iter() {
            
            for (idx, res) in mods.enumerate() {
                if let Ok( (position, m) ) = res {
                    let position2 = if m.strand == 0 {
                        position + 1
                    } else {
                        position - 1
                    };
                    let qual = qual_to_prob(m.qual as f32);
                    if qual < min_prob {
                        continue;
                    }
                    
                    
                    if new_sequence[position2 as usize] == b'G' 
                        || new_sequence[position2 as usize] == b'g' {
                        new_sequence[position as usize] = m.modified_base as u8;
                        remove_g_idx.push(position2 as usize);
                    }
                }
                
            }

            remove_g_idx.sort();
            for idx in remove_g_idx.iter().rev() {
                new_sequence.remove(*idx);
                new_quality.remove(*idx);
            }
            new_quality.iter_mut().for_each(|x| *x += 33);
            write!(wtr, "@{}\n{}\n+\n{}\n", String::from_utf8(read_id.to_vec()).unwrap(),
                    String::from_utf8(new_sequence).unwrap(), 
                    String::from_utf8(new_quality).unwrap()).unwrap();

            
        };
    }

    
    Ok(())
}

pub fn modify_fasta(fasta_path: &String, bed: &String, 
        min_score: u32, min_frac: f32, bedcpg: &String,
        output: &String) -> Result<(), Box<dyn std::error::Error>> {
    
    
    let fasta_path = Path::new(fasta_path);
    
    // Create index file
    let index_file = format!("{}.fai", fasta_path.display().to_string());
    let index_file = Path::new(&index_file);
    if !index_file.exists() {
        let index = fasta::index(fasta_path)?;
        let mut writer = fai::Writer::new(File::create(index_file)?);
        writer.write_index(&index)?;
    }
    let mut reader = fasta::indexed_reader::Builder::default().build_from_path(fasta_path).unwrap();
    let mut wtr = common_writer(output);
    let mut writer = fasta::Writer::new(wtr);
    let mut hash = HashMap::<String, Vec<ModRecord>>::new();
    let bed_iter = if bedcpg == "bedMethyl" {
        BedCpG::new(bed)
    } else {
        BedCpG::new(bed)
    };
    log::info!("Loading bed file");
    for record in bed_iter {
        if !hash.contains_key(&record.chrom) {
            let mut ranges = Vec::<ModRecord>::new();

            hash.insert(record.chrom.clone(), ranges);
        } 
        
        if record.score < min_score || record.frac < min_frac {
            continue;
        }
        hash.get_mut(&record.chrom).expect("REASON").push(record);

    }

    log::info!("Start to modify fasta file");

    for chrom in hash.keys() {
        let ranges = hash.get(chrom).unwrap();
        let region = chrom.parse().unwrap();
        let record = reader.query(&region).unwrap();
        let mut seq = record.sequence().as_ref().to_vec();
        let mut new_seq: Vec<&u8> = Vec::new();
        let mut replace_m_idx: HashSet<usize> = HashSet::new();
        let mut remove_g_idx: HashSet<usize> = HashSet::new();
        for r in ranges {
            replace_m_idx.insert(r.start);
            remove_g_idx.insert(r.end);
   
        }
        log::info!("Processing {}", chrom);
        
        for (idx, c) in seq.iter_mut().enumerate() {
            if replace_m_idx.contains(&idx) {
                *c = b'm';
            } else if remove_g_idx.contains(&idx) {
                // Remove element by setting it to a "null" value
                *c = b'\0';
            }
        }

        seq.retain(|&c| c != b'\0');
        

        let definition = fasta::record::Definition::new(record.name().to_string(), None);
        let new_seq = fasta::record::Sequence::from(seq);
        let new_record = fasta::Record::new(definition, new_seq);

        writer.write_record(&new_record).unwrap();
        
        log::info!("Finish to modify {}", chrom);
 
    }
    
    Ok(())
}