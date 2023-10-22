use anyhow::Result as anyResult;
// use crossbeam::queue::ArrayQueue;
use minimap2::*;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, record::Cigar,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use std::path::PathBuf;
use std::collections::HashMap;
use std::rc::Rc;
use std::str;
use std::mem;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

use crate::fastx::Fastx;
use crate::methy::qual_to_prob;

enum WorkQueue<T> {
    Work(T),
    Result(T),
}


#[derive(Clone, Debug)]
pub struct ReadUnit {
    data: Vec<Record>,
}

impl ReadUnit {
    fn new() -> Self {
        ReadUnit {
            data: Vec::new(),
        }
    }

    fn clear (&mut self) {
        self.data.clear();
    }
}

fn calculate_alignment_score_from_cigar(cigar: &Vec<(u32, u8)>) -> i32 {
    let mut score = 0;
    for op in cigar.iter() {
        match op.1 {
            0 => score += op.0 as i32 * 2,
            1 => score -= op.0 as i32 * 4,
            2 => score -= op.0 as i32 * 4,
            _ => {}
        }
    }

    score 
}


// calculate alignment score from CigarStringView
fn calculate_alignment_score_from_cigar_string_view(cigar: CigarStringView) -> i32 {
    let mut score = 0;
    for op in cigar.iter() {
        match op {
            Cigar::Match(len) => score += len * 2,
            Cigar::Diff(len) => score -= len * 4,
            Cigar::Del(len) => score -= len * 4,
            Cigar::Ins(len) => score -= len * 4,
            _ => {} 
        }
    }

    score.try_into().unwrap()
}

fn mod_seq_complement(seq: &mut Vec<u8>) {
    let mut flag = false;
    for i in 0..seq.len() {
        seq[i] = match seq[i] {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => {
                if flag {
                    flag = false;
                    b'm' 
                } else {
                    b'C'
                }
            }
            b'C' => b'G',
            b'm' => {
                flag = true;  
                b'G'
            },
            _ => seq[i],
        }
    }
    seq.reverse();
}


#[derive(Debug)]
pub struct AlignmentUnit {
    Primary: Vec<Record>,
    Secondary: Vec<Vec<Record>>,
}

impl AlignmentUnit {
    fn new() -> Self {
        AlignmentUnit {
            Primary: Vec::new(),
            Secondary: Vec::new(),
        }
    }

    fn clear (&mut self) {
        self.Primary.clear();
        self.Secondary.clear();
    }

    fn add_primary(&mut self, record: Record) {
        self.Primary.push(record);
        self.Secondary.push(Vec::new());
    }

    fn add_secondary(&mut self, record: Record, idx: usize) {
        // idx: index of primary alignment
        self.Secondary[idx].push(record);
    }

    fn is_empty(&self) -> bool {
        self.Primary.is_empty()
    }

    fn read_id(&self) -> String {
        String::from_utf8(self.Primary[0].qname().to_vec())
                                .unwrap().to_string()
    }

}

// align modified seq to modified reference by meth-minimap2 (rust)
fn map_seq_to_seqs(mod_seq: &[u8], targets: Vec<Vec<u8>>) -> anyResult<Vec::<Mapping>> {
    let mut align_results = Vec::new();

    for (target, ref_seq) in targets.iter().enumerate() {
        // println!(">{}\n{}", target.to_string(), String::from_utf8(seq.to_vec()).unwrap());
        let aligner = Aligner::builder()
        .with_seq_and_id(&ref_seq, target.to_string().as_bytes())
        .unwrap()
        .with_cigar();
        let res = aligner.map(&mod_seq, true, true, None, None);
        
        let res = &res.as_ref().unwrap();
        if res.len() > 0 {
            let best = &res[0];
            align_results.push(best.clone());
        }

    }

    Ok(align_results)
}

fn do_some(au: &AlignmentUnit, seqs: &HashMap<String, String>, 
    bam_header: &HeaderView, min_quality: u8, writer: &mut Writer) {
       
    if !au.is_empty() {
        let read_id = au.read_id();
        'outer: for (p, s) in au.Primary.iter().zip(au.Secondary.iter()) {
            let seq_len = p.seq_len() as u64;
            if p.mapq() >= min_quality || s.len() == 0 {
                writer.write(&p).unwrap();
                if s.len() > 0 {
                    for ss in s {
                        writer.write(&ss).unwrap();
                    }
                }
                continue 'outer;
            } else {
        
                let mut targets = Vec::new();
                let p_parser = RecordParser::parse(&p, &bam_header, seq_len);
                match sub_seq(&seqs, &p_parser.target, 
                                            p_parser.target_start as u64, 
                                            p_parser.target_end as u64) {
                    Some(seq) => {
                        if p_parser.is_reverse {
                            let mut seq = seq.as_bytes().to_vec();
                            mod_seq_complement(&mut seq);
                            targets.push(seq);
                         
                        } else {
                      
                            targets.push(seq.as_bytes().to_vec());
                        }
                    
                    },
                    None => log::debug!("{}: {:?}", read_id, p_parser),
                };
                let mut seq = p.seq().as_bytes().clone();
                let mut query_seq = &mut seq[p_parser.query_start as usize..p_parser.query_end as usize];
               
                if p_parser.is_reverse {
                    // query_seq.reverse();
                    for i in 0..query_seq.len() {
                        query_seq[i] = match query_seq[i] {
                            b'A' => b'T',
                            b'T' => b'A',
                            b'G' => b'C',
                            b'C' => b'G',
                            _ => query_seq[i],
                        }
                    }
                }
                
                let mut mod_seq = query_seq.to_vec();
                if let Ok(mods) = p.basemods_iter() {
                    let mut mod_flag = false; // to test whether change the sequence
                    'inner: for (idx, mod_res) in mods.enumerate() {
                        if let Ok( (position, m) ) = mod_res {
                            if position < p_parser.query_start as i32 || 
                                position >= p_parser.query_end as i32 {
                                continue 'inner;
                            }

                            let qual = qual_to_prob(m.qual as f32);
                            if qual >= 0.75 {
                                mod_flag = true;
                                let new_position = position - p_parser.query_start as i32;
                                mod_seq[new_position as usize] = m.modified_base as u8;
                            }
                        }
                    }
                    if !mod_flag {
                        writer.write(&p).unwrap();
                        if s.len() > 0 {
                            for ss in s {
                                writer.write(&ss).unwrap();
                            }
                        }
                        continue 'outer;
                    }

                } else {
                    writer.write(&p).unwrap();
                    if s.len() > 0 {
                        for ss in s {
                            writer.write(&ss).unwrap();
                        }
                    }
                    continue 'outer;
                }
    
                if p.is_reverse() {
                    mod_seq.reverse();
                }
                for ss in s {
                    let ss_parser = RecordParser::parse(&ss, &bam_header, seq_len);
                    match sub_seq(&seqs, &ss_parser.target, 
                                                ss_parser.target_start as u64, 
                                                ss_parser.target_end as u64) {
                        Some(seq) => {
                            if ss_parser.is_reverse {
                                let mut seq = seq.as_bytes().to_vec();
                                mod_seq_complement(&mut seq); // complementary
                                targets.push(seq);
                                // targets.insert(format!("{}_{}", ss_parser.target, &ss_parser.target_start.to_string()), seq);
                            } else {
                                // targets.insert(format!("{}_{}", ss_parser.target, &ss_parser.target_start.to_string()), seq.as_bytes().to_vec());
                                targets.push(seq.as_bytes().to_vec());
                            }
                        },
                        None => log::debug!("{}: {:?}", read_id, ss_parser),
                        }; 
                }
               
               
             

                let mut align_results = map_seq_to_seqs(&mod_seq, targets).unwrap();
                
                let align_result_len = align_results.len();
                if align_result_len > 0 {
                    align_results.sort_by(|a, b| a.alignment.as_ref().unwrap().nm.cmp(&b.alignment.as_ref().unwrap().nm));
                    
                    let mut a_mapq = match align_result_len == 1 {
                        true => align_results[0].mapq,
                        false => match align_results[0].alignment.as_ref().unwrap().nm == align_results[1].alignment.as_ref().unwrap().nm {
                            true => 0,
                            false => align_results[0].mapq,
                        }
                    };

                    for (a_idx, a) in align_results.iter().enumerate() {
                        if a_idx > 0 {
                            a_mapq = 0;

                        } 
                        // match p.aux(b"NM") {
                        // Ok(value) => {
                        //     println!("{} {} {:.?}", self.read_id(), p.mapq(), value);
                        // },
                        // Err(e) => {}
                        // }
                        // println!("{} {} {} {:?}", self.read_id(), 
                        //                         a.target_name.as_ref().unwrap(),
                        //                         a_mapq,
                        //                        a.alignment.as_ref().unwrap().nm);
                        let mut new_record = if a.target_name.as_ref().unwrap() == &"0".to_string() {
                            Record::from(p.clone())
                            
                        } else {
                            Record::from(s[a.target_name.as_ref().unwrap().parse::<usize>().unwrap() - 1 as usize].clone())
                        };
                        let new_flag = match a_idx {
                            0 => {
                                match new_record.is_reverse() {
                                    true => match new_record.is_supplementary() {
                                        true => 0x910,
                                        false => 0x10,
                                    },
                                    false => match new_record.is_supplementary() {
                                        true => 0x900,
                                        false => 0x0,
                                    }
                                }
                            },
                            _ => {
                                match new_record.is_reverse() {
                                    true => 0x110,
                                    false => 0x100,
                                }
                            }
                        };
                        new_record.set_flags(new_flag);
                        
                        match a_idx {
                            0 => {
                                new_record.set(p.qname(), None, &p.seq().as_bytes(), p.qual());
                            },
                            _ => {
                                new_record.set(p.qname(), None, &[], &[]);
                            }
                        }
                        
                        
                        let new_score = calculate_alignment_score_from_cigar(a.alignment.as_ref().unwrap().cigar.as_ref().unwrap());
                        new_record.set_mapq(a_mapq.try_into().unwrap());
                        new_record.remove_aux(b"NM");
                        new_record.remove_aux(b"AS");
                        
                        new_record.push_aux(b"NM", Aux::I32(a.alignment.as_ref().unwrap().nm));
                        new_record.push_aux(b"AS", Aux::I32(new_score));
                        
                        new_record.remove_aux(b"tp");
                        match a_idx {
                            0 => {
                                new_record.push_aux(b"tp", Aux::Char(b"P"[0]));
                            },
                            _ => {
                                new_record.push_aux(b"tp", Aux::Char(b"S"[0]));
                            }
                        }
                        writer.write(&new_record).unwrap();
                        
                    }
                } 

            }
        }
    }
}



#[derive(Debug)]
pub struct RecordParser {
    pub target: String,
    pub target_start: i64,
    pub target_end: i64,
    pub query_start: i64,
    pub query_end: i64,
    pub query_len: i64,
    pub is_reverse: bool,
}

impl RecordParser {
    fn parse(r: &Record, h: &HeaderView, seq_len: u64) -> RecordParser {
        let cigar = r.cigar();
        let target = str::from_utf8(
                h.tid2name(r.tid() as u32)).unwrap();
        let target_start = r.pos();
        let target_end = cigar.end_pos();
        let mut query_start = cigar.leading_hardclips() + cigar.leading_softclips();
        let mut query_end = cigar.leading_hardclips() + 1 + 
                            cigar.read_pos(
                                target_end as u32 - 1, 
                                false, false).unwrap().unwrap() as i64;
        let query_len = cigar.leading_hardclips() + 
                            i64::try_from(seq_len).unwrap() + 
                            cigar.trailing_hardclips();
        if r.is_reverse() {
            let temp = query_start;
            query_start = query_len - query_end;
            query_end = query_len - temp;
        }

        RecordParser {
            target: target.to_string(),
            target_start: target_start,
            target_end: target_end,
            query_start: query_start,
            query_end: query_end,
            query_len: query_len,
            is_reverse: r.is_reverse(),
        }
    }

}

pub fn sub_seq(seqs: &HashMap<String, String>, 
                chrom: &String, start: u64, end: u64) -> Option<String> {
    let seq = match seqs.get(chrom) {
        Some(seq) => seq,
        None => return None
    };
    let subseq = seq.get(start as usize..end as usize).unwrap();
    Some(subseq.to_string())
}


pub fn parse_read_unit(read_unit: &ReadUnit) -> AlignmentUnit {
    let mut idx: u64 = 0;
    let mut au = AlignmentUnit::new();
    for r in &read_unit.data {
        let read_id = String::from_utf8(r.qname().to_vec())
                                .unwrap().to_string();
        

        if !r.is_secondary() {
            au.add_primary(r.clone());
            idx += 1;
        } else {
            if idx == 0 {
                log::warn!("Secondary alignment `{:?}` could not found primary, \
                            skipped.\
                            The input bam should be sorted by read name.", read_id);
                continue;
            }
            let idx2 = idx - 1;
            au.add_secondary(r.clone(), idx2 as usize);
        }
        
    }
    au
}

pub fn read_bam(input_bam: &String, seqs: &HashMap<String, String>, min_quality: u8,
                output: &String) {
    let mut bam = Reader::from_path(input_bam).unwrap();
    let mut total_reads: u64 = 0;
    let mut total_alignments: u64 = 0;
    let mut total_unmapped: u64 = 0;
    let mut old_read_id = String::from("");
    let mut read_unit = ReadUnit::new();
    let bam_header = Header::from_template(bam.header());
    let bam_header = HeaderView::from_header(&bam_header);

    let header = Header::from_template(&bam_header);
    let mut writer = Writer::from_path(output, &header, bam::Format::Bam).unwrap();

    for r in bam.records() {
        let record = r.unwrap();
        total_alignments += 1;
        
        let read_id = String::from_utf8(record.qname().to_vec())
                                .unwrap().to_string();

        if record.is_unmapped() {
            total_unmapped += 1;
            total_reads += 1;
            continue;
        }
        
        if old_read_id != read_id {
            if old_read_id != "" {
                let au = parse_read_unit(&read_unit);
                do_some(&au, seqs, &bam_header, min_quality, &mut writer);
                
            }

            total_reads += 1;
            read_unit.clear();
            read_unit.data.push(record);
        } else {
            read_unit.data.push(record);
        }
        
        old_read_id = read_id;
    }
    log::info!("Alignments: {:?}", total_alignments);
    log::info!("Umapped: {:?}", total_unmapped);
    log::info!("Total reads: {:?}", total_reads);
    
}

pub fn mapping_unit(index_file: &String, record: &Record) {
    let aligner = Aligner::builder()
        .with_index(index_file, None)
        .unwrap()
        .with_cigar();
    let query = record.seq().as_bytes();
    let res = aligner.map(&query, true, true, None, None);
    
}

pub fn mapping(index_file: &String, input_bam: &String) {
    let mut bam = Reader::from_path(input_bam).unwrap();
    
    let aligner = Aligner::builder()
        .with_index(index_file, None)
        .unwrap()
        .with_cigar();

    for r in bam.records() {
        let record = r.unwrap();
        mapping_unit(index_file, &record);
    }
}
