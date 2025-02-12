#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as anyResult;
// use crossbeam::queue::ArrayQueue;
// use minimap2::*;
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;

use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use rayon::prelude::*;
use std::path::PathBuf;
use std::collections::HashMap;
use std::io::{Cursor, Read as StdRead, Write, Seek};
use std::process::{Command, Stdio};
use std::rc::Rc;
use std::str;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{ Duration, Instant };
use tempfile::NamedTempFile;

use crate::fastx::Fastx;
use crate::methy::qual_to_prob;

enum WorkQueue<T> {
    Work(T),
    Result(T),
}

#[derive(Clone, Debug)]
pub struct SeqRecord {
    pub name: String,
    pub seq: Vec<u8>,
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
fn calculate_alignment_score_from_cigar_string_view(cigar: CigarStringView) -> anyResult<i32> {
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
   
    let score: i32 = score.try_into().unwrap();
    Ok(score)
}

pub fn mod_seq_complement(seq: &mut Vec<u8>) {
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
// fn map_seq_to_seqs(mod_seq: &[u8], targets: Vec<Vec<u8>>) -> anyResult<Vec::<Mapping>> {
//     let mut align_results = Vec::new();

//     for (target, ref_seq) in targets.iter().enumerate() {
//         // println!(">{}\n{}", target.to_string(), String::from_utf8(seq.to_vec()).unwrap());
//         let aligner = Aligner::builder()
//         .with_seq_and_id(&ref_seq, target.to_string().as_bytes())
//         .unwrap()
//         .with_cigar();
//         let res = aligner.map(&mod_seq, true, true, None, None);
        
//         let res = &res.as_ref().unwrap();
//         if res.len() > 0 {
//             let best = &res[0];
//             align_results.push(best.clone());
//         }

//     }

//     Ok(align_results)
// }

// align modified seq to modified reference by meth-minimap2 (command)
fn map_seq_to_seqs_by_command(read_id: &String,mod_seq: &[u8], targets: Vec<SeqRecord>, 
                                ) -> anyResult<(Vec::<Record>, HeaderView)> {
    let mut align_results = Vec::new();
    let start_time = Instant::now();
    // let read_id_tmp = read_id.clone();
    // let mod_seq = Arc::new(mod_seq.to_vec());
    // let handle = thread::spawn(move || {
    //     let mut temp_file = NamedTempFile::new().unwrap();
    //     temp_file.write_all(format!(">{}\n", read_id_tmp).as_bytes()).unwrap();
    //     temp_file.write_all(&mod_seq).unwrap();
    //     temp_file
    // });
    let mut temp_file = NamedTempFile::new().unwrap();
    temp_file.write_all(format!(">{}\n", read_id).as_bytes()).unwrap();
    temp_file.write_all(mod_seq).unwrap();
   
    // temp_file.seek(std::io::SeekFrom::Start(0)).unwrap();
    // let mut contents = String::new();
    // temp_file.read_to_string(&mut contents)?;
    // println!("Temporary file contents: {}", contents);
    // ">i\nseq\n"
    let mut header = Header::new();
   
    let mut record = HeaderRecord::new(b"HD");
    record.push_tag(b"VN", &"1.6") ;
    record.push_tag(b"SO", &"coordinate");    
    header.push_record(&record);
    for i in 0..targets.len() {
      
        let mut record = HeaderRecord::new(b"SQ");
        record.push_tag(b"SN", &targets[i].name);
        record.push_tag(b"LN", &(targets[i].seq.len() as u64));   
        header.push_record(&record);
    }
    let mut record = HeaderRecord::new(b"RG");
    record.push_tag(b"ID", &read_id);
    header.push_record(&record);
    let header_view = HeaderView::from_header(&header);

    let ref_seqs = targets.iter().map(|x| {
        format!(">{}\n{}\n", x.name, String::from_utf8(x.seq.to_vec()).unwrap())
    }).collect::<Vec<String>>().join("");

    // let temp_file = handle.join().unwrap();

    let mut child = Command::new("minimap2")
                        .arg("-x")
                        .arg("map-ont")
                        .arg("-c")
                        .arg("-a")
                        .arg("-")
                        .arg(temp_file.path())
                        .stdin(Stdio::piped())
                        .stdout(Stdio::piped())
                        .stderr(Stdio::piped())
                        .spawn()
                        .unwrap();
    {
        let stdin = child.stdin.as_mut().expect("Failed to open stdin");
        stdin.write_all(ref_seqs.as_bytes()).expect("Failed to write to stdin");
    }    

    let output = child.wait_with_output().expect("Failed to read stdout");
    let output = String::from_utf8_lossy(&output.stdout).to_string();
   
    let end_time = Instant::now();
    let cpu_time = end_time - start_time;
    // println!("CPU time: {:.10}s", cpu_time.as_secs_f64());
    let lines = output.lines();
    
    for line in lines {
        if line.starts_with("@") {
            continue;
        }
        let record = Record::from_sam(&header_view, &line.as_bytes()).unwrap();
        if record.is_unmapped() {
            continue;
        }
        align_results.push(record);
    }

    Ok((align_results, header_view))
               
}

fn map_seq_to_seqs_by_edit_distance(read_id: &String, 
                mod_seq: &[u8], targets: Vec<SeqRecord>,

            ) -> anyResult<Vec::<(usize, u32)>> {
    
    use bio::alignment::distance::simd::*;

    
    let mut align_results: Vec::<(usize, u32)> = Vec::new();

    // only retain the C and m in mod seq, and convert to str
    let mod_seq = mod_seq.iter().filter(|x| {
        match x {
            b'C' | b'm' => true,
            _ => false,
        }
    }).collect::<Vec<&u8>>();

    // convert &u8 to u8 
    let mod_seq = mod_seq.iter().map(|x| {
        **x
    }).collect::<Vec<u8>>();


    for i in 0..targets.len() {
        let target = &targets[i];
        // only retain the C and m in target seq
        let target_seq = target.seq.iter().filter(|x| {
            match x {
                b'C' | b'm' => true,
                _ => false,
            }
        }).collect::<Vec<&u8>>();
        let target_seq = target_seq.iter().map(|x| {
            **x
        }).collect::<Vec<u8>>();


        let distance = levenshtein(&mod_seq, &target_seq);
        align_results.push((i, distance));
    }

    // sort align_results by distance
    align_results.sort_by(|a, b| a.1.cmp(&b.1));


    Ok(align_results)
}

fn map_seq_to_seqs_by_local(read_id: &String, 
                mod_seq: &[u8], targets: Vec<SeqRecord>,

            ) -> anyResult<Vec::<(usize, i32)>> {
    
    let mut align_results: Vec::<(usize, i32)> = Vec::new();
    let match_fn = |a: u8, b: u8| {
        if a == b'm' {
            if b == b'm' {
                4i32
            } else {
                -5i32
            }
        } else {
            if b == b'm' {
                -5i32
            } else {
                if a == b {
                    2i32
                } else {
                    -5i32
                }
            
            }
        }
    };
    let scoring = Scoring::new(-5, -2, &match_fn).xclip(-5).yclip(-5);


    // only retain the C and m in mod seq, and convert to str
    // let mut mod_seq = mod_seq.iter().filter(|x| {
    //     match x {
    //         b'C' | b'm' => true,
    //         _ => false,
    //     }
    // }).collect::<Vec<&u8>>();

    // convert &u8 to u8 
    let mod_seq = mod_seq.iter().map(|x| {
        *x
    }).collect::<Vec<u8>>();
    let mod_seq = String::from_utf8(mod_seq.to_vec()).unwrap();

    for i in 0..targets.len() {
        let target = &targets[i];
        // only retain the C and m in target seq
        // let mut target_seq = target.seq.iter().filter(|x| {
        //     match x {
        //         b'C' | b'm' => true,
        //         _ => false,
        //     }
        // }).collect::<Vec<&u8>>();
        let mut target_seq = target.seq.iter().map(|x| {
            *x
        }).collect::<Vec<u8>>();
        let target_seq = String::from_utf8(target_seq.to_vec()).unwrap();

        let mut aligner = Aligner::with_capacity_and_scoring(mod_seq.len(), target_seq.len(), scoring);
        let alignment = aligner.semiglobal(mod_seq.as_bytes(), target_seq.as_bytes());
        let score = alignment.score;

        align_results.push((i, score));
    }

    // sort align_results by score ascending = False
    align_results.sort_by(|a, b| a.1.cmp(&b.1).reverse());
  

    Ok(align_results)
}

fn map_seq_to_seqs_by_banded_semiglobal(read_id: &String, 
                mod_seq: &[u8], targets: Vec<SeqRecord>,

            ) -> anyResult<Vec::<(usize, i32)>> {
    use bio::alignment::pairwise::banded::*;
    use bio::alignment::sparse::hash_kmers;
    use std::iter::repeat;
    let mut align_results: Vec::<(usize, i32)> = Vec::new();
    let scoring = Scoring {
        match_fn: |a: u8, b: u8| if a == b { 2i32 } else { -5i32 },
        match_scores: Some((2, -5)),
        gap_open: -2,
        gap_extend: -1,
        xclip_prefix: -5,
        xclip_suffix: -5,
        yclip_prefix: -5,
        yclip_suffix: -5,
    };
    
    // only retain the C and m in mod seq, and convert to str
    // let mut mod_seq = mod_seq.iter().filter(|x| {
    //     match x {
    //         b'C' | b'm' => true,
    //         _ => false,
    //     }
    // }).collect::<Vec<&u8>>();

    // convert &u8 to u8 
    let mod_seq = mod_seq.iter().map(|x| {
        *x
    }).collect::<Vec<u8>>();
    let mod_seq_kmers_hash = hash_kmers(&mod_seq, 15);

    for i in 0..targets.len() {
        let target = &targets[i];
        // only retain the C and m in target seq
        // let mut target_seq = target.seq.iter().filter(|x| {
        //     match x {
        //         b'C' | b'm' => true,
        //         _ => false,
        //     }
        // }).collect::<Vec<&u8>>();
        let mut target_seq = target.seq.iter().map(|x| {
            *x
        }).collect::<Vec<u8>>();
        

        let mut aligner = Aligner::with_capacity_and_scoring(target_seq.len(), mod_seq.len(), scoring, 15, 5);
        let alignment = aligner.semiglobal_with_prehash(&target_seq, &mod_seq, &mod_seq_kmers_hash);
        let score = alignment.score;
  
        align_results.push((i, score));
    }

    // sort align_results by score ascending = False
    align_results.sort_by(|a, b| a.1.cmp(&b.1).reverse());


    Ok(align_results)
}


fn map_seq_to_seqs_by_block_aligner(read_id: &String, 
                mod_seq: &[u8], targets: Vec<SeqRecord>,

            ) -> anyResult<Vec::<(usize, i32)>> {

    use block_aligner::{cigar::*, scan_block::*, scores::*};
    let min_block_size = 32;
    let max_block_size = 256;
    let gaps = Gaps { open: -2, extend: -1 };

    let mut align_results: Vec::<(usize, i32)> = Vec::new();
    let matrix = ByteMatrix::new_simple(4, -2);
    
    // only retain the C and m in mod seq, and convert to str
    let mut mod_seq = mod_seq.iter().filter(|x| {
        match x {
            b'C' | b'm' => true,
            _ => false,
        }
    }).collect::<Vec<&u8>>();
    // convert m to G 
    for i in 0..mod_seq.len() {
        if mod_seq[i] == &b'm' {
            mod_seq[i] = &b'G';
        }
    }

    // convert &u8 to u8 
    let mut mod_seq = mod_seq.iter().map(|x| {
        **x
    }).collect::<Vec<u8>>();
    
    let r = PaddedBytes::from_bytes::<ByteMatrix>(&mod_seq, max_block_size);

    for i in 0..targets.len() {
        let target = &targets[i];
        // only retain the C and m in target seq
        let mut target_seq = target.seq.iter().filter(|x| {
            match x {
                b'C' | b'm' => true,
                _ => false,
            }
        }).collect::<Vec<&u8>>();

        // convert m to G
        for i in 0..target_seq.len() {
            if target_seq[i] == &b'm' {
                target_seq[i] = &b'G';
            }
        }
        let mut target_seq = target_seq.iter().map(|x| {
            **x
        }).collect::<Vec<u8>>();
        
        let q = PaddedBytes::from_bytes::<ByteMatrix>(&target_seq, max_block_size);

        let mut a = Block::<true, false>::new(q.len(), r.len(), max_block_size);
        a.align(&q, &r, &matrix, gaps, min_block_size..=max_block_size, 0);

        let res = a.res();

        if read_id == &String::from("000c27e0-d06a-45f1-b108-38b81c661ee1_0_9_10") {
            println!("{}: {:?}", read_id, res);
        }
        align_results.push((i, res.score));
    }

    // sort align_results by score ascending = False
    align_results.sort_by(|a, b| a.1.cmp(&b.1).reverse());


    Ok(align_results)
}



fn map_seq_to_seqs_by_sparse(read_id: &String, 
                mod_seq: &[u8], targets: Vec<SeqRecord>,

            ) -> anyResult<Vec::<(usize, u32)>> {
    use bio::alignment::sparse::*;

    let mut align_results: Vec::<(usize, u32)> = Vec::new();
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    // only retain the C and m in mod seq, and convert to str
    let mut mod_seq = mod_seq.iter().filter(|x| {
        match x {
            b'C' | b'm' => true,
            _ => false,
        }
    }).collect::<Vec<&u8>>();

    // convert &u8 to u8 
    let mut mod_seq = mod_seq.iter().map(|x| {
        **x
    }).collect::<Vec<u8>>();
    let mod_seq = String::from_utf8(mod_seq.to_vec()).unwrap();

    for i in 0..targets.len() {
        let target = &targets[i];
        // only retain the C and m in target seq
        let mut target_seq = target.seq.iter().filter(|x| {
            match x {
                b'C' | b'm' => true,
                _ => false,
            }
        }).collect::<Vec<&u8>>();
        let mut target_seq = target_seq.iter().map(|x| {
            **x
        }).collect::<Vec<u8>>();
        let target_seq = String::from_utf8(target_seq.to_vec()).unwrap();

        let matches = find_kmer_matches(mod_seq.as_bytes(), target_seq.as_bytes(), 15);
        let sparse_al = lcskpp(&matches, 15);
        let score = sparse_al.score;
  
        align_results.push((i, score));
    }

    // sort align_results by score ascending = False
    align_results.sort_by(|a, b| a.1.cmp(&b.1).reverse());


    Ok(align_results)
}

fn do_some(au: &AlignmentUnit, seqs: &HashMap<String, String>, 
    bam_header: &HeaderView, min_quality: u8, min_prob: f32,
    writer: &mut Writer) {
       
    if !au.is_empty() {
        let read_id = au.read_id();
        
        'outer: for (p, s) in au.Primary.iter().zip(au.Secondary.iter()) {
            let mut read_idx = 0;
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
        
                let mut targets: Vec::<SeqRecord> = Vec::new();
                let p_parser = RecordParser::parse(&p, &bam_header, seq_len);
                match sub_seq(&seqs, &p_parser.target, 
                                            p_parser.target_start as u64, 
                                            p_parser.target_end as u64) {
                    Some(seq) => {
                        if p_parser.is_reverse {
                            let mut seq = seq.as_bytes().to_vec();
                            mod_seq_complement(&mut seq);
                            targets.push(SeqRecord {
                                // name: p_parser.target.clone(),
                                name: read_idx.to_string(),
                                seq: seq,
                            });
                         
                        } else {
                      
                            targets.push(SeqRecord {
                                // name: p_parser.target.clone(), 
                                name: read_idx.to_string(),
                                seq: seq.as_bytes().to_vec(),
                        });
                        }
                    
                    },
                    None => log::debug!("{}: {:?}", read_id, p_parser),
                };
                read_idx += 1;
                
                let mut seq = p.seq().as_bytes().clone();
                let query_seq = &mut seq[p_parser.query_start as usize..p_parser.query_end as usize];
               
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
                            if qual >= min_prob {
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
                    mod_flag = false;

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
                for (ss_idx, ss) in s.iter().enumerate() {

                    let ss_parser = RecordParser::parse(&ss, &bam_header, seq_len);
                    match sub_seq(&seqs, &ss_parser.target, 
                                                ss_parser.target_start as u64, 
                                                ss_parser.target_end as u64) {
                        Some(seq) => {
                            if ss_parser.is_reverse {
                                let mut seq = seq.as_bytes().to_vec();
                                mod_seq_complement(&mut seq); // complementary
                                
                                targets.push(SeqRecord {
                                    // name: ss_parser.target.clone(),
                                    name: read_idx.to_string(),
                                    seq: seq,
                                });
                                // targets.insert(format!("{}_{}", ss_parser.target, &ss_parser.target_start.to_string()), seq);
                            } else {
                                // targets.insert(format!("{}_{}", ss_parser.target, &ss_parser.target_start.to_string()), seq.as_bytes().to_vec());
                                targets.push(SeqRecord {
                                    // name: ss_parser.target.clone(),
                                    name: read_idx.to_string(),
                                    seq: seq.as_bytes().to_vec(),
                                });
                            }
                        },
                        None => log::debug!("{}: {:?}", read_id, ss_parser),
                        }; 
                        read_idx += 1;
                    }
                
                let align_results =  map_seq_to_seqs_by_banded_semiglobal(&read_id, &mod_seq, targets.clone()).unwrap();
                let align_results_len = align_results.len();
                if align_results_len > 0 {
                    if align_results[0].1 == align_results[1].1 {
                        continue 'outer;
                    }
                    for (idx, a) in align_results.iter().enumerate() {
                        let a_idx = a.0;
                      
                        let mut new_record = if a_idx == 0 {
                            Record::from(p.clone())
                            
                        } else {
                            Record::from(s[a_idx - 1 ].clone())
                        };

                        if idx == 0 {
                            new_record.push_aux(b"RA", Aux::String(&"y")).unwrap();
                            // set mapq to 2 
                            new_record.set_mapq(2);
                            if new_record.is_reverse() {
                                if new_record.is_supplementary() {
                                    new_record.set_flags(0x910);
                                } else {
                                    new_record.set_flags(0x10);
                                }
                            } else {
                                if new_record.is_supplementary() {
                                    new_record.set_flags(0x900);
                                } else {
                                    new_record.set_flags(0x0);
                                }
                            }


                        } else {
                            // new_record.push_aux(b"RA", Aux::String(&"n")).unwrap();
                            new_record.set_mapq(0);
                            if new_record.is_reverse() {
                                if new_record.is_supplementary() {
                                    new_record.set_flags(0x910);
                                } else {
                                    new_record.set_flags(0x110);
                                }
                            } else {
                                if new_record.is_supplementary() {
                                    new_record.set_flags(0x1000);
                                } else {
                                    new_record.set_flags(0x100);
                                }
                                
                            }
                        }

                        new_record.remove_aux(b"NM").unwrap();
                        
                        new_record.remove_aux(b"tp").unwrap();
                        // new_record.push_aux(b"NM", Aux::U32(a.1.try_into().unwrap())).unwrap();
                      
                      
                        new_record.remove_aux(b"AS").unwrap();
                        new_record.push_aux(b"AS", Aux::I32(a.1)).unwrap();
                        match idx {
                            0 => {
                                new_record.push_aux(b"tp", Aux::Char(b"P"[0])).unwrap();
                            },
                            _ => {
                                new_record.push_aux(b"tp", Aux::Char(b"S"[0])).unwrap();
                            }
                        }
                        
                        writer.write(&new_record).unwrap();
                        
                    }


                // let (align_results, new_header_view) = map_seq_to_seqs_by_command(&read_id, &mod_seq, targets).unwrap();
                
                // let align_result_len = align_results.len();
                // if align_result_len > 0 {

                //     for (a_idx, a) in align_results.iter().enumerate() {
                //         let a_parser = RecordParser::parse(&a, &new_header_view, seq_len);

                     
                //         let mut new_record = if a_parser.target == "0" {
                //             Record::from(p.clone())
                            
                //         } else {
                //             Record::from(s[a_parser.target.parse::<usize>().unwrap() - 1 as usize].clone())
                //         };
                        
                //         if a_idx == 0 {
                //             new_record.push_aux(b"RA", Aux::String(&"y")).unwrap();
                //         }
                //         // add RA tag, RA:Z:
                        

                //         let cigar_string = new_record.cigar().take();
                //         let cigar_string: Option<&CigarString> = Some(&cigar_string);
                //         match !a.is_secondary() {
                //             true => {
                //                 new_record.set(p.qname(), cigar_string, &p.seq().as_bytes(), p.qual());
                //                 if new_record.is_reverse() {
                //                     if new_record.is_supplementary() {
                //                         new_record.set_flags(0x910);
                //                     } else {
                //                         new_record.set_flags(0x10);
                //                     }
                //                 } else {
                //                     if new_record.is_supplementary() {
                //                         new_record.set_flags(0x900);
                //                     } else {
                //                         new_record.set_flags(0x0);
                //                     }
                //                 }
                //             },
                //             false => {
                //                 if new_record.is_reverse() {
                //                     if new_record.is_supplementary() {
                //                         new_record.set_flags(0x910);
                //                     } else {
                //                         new_record.set_flags(0x110);
                //                     }
                //                 } else {
                //                     if new_record.is_supplementary() {
                //                         new_record.set_flags(0x1000);
                //                     } else {
                //                         new_record.set_flags(0x100);
                //                     }
                                    
                //                 }
                //                 new_record.set(p.qname(), cigar_string, &[], &[]);
                //             }
                //         }

                //         let new_score: u32 = match a.aux(b"AS") {
                //             Ok(value) => {
                //                 match value {
                //                     Aux::U8(v) => v.try_into().unwrap(),
                //                     Aux::U16(v) => v.try_into().unwrap(),
                //                     Aux::U32(v) => v.try_into().unwrap(),
                //                     _ => 0,
                //             }
                //         },
                //             Err(e) => calculate_alignment_score_from_cigar_string_view(a.cigar()).unwrap().try_into().unwrap_or(0),
                //         };
                //         let new_nm: u32 = match a.aux(b"NM") {
                //             Ok(value) => {
                //                 match value {
                //                     Aux::U8(v) => v.try_into().unwrap(),
                //                     Aux::U16(v) => v.try_into().unwrap(),
                //                     Aux::U32(v) => v.try_into().unwrap(),
                //                     _ => 0,
                //                 }
                //                 },
                //                 Err(e) => 0,
                //         };
                        
                //         new_record.set_mapq(a.mapq().try_into().unwrap());
                //         new_record.remove_aux(b"NM").unwrap();
                //         new_record.remove_aux(b"AS").unwrap();
                //         new_record.remove_aux(b"tp").unwrap();
                //         new_record.push_aux(b"NM", Aux::U32(new_nm)).unwrap();
                //         new_record.push_aux(b"AS", Aux::U32(new_score)).unwrap();
                //         match a_idx {
                //             0 => {
                //                 new_record.push_aux(b"tp", Aux::Char(b"P"[0])).unwrap();
                //             },
                //             _ => {
                //                 new_record.push_aux(b"tp", Aux::Char(b"S"[0])).unwrap();
                //             }
                //         }

                        
                    


                //         writer.write(&new_record).unwrap();
                //     }
                
                
                } else {
                    writer.write(&p).unwrap();
                    if s.len() > 0 {
                        for ss in s {
                            writer.write(&ss).unwrap();
                        }
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
    pub fn parse(r: &Record, h: &HeaderView, seq_len: u64) -> RecordParser {
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


fn parse_read_unit(read_unit: &ReadUnit) -> AlignmentUnit {
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
                min_prob: f32, output: &String) {
    let mut bam = Reader::from_path(input_bam).unwrap();
    let _ = bam.set_threads(8);
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
                do_some(&au, seqs, &bam_header, min_quality, 
                        min_prob, &mut writer);
                
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

// pub fn mapping_unit(index_file: &String, record: &Record) {
//     let aligner = Aligner::builder()
//         .with_index(index_file, None)
//         .unwrap()
//         .with_cigar();
//     let query = record.seq().as_bytes();
//     let res = aligner.map(&query, true, true, None, None);
    
// }

// pub fn mapping(index_file: &String, input_bam: &String) {
//     let mut bam = Reader::from_path(input_bam).unwrap();
    
//     let aligner = Aligner::builder()
//         .with_index(index_file, None)
//         .unwrap()
//         .with_cigar();

//     for r in bam.records() {
//         let record = r.unwrap();
//         mapping_unit(index_file, &record);
//     }
// }

// 