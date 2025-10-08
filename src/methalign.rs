use anyhow::Result as anyResult;
use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use noodles::bed::record;
use memchr::memchr_iter;
use indicatif::{ProgressBar, ProgressStyle, ProgressDrawTarget};
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use rayon::prelude::*;
use std::path::PathBuf;
use std::collections::{HashMap, HashSet};
use std::io::{Cursor, Read as StdRead, Write, Seek};
use std::process::{Command, Stdio};
use std::rc::Rc;
use std::str;
use std::sync::{Arc, Mutex};
use std::thread;
use std::fs::File;
use std::path::Path;
use std::io::BufRead;
use std::time::{ Duration, Instant };
use tempfile::NamedTempFile;
use polars::prelude::*;

use crate::core::BaseTable;
use crate::core::common_reader;
use crate::fastx::Fastx;
use crate::methy::qual_to_prob;


pub fn parse_bedgraph(bedgraph: &String, cov_cutoff: f64) -> anyResult<HashMap<String, HashSet<i64>>> {
    let mut fh = common_reader(bedgraph);
    let mut map: HashMap<String, HashSet<i64>> = HashMap::new();
    let mut line = String::new();

    loop {
        line.clear();
        let n = fh.read_line(&mut line)?;
        if n == 0 { break; } 
        if line.is_empty() { continue; }
      
        if line.as_bytes().first() == Some(&b'#') { continue; }

        let mut it = line.split_ascii_whitespace();
        let ctg = match it.next() { Some(x) => x, None => continue };
        let start_str = match it.next() { Some(x) => x, None => continue };
        let _end = it.next(); 
        let cov_str = match it.next() { Some(x) => x, None => continue };

        let cov: f64 = match cov_str.parse() { Ok(v) => v, Err(_) => continue };
        if cov < cov_cutoff { continue; }
        let s: i64 = match start_str.parse() { Ok(v) => v, Err(_) => continue };

        map.entry(ctg.to_owned()).or_default().insert(s);
    }

    log::info!(
        "Load {} methylation sites from bedgraph {} after filter by ref_prob_cutoff {}",
        map.values().map(|s| s.len()).sum::<usize>(),
        bedgraph,
        cov_cutoff
    );
    Ok(map)
}

// split records to primary and supplementary with multiple secondary alignments, respectively
pub fn split_records(records: Vec<Record>) -> Vec<Vec<Record>> {
    let mut result: Vec<Vec<Record>> = Vec::new();

    let mut flag = 0;
    let mut current_group: Vec<Record> = Vec::new();
    for record in records.into_iter() {
        let is_supplementary = record.is_supplementary();
        let is_secondary = record.is_secondary();
        let is_primary = !is_supplementary && !is_secondary;
        if is_primary || is_supplementary {
            if !current_group.is_empty() {
                result.push(std::mem::take(&mut current_group));
            }
            current_group.push(record);
        } else if is_secondary {
            current_group.push(record);
        }

    }

    if !current_group.is_empty() {
        result.push(current_group);
    }

    result
    

}

pub fn get_as(record: &Record) -> Option<i32> {
    match record.aux(b"AS") {
        Ok(Aux::I32(val)) => Some(val),
        Ok(Aux::U8(val)) => Some(val as i32),
        Ok(Aux::U16(val)) => Some(val as i32),
        Ok(Aux::U32(val)) => Some(val as i32),
        _ => None,
    }
}

pub fn is_contain_methylation(record: &Record) -> bool {
    record.aux(b"MM").is_ok() && record.aux(b"ML").is_ok()
}

fn revcomp(bs: &[u8]) -> Vec<u8> {
    fn comp(b: u8) -> u8 {
        match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        }
    }
    let mut v: Vec<u8> = bs.iter().map(|b| comp(*b)).collect();
    v.reverse();
    v
}

fn query_alignment_termini(rec: &bam::Record) -> (i64, i64) {
    let cig = rec.cigar();
    let mut qstart = 0i64;
    let mut qlen_aln = 0i64;
    let mut left = true;
    for op in cig.iter() {
        match op {
            Cigar::HardClip(len) if left => { qstart += *len as i64; }
            Cigar::SoftClip(len) if left => { qstart += *len as i64; }
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => { qlen_aln += *len as i64; left = false; }
            Cigar::Ins(len) => { qlen_aln += *len as i64; left = false; }
            Cigar::Del(_) | Cigar::RefSkip(_) => { left = false; }
            Cigar::SoftClip(_) => { /* trailing soft */ }
            Cigar::HardClip(len) => { /* trailing hard will add to end */ let _ = len; }
            _ => { left = false; }
        }
    }
    // trailing hard clip
    let mut right_h = 0i64;
    if let Some(last) = cig.last() { if let Cigar::HardClip(len) = last { right_h = *len as i64; } }
    (qstart, qstart + qlen_aln + right_h)
}

fn condense_cigar(cig: &bam::record::CigarStringView) -> String {
    let mut lsoft = 0i64;
    let mut rsoft = 0i64;
    let mut hard_l = 0i64;
    let mut hard_r = 0i64;
    let mut m: i64 = 0;
    let mut ins: i64 = 0;
    let mut del: i64 = 0;
    let ops = cig.iter().collect::<Vec<_>>();
    for (i, op) in ops.iter().enumerate() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => m += *l as i64,
            Cigar::Ins(l) => ins += *l as i64,
            Cigar::Del(l) => del += *l as i64,
            Cigar::SoftClip(l) => {
                if i == 0 { lsoft += *l as i64; }
                if i == ops.len()-1 { rsoft += *l as i64; }
            }
            Cigar::HardClip(l) => {
                if i == 0 { hard_l += *l as i64; }
                if i == ops.len()-1 { hard_r += *l as i64; }
            }
            _ => {}
        }
    }
    let mut s = String::new();
    if lsoft + hard_l > 0 { s.push_str(&format!("{}S", lsoft + hard_l)); }
    let match_len = if ins > del { m + del } else { m + ins };
    s.push_str(&format!("{}M", match_len.max(0)));
    let indel = del - ins;
    if indel > 0 {
        s.push_str(&format!("{}D", indel));
    } else if indel < 0 {
        s.push_str(&format!("{}I", -indel));
    }
    if rsoft + hard_r > 0 { s.push_str(&format!("{}S", rsoft + hard_r)); }
    s
}

fn reconstruct_sa(rec: &bam::Record, hdr: &HeaderView) -> anyResult<String> {
    let tid = rec.tid();
    let rname = std::str::from_utf8(hdr.tid2name(tid as u32))?.to_string();
    let pos1 = rec.pos() + 1;
    let strand = if rec.is_reverse() { "-" } else { "+" };
    let cigar = condense_cigar(&rec.cigar());
    let mapq = rec.mapq();
    let nm = match rec.aux(b"NM") {
        Ok(Aux::I8(v)) => v as i64,
        Ok(Aux::I16(v)) => v as i64,
        Ok(Aux::I32(v)) => v as i64,
        Ok(Aux::U8(v)) => v as i64,
        Ok(Aux::U16(v)) => v as i64,
        Ok(Aux::U32(v)) => v as i64,
        _ => 0,
    };
    Ok(format!("{},{},{},{},{},{};", rname, pos1, strand, cigar, mapq, nm))
}

fn aligned_pairs(rec: &bam::Record) -> Vec<(i64, i64)> {
    let mut v = Vec::new();
    let mut q = 0i64;

    for op in rec.cigar().iter() {
        match op {
            Cigar::SoftClip(l) | Cigar::HardClip(l) => q += *l as i64,
            _ => break,
        }
    }
    let mut r = rec.pos();
    for op in rec.cigar().iter() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                for _ in 0..(*l as i64) {
                    v.push((q, r));
                    q += 1;
                    r += 1;
                }
            }
            Cigar::Ins(l) => { q += *l as i64; }
            Cigar::Del(l) | Cigar::RefSkip(l) => { r += *l as i64; }
            Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
            _ => {}
        }
    }
    v
}

// fn get_5mc_sites_from_read(read_seq: &[u8], mm: Option<String>, ml: Option<Vec<u8>>, prob_cutoff: u8) -> (HashMap<i64, usize>, Vec<i16>) {
//     let mut c_positions: Vec<i64> = Vec::new();
//     for (i, b) in read_seq.iter().enumerate() {
//         if *b == b'C' { c_positions.push(i as i64); }
//     }
//     let mut index_map: HashMap<i64, usize> = HashMap::new();
//     for (idx, pos) in c_positions.iter().enumerate() {
//         index_map.insert(*pos, idx);
//     }
//     let mut meth: Vec<i16> = vec![0; c_positions.len()];

//     let mm_s = match mm { Some(s) => s, None => return (index_map, meth) };
//     let ml_v = match ml { Some(v) => v, None => return (index_map, meth) };

//     let parts: Vec<&str> = mm_s.trim_end_matches(';').split(';').collect();
//     let mut offset_items = 0usize;
//     let mut shifts: Option<Vec<i64>> = None;
//     for part in parts {
//         let items: Vec<&str> = part.split(',').collect();
//         if items.is_empty() { continue; }
//         if items[0].contains("C+m") {
//             let mut v: Vec<i64> = Vec::new();
//             for x in &items[1..] {
//                 if x.is_empty() { continue; }
//                 if let Ok(n) = x.parse::<i64>() { v.push(n); }
//             }
//             shifts = Some(v);
//             break;
//         } else {
//             offset_items += items.len();
//         }
//     }
//     let Some(shifts_v) = shifts else { return (index_map, meth) };
  
//     let slice = if offset_items < ml_v.len() { &ml_v[offset_items..] } else { &[] };
//     let ml_slice = &slice[..shifts_v.len().min(slice.len())];

//     let mut idx_sum: i64 = 0;
//     for (n, sh) in shifts_v.iter().enumerate() {
//         idx_sum += *sh;
//         let pos_idx = (n as i64 + idx_sum) as usize;
//         if pos_idx < meth.len() && n < ml_slice.len() {
//             if ml_slice[n] >= prob_cutoff {
//                 meth[pos_idx] = 1;
//             }
//         }
//     }
//     (index_map, meth)
// }

// fn get_5mc_sites_from_read(
//     read_seq: &[u8],
//     mm: Option<String>,
//     ml: Option<Vec<u8>>,
//     prob_cutoff: u8
// ) -> (Vec<i64>, Vec<i16>) {
//     // 读中 'C' 的位置（0-based，按升序）
//     // let mut c_positions: Vec<i64> = Vec::new();
//     // for (i, b) in read_seq.iter().enumerate() {
//     //     if *b == b'C' {
//     //         c_positions.push(i as i64);
//     //     }
//     // }
//     let c_positions: Vec<i64> = memchr_iter(b'C', read_seq).map(|i| i as i64).collect();
//     // filter C position which is not the CG context
//     // let mut c_positions_filt: Vec<i64> = Vec::new();
//     // for &pos in &c_positions {
//     //     if (pos + 1) < read_seq.len() && read_seq[(pos + 1) as usize] == b'G' {
//     //         c_positions_filt.push(pos);
//     //     }
//     // }
//     // let c_positions = c_positions_filt;
//     let mut meth: Vec<i16> = vec![0; c_positions.len()];

//     let mm_s = match mm { Some(s) => s, None => return (c_positions, meth) };
//     let ml_v = match ml { Some(v) => v, None => return (c_positions, meth) };

//     let parts: Vec<&str> = mm_s.trim_end_matches(';').split(';').collect();
//     let mut offset_items = 0usize;
//     let mut shifts: Option<Vec<i64>> = None;
//     for part in parts {
//         let items: Vec<&str> = part.split(',').collect();
//         if items.is_empty() { continue; }
//         if items[0].contains("C+m") {
//             let mut v: Vec<i64> = Vec::new();
//             for x in &items[1..] {
//                 if x.is_empty() { continue; }
//                 if let Ok(n) = x.parse::<i64>() { v.push(n); }
//             }
//             shifts = Some(v);
//             break;
//         } else {
//             offset_items += items.len();
//         }
//     }
//     let Some(shifts_v) = shifts else { return (c_positions, meth) };

//     let slice = if offset_items < ml_v.len() { &ml_v[offset_items..] } else { &[] };
//     let ml_slice = &slice[..shifts_v.len().min(slice.len())];

//     let mut idx_sum: i64 = 0;
//     for (n, sh) in shifts_v.iter().enumerate() {
//         idx_sum += *sh;
//         let pos_idx = (n as i64 + idx_sum) as usize;
//         if pos_idx < meth.len() && n < ml_slice.len() {
//             if ml_slice[n] >= prob_cutoff {
//                 meth[pos_idx] = 1;
//             }
//         }
//     }
//     (c_positions, meth)
// }

// fn get_5mc_sites_from_read(
//     read_seq: &[u8],
//     mm: Option<&str>,
//     ml: Option<&[u8]>,
//     prob_cutoff: u8,
// ) -> (Vec<i64>, Vec<i16>) {
//     println!("{:?}, {:?}", mm, ml);
//     let mm_s = match mm { Some(s) => s, None => return (Vec::new(), Vec::new()) };
//     let ml_v = ml.unwrap_or(&[]);
//     println!("MM: {:?}, ML: {:?}", mm_s, ml_v);

//     let pos_c: Vec<i64> = memchr_iter(b'C', read_seq).map(|i| i as i64).collect();
//     let pos_g: Vec<i64> = memchr_iter(b'G', read_seq).map(|i| i as i64).collect();
//     let mut flag_c: Vec<i16> = vec![0; pos_c.len()];
//     let mut flag_g: Vec<i16> = vec![0; pos_g.len()];

//     let mut ml_idx: usize = 0;

//     for part in mm_s.trim_end_matches(';').split(';') {
//         if part.is_empty() { continue; }
//         let items: Vec<&str> = part.split(',').collect();
//         if items.is_empty() { continue; }
//         let head = items[0];

//         enum Target<'a> { C(&'a [i64], &'a mut [i16]), G(&'a [i64], &'a mut [i16]) }
//         let mut target = if head.starts_with("C+m") || head.starts_with("c+m") {
//             Some(Target::C(&pos_c, &mut flag_c))
//         } else if head.starts_with("C-m") || head.starts_with("c-m") {
//             Some(Target::G(&pos_g, &mut flag_g))
//         } else {
//             None
//         };

//         if let Some(tgt) = target.as_mut() {
//             let mut idx_sum: i64 = 0;
//             let mut n_evt: usize = 0;

//             for tok in &items[1..] {
//                 if tok.is_empty() { continue; }
//                 let t = tok.trim_end_matches('?');
//                 if t.is_empty() { continue; }
//                 let Ok(sh) = t.parse::<i64>() else {
                    
//                     ml_idx = ml_idx.saturating_add(1);
//                     n_evt += 1;
//                     continue;
//                 };
//                 idx_sum += sh;
//                 let pos_idx = (n_evt as i64 + idx_sum) as usize;

//                 match tgt {
//                     Target::C(pos_list, flag_list) => {
//                         if pos_idx < pos_list.len() && ml_idx < ml_v.len() && ml_v[ml_idx] >= prob_cutoff {
//                             flag_list[pos_idx] = 1;
//                         }
//                     }
//                     Target::G(pos_list, flag_list) => {
//                         if pos_idx < pos_list.len() && ml_idx < ml_v.len() && ml_v[ml_idx] >= prob_cutoff {
//                             flag_list[pos_idx] = 1;
//                         }
//                     }
//                 }

//                 n_evt += 1;
//                 ml_idx = ml_idx.saturating_add(1);
//             }
//         } else {
//             let mut cnt = 0usize;
//             for tok in &items[1..] {
//                 if tok.is_empty() { continue; }
//                 let t = tok.trim_end_matches('?');
//                 if t.is_empty() { continue; }
//                 if t.parse::<i64>().is_ok() { cnt += 1; }
//             }
//             ml_idx = ml_idx.saturating_add(cnt);
//         }
//     }
    
//     let mut positions: Vec<i64> = Vec::with_capacity(pos_c.len() + pos_g.len());
//     let mut flags: Vec<i16> = Vec::with_capacity(positions.capacity());
//     let (mut i, mut j) = (0usize, 0usize);
//     while i < pos_c.len() || j < pos_g.len() {
//         if j >= pos_g.len() || (i < pos_c.len() && pos_c[i] <= pos_g[j]) {
//             positions.push(pos_c[i]);
//             flags.push(flag_c[i]);
//             i += 1;
//         } else {
//             positions.push(pos_g[j]);
//             flags.push(flag_g[j]);
//             j += 1;
//         }
//     }
//     println!("Positions: {:?}, Flags: {:?}", positions, flags);
//     (positions, flags)
// }

fn get_5mc_sites_from_read(
    read_seq: &[u8],
    mm: Option<&str>,
    ml: Option<&[u8]>,
    prob_cutoff: u8,
) -> (Vec<i64>, Vec<i16>) {

    let mm_s = match mm { Some(s) => s, None => return (Vec::new(), Vec::new()) };
    let ml_v = ml.unwrap_or(&[]);
  
    let c_positions: Vec<i64> = memchr_iter(b'C', read_seq).map(|i| i as i64).collect();

    let mut flags: Vec<i16> = vec![0; c_positions.len()];
    let mut ml_idx: usize = 0;

    for part in mm_s.trim_end_matches(';').split(';') {
        if part.is_empty() { continue; }
        let mut it = part.split(',');
        let head = it.next().unwrap_or("");
        let head_lc = head.to_ascii_lowercase();

        let is_c_5mc = head_lc.starts_with("c+m") || head_lc.starts_with("c-m");
        if !is_c_5mc {
            for tok in it {
                if tok.is_empty() { continue; }
                let has_q = tok.as_bytes().last() == Some(&b'?');
                let t = tok.trim_end_matches('?');
                if t.is_empty() { continue; }
                if !has_q && t.parse::<i64>().is_ok() {
                    ml_idx = ml_idx.saturating_add(1);
                }
            }
            continue;
        }

        let mut idx_sum: i64 = 0;
        let mut n_evt: usize = 0;

        for tok in it {
            if tok.is_empty() { continue; }
            let has_q = tok.as_bytes().last() == Some(&b'?');
            let t = tok.trim_end_matches('?');
            if t.is_empty() { continue; }
            let Ok(sh) = t.parse::<i64>() else { continue; };

            idx_sum += sh;
            let pos_idx = (n_evt as i64 + idx_sum) as usize;

            if pos_idx < flags.len() {
                if has_q {
         
                } else if ml_idx < ml_v.len() && ml_v[ml_idx] >= prob_cutoff {
                    flags[pos_idx] = 1;
                }
            }

            n_evt += 1;
            if !has_q {
                ml_idx = ml_idx.saturating_add(1);
            }
        }
    }

    (c_positions, flags)
}

fn count_5mc_consistency_split(
    cigar: &bam::record::CigarStringView,
    c_positions: &[i64],
    meth_flags: &[i16],
    is_fwd: bool,
    mut r: i64,
    ref_bytes: &[u8],
    ref_len: i64,
    site_set: Option<&HashSet<i64>>,
    cpg: bool,
) -> (i32, i32, i32) {
    let mut q: i64 = 0;
    let mut i: usize = 0;
    let mut matches = 0i32;
    let mut read_missing = 0i32; 
    let mut ref_missing = 0i32; 

    for op in cigar.iter() {
       
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                let len = *l as i64;
                let end_q = q + len;
    
                
                while i < c_positions.len() && c_positions[i] < q { i += 1; }
                
                while i < c_positions.len() && c_positions[i] < end_q {
                    let cpos = c_positions[i];
                    let offset = cpos - q;
                    let rpos = r + offset;
                    
                    if rpos >= 0 && rpos < ref_len {
                        let rb = ref_bytes[rpos as usize];
                        let is_meth = meth_flags[i] != 0;
                        let hit = if is_fwd {
                            if cpg {
                                let r_after = if rpos + 1 < ref_len { ref_bytes[(rpos + 1) as usize] } else { b'N' };
                                r_after == b'G' && rb == b'C' && site_set.map(|s| s.contains(&rpos)).unwrap_or(false)
                            } else {
                                rb == b'C' && site_set.map(|s| s.contains(&rpos)).unwrap_or(false)
                            }
                        } else {
                          
                            if cpg {
                                let r_prev = if rpos > 0 { ref_bytes[(rpos - 1) as usize] } else { b'N' };
                                r_prev == b'C' && rb == b'G' && site_set.map(|s| s.contains(&(rpos - 1))).unwrap_or(false)
                            } else {
                                rb == b'G' && site_set.map(|s| s.contains(&(rpos - 1))).unwrap_or(false)
                            }  
                            
                            
                        };
                       
                        if hit == is_meth {
                            matches += 1;
                        } else {
                            if hit && !is_meth {
                                read_missing += 1;
                            } else if !hit && is_meth {
                                ref_missing += 1;
                            }
                        }
                        // println!("is_meth: {}, hit: {}, cpos: {}, q: {}, r: {}, offset: {}, rpos: {}", is_meth, hit, cpos, q, r, offset, rpos);
                    }
                    i += 1;
                }

                q = end_q;
                r += len;
            }
            Cigar::Ins(l) => {
                q += *l as i64;
                while i < c_positions.len() && c_positions[i] < q { i += 1; }
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => { r += *l as i64; }
            Cigar::SoftClip(l) => {
                q += *l as i64;
                while i < c_positions.len() && c_positions[i] < q { i += 1; }
            }
            Cigar::HardClip(l) | Cigar::Pad(l) => {
                q += *l as i64;
                while i < c_positions.len() && c_positions[i] < q { i += 1; }
            }
            _ => {}
        }
    }

    (matches, read_missing, ref_missing)
}

fn get_mm_ml_from_rec(
    rec: &bam::Record
) -> (Option<String>, Option<Vec<u8>>) {
    let mm_owned = match rec.aux(b"MM") {
        Ok(Aux::String(s)) => Some(String::from_utf8_lossy(s.as_bytes()).to_string()),
        _ => None,
    };
    let ml_owned: Option<Vec<u8>> = match rec.aux(b"ML") {
        Ok(Aux::ArrayU8(v))  => Some(v.iter().collect()),
        Ok(Aux::ArrayI32(v)) => Some(v.iter().map(|x| x.clamp(0, 255) as u8).collect()),
        Ok(Aux::ArrayI16(v)) => Some(v.iter().map(|x| (x as i32).clamp(0, 255) as u8).collect()),
        Ok(Aux::ArrayU16(v)) => Some(v.iter().map(|x| (x as i32).clamp(0, 255) as u8).collect()),
        _ => None,
    };
    (mm_owned, ml_owned)
}

fn recalc_score_and_update(
    rec: &mut bam::Record,
    tid_names: &[String],
    match_score: i32,
    ref_miss_penalty: i32,
    read_miss_penalty: i32,
    prob_cutoff: u8,
    fa: &HashMap<String, String>,
    met_sites: &HashMap<String, HashSet<i64>>,
    cpg: bool
) -> anyResult<()> {

    // let read_seq_raw = rec.seq().as_bytes();
    // let read_seq = if rec.is_reverse() { revcomp(&read_seq_raw) } else { read_seq_raw.to_vec() };
    
    let read_seq = rec.seq().as_bytes().to_vec();
   
    // let mm_owned = match rec.aux(b"MM") {
    //     Ok(Aux::String(s)) => Some(String::from_utf8_lossy(s.as_bytes()).to_string()),
    //     _ => None,
    // };
    // let ml_owned: Option<Vec<u8>> = match rec.aux(b"ML") {
    //     Ok(Aux::ArrayU8(v))  => Some(v.iter().collect()),
    //     Ok(Aux::ArrayI32(v)) => Some(v.iter().map(|x| x.clamp(0, 255) as u8).collect()),
    //     Ok(Aux::ArrayI16(v)) => Some(v.iter().map(|x| (x as i32).clamp(0, 255) as u8).collect()),
    //     Ok(Aux::ArrayU16(v)) => Some(v.iter().map(|x| (x as i32).clamp(0, 255) as u8).collect()),
    //     _ => None,
    // };
    let (mm_owned, ml_owned) = get_mm_ml_from_rec(rec);
    let (pos_read, flags_read) = get_5mc_sites_from_read(
        &read_seq,
        mm_owned.as_deref(),
        ml_owned.as_deref(),
        prob_cutoff
    );


    let mut pairs: Vec<(i64, i16)> = pos_read.into_iter().zip(flags_read.into_iter()).collect();
    if rec.is_reverse() {
        let l = read_seq.len() as i64;
        for (p, _) in pairs.iter_mut() {
           
            *p = l - 1 - *p;
        }
    }
    pairs.sort_unstable_by_key(|x| x.0);
    let (c_positions_q, meth_flags_q): (Vec<i64>, Vec<i16>) = pairs.into_iter().unzip();

    let rname = match tid_names.get(rec.tid() as usize) {
        Some(s) => s.as_str(),
        None => return Ok(()),
    };
    let Some(ref_seq) = fa.get(rname) else { return Ok(()); };
    let ref_bytes = ref_seq.as_bytes();
    let ref_len = ref_bytes.len() as i64;
    let is_fwd = !rec.is_reverse();
    let site_set = met_sites.get(rname);

    let (mc, mr, mf) = count_5mc_consistency_split(
        &rec.cigar(),
        &c_positions_q,
        &meth_flags_q,
        is_fwd,
        rec.pos(),
        ref_bytes,
        ref_len,
        site_set,
        cpg
    );
    
    let mismatches = mr + mf;
    
    let orig_as = match rec.aux(b"AS") {
        Ok(Aux::I8(v)) => v as i32,
        Ok(Aux::I16(v)) => v as i32,
        Ok(Aux::I32(v)) => v,
        Ok(Aux::U8(v)) => v as i32,
        Ok(Aux::U16(v)) => v as i32,
        Ok(Aux::U32(v)) => v as i32,
        _ => 0,
    };
    // let new_as = orig_as - penalty * mismatches;
    let new_as = orig_as + match_score * mc - read_miss_penalty * mr - ref_miss_penalty * mf;
    rec.push_aux(b"MA", Aux::I32(mismatches as i32)).ok();
    rec.remove_aux(b"AS").ok();
    rec.push_aux(b"AS", Aux::I32(new_as)).ok();

    Ok(())
}

fn recalc_score_and_update2(
    recs: &mut Vec<bam::Record>,
    read_seq: &[u8],
    c_positions_base: &[i64],
    meth_flags_base: &[i16],
    tid_names: &[String],
    match_score: i32,
    ref_miss_penalty: i32,
    read_miss_penalty: i32,
    prob_cutoff: u8,
    fa: &HashMap<String, String>,
    met_sites: &HashMap<String, HashSet<i64>>,
    cpg: bool,
) -> anyResult<()> {
    for record in recs.iter_mut() {
        // println!("c_positions_base: {:?}", c_positions_base);
        let mut pairs: Vec<(i64, i16)> = c_positions_base
                    .iter()
                    .copied()
                    .zip(meth_flags_base.iter().copied())
                    .collect();
        if record.is_reverse() {
            let l = read_seq.len() as i64;
            for (p, _) in pairs.iter_mut() {
                *p = l - 1 - *p;
            }
        }
        pairs.sort_unstable_by_key(|x| x.0);
        let (c_positions_q, meth_flags_q): (Vec<i64>, Vec<i16>) = pairs.into_iter().unzip();
        
        let is_fwd: bool = !record.is_reverse();
        let rname = match tid_names.get(record.tid() as usize) {
            Some(s) => s.as_str(),
            None => return Ok(()),
        };
        let Some(ref_seq) = fa.get(rname) else { return Ok(()); };
        let ref_bytes = ref_seq.as_bytes();
        let ref_len = ref_bytes.len() as i64;
        let site_set = met_sites.get(rname);
        // println!("rname : {}", rname);
        // println!("{:?}", site_set);
        // println!("site_set_len: {}", site_set.map(|s| s.len()).unwrap_or(0));
        // println!("c_positions_q: {:?}", c_positions_q);
        // println!("read: {}, is_fwd: {}", std::str::from_utf8(record.qname()).unwrap_or("N/A"), is_fwd);
        // println!("pos_len: {}, meth_len: {}", c_positions_q.len(), meth_flags_q.len());
        // let filtered_flags = meth_flags_q.iter().filter(|&&f| f != 0).count();
        // println!("filtered flags (meth): {}", filtered_flags);
        let (mc, mr, mf) = count_5mc_consistency_split(
            &record.cigar(),
            &c_positions_q,
            &meth_flags_q,
            !record.is_reverse(),
            record.pos(),
            ref_bytes,
            ref_len,
            site_set,
            cpg
        );
    
       
        
        // println!("mc: {}, mr: {}, mf: {}", mc, mr, mf);
        let mismatches = mr + mf;
        
        let orig_as = match record.aux(b"AS") {
            Ok(Aux::I8(v)) => v as i32,
            Ok(Aux::I16(v)) => v as i32,
            Ok(Aux::I32(v)) => v,
            Ok(Aux::U8(v)) => v as i32,
            Ok(Aux::U16(v)) => v as i32,
            Ok(Aux::U32(v)) => v as i32,
            _ => 0,
        };
        // let new_as = orig_as - penalty * mismatches;
        let new_as = orig_as + match_score * mc - read_miss_penalty * mr - ref_miss_penalty * mf;
        record.push_aux(b"MA", Aux::I32(mismatches as i32)).ok();
        record.remove_aux(b"AS").ok();
        record.push_aux(b"AS", Aux::I32(new_as)).ok();
    
    }
    
    Ok(())
}

// pub fn parse_bam(input_bam: &String, fa: &Option<String>, bg: &Option<String>, 
//                     output_bam: &String, threads: usize) {
//     let mut bam = if input_bam == &String::from("-") {
//         Reader::from_stdin().expect("Failed to read from stdin")
//     } else {
//         Reader::from_path(input_bam).expect("Failed to read from the provided path")
//     };

//     let fa_hash: HashMap<String, String> = if let Some(fa) = fa {
//         let fasta = Fastx::new(fa);
//         fasta.get_chrom_seqs().unwrap()
//     } else {
//         HashMap::new()
//     };
//     let bg_map: HashMap<String, HashSet<i64>> = if let Some(bg) = bg {
//         parse_bedgraph(bg, 0.5).expect("Failed to parse bedgraph")
//     } else {
//         HashMap::new()
//     };


//     let header = Header::from_template(bam.header());
//     let headerview = HeaderView::from_header(&header);
//     let tid_names: Vec<String> = (0..headerview.target_count())
//         .map(|i| String::from_utf8_lossy(headerview.tid2name(i as u32)).to_string())
//         .collect();
//     let _ = bam.set_threads(threads);

//     let (sender, receiver) = bounded::<Vec<(usize, Vec<Record>)>>(1000);
//     let mut handles = vec![];
//     // let mut writer = Writer::from_path(output_bam, &header, bam::Format::Bam).expect("Failed to create BAM writer");
//     // let _ = writer.set_threads(threads);
//     // let writer = Arc::new(Mutex::new(writer));
//     let (wtx, wrx) = crossbeam_channel::bounded::<Vec<Record>>(10_000);

//     let out_path = output_bam.clone();
//     let header_for_writer = header; 
//     let writer_threads = threads as i32;
//     let writer_handle = std::thread::spawn(move || {
//          let mut writer = Writer::from_path(&out_path, &header_for_writer, bam::Format::Bam)
//              .expect("Failed to create BAM writer");
        
//         let _ = writer.set_threads(writer_threads.max(1).try_into().unwrap());
//         while let Ok(group) = wrx.recv() {
//             for rec in group {
//                 writer.write(&rec).expect("Failed to write record");
//             }
//         }
//     });

//     for _ in 0..8 {
//         let receiver: Receiver<_> = receiver.clone();
//         let tid_names = tid_names.clone();
//         let fa_hash = fa_hash.clone();
//         let bg_map = bg_map.clone();
//         let fa_len: usize = fa_hash.len();
//         let bg_len: usize = bg_map.len();
//         // let writer = Arc::clone(&writer);
//         let wtx = wtx.clone();
//         handles.push(thread::spawn(move || {
//             while let Ok(batch) = receiver.recv() {
//                 // let mut output_batch: Vec<_> = vec![];
//                 let mut flat_batch: Vec<Record> = Vec::new();
//                 for (idx, records) in batch {
//                     // Process each record
//                     // For example, print the read name and index
//                     let split_recs = split_records(records);
                    
//                     for mut group in split_recs.into_iter() {
//                         match group.len() {
//                             0 => continue, 
//                             1 => {
//                                 // output_batch.push(group);
//                                 flat_batch.extend(group.into_iter());
//                             },
//                             _ => {
                                
//                                 let primary_mapq = group[0].mapq();
//                                 if primary_mapq >= 2 {
//                                     // output_batch.push(group);
//                                     flat_batch.extend(group.into_iter());
//                                 } else {
//                                     let as_vec = group.iter().filter_map(|r| get_as(r)).collect::<Vec<i32>>();
//                                     let primary_as = as_vec[0];
//                                     let max_as = *as_vec.iter().max().unwrap();
//                                     let secondary_as = as_vec[1];
//                                     if primary_as > secondary_as {
//                                         group[0].set_mapq(2);
                                        
//                                         // output_batch.push(group);
//                                         flat_batch.extend(group.into_iter());
//                                     } else {
//                                         if (fa_len > 0) && (bg_len > 0) {
//                                             if !is_contain_methylation(&group[0]) {
//                                                 // output_batch.push(group);
//                                                 flat_batch.extend(group.into_iter());
                                                
//                                             } else {
//                                                 for record in group.iter_mut() {
//                                                     recalc_score_and_update( record, &tid_names, 2, 125, &fa_hash, &bg_map, 2).ok();
//                                                 }
//                                                 let as_vec = group.iter().filter_map(|r| get_as(r)).collect::<Vec<i32>>();
//                                                 let max_as = *as_vec.iter().max().unwrap();
//                                                 let best_recs_idx = as_vec.iter().position(|&x| x == max_as).unwrap_or(0);
//                                                 let secondary_as = as_vec.iter().enumerate().filter(|(i, _)| *i != best_recs_idx).map(|(_, &v)| v).max().unwrap_or(i32::MIN);

//                                                 if max_as > secondary_as {
//                                                     group[best_recs_idx].set_mapq(2);
//                                                     if best_recs_idx != 0 {
//                                                         group[best_recs_idx].push_aux(b"RF", Aux::String("Y")).ok();
//                                                         let _primary_flag = group[0].flags();
//                                                         let _primary_flag0 = if _primary_flag & 0x0 == 0 {
//                                                             0x0
//                                                         } else {
//                                                             0x800
//                                                         };
//                                                         let _flag = group[best_recs_idx].flags();
                                            
//                                                         // clear secondary flag and set primary flag
//                                                         group[best_recs_idx].set_flags( (_flag & !0x100) | _primary_flag0);
//                                                         group[0].set_flags(_primary_flag & !_primary_flag0 | 0x100);
//                                                         group[0].set_mapq(0);

//                                                         // set tp:A:P to primary, tp:A:S to secondary
//                                                         group[best_recs_idx].push_aux(b"tp", Aux::String("P")).ok();
//                                                         group[0].push_aux(b"tp", Aux::String("S")).ok();

//                                                         // swap 0 and best_recs_idx
//                                                         group.swap(0, best_recs_idx);
//                                                     }
                                                    
//                                                 }
                                                
//                                                 // output_batch.push(group);
//                                                 flat_batch.extend(group.into_iter());
//                                             }
//                                         } else {
//                                             // output_batch.push(group);
//                                             flat_batch.extend(group.into_iter());
//                                         }
                                        
//                                     }
//                                 }
                                
//                             }
//                         }
                        
//                     }
               
//                 }
                
//                 // let mut writer = writer.lock().unwrap();
//                 // for out_recs in &output_batch {
//                 //     for rec in out_recs {
//                 //         writer.write(&rec).expect("Failed to write record");
//                 //     }
//                 // }
//                 if !flat_batch.is_empty() {
//                     wtx.send(flat_batch).unwrap();
//                 }
//             }
//         }));
//     }
//     drop(wtx);
//     let mut previous_read = Vec::new();
//     let mut idx: usize = 0;
//     let mut batch = Vec::with_capacity(1000);
//     let mut records_to_process = vec![];
//     for record in bam.records() {
//         let record = record.expect("Failed to read record");
//         if record.is_unmapped() {
//             continue;
//         }
//         let read = record.qname();
        
//         if previous_read.is_empty() {
//             // Process previous_read's records
//             previous_read = read.to_vec();
//             records_to_process.push(record);
//             continue
//         }
//             // New read, process the previous one
        
//         if previous_read != read {
//             batch.push((idx, std::mem::take(&mut records_to_process)));

//             if batch.len() >= 1000 {
//                 sender.send(std::mem::take(&mut batch)).unwrap();
//             }
          
//             idx += 1;
//             previous_read = read.to_vec();
//         }

//         records_to_process.push(record);
        
//     }

//     if !records_to_process.is_empty() {
//         batch.push((idx, std::mem::take(&mut records_to_process)));
//     }

//     if !batch.is_empty() {
//         sender.send(batch).unwrap();
//     }

//     drop(sender);
//     for handle in handles {
//         handle.join().expect("Thread panicked");
//     }

//     writer_handle.join().expect("Writer thread panicked");

// }

fn process_batch_with_rayon(
    batch: Vec<(usize, Vec<Record>)>,
    tid_names: &[String],
    fa_hash: &HashMap<String, String>,
    bg_map: &HashMap<String, HashSet<i64>>,
    match_score: i32,
    ref_miss_penalty: i32,
    read_miss_penalty: i32,
    prob_cutoff: u8,
    designate_mapq: u8,
    is_set_y: bool,
    cpg: bool,
) -> Vec<Record> {
    let fa_len = fa_hash.len();
    let bg_len = bg_map.len();

    // let groups: Vec<Vec<Record>> = batch
    //     .into_iter()
    //     .flat_map(|(_, recs)| split_records(recs))
    //     .collect();

    batch
        .into_par_iter()
        .map(|(_, mut records)| {
            
            if records.is_empty() {
                return Vec::<Record>::new();
            }
            let (mm_owned, ml_owned, read_seq_raw) =
            if fa_len > 0 && bg_len > 0 && !is_set_y && is_contain_methylation(&records[0]) {
                let (mm, ml) = get_mm_ml_from_rec(&records[0]);
                let seq = records[0].seq().as_bytes().to_vec();
                if records[0].is_reverse() {
                    (mm, ml, revcomp(&seq))
                } else {
                    (mm, ml, seq)
                }
            } else {
                (None, None, Vec::new())
            };
            let primary_record_is_reverse = records[0].is_reverse();
            let mut groups = split_records(records);

            if fa_len > 0 && bg_len > 0  {
          
                for group in groups.iter_mut() {

                    if group.len() == 1 {
                        continue;
                    }

                    let primary_mapq = group[0].mapq();
                    if primary_mapq >= designate_mapq {
                        continue;
                    }
                    
                    // let as_vec: Vec<i32> = group.iter().filter_map(|r| get_as(r)).collect();
                    // if as_vec.len() < 2 {
                    //     continue;
                    // }

                    // let primary_as = as_vec[0];
                    // let (best_idx, max_as) = as_vec
                    //     .iter()
                    //     .copied()
                    //     .enumerate()
                    //     .max_by_key(|&(_, v)| v)
                    //     .unwrap_or((0, i32::MIN));
                    // let secondary_as = as_vec
                    //     .iter()
                    //     .enumerate()
                    //     .filter(|(i, _)| *i != best_idx)
                    //     .map(|(_, &v)| v)
                    //     .max()
                    //     .unwrap_or(i32::MIN);

                    // if (primary_as > secondary_as) && (best_idx == 0) {
                    //     group[0].set_mapq(designate_mapq);
                    //     group[0].push_aux(b"RF", Aux::String("Y")).ok();
                    //     continue;
                    // } 

                    if is_set_y {
                        for rec in group.iter_mut() {
                            let _ = recalc_score_and_update(rec, 
                                                            tid_names, 
                                                            match_score,
                                                            ref_miss_penalty,
                                                            read_miss_penalty,
                                                            prob_cutoff, 
                                                            fa_hash, 
                                                            bg_map, 
                                                            cpg);
                        }
                    } else {
                        let (c_positions_base, meth_flags_base) = get_5mc_sites_from_read(
                            &read_seq_raw,
                            mm_owned.as_deref(), 
                            ml_owned.as_deref(), 
                            prob_cutoff
                        );
                        if c_positions_base.is_empty() {
                            continue;
                        }
                        let _ = recalc_score_and_update2(group, 
                                                        &read_seq_raw,
                                                        &c_positions_base,
                                                        &meth_flags_base,
                                                        tid_names, 
                                                        match_score,
                                                        ref_miss_penalty,
                                                        read_miss_penalty,
                                                        prob_cutoff, 
                                                        fa_hash, 
                                                        bg_map,
                                                        cpg);
                    }
            
                    let as_vec2: Vec<i32> = group.iter().filter_map(|r| get_as(r)).collect();
                    let (best_recs_idx, max_as2) = as_vec2
                        .iter()
                        .copied()
                        .enumerate()
                        .max_by_key(|&(_, v)| v)
                        .unwrap_or((0, i32::MIN));
                    let secondary_as2 = as_vec2
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != best_recs_idx)
                        .map(|(_, &v)| v)
                        .max()
                        .unwrap_or(i32::MIN);

                    if max_as2 > secondary_as2 {
                        group[best_recs_idx].set_mapq(designate_mapq);
                        if best_recs_idx != 0 {
                            group[best_recs_idx].push_aux(b"RF", Aux::String("Y")).ok();
                            let _primary_flag = group[0].flags();
                            let _primary_flag0 = if _primary_flag & 0x0 == 0 { 0x0 } else { 0x800 };
                            let _flag = group[best_recs_idx].flags();
                            
                            group[best_recs_idx].set_flags((_flag & !0x100) | _primary_flag0);
                            group[0].set_flags(_primary_flag & !_primary_flag0 | 0x100);
                            group[0].set_mapq(0);
                            group[best_recs_idx].push_aux(b"tp", Aux::String("P")).ok();
                            group[0].push_aux(b"tp", Aux::String("S")).ok();
                            group.swap(0, best_recs_idx);
                        }
                    }
                }
            } else {
                for group in groups.iter_mut() {
                    if group.len() == 1 {
                        continue;
                    }

                    let primary_mapq = group[0].mapq();
                    if primary_mapq >= designate_mapq {
                        continue;
                    }

                    let as_vec: Vec<i32> = group.iter().filter_map(|r| get_as(r)).collect();
                    if as_vec.len() < 2 {
                        continue;
                    }

                    let primary_as = as_vec[0];
                    let (best_idx, max_as) = as_vec
                        .iter()
                        .copied()
                        .enumerate()
                        .max_by_key(|&(_, v)| v)
                        .unwrap_or((0, i32::MIN));
                    let secondary_as = as_vec
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != best_idx)
                        .map(|(_, &v)| v)
                        .max()
                        .unwrap_or(i32::MIN);

                    if (primary_as > secondary_as) && (best_idx == 0) {
                        group[0].set_mapq(designate_mapq);
                        group[0].push_aux(b"RF", Aux::String("Y")).ok();
                    } 
                }
            }
            
            groups.into_iter().flatten().collect::<Vec<Record>>()
        })
        .flatten()
        .collect()
}


pub fn parse_bam(input_bam: &String, 
                    fa: &Option<String>, 
                    bg: &Option<String>, 
                    match_score: i32,
                    ref_miss_penalty: i32,
                    read_miss_penalty: i32,
                    ref_prob_cutoff: f64,
                    prob_cutoff: u8,
                    designate_mapq: u8,
                    is_set_y: bool,
                    cpg: bool,
                    output_secondary: bool,
                    output_bam: &String, 
                    threads: usize) {
    let chunksize: usize = 50_000;
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };

    let fa_hash: HashMap<String, String> = if let Some(fa) = fa {
        let fasta = Fastx::new(fa);
        fasta.get_chrom_seqs().unwrap()
    } else {
        HashMap::new()
    };
    let bg_map: HashMap<String, HashSet<i64>> = if let Some(bg) = bg {
        parse_bedgraph(bg, ref_prob_cutoff).expect("Failed to parse bedgraph")
    } else {
        HashMap::new()
    };

    let header = Header::from_template(bam.header());
    let headerview = HeaderView::from_header(&header);
    let tid_names: Vec<String> = (0..headerview.target_count())
        .map(|i| String::from_utf8_lossy(headerview.tid2name(i as u32)).to_string())
        .collect();
    let _ = bam.set_threads(threads);

    let (wtx, wrx) = crossbeam_channel::bounded::<Vec<Record>>(10_000);

    let out_path = output_bam.clone();
    let header_for_writer = header;
    let writer_threads = threads as i32;
    let writer_handle = std::thread::spawn(move || {
        let mut writer = Writer::from_path(&out_path, &header_for_writer, bam::Format::Bam)
            .expect("Failed to create BAM writer");
        let _ = writer.set_threads(writer_threads.max(1).try_into().unwrap());
        while let Ok(group) = wrx.recv() {
            for rec in group {
                if !output_secondary && rec.is_secondary() {
                    continue;
                }
                writer.write(&rec).expect("Failed to write record");
            }
        }
    });

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} Loaded {pos} records | speed {per_sec} | elapsed {elapsed}")
            .unwrap()
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100));
    pb.set_draw_target(ProgressDrawTarget::stderr_with_hz(10)); 
    let mut seen: u64 = 0;

    let mut previous_read: Vec<u8> = Vec::new();
    let mut idx: usize = 0;
    let mut batch: Vec<(usize, Vec<Record>)> = Vec::with_capacity(chunksize );
    let mut records_to_process: Vec<Record> = Vec::new();

    for rec in bam.records() {
        seen += 1;
        if seen % 10_000 == 0 { pb.set_position(seen); } 
        let record = rec.expect("Failed to read record");
        if record.is_unmapped() {
            continue;
        }
        let qname = record.qname();
        if previous_read.is_empty() {
            previous_read.extend_from_slice(qname);
            records_to_process.push(record);
            continue;
        }
        if qname != &previous_read[..] {
            batch.push((idx, std::mem::take(&mut records_to_process)));
            if batch.len() >= chunksize {
                let flat = process_batch_with_rayon(std::mem::take(&mut batch), &tid_names, &fa_hash, &bg_map,
                    match_score, ref_miss_penalty, read_miss_penalty,
                    prob_cutoff, designate_mapq, is_set_y, cpg);  
                if !flat.is_empty() {
                    wtx.send(flat).unwrap();
                }
            }
            idx += 1;
            previous_read.clear();
            previous_read.extend_from_slice(qname);
           
        }
        records_to_process.push(record);
    }

    if !records_to_process.is_empty() {
        batch.push((idx, std::mem::take(&mut records_to_process)));
    }
    if !batch.is_empty() {
        let flat = process_batch_with_rayon(std::mem::take(&mut batch), &tid_names, &fa_hash, &bg_map,
            match_score, ref_miss_penalty, read_miss_penalty, prob_cutoff, designate_mapq, is_set_y, cpg);
        if !flat.is_empty() {
            wtx.send(flat).unwrap();
        }
    }

    drop(wtx);
    writer_handle.join().expect("Writer thread panicked");
    pb.set_position(seen);
    pb.finish_with_message("Finished processing BAM file.");
    
    log::info!("Successfully refined alignments and wrote to {}", output_bam);
    
}