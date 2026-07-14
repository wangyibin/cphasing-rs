
#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
use std::collections::HashMap;
use bio::io::fastq::{Reader as FastqReader, Record as FastqRecord, Writer as FastqWriter};
use bio::io::fasta::{Reader as FastaReader, Record as FastaRecord, Writer as FastaWriter};
use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer, ext::BamRecordExtensions};
use rayon::prelude::*;
use std::path::{ Path, PathBuf };
use std::thread;
use std::sync::{ Arc, Mutex};

use crate::core::{ 
        ChromSizeRecord,
        common_writer};
use crate::pairs::{ Pairs, PairHeader };


// split bam by record number and write to different files 
pub fn split_bam(input_bam: &String, output_prefix: &String, 
             record_num: usize) -> Result<(), Box<dyn std::error::Error>> {
    
    let mut output_prefix = output_prefix.clone();
    if output_prefix.ends_with("/") {
        output_prefix = format!("{}/{}", output_prefix, "split");
    }
    let parent = Path::new(&output_prefix).parent().unwrap().to_path_buf();
    match parent.exists() {
        true => {},
        false => {
            std::fs::create_dir_all(&parent).unwrap();
        }
    }

    
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(14);
    let header = Header::from_template(bam.header());
    let mut i = 0;
    let mut j = 0;
    // let mut wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam).unwrap();
    // let _ = wtr.set_threads(8);
    // log::info!("write {} records to {}", record_num, format!("{}_{}.bam", output_prefix, j));

    let mut read_name: Vec<u8> = Vec::new();
    let mut previous_read_name: Vec<u8> = Vec::new();
    // let mut line_num = 0;
    // for r in bam.records() {
    //     let record = r?;
    //     read_name.clear();
    //     read_name.extend_from_slice(record.qname());
        
    //     i += 1;
        
    //     if i == record_num {
    //         if read_name != previous_read_name {

    //             j += 1;
    //             wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam)?;
    //             let _ = wtr.set_threads(8);
    //             log::info!("write {} records to {}", line_num, format!("{}_{}.bam", output_prefix, j));
                
    //             i = 0;
    //             line_num = 0;
    //         } else {
    //             i -= 1;
    //         }
    //     }
    //     line_num += 1;
    //     wtr.write(&record)?;
    //     std::mem::swap(&mut read_name, &mut previous_read_name);
    // } 

    let (sender, receiver) = bounded::<(usize, Vec<Record>)>(10);

    let mut handles = Vec::new();
    for _ in 0..8 {
        let receiver = receiver.clone();
        let output_prefix = output_prefix.clone();
        let header = header.clone();
        let handle = thread::spawn(move || {
            while let Ok((chunk_id, records)) = receiver.recv() {
                let mut wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, chunk_id), &header, bam::Format::Bam).unwrap();
                let _ = wtr.set_threads(8);
                let length = records.len();
                for record in records {
                    wtr.write(&record).unwrap();
                }
                log::info!("write {} records to {}", length, format!("{}_{}.bam", output_prefix, chunk_id));
            }
        });
        handles.push(handle);
    }

    let mut batch = Vec::with_capacity(record_num + 100);
    let mut chunk_id: usize = 0;
    for r in bam.records() {
        let record = r?;
        read_name.clear();
        read_name.extend_from_slice(record.qname());
        
        i += 1;
        
        if i == record_num {
            if read_name != previous_read_name {

                j += 1;
                
                sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                chunk_id += 1;

                i = 0;
            } else {
                i -= 1;
            }
        }

        batch.push(record);
        std::mem::swap(&mut read_name, &mut previous_read_name);
    }
    if batch.len() > 0 {
        sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
    }
    
    drop(sender);
    
    for handle in handles {
        handle.join().unwrap();
    }

    Ok(())
}

pub fn slide2raw(input_bam: &String, output: &String, threads: usize) {
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(threads);

    let header = Header::from_template(bam.header());
    

    let mut wtr = Writer::from_path(output, &header, bam::Format::Bam).unwrap();
    let _ = wtr.set_threads(threads);
    while let Some(r) = bam.records().next() {
        let record = r.unwrap();
        let mut new_record = record.clone();
        let read_id = std::str::from_utf8(record.qname()).unwrap();
        let (read_id, suffix) = read_id.rsplit_once("_").unwrap();
        let mut flag = record.flags();
        
        if suffix != "0" {
            flag  += 2048;
        }

        new_record.set_qname(read_id.as_bytes());
        new_record.set_flags(flag);

        wtr.write(&new_record).unwrap();
        

    }

}

fn parse_cigar_for_paf(record: &Record) -> (i64, i64, i64, i64, u32, u32, u32) {
    let mut match_len: i64 = 0;
    let mut aln_len: i64 = 0;
    let mut ins_len: i64 = 0;
    let mut del_len: i64 = 0;
    
    let mut qstart: u32 = 0;
    let mut qpos: u32 = 0;
    let mut qlen: u32 = 0;
    let mut qend_consumed: u32 = 0;
    let mut found_start = false;

    for op in record.cigar().iter() {
        let l = op.len();
        match op {
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                let len = l as i64;
                match_len += len;
                aln_len += len;
                if !found_start { qstart = qpos; found_start = true; }
                qpos += l;
                qlen += l;
                qend_consumed += l;
            }
            Cigar::Del(_) => {
                let len = l as i64;
                del_len += len;
                aln_len += len;
            }
            Cigar::Ins(_) => {
                let len = l as i64;
                ins_len += len;
                aln_len += len;
                if !found_start { qstart = qpos; found_start = true; }
                qpos += l;
                qlen += l;
                qend_consumed += l;
            }
            Cigar::RefSkip(_) | Cigar::Pad(_) => {
                aln_len += l as i64;
            }
            Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                if !found_start { qstart = qpos + l; found_start = true; }
                qpos += l;
                qlen += l;
            }
        }
    }
    let qend = qstart + qend_consumed;
    (match_len, aln_len, ins_len, del_len, qstart, qend, qlen)
}


fn get_query_start_end(record: &Record) -> (u32, u32, u32) {
    let mut qlen: u32 = 0;

    for op in record.cigar().iter() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | 
            Cigar::Ins(l) | Cigar::SoftClip(l) | Cigar::HardClip(l) => {
                qlen += *l as u32;
            }
            _ => {}
        }
    }

    let mut bam_start: u32 = 0;
    let mut bam_end: u32 = 0;
    let mut current_pos_on_read: u32 = 0;
    let mut is_first_match = true;

    for op in record.cigar().iter() {
        let l = op.len() as u32;
        match op {
            Cigar::HardClip(_) => {
                current_pos_on_read += l;
            }
            Cigar::SoftClip(_) => {
                current_pos_on_read += l;
            }
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::Ins(_) => {
                if is_first_match {
                    bam_start = current_pos_on_read;
                    is_first_match = false;
                }
                current_pos_on_read += l;
                bam_end = current_pos_on_read; 
            }
            _ => {}
        }
    }

    if record.is_reverse() {
        (qlen - bam_end, qlen - bam_start, qlen)
    } else {
        (bam_start, bam_end, qlen)
    }
}

pub fn bam2paf(input_bam: &String, output: &String, threads: usize, is_secondary: bool) {
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };

    let header = Header::from_template(bam.header());
    let hv = HeaderView::from_header(&header);

    let tnames: Vec<&str> = (0..hv.target_count())
        .map(|tid| std::str::from_utf8(hv.tid2name(tid)).unwrap())
        .collect();
    let tlens: Vec<usize> = (0..hv.target_count())
        .map(|tid| hv.target_len(tid).unwrap() as usize)
        .collect();

    let _ = bam.set_threads(threads);

    let mut writer = common_writer(output);

    for r in bam.records() {
        let record = r.unwrap();
        if record.is_unmapped() {
            continue;
        }
        if record.is_secondary() && !is_secondary {
            continue;
        }

        let flag = if record.is_secondary() { "tp:A:S" } else { "tp:A:P" };
        let strand = if record.is_reverse() { "-" } else { "+" };

        let mut match_len: i64 = 0;
        let mut del_len: i64 = 0;
        let mut ins_len: i64 = 0;
        let mut aln_len: i64 = 0;
        for op in record.cigar().iter() {
            match op {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                    let len = *l as i64;
                    match_len += len;
                    aln_len += len;
                }
                Cigar::Del(l) => {
                    let len = *l as i64;
                    del_len += len;
                    aln_len += len;
                }
                Cigar::Ins(l) => {
                    let len = *l as i64;
                    ins_len += len;
                    aln_len += len;
                }
                Cigar::RefSkip(l) | Cigar::Pad(l) => {
                    aln_len += *l as i64;
                }
                _ => {}
            }
        }

        let nm: i64 = match record.aux(b"NM") {
            Ok(Aux::U8(v)) => v as i64,
            Ok(Aux::U16(v)) => v as i64,
            Ok(Aux::U32(v)) => v as i64,
            Ok(Aux::I32(v)) => v as i64,
            _ => 0,
        };
        let ascore: i64 = match record.aux(b"AS") {
            Ok(Aux::U8(v)) => v as i64,
            Ok(Aux::U16(v)) => v as i64,
            Ok(Aux::U32(v)) => v as i64,
            Ok(Aux::I32(v)) => v as i64,
            _ => 0,
        };

        match_len = match_len - (nm - (ins_len + del_len));

        let (qstart, qend, qlen) = get_query_start_end(&record);

        let qname = std::str::from_utf8(record.qname()).unwrap();
        let tid = record.tid() as usize;
        let tname = tnames[tid];
        let tlen = tlens[tid];
        let tstart = record.reference_start();
        let tend = record.reference_end();
        let mapq = record.mapq();

        use std::io::Write;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\tAS:i:{}\t{}",
            qname, qlen, qstart, qend, strand,
            tname, tlen, tstart, tend,
            match_len, aln_len, mapq, nm, ascore, flag
        ).unwrap();
    }

    log::info!("Successfully output paf to {}", output);
}

// pub fn bam2paf(input_bam: &String, output: &String, threads: usize, is_secondary: bool) {
//     let mut bam = if input_bam == &String::from("-") {
//         Reader::from_stdin().expect("Failed to read from stdin")
//     } else {
//         Reader::from_path(input_bam).expect("Failed to read from the provided path")
//     };
    
//     let header = Header::from_template(bam.header());
//     let header = HeaderView::from_header(&header);
//     let mut chromsizes = Vec::new();
    

//     for tid in 0..header.target_count() {
//         let name = header.tid2name(tid);
//         let len = header.target_len(tid).unwrap();
//         let csr: ChromSizeRecord = ChromSizeRecord {
//             chrom: std::str::from_utf8(name).unwrap().to_string(), 
//             size: len
//         };
//         chromsizes.push(csr);
//     }

//     let chromsizes_db: HashMap<String, u64> = chromsizes.iter().map(|x| (x.chrom.clone(), x.size)).collect();

//     let tid2name: HashMap<usize, String> = (0..header.target_count()).map(|tid| {
//         let name = header.tid2name(tid);
//         (tid as usize, std::str::from_utf8(name).unwrap().to_string())
//     }).collect();

//     let _ = bam.set_threads(threads);
//     let mut idx = 0;

//     let mut writer = common_writer(output);
    
//     while let Some(r) = bam.records().next() {
//         let record = r.unwrap();

//         if record.is_unmapped() { 
//             continue;
//         }

//         if record.is_secondary() & !is_secondary {
//             continue
//         }
//         let flag = if record.is_secondary() {
//             "tp:A:S"
//         } else {
//             "tp:A:P"
//         };

//         let chrom = std::str::from_utf8(header.tid2name(record.tid().try_into().unwrap())).unwrap().to_string();
//         let strand = if record.is_reverse() { "-" } else { "+" };
        
//         let (mut match_length, mut deletion_length, mut insertion_length, mut mm, mut alignment_length) = (0, 0, 0, 0, 0);

//         for cigar in record.cigar().iter() {
//             let len = cigar.len() as i64;
//             match cigar.char() {
//                 'M' | '=' | 'X' => {
//                     match_length += len;
//                     alignment_length += len;
//                     if cigar.char() == 'X' {
//                         mm += len;
//                     }
//                 },
//                 'D' => {
//                     deletion_length += len;
//                     alignment_length += len;
//                 },
//                 'I' => {
//                     insertion_length += len;
//                     alignment_length += len;
//                 },
//                 'N' | 'P' => {
//                     alignment_length += len;
//                 },
//                 _ => {}
//             }
//         }


//         let nm: i64 = match record.aux(b"NM") {
//             Ok(value) => {
//                 match value {
//                     Aux::U8(v) => v.try_into().unwrap(),
//                     Aux::U16(v) => v.try_into().unwrap(),
//                     Aux::U32(v) => v.try_into().unwrap(),
//                     Aux::I32(v) => v.try_into().unwrap(),
//                     _ => 0, 
//                 }
                
//             },
//             Err(e) => 0
//         };

//         let _as: i64 =  match record.aux(b"AS") {
//             Ok(value) => {
//                 match value {
//                     Aux::U8(v) => v.try_into().unwrap(),
//                     Aux::U16(v) => v.try_into().unwrap(),
//                     Aux::U32(v) => v.try_into().unwrap(),
//                     Aux::I32(v) => v.try_into().unwrap(),
//                     _ => 0, 
//                 }
                
//             },
//             Err(e) => 0
//         };

//         match_length = match_length - (nm - (insertion_length + deletion_length));

//         let as_string = format!("AS:i:{}", _as);
//         let nm_string = format!("NM:i:{}", nm);

//         let (qstart, qend, qlen) = get_query_start_end(&record);
    
//         let qname = std::str::from_utf8(record.qname()).unwrap();
//         let tname = std::str::from_utf8(header.tid2name(record.tid().try_into().unwrap())).unwrap();
//         let tlen =  header.target_len(record.tid().try_into().unwrap()).unwrap() as usize; 
//         let tstart = record.reference_start();
//         let tend = record.reference_end();
//         let mapq = record.mapq();
//         writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
//                 qname, qlen, qstart, qend, strand, 
//                 tname, tlen, tstart, tend, match_length, alignment_length, mapq,
//                 nm_string, as_string, flag).unwrap();

//         idx += 1;

//     }

//     log::info!("Successfully output paf to {}", output);

// }

pub fn bam2pairs_pairend(input_bam: &String, min_mapq: u8, output: &String, threads: usize) {
    
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
   
    let header = Header::from_template(bam.header());
    let header = HeaderView::from_header(&header);
    let mut chromsizes = Vec::new();
    

    for tid in 0..header.target_count() {
        let name = header.tid2name(tid);
        let len = header.target_len(tid).unwrap();
        let csr: ChromSizeRecord = ChromSizeRecord {
            chrom: std::str::from_utf8(name).unwrap().to_string(), 
            size: len
        };
        chromsizes.push(csr);
    }

    let _ = bam.set_threads(threads);

    let mut ph: PairHeader = PairHeader::new();
    ph.from_chromsizes(chromsizes);

    let mut writer = common_writer(output);
    writer.write_all(ph.to_string().as_bytes()).unwrap();
  
    let mut idx = 0;

    while let Some(r) = bam.records().next() {
        let record = r.unwrap();
        
        if !record.is_paired() {
            continue
        }
        let Some(r2) = bam.records().next() else {
            continue
        };
        let record2 = r2.unwrap();

        if record.is_unmapped() || record2.is_unmapped(){
            continue 
        }

        if record.is_secondary() || record.is_supplementary() || record.is_duplicate() || record.is_quality_check_failed() {
            continue;
        }
        if record2.is_secondary() || record2.is_supplementary() || record2.is_duplicate() || record2.is_quality_check_failed() {
            continue;
        }


        if record.mapq() < min_mapq || record2.mapq() < min_mapq {
            continue
        }

        idx += 1;
        
        let mut chrom1 = std::str::from_utf8(header.tid2name(record.tid().try_into().unwrap())).unwrap().to_string();
        let mut chrom2 = std::str::from_utf8(header.tid2name(record2.tid().try_into().unwrap())).unwrap().to_string();
        let mut pos1 = record.pos() + 1;
        let mut pos2 = record2.pos() + 1;
        let mut strand1 = if record.is_reverse() { "-" } else { "+" };
        let mut strand2 = if record2.is_reverse() { "-" } else { "+" };

        if chrom1 > chrom2 || (chrom1 == chrom2 && pos1 > pos2) {
            std::mem::swap(&mut chrom1, &mut chrom2);
            std::mem::swap(&mut pos1, &mut pos2);
            std::mem::swap(&mut strand1, &mut strand2);
        }

        let mapq = std::cmp::min(record.mapq(), record2.mapq());
        
        writer.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            idx, chrom1, pos1,  chrom2, pos2, strand1, strand2, mapq,).as_bytes()).unwrap();

    }
    
} 

pub fn bam2pqs_pairend(
    input_bam: &String,
    min_mapq: u8,
    chunksize: usize,
    output: &String,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    use polars::prelude::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::Arc;
    use std::io::Write;

    polars::enable_string_cache();
    struct Contact {
        read_idx: u32,
        tid1: u32,
        pos1: u32,
        tid2: u32,
        pos2: u32,
        is_reverse1: bool,
        is_reverse2: bool,
        mapq: u8,
    }

    std::fs::create_dir_all(output)?;
    std::fs::create_dir_all(format!("{}/q0", output))?;
    std::fs::create_dir_all(format!("{}/q1", output))?;

    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(8);

    let header = Header::from_template(bam.header());
    let header_view = HeaderView::from_header(&header);

    let chrom_names: Arc<Vec<String>> = Arc::new(
        (0..header_view.target_count())
            .map(|tid| std::str::from_utf8(header_view.tid2name(tid)).unwrap().to_string())
            .collect()
    );

    {
        let mut cs_writer = common_writer(&format!("{}/_contigsizes", output));
        for tid in 0..header_view.target_count() {
            writeln!(cs_writer, "{}\t{}", chrom_names[tid as usize], header_view.target_len(tid).unwrap())?;
        }
    }

    
    {
        let mut md_writer = common_writer(&format!("{}/_metadata", output));
        writeln!(md_writer, "{{\"type\": \"pqs\"}}")?;
        
        let mut readme_writer = common_writer(&format!("{}/_readme", output));
        writeln!(readme_writer, "PQS format generated directly from BAM by cphasing-rs")?;
    }

    let (tx, rx) = bounded::<(usize, Vec<Contact>)>(10);
    let q0_total = Arc::new(AtomicU64::new(0));
    let q1_total = Arc::new(AtomicU64::new(0));
    let mut writer_handles = vec![];
    let num_writer_threads = (threads / 8).max(1).min(2);
    for _ in 0..num_writer_threads {
        let rx = rx.clone();
        let output = output.clone();
        let q0_total = Arc::clone(&q0_total);
        let q1_total = Arc::clone(&q1_total);
        let chrom_names = Arc::clone(&chrom_names);
        writer_handles.push(std::thread::spawn(move || {
            while let Ok((chunk_idx, contacts)) = rx.recv() {
                let n = contacts.len();
                let mut q0_read_idx = Vec::with_capacity(n);
                let mut q0_chrom1 = Vec::with_capacity(n);
                let mut q0_pos1 = Vec::with_capacity(n);
                let mut q0_chrom2 = Vec::with_capacity(n);
                let mut q0_pos2 = Vec::with_capacity(n);
                let mut q0_strand1 = Vec::with_capacity(n);
                let mut q0_strand2 = Vec::with_capacity(n);
                let mut q0_mapq = Vec::with_capacity(n);

                let mut q1_read_idx = Vec::with_capacity(n);
                let mut q1_chrom1 = Vec::with_capacity(n);
                let mut q1_pos1 = Vec::with_capacity(n);
                let mut q1_chrom2 = Vec::with_capacity(n);
                let mut q1_pos2 = Vec::with_capacity(n);
                let mut q1_strand1 = Vec::with_capacity(n);
                let mut q1_strand2 = Vec::with_capacity(n);
                let mut q1_mapq = Vec::with_capacity(n);

                for c in contacts {
                    let chr1 = chrom_names[c.tid1 as usize].as_str();
                    let chr2 = chrom_names[c.tid2 as usize].as_str();
                    let str1 = if c.is_reverse1 { "-" } else { "+" };
                    let str2 = if c.is_reverse2 { "-" } else { "+" };

                    q0_read_idx.push(c.read_idx);
                    q0_chrom1.push(chr1);
                    q0_pos1.push(c.pos1);
                    q0_chrom2.push(chr2);
                    q0_pos2.push(c.pos2);
                    q0_strand1.push(str1);
                    q0_strand2.push(str2);
                    q0_mapq.push(c.mapq);

                    if c.mapq >= 1 {
                        q1_read_idx.push(c.read_idx);
                        q1_chrom1.push(chr1);
                        q1_pos1.push(c.pos1);
                        q1_chrom2.push(chr2);
                        q1_pos2.push(c.pos2);
                        q1_strand1.push(str1);
                        q1_strand2.push(str2);
                        q1_mapq.push(c.mapq);
                    }
                }

                let mut df_q0 = DataFrame::new(vec![
                    Series::new("read_idx".into(), &q0_read_idx).into(),
                    Series::new("chrom1".into(), &q0_chrom1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos1".into(), &q0_pos1).into(),
                    Series::new("chrom2".into(), &q0_chrom2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos2".into(), &q0_pos2).into(),
                    Series::new("strand1".into(), &q0_strand1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("strand2".into(), &q0_strand2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("mapq".into(), &q0_mapq).into(),
                ]).unwrap();

                let out_q0_path = format!("{}/q0/{}.parquet", output, chunk_idx);
                let f0 = std::fs::File::create(&out_q0_path).unwrap();
                let _ = ParquetWriter::new(f0).finish(&mut df_q0);

                let mut df_q1 = DataFrame::new(vec![
                    Series::new("read_idx".into(), &q1_read_idx).into(),
                    Series::new("chrom1".into(), &q1_chrom1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos1".into(), &q1_pos1).into(),
                    Series::new("chrom2".into(), &q1_chrom2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos2".into(), &q1_pos2).into(),
                    Series::new("strand1".into(), &q1_strand1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("strand2".into(), &q1_strand2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("mapq".into(), &q1_mapq).into(),
                ]).unwrap();
                
                let out_q1_path = format!("{}/q1/{}.parquet", output, chunk_idx);
                let f1 = std::fs::File::create(&out_q1_path).unwrap();
                let _ = ParquetWriter::new(f1).finish(&mut df_q1);

                q0_total.fetch_add(df_q0.height() as u64, Ordering::Relaxed);
                q1_total.fetch_add(df_q1.height() as u64, Ordering::Relaxed);
               
            }
        }));
    }

    let mut idx = 0;
    let mut chunk_idx = 0;
    let mut batch = Vec::with_capacity(chunksize);
    let mut record = Record::new();
    let mut record2 = Record::new();
    loop {
        match bam.read(&mut record) {
            Some(Ok(())) => {},
            _ => break,
        }
        
        if !record.is_paired() {
            continue;
        }

        match bam.read(&mut record2) {
            Some(Ok(())) => {},
            _ => break, 
        }

        if record.is_unmapped() || record2.is_unmapped() {
            continue;
        }

        if record.is_secondary() || record.is_supplementary() || record.is_duplicate() || record.is_quality_check_failed() {
            continue;
        }
        if record2.is_secondary() || record2.is_supplementary() || record2.is_duplicate() || record2.is_quality_check_failed() {
            continue;
        }

        if record.mapq() < min_mapq || record2.mapq() < min_mapq {
            continue;
        }

        idx += 1;
        let mut tid1 = record.tid() as u32;
        let mut tid2 = record2.tid() as u32;
        let mut pos1 = record.pos() as u32 + 1;
        let mut pos2 = record2.pos() as u32 + 1;
        let mut is_reverse1 = record.is_reverse();
        let mut is_reverse2 = record2.is_reverse();

        if tid1 > tid2 || (tid1 == tid2 && pos1 > pos2) {
            std::mem::swap(&mut tid1, &mut tid2);
            std::mem::swap(&mut pos1, &mut pos2);
            std::mem::swap(&mut is_reverse1, &mut is_reverse2);
        }

        let mapq = std::cmp::min(record.mapq(), record2.mapq());
        
        batch.push(Contact {
            read_idx: idx as u32,
            tid1,
            pos1,
            tid2,
            pos2,
            is_reverse1,
            is_reverse2,
            mapq,
        });

        if batch.len() >= chunksize {
            tx.send((chunk_idx, batch)).unwrap();
            batch = Vec::with_capacity(chunksize);
            chunk_idx += 1;
        }
    }

    

    if !batch.is_empty() {
        tx.send((chunk_idx, batch)).unwrap();
    }

    drop(tx);
    for handle in writer_handles {
        handle.join().unwrap();
    }

    {
        let mut f = std::fs::File::create(format!("{}/_metadata_counts", output))?;
        writeln!(f, "q0\t{}", q0_total.load(Ordering::Relaxed))?;
        writeln!(f, "q1\t{}", q1_total.load(Ordering::Relaxed))?;
    }
   
    log::info!("Successfully converted BAM directly to PQS format: {}", output);
    Ok(())
}

pub fn bam2pairs(input_bam: &String, min_mapq: u8, output: &String, threads: usize) {
    
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
   
    let header = Header::from_template(bam.header());
    let header = HeaderView::from_header(&header);
    let mut chromsizes = Vec::new();
    

    for tid in 0..header.target_count() {
        let name = header.tid2name(tid);
        let len = header.target_len(tid).unwrap();
        let csr: ChromSizeRecord = ChromSizeRecord {
            chrom: std::str::from_utf8(name).unwrap().to_string(), 
            size: len
        };
        chromsizes.push(csr);
    }

    let _ = bam.set_threads(threads);

    let mut ph: PairHeader = PairHeader::new();
    ph.from_chromsizes(chromsizes);

    let mut writer = common_writer(output);
    writer.write_all(ph.to_string().as_bytes()).unwrap();
  
    let mut idx = 0;

    let mut current_qname: Vec<u8> = Vec::new();
    let mut group: Vec<Record> = Vec::new();

    let mut process_group = |records: &[Record], w: &mut dyn std::io::Write, idx_ref: &mut usize| {
        let mut valid: Vec<&Record> = records.iter().filter(|r| {
            !r.is_unmapped() &&
            !r.is_secondary() &&
            !r.is_duplicate() &&
            !r.is_quality_check_failed() &&
            r.mapq() >= min_mapq
        }).collect();

        if valid.len() < 2 {
            return;
        }

        // 按 Read 上的匹配起点顺序进行 5' 到 3' 排序
        valid.sort_by_key(|r| {
            let (qstart, _, _) = get_query_start_end(r);
            qstart
        });

        // 组内所有片段两两配对 (Pore-C Concatemer 展开)
        for i in 0..valid.len() {
            for j in (i + 1)..valid.len() {
                let r1 = valid[i];
                let r2 = valid[j];

                *idx_ref += 1;
                let mut chrom1 = std::str::from_utf8(header.tid2name(r1.tid().try_into().unwrap())).unwrap().to_string();
                let mut chrom2 = std::str::from_utf8(header.tid2name(r2.tid().try_into().unwrap())).unwrap().to_string();
                let mut pos1 = r1.pos() + 1;
                let mut pos2 = r2.pos() + 1;
                let mut strand1 = if r1.is_reverse() { "-" } else { "+" };
                let mut strand2 = if r2.is_reverse() { "-" } else { "+" };

                if chrom1 > chrom2 || (chrom1 == chrom2 && pos1 > pos2) {
                    std::mem::swap(&mut chrom1, &mut chrom2);
                    std::mem::swap(&mut pos1, &mut pos2);
                    std::mem::swap(&mut strand1, &mut strand2);
                }

                let mapq = std::cmp::min(r1.mapq(), r2.mapq());
                w.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                    *idx_ref, chrom1, pos1, chrom2, pos2, strand1, strand2, mapq).as_bytes()).unwrap();
            }
        }
    };

    for r in bam.records() {
        let record = r.unwrap();
        let qname = record.qname().to_vec();

        if current_qname.is_empty() {
            current_qname = qname;
            group.push(record);
        } else if qname == current_qname {
            group.push(record);
        } else {
            process_group(&group, &mut writer, &mut idx);
            current_qname = qname;
            group.clear();
            group.push(record);
        }
    }

    if !group.is_empty() {
        process_group(&group, &mut writer, &mut idx);
    }
} 

pub fn bam2pqs(
    input_bam: &String,
    min_mapq: u8,
    chunksize: usize,
    output: &String,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    use polars::prelude::*;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::Arc;
    use std::io::Write;

    polars::enable_string_cache();
    struct Contact {
        read_idx: u32,
        tid1: u32,
        pos1: u32,
        tid2: u32,
        pos2: u32,
        is_reverse1: bool,
        is_reverse2: bool,
        mapq: u8,
    }

    std::fs::create_dir_all(output)?;
    std::fs::create_dir_all(format!("{}/q0", output))?;
    std::fs::create_dir_all(format!("{}/q1", output))?;

    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(8);

    let header = Header::from_template(bam.header());
    let header_view = HeaderView::from_header(&header);

    let chrom_names: Arc<Vec<String>> = Arc::new(
        (0..header_view.target_count())
            .map(|tid| std::str::from_utf8(header_view.tid2name(tid)).unwrap().to_string())
            .collect()
    );

    {
        let mut cs_writer = common_writer(&format!("{}/_contigsizes", output));
        for tid in 0..header_view.target_count() {
            writeln!(cs_writer, "{}\t{}", chrom_names[tid as usize], header_view.target_len(tid).unwrap())?;
        }
    }

    
    {
        let mut md_writer = common_writer(&format!("{}/_metadata", output));
        writeln!(md_writer, "{{\"type\": \"pqs\"}}")?;
        
        let mut readme_writer = common_writer(&format!("{}/_readme", output));
        writeln!(readme_writer, "PQS format generated directly from BAM by cphasing-rs")?;
    }

    let (tx, rx) = bounded::<(usize, Vec<Contact>)>(10);
    let q0_total = Arc::new(AtomicU64::new(0));
    let q1_total = Arc::new(AtomicU64::new(0));
    let mut writer_handles = vec![];
    let num_writer_threads = (threads / 8).max(1).min(2);
    for _ in 0..num_writer_threads {
        let rx = rx.clone();
        let output = output.clone();
        let q0_total = Arc::clone(&q0_total);
        let q1_total = Arc::clone(&q1_total);
        let chrom_names = Arc::clone(&chrom_names);
        writer_handles.push(std::thread::spawn(move || {
            while let Ok((chunk_idx, contacts)) = rx.recv() {
                let n = contacts.len();
                let mut q0_read_idx = Vec::with_capacity(n);
                let mut q0_chrom1 = Vec::with_capacity(n);
                let mut q0_pos1 = Vec::with_capacity(n);
                let mut q0_chrom2 = Vec::with_capacity(n);
                let mut q0_pos2 = Vec::with_capacity(n);
                let mut q0_strand1 = Vec::with_capacity(n);
                let mut q0_strand2 = Vec::with_capacity(n);
                let mut q0_mapq = Vec::with_capacity(n);

                let mut q1_read_idx = Vec::with_capacity(n);
                let mut q1_chrom1 = Vec::with_capacity(n);
                let mut q1_pos1 = Vec::with_capacity(n);
                let mut q1_chrom2 = Vec::with_capacity(n);
                let mut q1_pos2 = Vec::with_capacity(n);
                let mut q1_strand1 = Vec::with_capacity(n);
                let mut q1_strand2 = Vec::with_capacity(n);
                let mut q1_mapq = Vec::with_capacity(n);

                for c in contacts {
                    let chr1 = chrom_names[c.tid1 as usize].as_str();
                    let chr2 = chrom_names[c.tid2 as usize].as_str();
                    let str1 = if c.is_reverse1 { "-" } else { "+" };
                    let str2 = if c.is_reverse2 { "-" } else { "+" };

                    q0_read_idx.push(c.read_idx);
                    q0_chrom1.push(chr1);
                    q0_pos1.push(c.pos1);
                    q0_chrom2.push(chr2);
                    q0_pos2.push(c.pos2);
                    q0_strand1.push(str1);
                    q0_strand2.push(str2);
                    q0_mapq.push(c.mapq);

                    if c.mapq >= 1 {
                        q1_read_idx.push(c.read_idx);
                        q1_chrom1.push(chr1);
                        q1_pos1.push(c.pos1);
                        q1_chrom2.push(chr2);
                        q1_pos2.push(c.pos2);
                        q1_strand1.push(str1);
                        q1_strand2.push(str2);
                        q1_mapq.push(c.mapq);
                    }
                }

                let mut df_q0 = DataFrame::new(vec![
                    Series::new("read_idx".into(), &q0_read_idx).into(),
                    Series::new("chrom1".into(), &q0_chrom1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos1".into(), &q0_pos1).into(),
                    Series::new("chrom2".into(), &q0_chrom2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos2".into(), &q0_pos2).into(),
                    Series::new("strand1".into(), &q0_strand1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("strand2".into(), &q0_strand2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("mapq".into(), &q0_mapq).into(),
                ]).unwrap();

                let out_q0_path = format!("{}/q0/{}.parquet", output, chunk_idx);
                let f0 = std::fs::File::create(&out_q0_path).unwrap();
                let _ = ParquetWriter::new(f0).finish(&mut df_q0);

                let mut df_q1 = DataFrame::new(vec![
                    Series::new("read_idx".into(), &q1_read_idx).into(),
                    Series::new("chrom1".into(), &q1_chrom1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos1".into(), &q1_pos1).into(),
                    Series::new("chrom2".into(), &q1_chrom2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("pos2".into(), &q1_pos2).into(),
                    Series::new("strand1".into(), &q1_strand1).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("strand2".into(), &q1_strand2).cast(&DataType::Categorical(None, CategoricalOrdering::Physical)).unwrap().into(),
                    Series::new("mapq".into(), &q1_mapq).into(),
                ]).unwrap();
                
                let out_q1_path = format!("{}/q1/{}.parquet", output, chunk_idx);
                let f1 = std::fs::File::create(&out_q1_path).unwrap();
                let _ = ParquetWriter::new(f1).finish(&mut df_q1);

                q0_total.fetch_add(df_q0.height() as u64, Ordering::Relaxed);
                q1_total.fetch_add(df_q1.height() as u64, Ordering::Relaxed);
               
            }
        }));
    }

    let mut idx = 0;
    let mut chunk_idx = 0;
    let mut batch = Vec::with_capacity(chunksize);
   
    let mut current_qname: Vec<u8> = Vec::new();
    let mut group: Vec<Record> = Vec::new();

    let mut process_pqs_group = |records: &[Record], batch_ref: &mut Vec<Contact>, idx_ref: &mut usize, chunk_idx_ref: &mut usize| {
        let mut valid: Vec<&Record> = records.iter().filter(|r| {
            !r.is_unmapped() &&
            !r.is_secondary() &&
            !r.is_duplicate() &&
            !r.is_quality_check_failed() &&
            r.mapq() >= min_mapq
        }).collect();

        if valid.len() < 2 {
            return;
        }

        valid.sort_by_key(|r| {
            let (qstart, _, _) = get_query_start_end(r);
            qstart
        });

        for i in 0..valid.len() {
            for j in (i + 1)..valid.len() {
                let r1 = valid[i];
                let r2 = valid[j];

                *idx_ref += 1;
                let mut tid1 = r1.tid() as u32;
                let mut tid2 = r2.tid() as u32;
                let mut pos1 = r1.pos() as u32 + 1;
                let mut pos2 = r2.pos() as u32 + 1;
                let mut is_reverse1 = r1.is_reverse();
                let mut is_reverse2 = r2.is_reverse();

                if tid1 > tid2 || (tid1 == tid2 && pos1 > pos2) {
                    std::mem::swap(&mut tid1, &mut tid2);
                    std::mem::swap(&mut pos1, &mut pos2);
                    std::mem::swap(&mut is_reverse1, &mut is_reverse2);
                }

                let mapq = std::cmp::min(r1.mapq(), r2.mapq());

                batch_ref.push(Contact {
                    read_idx: *idx_ref as u32,
                    tid1,
                    pos1,
                    tid2,
                    pos2,
                    is_reverse1,
                    is_reverse2,
                    mapq,
                });

                if batch_ref.len() >= chunksize {
                    tx.send((*chunk_idx_ref, std::mem::take(batch_ref))).unwrap();
                    *batch_ref = Vec::with_capacity(chunksize);
                    *chunk_idx_ref += 1;
                }
            }
        }
    };

    let mut record = Record::new();
    loop {
        match bam.read(&mut record) {
            Some(Ok(())) => {
                let qname = record.qname().to_vec();
                if current_qname.is_empty() {
                    current_qname = qname;
                    group.push(record.clone());
                } else if qname == current_qname {
                    group.push(record.clone());
                } else {
                    process_pqs_group(&group, &mut batch, &mut idx, &mut chunk_idx);
                    current_qname = qname;
                    group.clear();
                    group.push(record.clone());
                }
            },
            _ => break,
        }
    }

    if !group.is_empty() {
        process_pqs_group(&group, &mut batch, &mut idx, &mut chunk_idx);
    }

    if !batch.is_empty() {
        tx.send((chunk_idx, batch)).unwrap();
    }

    drop(tx);
    for handle in writer_handles {
        handle.join().unwrap();
    }

    {
        let mut f = std::fs::File::create(format!("{}/_metadata_counts", output))?;
        writeln!(f, "q0\t{}", q0_total.load(Ordering::Relaxed))?;
        writeln!(f, "q1\t{}", q1_total.load(Ordering::Relaxed))?;
    }
   
    log::info!("Successfully converted BAM directly to PQS format: {}", output);
    Ok(())
}

pub fn bam2fastq(input_bams: &Vec<&String>, output: &String, threads: usize) {
    let writer = common_writer(output);
    let mut wtr = FastqWriter::new(writer);

    for input_bam in input_bams {
        let mut bam = if *input_bam == &String::from("-") {
            Reader::from_stdin().expect("Failed to read from stdin")
        } else {
            Reader::from_path(input_bam).expect("Failed to read from the provided path")
        };
        
        let _ = bam.set_threads(threads);

        while let Some(r) = bam.records().next() {
            let record = r.unwrap();
            
            let id = String::from_utf8(record.qname().to_vec()).unwrap();
            let seq = record.seq().as_bytes();
            let qual = record.qual();
            
            let mut ascii_qual = qual.to_vec();
            for q in &mut ascii_qual {
                *q += 33;
            }

            if let Err(err) = wtr.write(&id, None, &seq, &ascii_qual) {
                if err.kind() == std::io::ErrorKind::BrokenPipe {
                    std::process::exit(0);
                } else {
                    panic!("Failed to write fastq: {:?}", err);
                }
            }
        }
    }
}

pub fn bam2fasta(input_bam: &String, output:&String, threads: usize) {
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    
    let header = Header::from_template(bam.header());
    let header = HeaderView::from_header(&header);

    let _ = bam.set_threads(threads);

    let writer = common_writer(output);
    let mut wtr = FastaWriter::new(writer);
    while let Some(r) = bam.records().next() {
        let record = r.unwrap();
        
        let id = String::from_utf8(record.qname().to_vec()).unwrap();
        let seq = record.seq().as_bytes();
        
        let seq_record = FastaRecord::with_attrs(&id, None, &seq,);
        let _ = wtr.write_record(&seq_record);

    }

}


fn calculate_quartiles(data: &Vec<usize>) -> (usize, usize, usize) {
   
    let n = data.len();
    let q1_index = (n as f64 * 0.25).ceil() as usize - 1;
    let q2_index = (n as f64 * 0.5).ceil() as usize - 1;
    let q3_index = (n as f64 * 0.75).ceil() as usize - 1;

    let q1 = data[q1_index];
    let q2 = data[q2_index];
    let q3 = data[q3_index];

    (q1, q2, q3)
}

pub fn bamstat_hic(input_bams: &Vec<&String>, output: &String, threads: usize) {
    let mut writer = common_writer(output);
    let res = input_bams.par_iter().map(|input_bam| {
        let mut bam = if *input_bam == &String::from("-") {
            Reader::from_stdin().expect("Failed to read from stdin")
        } else {
            Reader::from_path(input_bam).expect("Failed to read from the provided path")
        };

        let header = Header::from_template(bam.header());
        let header = HeaderView::from_header(&header);
    
        let _ = bam.set_threads(threads);
        let mut unmap_counts = 0;
        let mut multiple_counts = 0;
        let mut unique_counts = 0;
        let mut singleton_counts = 0;
        while let Some(r) = bam.records().next() {
            let record = r.unwrap();
            
            let Some(r2) = bam.records().next() else {
                continue
            };
            let record2 = r2.unwrap();
    
            if record.is_unmapped() || record2.is_unmapped(){
                unmap_counts += 1;
                continue 
            }


            let mapq1 = record.mapq();
            let mapq2 = record2.mapq();
            
            if (mapq1 == 0) && (mapq2 == 0) {
                multiple_counts += 1;
            } else if (mapq1 == 0) | (mapq2 == 0) {
                singleton_counts += 1;
            } else {
                unique_counts += 1;
            }
        }

        (input_bam, unmap_counts, unique_counts, multiple_counts, singleton_counts)
    }).collect::<Vec<_>>();

    
    writeln!(writer, "file\tunmapped\tunique\tmultiple\tsingleton").unwrap();
    for record in res {
        
        let (input_bam, unmap_counts, unique_counts, multiple_counts, singleton_counts) = record;
        // basename of input_bam 
        let input_bam = std::path::Path::new(input_bam);
        let input_bam = input_bam.file_name().unwrap().to_str().unwrap();
        writeln!(writer, "{}\t{}\t{}\t{}\t{}", 
                    input_bam, unmap_counts, unique_counts, 
                    multiple_counts, singleton_counts).unwrap();
    }
    

}


pub fn bamstat_porec(input_bams: &Vec<&String>, output: &String, threads: usize) {
    let mut writer = common_writer(output);
    let res = input_bams.par_iter().map(|input_bam| {
        let mut bam = if *input_bam == &String::from("-") {
            Reader::from_stdin().expect("Failed to read from stdin")
        } else {
            Reader::from_path(input_bam).expect("Failed to read from the provided path")
        };

        let header = Header::from_template(bam.header());
        let header = HeaderView::from_header(&header);
    
        let _ = bam.set_threads(threads);
        let mut unmap_counts = 0;
        let mut multiple_counts = 0;
        let mut unique_counts = 0;
        let mut unique_q2_counts = 0;
        while let Some(r) = bam.records().next() {
            let record = r.unwrap();
            
            let Some(r2) = bam.records().next() else {
                continue
            };

            if record.is_secondary() {
                continue;
            }
        
            if record.is_unmapped() {
                unmap_counts += 1;
                continue 
            }

            let mapq = record.mapq();
            
            if mapq == 0 {
                multiple_counts += 1;
            } else if mapq >= 1 {
                unique_counts += 1;
                if mapq >= 2 {
                    unique_q2_counts += 1;
                }
            }

            
        }

        (input_bam, unmap_counts, unique_counts, multiple_counts, unique_q2_counts)
    }).collect::<Vec<_>>();

    
    writeln!(writer, "file\tunmapped\tunique\tmultiple\tunique_q2").unwrap();
    for record in res {
        
        let (input_bam, unmap_counts, unique_counts, multiple_counts, unique_q2_counts) = record;
        // basename of input_bam 
        let input_bam = std::path::Path::new(input_bam);
        let input_bam = input_bam.file_name().unwrap().to_str().unwrap();
        writeln!(writer, "{}\t{}\t{}\t{}\t{}", 
                    input_bam, unmap_counts, unique_counts, 
                    multiple_counts, unique_q2_counts).unwrap();
    }
    


}


pub fn bamstat(input_bams: &Vec<&String>, output: &String, threads: usize) {

    let mut writer = common_writer(output);
        
    writeln!(writer, "file\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tN50\tGC(%)").unwrap();
    let res = input_bams.par_iter().map(|input_bam| {
        let mut bam = if *input_bam == &String::from("-") {
            Reader::from_stdin().expect("Failed to read from stdin")
        } else {
            Reader::from_path(input_bam).expect("Failed to read from the provided path")
        };
        
        let header = Header::from_template(bam.header());
        let header = HeaderView::from_header(&header);
    
        let _ = bam.set_threads(threads);
    
        let mut seq_len_vec: Vec<usize> = Vec::new();
        let mut gc_count: u64 = 0;
        let total_qual: u64 = 0;
        while let Some(r) = bam.records().next() {
            let record = r.unwrap();
            
            let seq_len = record.seq().len();
            
            let seq = record.seq();
            for base in seq.as_bytes() {
                if base == b'G' || base == b'C' {
                    gc_count += 1;
                }
            }
    
            // let qual = record.qual();
            // for &q in qual {
            //     total_qual += q as u64;
            // }
            
            seq_len_vec.push(seq_len);
    
        } 
    
        let total_count = seq_len_vec.len();
    
        seq_len_vec.par_sort();
        let total_len: usize = seq_len_vec.par_iter().sum();
        let gc_content = if total_len > 0 {
            (gc_count as f64 / total_len as f64) * 100.0
        } else {
            0.0
        };
    
        let mut n50 = 0;
        let mut n50_len = 0;
        let mut n50_count = 0;
        for len in seq_len_vec.iter().rev() {
            n50_len += len;
            n50_count += 1;
            if n50_len >= total_len / 2 {
                n50 = *len;
                break;
            }
        }
    
        // calculate Q1 Q2 Q3 
        let (q1, q2, q3) = calculate_quartiles(&seq_len_vec);
    
        let min_len = seq_len_vec[0].clone();
        let max_len = seq_len_vec.last().unwrap().clone();
    
        let avg_len = total_len as f64 / total_count as f64;
    
        
        (input_bam, total_count, total_len, min_len, 
            avg_len, max_len, q1, q2, q3, n50, gc_content)
        
    }).collect::<Vec<_>>();
    

    
    for r in res {
        let (input_bam, total_count, total_len, min_len, 
            avg_len, max_len, q1, q2, q3,
            n50, gc_content) = r;
        writeln!(writer, "{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{:.2}%", 
                    input_bam, total_count, total_len, min_len, 
                    avg_len, max_len, q1, q2, q3,
                    n50, gc_content).unwrap();
    }
    
}

pub fn phase_reads(
    input_bam: &String,
    contigs_group_file: &String,
    output: &String,
    min_mapq: u8,
    threads: usize,
) {
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};
    use std::fs::File;

    let mut contig_to_group: HashMap<String, String> = HashMap::new();
    let file = File::open(contigs_group_file).expect("Failed to open contigs group file");
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.unwrap();
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.len() >= 2 {
            contig_to_group.insert(parts[0].to_string(), parts[1].to_string());
        }
    }
    log::info!("Loaded {} contig-to-group mappings.", contig_to_group.len());

    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(threads);

    let header = Header::from_template(bam.header());
    let hv = HeaderView::from_header(&header);

    let tnames: Vec<String> = (0..hv.target_count())
        .map(|tid| std::str::from_utf8(hv.tid2name(tid)).unwrap().to_string())
        .collect();

    let mut wtr = common_writer(output);
    let mut prev_q: Option<String> = None;
    let mut group_scores: HashMap<String, u32> = HashMap::new();

    let mut emit_read = |qname: &str, scores: &mut HashMap<String, u32>, w: &mut dyn Write| {
        if let Some((best_group, _)) = scores.iter().max_by_key(|entry| entry.1) {
            writeln!(w, "{}\t{}", qname, best_group).unwrap();
        } else {
            writeln!(w, "{}\tUnassigned", qname).unwrap();
        }
        scores.clear();
    };

    log::info!("Scanning BAM file to phase reads (requires name-sorted BAM)...");

    let mut total_reads = 0;
    let mut phased_reads = 0;

    for r in bam.records() {
        let record = r.unwrap();
        
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        
        let is_switch = match prev_q.as_ref() {
            Some(pq) => pq != &qname,
            None => false,
        };

        if is_switch {
            emit_read(prev_q.as_ref().unwrap(), &mut group_scores, &mut wtr);
            total_reads += 1;
            if !group_scores.is_empty() { phased_reads += 1; }
            prev_q = Some(qname.clone());
        } else if prev_q.is_none() {
            prev_q = Some(qname.clone());
        }

        if !record.is_unmapped() && record.mapq() >= min_mapq && record.tid() >= 0 {
            let tname = &tnames[record.tid() as usize];
            if let Some(group) = contig_to_group.get(tname) {
                *group_scores.entry(group.clone()).or_insert(0) += 1;
            }
        }
    }

    if let Some(pq) = prev_q {
        emit_read(&pq, &mut group_scores, &mut wtr);
        total_reads += 1;
        if !group_scores.is_empty() { phased_reads += 1; }
    }

    log::info!(
        "Phasing finished. Total reads processed: {}, successfully assigned: {}",
        total_reads, phased_reads
    );
    log::info!("Output written to: {}", output);
}

pub fn prune_bam(
    input_bam: &String,
    prune_table: &String,
    output: &String,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::HashSet;
    use std::io::{BufRead, BufReader};
    use crate::core::common_reader;

    let mut blacklist_names = HashSet::new();
    let f = common_reader(prune_table);
    let rdr = BufReader::new(f);
    for line in rdr.lines().flatten() {
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') { continue; }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 2 { continue; }
        let c1 = fields[0].to_string();
        let c2 = fields[1].to_string();
        if c1 <= c2 {
            blacklist_names.insert((c1, c2));
        } else {
            blacklist_names.insert((c2, c1));
        }
    }
    log::info!("Loaded {} blacklist contig pairs from {}", blacklist_names.len(), prune_table);

    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(threads);

    let header = Header::from_template(bam.header());
    let hv = HeaderView::from_header(&header);

    let mut name_to_tid = HashMap::new();
    for tid in 0..hv.target_count() {
        let name_bytes = hv.tid2name(tid);
        if let Ok(name_str) = std::str::from_utf8(name_bytes) {
            name_to_tid.insert(name_str.to_string(), tid as i32);
        }
    }


    let mut blacklist_tids = HashSet::new();
    for (c1, c2) in blacklist_names {
        if let (Some(&tid1), Some(&tid2)) = (name_to_tid.get(&c1), name_to_tid.get(&c2)) {
            let pair = if tid1 <= tid2 { (tid1, tid2) } else { (tid2, tid1) };
            blacklist_tids.insert(pair);
        }
    }


    let mut wtr = Writer::from_path(output, &header, bam::Format::Bam).unwrap();
    let _ = wtr.set_threads(threads);

    let mut total_records = 0;
    let mut pruned_records = 0;

    for r in bam.records() {
        let record = r?;
        total_records += 1;

        if !record.is_unmapped() {
            let tid = record.tid();
            let mtid = record.mtid();

            if tid >= 0 && mtid >= 0 {
                let pair = if tid <= mtid { (tid, mtid) } else { (mtid, tid) };
                if blacklist_tids.contains(&pair) {
                    pruned_records += 1;
                    continue; 
                }
            }
        }

        wtr.write(&record)?;
    }

    log::info!(
        "BAM pruning finished. Total records: {}, Pruned (removed): {}, Output written to: {}",
        total_records, pruned_records, output
    );

    Ok(())
}


struct ContigInterval {
    start: u32,
    end: u32,
    name: String,
}

fn load_contig_map(
    contig_bed: &str,
) -> Result<(HashMap<String, Vec<ContigInterval>>, HashMap<String, u64>), Box<dyn std::error::Error>> {
    use std::io::{BufRead, BufReader};
    use std::fs::File;

    let mut map: HashMap<String, Vec<ContigInterval>> = HashMap::new();
    let mut contig_sizes: HashMap<String, u64> = HashMap::new();
    let file = File::open(contig_bed)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.trim().split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let chrom = parts[0].to_string();
        let start = parts[1].parse::<u32>()?;
        let end = parts[2].parse::<u32>()?;
        let contig = parts[3].to_string();

        contig_sizes.insert(contig.clone(), (end - start) as u64);
        map.entry(chrom).or_default().push(ContigInterval {
            start,
            end,
            name: contig,
        });
    }

    for intervals in map.values_mut() {
        intervals.sort_by_key(|x| x.start);
    }

    Ok((map, contig_sizes))
}

fn map_chr_to_ctg(
    chrom: &str,
    pos: u32,
    map: &HashMap<String, Vec<ContigInterval>>,
) -> Option<(String, u32)> {
    let intervals = map.get(chrom)?;
    let idx = intervals.binary_search_by(|val| {
        if pos < val.start {
            std::cmp::Ordering::Greater
        } else if pos >= val.end {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Equal
        }
    });

    if let Ok(idx) = idx {
        let interval = &intervals[idx];
        Some((interval.name.clone(), pos - interval.start))
    } else {
        None
    }
}

pub fn chr_to_ctg(
    input_bam: &String,
    contig_bed: &String,
    output_bam: &String,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let (contig_map, contig_sizes) = load_contig_map(contig_bed)?;
    log::info!("Loaded {} chromosome maps from BED file.", contig_map.len());

    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(threads);

    let old_header = Header::from_template(bam.header());
    let old_header_view = HeaderView::from_header(&old_header);

    let mut new_header = Header::new();
    let mut names: Vec<_> = contig_sizes.keys().cloned().collect();
    names.sort();
    for name in &names {
        let size = contig_sizes.get(name).unwrap();
        let mut header_rec = HeaderRecord::new(b"SQ");
        header_rec.push_tag(b"SN", name);
        header_rec.push_tag(b"LN", &size.to_string());
        new_header.push_record(&header_rec);
    }
    let new_header_view = HeaderView::from_header(&new_header);

    let mut wtr = Writer::from_path(output_bam, &new_header, bam::Format::Bam)?;
    let _ = wtr.set_threads(threads);

    let mut total_records = 0;
    let mut mapped_records = 0;

    for r in bam.records() {
        let record = r?;
        total_records += 1;

        if record.is_unmapped() {
            let mut new_record = record.clone();
            let mtid = record.mtid();
            if mtid >= 0 {
                let mchrom_name = std::str::from_utf8(old_header_view.tid2name(mtid as u32))?;
                let mpos = record.mpos() as u32;
                if let Some((mcontig_name, new_mpos)) = map_chr_to_ctg(mchrom_name, mpos, &contig_map) {
                    if let Some(new_mtid) = new_header_view.tid(mcontig_name.as_bytes()) {
                        new_record.set_mtid(new_mtid as i32);
                        new_record.set_mpos(new_mpos as i64);
                    } else {
                        new_record.set_mtid(-1);
                        new_record.set_mpos(-1);
                    }
                } else {
                    new_record.set_mtid(-1);
                    new_record.set_mpos(-1);
                }
            }
            wtr.write(&new_record)?;
            continue;
        }

        let tid = record.tid();
        let pos = record.pos() as u32;
        let chrom_name = std::str::from_utf8(old_header_view.tid2name(tid as u32))?;

        if let Some((contig_name, new_pos)) = map_chr_to_ctg(chrom_name, pos, &contig_map) {
            let mut new_record = record.clone();
            let new_tid = new_header_view.tid(contig_name.as_bytes()).ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::NotFound, format!("Contig {} not found in new header", contig_name))
            })? as i32;

            new_record.set_tid(new_tid);
            new_record.set_pos(new_pos as i64);

            let mtid = record.mtid();
            if mtid >= 0 {
                let mchrom_name = std::str::from_utf8(old_header_view.tid2name(mtid as u32))?;
                let mpos = record.mpos() as u32;
                if let Some((mcontig_name, new_mpos)) = map_chr_to_ctg(mchrom_name, mpos, &contig_map) {
                    if let Some(new_mtid_val) = new_header_view.tid(mcontig_name.as_bytes()) {
                        new_record.set_mtid(new_mtid_val as i32);
                        new_record.set_mpos(new_mpos as i64);
                    } else {
                        new_record.set_mtid(-1);
                        new_record.set_mpos(-1);
                    }
                } else {
                    new_record.set_mtid(-1);
                    new_record.set_mpos(-1);
                }
            }

            wtr.write(&new_record)?;
            mapped_records += 1;
        }
    }

    log::info!(
        "BAM coordinates conversion finished. Total records: {}, Mapped to contigs: {}, Written to: {}",
        total_records, mapped_records, output_bam
    );

    Ok(())
}