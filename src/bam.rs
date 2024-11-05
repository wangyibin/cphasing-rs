
use std::collections::HashMap;
use bio::io::fastq::{Reader as FastqReader, Record as FastqRecord, Writer as FastqWriter};
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer, ext::BamRecordExtensions};


    use crate::core::{ 
        ChromSizeRecord,
        common_writer};
use crate::pairs::{ Pairs, PairHeader };


// split bam by record number and write to different files 
pub fn split_bam(input_bam: &String, output_prefix: &String, 
             record_num: usize) -> Result<(), Box<dyn std::error::Error>> {

                
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    let _ = bam.set_threads(8);
    let header = Header::from_template(bam.header());
    let mut i = 0;
    let mut j = 0;
    let mut wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam).unwrap();
    let _ = wtr.set_threads(8);
    log::info!("write {} records to {}", record_num, format!("{}_{}.bam", output_prefix, j));

    let mut read_name: Vec<u8> = Vec::new();
    let mut previous_read_name: Vec<u8> = Vec::new();
    let mut line_num = 0;
    for r in bam.records() {
        let record = r?;
        read_name.clear();
        read_name.extend_from_slice(record.qname());
        
        i += 1;
        
        if i == record_num {
            if read_name != previous_read_name {

                j += 1;
                wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam)?;
                let _ = wtr.set_threads(8);
                log::info!("write {} records to {}", line_num, format!("{}_{}.bam", output_prefix, j));
                
                i = 0;
                line_num = 0;
            } else {
                i -= 1;
            }
        }
        line_num += 1;
        wtr.write(&record)?;
        std::mem::swap(&mut read_name, &mut previous_read_name);
    } 

    log::info!("write {} records to {}", i, format!("{}_{}.bam", output_prefix, j));
    

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
    wtr.set_threads(threads);
    while let Some(r) = bam.records().next() {
        let record = r.unwrap();
        let mut new_record = record.clone();
        let mut read_id = std::str::from_utf8(record.qname()).unwrap();
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

fn get_query_start_end(record: &Record) -> (u32, u32, u32) {
    let mut query_start = 0;
    let mut query_end = 0;
    let mut query_pos = 0;
    let mut query_length = 0;
    for cigar in record.cigar().iter() {
        match cigar.char() {
            'M' | '=' | 'X' | 'I' => {
                if query_start == 0 {
                    query_start = query_pos;
                }
                query_pos += cigar.len();
                query_length += cigar.len();
                query_end += cigar.len()
            }
            'S' | 'H' => {
                if query_start == 0 {
                    query_start = query_pos + cigar.len();
                }
                query_pos += cigar.len();
                query_length += cigar.len();
            }
            'D' | 'N' | 'P' => {

            }
            _ => {}
        }
    }
    // query_end = query_pos;
    query_end = query_start + query_end;

    (query_start, query_end, query_length)
}

pub fn bam2paf(input_bam: &String, output: &String, threads: usize, is_secondary: bool) {
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
    let mut idx = 0;

    let mut writer = common_writer(output);
    
    while let Some(r) = bam.records().next() {
        let record = r.unwrap();

        if record.is_unmapped() { 
            continue;
        }

        if record.is_secondary() & !is_secondary {
            continue
        }

        let chrom = std::str::from_utf8(header.tid2name(record.tid().try_into().unwrap())).unwrap().to_string();
        let strand = if record.is_reverse() { "-" } else { "+" };
        
        let (mut match_length, mut deletion_length, mut insertion_length, mut mm, mut alignment_length) = (0, 0, 0, 0, 0);

        for cigar in record.cigar().iter() {
            let len = cigar.len() as i64;
            match cigar.char() {
                'M' | '=' | 'X' => {
                    match_length += len;
                    alignment_length += len;
                    if cigar.char() == 'X' {
                        mm += len;
                    }
                },
                'D' => {
                    deletion_length += len;
                    alignment_length += len;
                },
                'I' => {
                    insertion_length += len;
                    alignment_length += len;
                },
                'N' | 'P' => {
                    alignment_length += len;
                },
                _ => {}
            }
        }


        let nm: i64 = match record.aux(b"NM") {
            Ok(value) => {
                match value {
                    Aux::U8(v) => v.try_into().unwrap(),
                    Aux::U16(v) => v.try_into().unwrap(),
                    Aux::U32(v) => v.try_into().unwrap(),
                    Aux::I32(v) => v.try_into().unwrap(),
                    _ => 0, 
                }
                
            },
            Err(e) => 0
        };

        match_length = match_length - (nm - (insertion_length + deletion_length));


        let (qstart, qend, qlen) = get_query_start_end(&record);
    
        let qname = std::str::from_utf8(record.qname()).unwrap();
        let tname = std::str::from_utf8(header.tid2name(record.tid().try_into().unwrap())).unwrap();
        let tlen =  header.target_len(record.tid().try_into().unwrap()).unwrap() as usize; 
        let tstart = record.reference_start();
        let tend = record.reference_end();
        let mapq = record.mapq();
        writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                qname, qlen, qstart, qend, strand, 
                tname, tlen, tstart, tend, match_length, alignment_length, mapq);

        idx += 1;

    }
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

pub fn bam2fastq(input_bam: &String, output:&String, threads: usize) {
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    
    let header = Header::from_template(bam.header());
    let header = HeaderView::from_header(&header);

    let _ = bam.set_threads(threads);

    let mut writer = common_writer(output);
    let mut wtr = FastqWriter::new(writer);
    while let Some(r) = bam.records().next() {
        let record = r.unwrap();
        
        let id = String::from_utf8(record.qname().to_vec()).unwrap();
        let seq = record.seq().as_bytes();
        let qual = record.qual();
        let seq_record = FastqRecord::with_attrs(&id, None, &seq, &qual);
        wtr.write_record(&seq_record);

    }
    


}

