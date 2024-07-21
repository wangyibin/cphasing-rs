
use std::collections::HashMap;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};

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
    bam.set_threads(8);
    let header = Header::from_template(bam.header());
    let mut i = 0;
    let mut j = 0;
    let mut wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam).unwrap();
    log::info!("write {} records to {}", i, format!("{}_{}.bam", output_prefix, j));
    for r in bam.records() {
        let record = r?;
        wtr.write(&record)?;
        i += 1;
        if i == record_num {
            j += 1;
            wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam)?;
            log::info!("write {} records to {}", i, format!("{}_{}.bam", output_prefix, j));
            
            i = 0;
        }
    }
    
    Ok(())
}

pub fn slide2raw(input_bam: &String, output: &String, threads: usize) {
    let mut bam = if input_bam == &String::from("-") {
        Reader::from_stdin().expect("Failed to read from stdin")
    } else {
        Reader::from_path(input_bam).expect("Failed to read from the provided path")
    };
    bam.set_threads(threads);

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

    bam.set_threads(threads);

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

