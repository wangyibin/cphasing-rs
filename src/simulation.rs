use anyhow::Result as anyResult;
use noodles::{core::Region,
              core::Position,
              vcf};
use rust_htslib::bam::{
    self, 
    record::Aux, record::AuxArray, 
    record::Cigar, record::CigarStringView,
    Header, HeaderView, header::HeaderRecord,
    Read, Reader, Record,
    Writer, 
};

use rand::distributions::{Distribution, Uniform};
use rand::{ Rng, SeedableRng };
use rand::rngs::StdRng;

use std::path::Path;

use crate::core::{ BaseTable, common_writer };
use crate::fastx::Fastx;


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
        let target =std::str::from_utf8(
                h.tid2name(r.tid() as u32)).unwrap();
        let target_start = r.pos();
        let target_end = cigar.end_pos();
        let mut query_start = 0;
        let mut query_end =  cigar.read_pos(
                                target_end as u32 - 1, 
                                false, false).unwrap().unwrap() as i64;
        let query_len = i64::try_from(seq_len).unwrap();
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


// simulation a read bam through split read in a input bam file 
// input: bam file, output: bam file
pub fn simulation_from_split_read(input_bam: &String, output_bam: &String, min_quality: u8) {
    let mut bam = Reader::from_path(&Path::new(input_bam)).unwrap();
    let header = bam::Header::from_template(bam.header());
    let header_view = HeaderView::from_header(&header);

    let mut new_header = Header::new();
    let mut header_records = HeaderRecord::new(b"HD");
    header_records.push_tag(b"VN", &"1.6");
    header_records.push_tag(b"SO", &"coordinate");
    new_header.push_record(&header_records);
    let mut header_records = HeaderRecord::new(b"RG");
    header_records.push_tag(b"ID", &"1");
    new_header.push_record(&header_records);


    let mut writer = Writer::from_path(output_bam, &new_header, bam::Format::Bam).unwrap();

    let mut old_read_id = Vec::new();
    let mut alignment_idx = 0;
    'outer: for r in bam.records() {
        let record = r.unwrap();
        // if String::from_utf8(record.qname().to_vec()).unwrap() != "0a23f67c-3d91-4f24-8561-9e211abd2021" {
        //     continue
        // }
        if record.mapq() < min_quality {
            continue 'outer;
        }

        if record.is_unmapped() {
            continue 'outer;
        }
        if record.is_secondary() {
            continue 'outer;
        }
        
        let qname_raw = record.qname().to_vec();
        let mut qname = record.qname().to_vec();

        if old_read_id == qname_raw {
                alignment_idx += 1;
        } else {
                alignment_idx = 0;
        }
        qname.push(b'_');
        qname.push(alignment_idx.to_string().as_bytes()[0]);

        

        let cigar = record.cigar();
        let seq_len = record.seq_len() as u64;
        let r_parser = RecordParser::parse(&record, &header_view, seq_len);
        
        let mut query_seq = record.seq().as_bytes().to_vec();
        let query_len = query_seq.len() as i64;
        if record.is_reverse() {
            query_seq.reverse();

            for i in 0..query_seq.len() {
                query_seq[i] = match query_seq[i] {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'G' => b'C',
                    b'C' => b'G',
                    _ => query_seq[i],
                };
               
            }
            
        }
        
        
        let c_counts = query_seq.iter().filter(|&c| *c == b'C').count();
        // println!("{}", record.qname().iter().map(|x| *x as char).collect::<String>());
      
        // let mut mod_pos = Vec::new();
        // if let Ok(mods) = record.basemods_iter() {
        //     for (idx, mod_res) in mods.enumerate() {
        //         if let Ok((position, m)) = mod_res {
        //             // println!("{}: {:.?}", position, record.seq().as_bytes().to_vec()[position as usize..=(position + 1) as usize].iter().map(|x| *x as char).collect::<String>());
        //             mod_pos.push(position);
        //         }
        //     }
        // };

        if query_len < 1000 {
            // let mut new_record = Record::from(record.clone());
            // new_record.set(&qname, None, &seq, record.qual());
            // new_record.set_pos(0);
            // new_record.set_flags(0x4);
            // new_record.set_mapq(0);
            // new_record.set_tid(-1);
            // new_record.remove_aux(b"MD");
            // new_record.remove_aux(b"AS");
            // new_record.remove_aux(b"NM");
            // new_record.remove_aux(b"XA");
            // new_record.remove_aux(b"SA");
            // writer.write(&new_record).unwrap();
            continue
        }

        let mut query_qual = record.qual().to_vec();
        // let mut query_qual = &mut qual[r_parser.query_start as usize..r_parser.query_end as usize];
        if record.is_reverse() {
            query_qual.reverse();
        }

        let range = Uniform::new(1, query_len/500);
        let seed: u64 = 345;
    
        let mut rng: StdRng = StdRng::seed_from_u64(seed);
        let count = range.sample(&mut rng);
       
        let range = Uniform::new(0, query_len);
        let mut rng = rand::thread_rng();
        let mut pos_list = Vec::new(); 
        for _ in 0..count * 2 {
            let random_pos = range.sample(&mut rng);
            if pos_list.contains(&random_pos) {
                continue
            }
            if random_pos < &r_parser.query_start as &i64 +500 || random_pos > &r_parser.query_end as &i64 - 500 {
                continue
            }
            pos_list.push(random_pos);

        }
        pos_list.sort();

        // remove too short 
        let mut new_pos_list = Vec::new();
        for (idx, pos) in pos_list.iter().enumerate() {
            if idx == 0 {
                new_pos_list.push(*pos);
            } else {
                if pos - pos_list[idx - 1] > 500 {
                    new_pos_list.push(*pos);
                }
            }
        }
        let mut pos_list_1 = new_pos_list.clone();
        let mut pos_list_2 = new_pos_list.clone();
        pos_list_1.insert(0, 0);
        pos_list_2.push((query_len as usize).try_into().unwrap());
        
        let mut mm_list = match record.aux(b"MM").unwrap() {
            Aux::String(v) => {
                v.strip_suffix(";").unwrap().split(",").collect::<Vec<_>>()
            },
            _ => panic!("MD tag not found"),
        };
        let mm_cap = mm_list[0];
        mm_list.remove(0);
        let mm_list = mm_list.iter().map(|x| x.parse::<u32>().unwrap()).collect::<Vec<u32>>();
        // println!("{:?} {}", mm_list, mm_list.len());
        let ml_list = match record.aux(b"ML").unwrap() {
            Aux::ArrayU8(v) => {
                v.iter().map(|x| x as u16).collect::<Vec<_>>()
            },
            _ => panic!("ML tag not found"),
        };

        let mut mm_boundary: usize = 0;
        let mut mm_init: usize = 0;
        let mut least_c_counts = 0;

        for (i, (a, b)) in pos_list_1.iter().zip(pos_list_2.iter()).enumerate() {

            let mut new_record = Record::from(record.clone());
            let mut new_qname = qname.clone();
            let new_seq = &query_seq[*a as usize..*b as usize];
            let c_counts = new_seq.iter().filter(|&c| *c == b'C').count();

            // println!("{} {}", c_scounts, new_seq.len());
            let new_qual = &query_qual[*a as usize..*b as usize];

            new_qname.push(b'_');
            i.to_string().as_bytes().iter().for_each(|x| new_qname.push(*x));
            new_record.set(&new_qname, None, &new_seq, &new_qual);
            new_record.set_pos(0);
            new_record.set_flags(0x4);
            new_record.set_mapq(0);
            new_record.set_tid(-1);
            new_record.remove_aux(b"MD").unwrap();
            new_record.remove_aux(b"AS").unwrap();
            new_record.remove_aux(b"NM").unwrap();
            new_record.remove_aux(b"XA").unwrap();
            new_record.remove_aux(b"SA").unwrap();
            new_record.remove_aux(b"tp").unwrap();

            let mut new_mm_list = Vec::new();
            let mut new_ml_list: Vec<u8> = Vec::new();
            let mut mm_cum = 0;
            
            'inner2: for (mm_idx, mm) in mm_list[mm_init..].iter().enumerate() {
    
                let mut mm_1 = mm.to_owned();
                               
                if mm_idx == 0 {
                    mm_1 =  mm_1 - least_c_counts;
                }
                
                mm_cum += mm_1 + 1;
                if mm_cum <= c_counts as u32{
                    // println!("{} {} {} {} {}", mm_idx, mm_init, mm_boundary,  mm_cum, c_counts);
                    mm_boundary += 1;
                    new_mm_list.push(mm_1);
                    new_ml_list.push(ml_list[mm_init + mm_idx] as u8);
                } else {
                    // println!("{} {} {} {} {}", mm_idx, mm_init, mm_boundary,  mm_cum, c_counts);
                    mm_cum -= mm_1 + 1;
                    // println!("{} {} {} {} {}", mm_idx, mm_init, mm_boundary,  mm_cum, c_counts);
                    least_c_counts = c_counts as u32 - mm_cum as u32;
                    break 'inner2;
                }
            }
            
            if new_mm_list.len() != 0 {

                new_record.remove_aux(b"MM").unwrap();
                new_record.remove_aux(b"ML").unwrap();
                // println!("{:?} {} {}", new_mm_list, new_mm_list.iter().sum::<u32>() + new_mm_list.len() as u32, new_mm_list.len());
                // push MM tag
                let mut mm_str = String::new();
                mm_str.push_str(mm_cap);
                mm_str.push(',');
                for (mm_idx2, mm) in new_mm_list.iter().enumerate() {
                    mm_str.push_str(&mm.to_string());
                    if mm_idx2 != new_mm_list.len() - 1 {
                        mm_str.push(',');
                    }
                
                }
                mm_str.push(';');
                new_record.push_aux(b"MM", Aux::String(&mm_str)).unwrap();

                // push ML tag 
                let ml_aux_array: AuxArray<u8> = AuxArray::from(&new_ml_list);
                new_record.push_aux(b"ML", Aux::ArrayU8(ml_aux_array)).unwrap();
            }
            
            writer.write(&new_record).unwrap();
            mm_init = mm_boundary;
            // println!("{}", new_qname.iter().map(|x| *x as char).collect::<String>());
            // println!("{}", String::from_utf8(new_seq.to_vec()).unwrap());
            let mut mod_pos = Vec::new();
            if let Ok(mods) = new_record.basemods_iter() {
                for (idx, mod_res) in mods.enumerate() {
                    if let Ok((position, m)) = mod_res {
                        // println!("{}: {:.?}", position, 
                        //             new_record.seq().as_bytes().to_vec()[position as usize..=(position + 1) as usize].iter().map(|x| *x as char).collect::<String>());
                        mod_pos.push(position);
                    }
                }
            };
        }

        old_read_id = qname_raw;
    }

}


// simulation pore-c reads through consensus sequence from a vcf file 
// input: vcf file, output: fasta file
pub fn simulate_porec(fasta: &String, vcf: &String, bed: &String, output: &String) {
    let fa = Fastx::new(fasta);
    let vcf_prefix = vcf.split(".").collect::<Vec<&str>>()[0];
    let seqs = fa.get_chrom_seqs().unwrap();
    let mut writer = common_writer(&output);
    let mut bed_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(bed).unwrap();

    let mut prev_read_idx = String::new();
    let mut read = Vec::new();
    for record in bed_reader.records() {
        let record = record.unwrap();
        let chrom = record.get(0).unwrap();
        let start = record.get(1).unwrap().parse::<usize>().unwrap();
        let end = record.get(2).unwrap().parse::<usize>().unwrap();
        let seq = seqs.get(chrom).unwrap().get(start..end).unwrap();
        let read_idx = record.get(3).unwrap();
        let region = format!("{}:{}-{}", chrom, start + 1, end);
        let region = region.parse().unwrap();
        let mut vcf_reader = vcf::indexed_reader::Builder::default().build_from_path(vcf).unwrap();
        let vcf_header = vcf_reader.read_header().unwrap();
        
        let vcf_records = vcf_reader.query(&vcf_header, &region).unwrap();
        let seq = seq.chars();
        let seq_vec = seq.collect::<Vec<char>>();

        let mut pos_vec = Vec::new();
        let mut ref_vec = Vec::new();
        let mut alt_vec = Vec::new();
        'inner: for result in vcf_records {
            let r = result.unwrap();
            let ref_len = r.reference_bases().len();
            let alt_len = r.alternate_bases().len();
            let ref_seq = r.reference_bases().to_vec();
            let alt_seq = r.alternate_bases().to_string().chars().collect::<Vec<char>>();
            let pos = usize::from(r.position()) as i64 - start as i64 - 1;
            if pos < 0 {
                continue 'inner;
            }
            pos_vec.push(pos);
            ref_vec.push(ref_len);
            alt_vec.push(alt_seq);
     
        }

        if pos_vec.len() == 0 {
            read.push(seq_vec.iter().collect::<String>());
        } else {
            let mut new_seq_vec = Vec::new();
            let mut flag = 0;
            'inner: for (i, c) in seq_vec.iter().enumerate() {
                while flag != 0 {
                    flag -= 1;
                    continue 'inner;
                }
                if pos_vec.contains(&(i as i64)) {
                    let idx = pos_vec.iter().position(|&x| x == i as i64).unwrap();
                    let ref_len = ref_vec[idx];
                    let alt_len = alt_vec[idx].len();
                    if ref_len == 1 {
                        if alt_len == 1 {
                            new_seq_vec.push(alt_vec[idx][0]);

                        } else if alt_len > 1 {
                            for s in alt_vec[idx].iter() {
                                new_seq_vec.push(*s);
                            }
                        } else {
                            continue 'inner;
                        }
                    } else if ref_len > 1 {
                        if alt_len == 1 {
                            new_seq_vec.push(alt_vec[idx][0]);

                        } else if alt_len > 1 {
                            for s in alt_vec[idx].iter() {
                                new_seq_vec.push(*s);
                            }
                        } else {
                            continue 'inner;
                        }
                        flag = ref_len - 1
                    } else {
                        if alt_len == 1 {
                            new_seq_vec.push(alt_vec[idx][0]);

                        } else if alt_len > 1 {
                            for s in alt_vec[idx].iter() {
                                new_seq_vec.push(*s);
                            }
                        } else {
                            continue 'inner;
                        }
                    }
                } else {
                    new_seq_vec.push(*c);
                }
            }
            let new_seq = new_seq_vec.iter().collect::<String>();
            read.push(new_seq);
        }
        if prev_read_idx != read_idx.to_string() && prev_read_idx != "" {
            writer.write_all(format!(">{}_{}\n{}\n", prev_read_idx, vcf_prefix, read.join("")).as_bytes()).unwrap();
           
            read.clear();
        }
     

        prev_read_idx = read_idx.to_string();

    }


}