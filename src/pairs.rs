use anyhow::Result as anyResult;
use log::LevelFilter;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::sync::{ Arc, Mutex };
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};

use crate::bed::Bed3;
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ChromSizeRecord, ContigPair};
use crate::mnd::MndRecord;
use crate::porec::PoreCRecord;
use crate::contacts::{ Contacts, ContactRecord };



#[derive(Debug, Clone)]
pub struct PairHeader {
    pub chromsizes: Vec<ChromSizeRecord>,
    pub header: Vec<String>,
}

impl PairHeader {
    pub fn new() -> PairHeader {
        let chromsizes: Vec<ChromSizeRecord> = Vec::new();
        let header: Vec<String> = Vec::new();

        PairHeader { chromsizes: chromsizes, header: header }
    }

    pub fn from_pairs(&mut self, name: &String) {
        log::set_max_level(LevelFilter::Off);
        let reader = common_reader(name);
        log::set_max_level(LevelFilter::Info);
        for line in reader.lines() {
            if let Ok(line) = line {
                if line.starts_with("#") {
                    if &line[0..10] == "#chromsize" {
                        
                        let s: String = line.replace("#chromsize: ", "");
                        let s: Vec<&str> = s.split(" ").collect();
                        
                        let size: u64 = s[1].clone().parse::<u64>().unwrap();
                        let chrom: String = s[0].to_string();
                        let c: ChromSizeRecord = ChromSizeRecord { chrom: chrom, size: size };

                        self.chromsizes.push(c);

                    } else if &line[0..8] == "#columns" {
                        
                        let s: String = line.replace("#columns: ", "");
                        let s: Vec<&str> = s.split(" ").collect();
                        
                        let s: Vec<String> = s.iter().map(|&x| x.to_owned()).collect();
                        self.header = s;

                    } else {
                        continue
                    }
                } else {
                    break
                }
            }
        }
    }

    pub fn from_chromsizes(&mut self, chromsizes: Vec<ChromSizeRecord>) {
        let pair_header: Vec<String> = vec!["readID".to_string(), 
                                    "chrom1".to_string(), 
                                    "pos1".to_string(), 
                                    "chrom2".to_string(), 
                                    "pos2".to_string(), 
                                    "strand1".to_string(),
                                    "strand2".to_string()
                                    ];
        self.chromsizes = chromsizes;
        self.header = pair_header;
    }

    

    pub fn to_string(&self) -> String {
        let mut result = String::new();
        
        result.push_str("## pairs format 1.0\n");
        result.push_str("#shape: upper triangle\n");

        for chromsize in &self.chromsizes {
            result.push_str(&format!("#chromsize: {} {}\n", 
                                chromsize.chrom, chromsize.size));
        }

        result.push_str(&format!("#columns: {}\n", &self.header.join(" ")));

        result
    }

}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PairRecord {
    pub readID: u64,
    pub chrom1: String,
    pub pos1: u64,
    pub chrom2: String,
    pub pos2: u64,
    pub strand1: char,
    pub strand2: char,
}

impl PairRecord {
    pub fn from_pore_c_pair(pair: Vec<PoreCRecord>, readID: u64) -> PairRecord {
        let (pair1, pair2) = (pair[0].clone(), pair[1].clone());
        
        let pos1: u64 = (pair1.target_end + pair1.target_start) / 2;
        let pos2: u64 = (pair2.target_end + pair2.target_start) / 2;

        PairRecord {
            readID: readID,
            chrom1: pair1.target,
            pos1: pos1, 
            chrom2: pair2.target,
            pos2: pos2,
            strand1: pair1.query_strand,
            strand2: pair2.query_strand,
        }
    }

    // pub fn is_in_regions(&self, interval_hash: &HashMap<String, Lapper<usize, u8>>) -> bool {
    //     let chrom1 = &self.chrom1;
    //     let chrom2 = &self.chrom2;
    //     let pos1 = self.pos1 as usize;
    //     let pos2 = self.pos2 as usize;

    //     let mut is_in_regions = false;
    //     if let Some(interval1) = interval_hash.get(&chrom1) {
    //         if let Some(interval2) = interval_hash.get(&chrom2) {
    //             let iv1 = interval1.count(pos1, pos1+1);
    //             let iv2 = interval2.count(pos2, pos2+1);

    //             if iv1 > 0 && iv2 > 0 {
    //                 is_in_regions = true;
    //             }
    //         }
    //     }
    //     is_in_regions
    // }

}

#[derive(Debug, Clone)]
pub struct Pairs {
    file: String,
    header: PairHeader,
}

impl BaseTable for Pairs {
    fn new(name: &String) -> Pairs {
        Pairs { file: name.clone(), header: PairHeader::new() }
    }

    fn file_name(&self) -> Cow<'_, str> {
        let path = Path::new(&self.file);
        path.file_name().expect("REASON").to_string_lossy()
    }

    fn prefix(&self) -> String {
        let binding = self.file_name().to_string();
        let file_path = Path::new(&binding);
        let file_prefix = file_path.file_stem().unwrap().to_str().unwrap();

        (*file_prefix).to_string()
    }
}

impl Pairs {
    pub fn parse(&mut self) -> anyResult<csv::Reader<Box<dyn BufRead + Send>>> {
        let input = common_reader(&self.file);
        self.header = PairHeader::new();
        self.header.from_pairs(&self.file);

        let rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .from_reader(input);

        Ok(rdr)
    }

    pub fn remove_by_contig_pairs(&mut self, contigs: HashSet<ContigPair>, output: &String) -> anyResult<()> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let ph = self.header.clone();
        let mut writer = common_writer(output);
        writer.write_all(ph.to_string().as_bytes()).unwrap();
        let mut wtr = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(writer);

        for result in rdr.deserialize() {
            if result.is_err() {
                println!("{:?}", result);
                continue
            }
            let record: PairRecord = result?;
            let cp = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
            if !contigs.contains(&cp) {
                wtr.serialize(record)?;
            }
        }

        Ok(())
    }

    pub fn to_mnd(&mut self, output: &String) -> anyResult<()>{
        let parse_result = self.parse();
        let mnd_defualt = MndRecord::default();

        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let wtr = common_writer(output);
        let mut wtr = csv::WriterBuilder::new()
                .has_headers(false)
                .delimiter(b' ')
                .from_writer(wtr);
        
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    
                    let record = PairRecord {
                        readID: idx as u64,
                        chrom1: record[1].to_string(),
                        pos1: record[2].parse::<u64>().unwrap(),
                        chrom2: record[3].to_string(),
                        pos2: record[4].parse::<u64>().unwrap(),
                        strand1: record[5].parse::<char>().unwrap(),
                        strand2: record[6].parse::<char>().unwrap(),
                    };

                    let strand1 = if record.strand1 == '+' { 0 } else { -1 };
                    let strand2 = if record.strand2 == '+' { 0 } else { -1 };
                    let mnd_record = MndRecord {
                        strand1: strand1,
                        chrom1: record.chrom1,
                        pos1: record.pos1,
                        frag1: mnd_defualt.frag1,
                        strand2: strand2,
                        chrom2: record.chrom2,
                        pos2: record.pos2,
                        frag2: mnd_defualt.frag2,
                        mapq1: mnd_defualt.mapq1,
                        cigar1: mnd_defualt.cigar1,
                        sequence1: mnd_defualt.sequence1,
                        mapq2: mnd_defualt.mapq2,
                        cigar2: mnd_defualt.cigar2,
                        sequence2: mnd_defualt.sequence2,
                        readname1: mnd_defualt.readname1,
                        readname2: mnd_defualt.readname2,
                    };
                    wtr.serialize(mnd_record)?;
                }
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
            }
        }

            // wtr.serialize(mnd_record)?;
        }

        wtr.flush()?;

        log::info!("Successful output mnd file `{}`", output);
        
        Ok(())
    }

    // convert pairs to pseudo bam file 
    pub fn to_bam(&mut self, output: &String) {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };
        let pair_header = self.header.clone();
        let mut bam_header = Header::new();
        let mut record = HeaderRecord::new(b"HD");
        record.push_tag(b"VN", &"1.6") ;
        record.push_tag(b"SO", &"coordinate");    
        bam_header.push_record(&record);

        for chromsize in &pair_header.chromsizes {
            let mut record = HeaderRecord::new(b"SQ");
            record.push_tag(b"SN", &chromsize.chrom);
            record.push_tag(b"LN", &chromsize.size);
            bam_header.push_record(&record);
        }
        let mut record = HeaderRecord::new(b"PG");
        record.push_tag(b"ID", &"pairs2bam");
        record.push_tag(b"PN", &"cphasing-rs");
        let version = env!("CARGO_PKG_VERSION");
        record.push_tag(b"VN", &version);
        record.push_tag(b"CL", &"cphasing-rs pairs2bam");
        bam_header.push_record(&record);
        
        let bam_header_view = HeaderView::from_header(&bam_header);

        let mut wtr = Writer::from_path(output, &bam_header, bam::Format::Bam).unwrap();

        for (id, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    // let flag1 = match record[5].parse::<char>().unwrap() {
                    //     '+' => 0,
                    //     '-' => 16,
                    //     _ => 0 
                    // };
                    let flag1 = 65;

                    // let flag2 = match record[6].parse::<char>().unwrap() {
                    //     '+' => 0,
                    //     '-' => 16,
                    //     _ => 0 
                    // };
                    let flag2 = 145;
                    {
                        let mut bam_record1 = Record::new();
                        let cigar = CigarString(vec![Cigar::Match(150)]);
                        let cigar: Option<&CigarString> = Some(&cigar);
                        bam_record1.set(record[0].as_bytes(), cigar, &[], &[]);
                        bam_record1.set_tid(bam_header_view.tid(record[1].as_bytes()).unwrap().try_into().unwrap());
                        bam_record1.set_pos(record[2].parse::<i64>().unwrap());
                        bam_record1.set_mapq(60);
                        bam_record1.set_mtid(bam_header_view.tid(record[3].as_bytes()).unwrap().try_into().unwrap());
                        bam_record1.set_mpos(record[4].parse::<i64>().unwrap());
                        bam_record1.set_insert_size(0);
                        // bam_record1.set_qname(record[0].as_bytes());
                        bam_record1.set_flags(flag1);
                        bam_record1.push_aux(b"RG", Aux::String("1")).unwrap();
                        bam_record1.push_aux(b"AS", Aux::U8(0)).unwrap();
                        bam_record1.push_aux(b"NM", Aux::U8(0)).unwrap();
                        bam_record1.push_aux(b"MQ", Aux::U8(60)).unwrap();
                    
                        wtr.write(&bam_record1).unwrap();
                    }
                    
                    {
                        let mut bam_record2 = Record::new();
                        let cigar = CigarString(vec![Cigar::Match(150)]);
                        let cigar: Option<&CigarString> = Some(&cigar);
                        bam_record2.set(record[0].as_bytes(), cigar, &[], &[]);
                        bam_record2.set_tid(bam_header_view.tid(record[3].as_bytes()).unwrap().try_into().unwrap());
                        bam_record2.set_pos(record[4].parse::<i64>().unwrap());
                        bam_record2.set_mapq(60);
                        bam_record2.set_mtid(bam_header_view.tid(record[1].as_bytes()).unwrap().try_into().unwrap());
                        bam_record2.set_mpos(record[2].parse::<i64>().unwrap());
                        bam_record2.set_insert_size(0);
                        // bam_record2.set_qname(record[0].as_bytes());
                        bam_record2.set_flags(flag2);
                        bam_record2.push_aux(b"RG", Aux::String("1")).unwrap();
                        bam_record2.push_aux(b"AS", Aux::U8(0)).unwrap();
                        bam_record2.push_aux(b"NM", Aux::U8(0)).unwrap();
                        bam_record2.push_aux(b"MQ", Aux::U8(60)).unwrap();

                        wtr.write(&bam_record2).unwrap();
                    }

                }
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
            }
                
            }
            log::info!("Successful output bam file `{}`", output);

    }  
    // split contig to nparts by split_num and calculate the number of contacts
    pub fn to_split_contacts(&mut self, min_contacts: u32, split_num: u32) -> anyResult<Contacts> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let contigsizes = self.header.chromsizes.clone();
        let mut contacts = Contacts::new(&format!("{}.pixels", self.prefix()).to_string());       

        let split_contigsizes: HashMap<String, u64> = contigsizes
            .par_iter()
            .map(|x| (x.chrom.clone(), x.size / split_num as u64))
            .collect();

        let mut contact_hash: HashMap<ContigPair, u32>  = HashMap::new();
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
               
                    let split_index1 = record[2].parse::<u64>().unwrap() / split_contigsizes.get(&record[1]).unwrap();
                    let split_index2 = record[4].parse::<u64>().unwrap() / split_contigsizes.get(&record[3]).unwrap();
                    let mut cp = ContigPair::new(
                                    format!("{}_{}", record[1].to_string(), split_index1), 
                                    format!("{}_{}", record[3].to_string(), split_index2));
                    cp.order();

                    *contact_hash.entry(cp).or_insert(0) += 1;
                },
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
               
            }
            
        }
        
       // convert count to f64
       let contact_hash: HashMap<ContigPair, f64> = contact_hash
            .into_par_iter()
            .map(|(cp, count)| {
                let count = count as f64;
                (cp, count)
            })
            .collect();

        let mut contact_records: Vec<ContactRecord> = contact_hash.into_par_iter(
                ).filter_map(|(cp, count)| {
                if count >= (min_contacts as f64){
                let record = ContactRecord {
                chrom1: cp.Contig1,
                chrom2: cp.Contig2,
                count: count,
            };
            Some(record)
            } else {
            None
            }
        }).collect();

        contact_records.retain(|x| x.is_some());
        contacts.records = contact_records;
        Ok(contacts)
    }

    pub fn to_contacts(&mut self, min_contacts: u32) -> anyResult<Contacts> {
        use hashbrown::HashSet;
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let mut contacts = Contacts::new(&format!("{}.pixels", self.prefix()).to_string());
        
        let mut contact_hash: HashMap<ContigPair, f64>  = HashMap::new();
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                //    let record =  PairRecord {
                //         readID: idx as u64,
                //         chrom1: record[1].to_string(),
                //         pos1: record[2].parse::<u64>().unwrap(),
                //         chrom2: record[3].to_string(),
                //         pos2: record[4].parse::<u64>().unwrap(),
                //         strand1: record[5].parse::<char>().unwrap(),
                //         strand2: record[6].parse::<char>().unwrap(),
                //     };
                    // let mut cp = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
                    let mut cp = ContigPair::new(record[1].to_string(), record[3].to_string());
                    cp.order();
                    *contact_hash.entry(cp).or_insert(1.0) += 1.0;
                },
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
            }
            
        }
        
    
        let mut contact_records: Vec<ContactRecord> = contact_hash.par_iter(
            ).map(|(cp, count)| {
                if count >= &(min_contacts as f64){
                    let record = ContactRecord {
                        chrom1: cp.Contig1.to_owned(),
                        chrom2: cp.Contig2.to_owned(),
                        count: *count,
                    };
                    record
                } else {
                    let record = ContactRecord {
                        chrom1: "".to_string(),
                        chrom2: "".to_string(),
                        count: 0.0,
                    };
                    record
                }
                
            }).collect();
        
        contact_records.retain(|x| x.is_some());
        
        contacts.records = contact_records;
        Ok(contacts)

    }   

    pub fn to_clm(&mut self, min_contacts: u32, output: &String) {
        use hashbrown::HashMap;
        let parse_result = self.parse();

        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        // get output prefix
        let output_prefix = Path::new(&output).file_stem().unwrap().to_str().unwrap();

        let pair_header = self.header.clone();
        let chromsizes: Vec<ChromSizeRecord> = pair_header.chromsizes.clone();
        // to hashmap
        let chromsizes: HashMap<String, u32> = chromsizes
            .iter()
            .map(|x| (x.chrom.clone(), x.size.try_into().unwrap()))
            .collect();
        match chromsizes.len() {
            0 => log::warn!("chromsizes is empty !!!, please check it."),
            _ => {}
        }
        let mut contacts = Contacts::new(&format!("{}.pixels", self.prefix()).to_string());
        let split_contigsizes: HashMap<&String, u32> = chromsizes
            .iter()
            .map(|(k, v)| (k, *v / 2))
            .collect();

        let mut contact_hash: HashMap<ContigPair, u32>  = HashMap::new();
        let mut data: HashMap<ContigPair, Vec<Vec<u32>>> = HashMap::new();
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    let contig1 = record[1].to_string();
                    let contig2 = record[3].to_string();
                    let pos1 = record[2].parse::<u32>().unwrap_or_default();
                    let pos2 = record[4].parse::<u32>().unwrap_or_default();
                    
                    let (contig1, contig2, pos1, pos2) = if contig1 > contig2 {
                        (contig2, contig1, pos2, pos1)
                    } else {
                        (contig1, contig2, pos1, pos2)
                    };
                    let split_index1 = pos1 / split_contigsizes.get(&contig1).unwrap();
                    let split_index2 = pos2 / split_contigsizes.get(&contig2).unwrap();
                    let cp = ContigPair::new(
                        format!("{}_{}", contig1, split_index1), 
                        format!("{}_{}", contig2, split_index2));

                    let mut contig_pair = ContigPair::new(contig1, contig2);
                    contig_pair.order();
                    *contact_hash.entry(cp).or_insert(0) += 1;

                    data.entry(contig_pair)
                        .or_insert_with(Vec::new)
                        .push(vec![pos1, pos2]);
                },
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
            }
        }

        // // convert count to f64
        // let contact_hash: HashMap<ContigPair, f64> = contact_hash
        //                                             .into_par_iter()
        //                                             .map(|(cp, count)| {
        //                                                 let count = count as f64;
        //                                                 (cp, count)
        //                                             })
        //                                             .collect();

        let mut contact_records: Vec<ContactRecord> = contact_hash.into_par_iter(
            ).filter_map(|(cp, count)| {
                if count >= min_contacts{
                    let record = ContactRecord {
                        chrom1: cp.Contig1,
                        chrom2: cp.Contig2,
                        count: count as f64,
                    };
                    Some(record)
                } else {
                    None
                }
         }).collect();
    
        contact_records.retain(|x| x.is_some());
        contacts.records = contact_records;

        contacts.write(&format!("{}.split.contacts", output_prefix.to_string()));
        log::info!("Successful output split contacts file `{}`", &format!("{}.split.contacts", output_prefix.to_string()));


        log::info!("Calculating the distance between contigs");
        let result = data.par_iter(
            ).filter_map(
                |(cp, vec)| {
                    if vec.len() >= min_contacts as usize {
                        Some((cp, vec))
                    } else {
                        None
                    }
                }
            ).collect::<HashMap<_, _>>();

        let result = result.into_par_iter().map(|(cp, vec)|{
                let length1 = chromsizes.get(&cp.Contig1).unwrap();
                let length2 = chromsizes.get(&cp.Contig2).unwrap();
                
                let res = vec.par_iter().map(|x| {
                    let pos1 = x[0];
                    let pos2 = x[1];
                   
                    vec![length1 - pos1 + pos2, // ctg1+ ctg2+
                        length1 - pos1 + length2 - pos2, // ctg1+ ctg2-
                        pos1 + pos2, // ctg1- ctg2+
                        pos1 + length2 - pos2 // ctg1- ctg2-
                        ]
                }).collect::<Vec<_>>();

                //zip 
                let zipped: Vec<Vec<u32>> = res[0].iter().enumerate().map(|(i, _)| {
                    res.iter().map(|x| x[i]).collect::<Vec<_>>()
                }).collect::<Vec<_>>();
              
                (cp, zipped)
                
            }).collect::<HashMap<_, Vec<Vec<u32>>>>();
    
        
        let mut wtr = common_writer(output);
        
        for (cp, res) in result {
            let count = res[0].len();
            
            if count < min_contacts as usize {
                continue
            } 
            let mut buffer = Vec::new();
            for (i, res1) in res.iter().enumerate() {
                let res1 = res1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ");
                match i {
                    0 => {
                        write!(buffer, "{}+ {}+\t{}\t{}\n", cp.Contig1, cp.Contig2, count, res1).unwrap();
                    },
                    1 => {
                        write!(buffer, "{}+ {}-\t{}\t{}\n", cp.Contig1, cp.Contig2, count, res1).unwrap();
                    },
                    2 => {
                        write!(buffer, "{}- {}+\t{}\t{}\n", cp.Contig1, cp.Contig2, count, res1).unwrap();
                    },
                    3 => {
                        write!(buffer, "{}- {}-\t{}\t{}\n", cp.Contig1, cp.Contig2, count, res1).unwrap();
                    },
                    _ => panic!("Error: Invalid index"),
                };

            }
            wtr.write_all(&buffer).unwrap();

        }
        log::info!("Successful output clm file `{}`", output);
    }

    pub fn intersect(&mut self, hcr_bed: &String, invert: bool, output: &String) {
        type IvU8 = Interval<usize, u8>;
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        let mut wtr = common_writer(output);
        self.header = PairHeader::new();
        self.header.from_pairs(&self.file);
        wtr.write_all(self.header.to_string().as_bytes()).unwrap();
        
        for (idx, record) in self.parse().unwrap().records().enumerate() {
            let record = record.unwrap();
            let chrom1 = record[1].to_string();
            let pos1 = record[2].parse::<usize>().unwrap();
            let chrom2 = record[3].to_string();
            let pos2 = record[4].parse::<usize>().unwrap();

            let is_in_regions = if let Some(interval1) = interval_hash.get(&chrom1) {
                if let Some(interval2) = interval_hash.get(&chrom2) {
                  
                    if interval1.count(pos1, pos1+1) > 0 && interval2.count(pos2, pos2+1) > 0 {
                        true
                    } else {
                        false
                    }
                } else {
                    false
                }
            } else {
                false
            };

            if invert {
                match is_in_regions {
                    true => continue,
                    false => {
                        wtr.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                                                idx, chrom1, pos1, chrom2, pos2, 
                                                record[5].parse::<char>().unwrap(), 
                                                record[6].parse::<char>().unwrap()
                                            ).as_bytes()).unwrap();
                    }
                }
            } else {

                match is_in_regions {
                    true => {
                        wtr.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                                                idx, chrom1, pos1, chrom2, pos2, 
                                                record[5].parse::<char>().unwrap(), 
                                                record[6].parse::<char>().unwrap()
                                            ).as_bytes()).unwrap();
                    },
                    false => continue
                }
            }

        }

        log::info!("Successful output new pairs file into `{}`", output);
    }
}

