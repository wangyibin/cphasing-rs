#![allow(unused)]
use anyhow::Result as anyResult;
use bytecount;
use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use csv::StringRecord;
use log::LevelFilter;
use rand::prelude::*;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use std::borrow::Cow;
use std::collections::{ BTreeMap, HashMap, HashSet };
use std::hash::BuildHasherDefault;
use std::error::Error;
use std::fs::File;
use std::hash::{ Hash, Hasher };
use std::path::Path;
use std::thread;
use std::sync::{ Arc, Mutex , mpsc};
use std::io::{ Write, BufReader, BufRead, Read };
use serde::{ Deserialize, Serialize};
use smallvec::{ smallvec, SmallVec };
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use twox_hash::XxHash64;
// use tokio::fs::File;
// use tokio::io::AsyncWriteExt;

use crate::bed::{ Bed3, Bed4 };
use crate::core::{ common_reader, common_writer };
use crate::core::{ 
    BaseTable, 
    ChromSize,
    ChromSizeRecord, 
    ContigPair, ContigPair2,
    binify
};
use crate::mnd::MndRecord;
use crate::porec::{ PoreCRecord, PoreCRecordPlus };
use crate::contacts::{ Contacts, ContactRecord };

type SmallIntVec = SmallVec<[u32; 2]>;
type SmallIntVec4 = SmallVec<[u32; 4]>;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SplitIdx {
    pub contig_idx: u32, 
    pub split_idx: u8,
}

impl Hash for SplitIdx {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.contig_idx.hash(state);
        self.split_idx.hash(state);
    }
}



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
        
        log::info!("Load contigsizes from pairs file.")
    }

    pub fn from_chromsizes(&mut self, chromsizes: Vec<ChromSizeRecord>) {
        let pair_header: Vec<String> = vec!["readID".to_string(), 
                                    "chrom1".to_string(), 
                                    "pos1".to_string(), 
                                    "chrom2".to_string(), 
                                    "pos2".to_string(), 
                                    "strand1".to_string(),
                                    "strand2".to_string(),
                                    "mapq".to_string()
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
    pub mapq: u8,
}

impl PairRecord {
    pub fn from_pore_c_pair(pair: Vec<&PoreCRecord>, readID: u64) -> PairRecord {
        // let (pair1, pair2) = (pair[0].clone(), pair[1].clone());
        
        let pos1: u64 = (pair[0].target_end + pair[0].target_start) / 2;
        let pos2: u64 = (pair[1].target_end + pair[1].target_start) / 2;
        // let mapq: f64 = (pair1.mapq as f64 * pair2.mapq as f64).sqrt();
        // let mapq: u8 = mapq.round() as u8;
        let mapq = std::cmp::min(pair[0].mapq, pair[1].mapq);

        PairRecord {
            readID: readID,
            chrom1: pair[0].target.clone(),
            pos1: pos1, 
            chrom2: pair[1].target.clone(),
            pos2: pos2,
            strand1: pair[0].query_strand,
            strand2: pair[1].query_strand,
            mapq: mapq,
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

    pub fn to_mnd(&mut self, min_quality: u8, output: &String) -> anyResult<()>{
        let parse_result = self.parse();
        let mnd_defualt = MndRecord::default();

        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };
        let filter_mapq = min_quality > 0;
        let (sender, receiver) = bounded::<Vec<(usize, StringRecord)>>(100);
        let wtr = Arc::new(Mutex::new(common_writer(output)));
        let mut handles = vec![];
        for _ in 0..4 {
            let receiver: Receiver<_> = receiver.clone();
            let mut wtr = Arc::clone(&wtr);
            handles.push(thread::spawn(move || {
                
                while let Ok(records) = receiver.recv() {
                    let mut data = vec![];
                    let filter_needed = filter_mapq &&  records.get(0).map_or(false, |(_, record)| record.len() >= 8);

                    for (idx, record) in records.iter().filter(|(idx, record)| {
                        !filter_needed || record.get(7)
                            .and_then(|mapq_str| mapq_str.parse::<u8>().ok())
                            .map_or(true, |mapq| mapq > min_quality)
                    }) {
                        
                        let strand1 = match record[5].parse::<char>() {
                            Ok('+') => 0,
                            _ => -1,
                        };
                        let strand2 = match record[6].parse::<char>() {
                            Ok('+') => 0,
                            _ => -1,
                        };
    
                        let pos1 = record[2].parse::<u64>().unwrap_or_default();
                        let pos2 = record[4].parse::<u64>().unwrap_or_default();
                        
                        let line = format!("{} {} {} {} \
                                            {} {} {} {} \
                                            {} {} {} {} \
                                            {} {} {} {}", 
                                            strand1, record[1].to_string(), pos1, mnd_defualt.frag1,
                                            strand2, record[3].to_string(), pos2, mnd_defualt.frag2,
                                            mnd_defualt.mapq1, mnd_defualt.cigar1, mnd_defualt.sequence1,
                                            mnd_defualt.mapq2, mnd_defualt.cigar2, mnd_defualt.sequence2,
                                            mnd_defualt.readname1, mnd_defualt.readname2);

                        data.push(line);
                    }

                    let mut wtr = wtr.lock().unwrap();
                    for line in data.drain(..) {
                        writeln!(wtr, "{}", line).unwrap();
                    }
                    
                }

               
            }));
        }

        let mut batch = Vec::with_capacity(1000);
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    batch.push((idx, record));
                    if batch.len() >= 1000 {
                        // sender.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).unwrap();
                        sender.send(std::mem::take(&mut batch)).unwrap();
                    }
                },
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
            }
        }

        if !batch.is_empty() {
            sender.send(batch).unwrap();
        }

        drop(sender);
        for handle in handles {
            handle.join().unwrap();
        }



        // let mut wtr = common_writer(output);
        
        // for (idx, record) in rdr.records().enumerate() {
        //     match record {
        //         Ok(record) => {
                   
        //             if  filter_mapq {
        //                 if record.len() >= 8{

        //                     let mapq = match record.get(7).unwrap_or("60").parse::<u8>() {
        //                         Ok(mapq) => mapq,
        //                         Err(_) => 60,
        //                     };
        //                     if mapq < min_quality {
        //                         continue
        //                     }
        //                 }

        //             }
                    
        //             let strand1 = match record[5].parse::<char>() {
        //                 Ok('+') => 0,
        //                 _ => -1,
        //             };
        //             let strand2 = match record[6].parse::<char>() {
        //                 Ok('+') => 0,
        //                 _ => -1,
        //             };

        //             let pos1 = record[2].parse::<u64>().unwrap_or_default();
        //             let pos2 = record[4].parse::<u64>().unwrap_or_default();
                    
        //             let line = format!("{} {} {} {} \
        //                                 {} {} {} {} \
        //                                 {} {} {} {} \
        //                                 {} {} {} {}", 
        //                                 strand1, record[1].to_string(), pos1, mnd_defualt.frag1,
        //                                 strand2, record[3].to_string(), pos2, mnd_defualt.frag2,
        //                                 mnd_defualt.mapq1, mnd_defualt.cigar1, mnd_defualt.sequence1,
        //                                 mnd_defualt.mapq2, mnd_defualt.cigar2, mnd_defualt.sequence2,
        //                                 mnd_defualt.readname1, mnd_defualt.readname2);
                    
        //             writeln!(wtr, "{}\n", line)?;
                                    
        //         }
        //         Err(e) => {
        //             eprintln!("{:?}", e);
        //             continue
        //     }
        // }
        //     wtr.flush()?;

        // }

        // wtr.flush()?;

        log::info!("Successful output mnd file `{}`", output);
        
        Ok(())
    }

    // convert pairs to pseudo bam file 
    pub fn to_bam(&mut self, min_quality: u8, output: &String, threads: usize ) {
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
        wtr.set_threads(threads);
        let filter_mapq = min_quality > 0;
        for (id, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    // let flag1 = match record[5].parse::<char>().unwrap() {
                    //     '+' => 0,
                    //     '-' => 16,
                    //     _ => 0 

                    // };

                    if filter_mapq {
                        if record.len() >= 8 {
                            let mapq = record.get(7).unwrap_or(&"60").parse::<u8>().unwrap_or_default();
                            if mapq <= min_quality {
                                continue
                            }
                        }
                    }

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
    pub fn to_split_contacts(&mut self, min_contacts: u32, split_num: u32, min_quality: u8) -> anyResult<Contacts> {
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
        let filter_mapq = min_quality > 0;
        let mut contact_hash: HashMap<ContigPair, u32>  = HashMap::new();
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    if filter_mapq {
                        if record.len() >= 8 {
                            let mapq = record.get(7).unwrap_or(&"60").parse::<u8>().unwrap_or_default();
                            if mapq <= min_quality {
                                continue
                            }
                        }
                    }
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

    pub fn to_contacts(&mut self, min_contacts: u32, min_quality: u8) -> anyResult<Contacts> {
        use hashbrown::HashSet;
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let pair_header = self.header.clone();
        let chromsizes: Vec<ChromSizeRecord> = pair_header.chromsizes.clone();
        let chromsizes: HashMap<String, u32> = chromsizes
            .iter()
            .map(|x| {
                let size = if x.size >  u32::MAX as u64 {
                    u32::MAX
                } else {
                    x.size.try_into().unwrap()
                };
                (x.chrom.clone(), size)
            })
            .collect();

        let contig_idx: HashMap<String, usize> = chromsizes.keys()
            .enumerate()
            .map(|(index, key)| (key.clone(), index))
            .collect();

        let idx_sizes: HashMap<usize, u32> = chromsizes.iter()
            .filter_map(|(key, size)| contig_idx.get(key).map(|&index| (index, *size)))
            .collect();

        let idx_contig: HashMap<usize, String> = contig_idx.clone()
            .into_iter()
            .map(|(key, index)| (index, key))
            .collect();
        match chromsizes.len() {
            0 => log::warn!("chromsizes is empty !!!, please check it."),
            _ => {}
        }

        let filter_mapq = min_quality > 0;
        let mut contacts = Contacts::new(&format!("{}.pixels", self.prefix()).to_string());
        
        let mut contact_hash: HashMap<(&usize, &usize), f64>  = HashMap::new();
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    
                    if filter_mapq {
                        if record.len() >= 8 {
                            let mapq = record.get(7).unwrap_or(&"60").parse::<u8>().unwrap_or_default();
                            if mapq <= min_quality {
                                continue
                            }
                        }
                    }

                    let contig1_str = record[1].to_string();
                    let contig2_str = record[3].to_string();
                    
                    let (contig1, contig2) = if contig1_str > contig2_str {
                        (contig2_str, contig1_str)
                    } else {
                        (contig1_str, contig2_str)
                    };

                    if let (Some(idx1), Some(idx2)) = (contig_idx.get(&contig1), contig_idx.get(&contig2)) {
                        *contact_hash.entry((idx1, idx2)).or_insert(1.0) += 1.0;
                    }

                    
                    // let mut cp = ContigPair::new(record[1].to_string(), record[3].to_string());
                    // cp.order();
                    // *contact_hash.entry(cp).or_insert(1.0) += 1.0;
                },
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
            }
            
        }
        
    
        let mut contact_records: Vec<ContactRecord> = contact_hash.par_iter(
            ).map(|(cp_idx, count)| {
                let contig1 = idx_contig.get(cp_idx.0).unwrap();
                let contig2 = idx_contig.get(cp_idx.1).unwrap();
                if count >= &(min_contacts as f64){
                    let record = ContactRecord {
                        chrom1: contig1.to_owned(),
                        chrom2: contig2.to_owned(),
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

    pub fn to_clm(&mut self, min_contacts: u32, 
                // binsize: u32,
                min_quality: u8,  output: &String,
                output_split_contacts: bool, 
                threads: usize,
                low_memory: bool) {
        use hashbrown::HashMap;
        let parse_result = self.parse();

        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        // get output prefix
        let output_prefix =  if output.ends_with(".gz") {
            Path::new(output.trim_end_matches(".gz")).file_stem().unwrap().to_str().unwrap()
        } else {
            Path::new(&output).file_stem().unwrap().to_str().unwrap()
        };

        let pair_header = self.header.clone();
        let chromsizes: &Vec<ChromSizeRecord> = &pair_header.chromsizes;
        let chromsizes: HashMap<String, u32> = chromsizes
            .iter()
            .map(|x| (x.chrom.clone(), x.size.try_into().unwrap()))
            .collect();

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = chromsizes.keys()
            .enumerate()
            .map(|(index, key)| (key.clone(), index as u32))
            .collect();

        let idx_sizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = chromsizes.par_iter()
            .filter_map(|(key, size)| contig_idx.get(key).map(|&index| (index, *size)))
            .collect();

        let idx_contig: HashMap<u32, &String, BuildHasherDefault<XxHash64>> = contig_idx.par_iter()
            .map(|(key, index)| (*index, key))
            .collect();

        if chromsizes.is_empty(){
            log::warn!("chromsizes is empty !!!, please check it.");
        }
        
        let mut contacts = Contacts::new(&format!("{}.pixels", self.prefix()).to_string());
        let split_contigsizes: HashMap<u32, u32> = idx_sizes
            .par_iter()
            .map(|(k, v)| (*k, *v / 2))
            .collect(); 

        log::info!("Parsing pairs file ...");
        let mut contact_hash: HashMap<(SplitIdx, SplitIdx), u32>  = HashMap::new();
        let filter_mapq = min_quality > 0;
    
        
        let mut data: Arc<Mutex<HashMap<(u32, u32), Vec<SmallIntVec>>>> = Arc::new(Mutex::new(HashMap::new()));
        let (sender, receiver) = bounded::<Vec<(usize, StringRecord)>>(1000);

        let mut handles = vec![];
        for _ in 0..4 {
            let receiver: Receiver<_> = receiver.clone();
            
            let data = Arc::clone(&data);
            let contig_idx = contig_idx.clone();

            handles.push(thread::spawn(move || {
                while let Ok((records)) = receiver.recv() {

                    let mut updates: HashMap<(u32, u32), Vec<SmallIntVec>> = HashMap::with_capacity(1000);
                    let filter_needed = filter_mapq &&  records.get(0).map_or(false, |(_, record)| record.len() >= 8);

                    for (idx, record) in records.iter().filter(|(idx, record)| {
                        !filter_needed || record.get(7)
                            .and_then(|mapq_str| mapq_str.parse::<u8>().ok())
                            .map_or(true, |mapq| mapq > min_quality)
                    }) {
                        let pos1 = record[2].parse::<u32>().unwrap_or_default();
                        let pos2 = record[4].parse::<u32>().unwrap_or_default();
                        
                        let (contig1, contig2, pos1, pos2) = if record[1] > record[3] {
                            (&record[3], &record[1],  pos2, pos1)
                        } else {
                            (&record[1], &record[3], pos1, pos2)
                        };
                        if let (Some(&idx1), Some(&idx2)) = (contig_idx.get(contig1), contig_idx.get(contig2)) {
                            updates.entry((idx1, idx2))
                                .and_modify(|e| e.push(smallvec![pos1, pos2]))
                                .or_insert_with(|| vec![smallvec![pos1, pos2]]);
                        }
                    }
                    {
                        let mut data = data.lock().unwrap();
                        for (key, mut local_values) in updates {
                            data.entry(key).or_insert_with(Vec::new).append(&mut local_values);
                        }
                    }
                    
                }
              
            }));
        }

        let mut batch = Vec::with_capacity(1000);
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    batch.push((idx, record));
                    if batch.len() >= 1000 {
                        // sender.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).unwrap();  
                        sender.send(std::mem::take(&mut batch)).unwrap();

                    }
                },
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
            }
            
        }

        if !batch.is_empty() {
            sender.send(batch).unwrap();
        }

        drop(sender);
        for handle in handles {
            handle.join().unwrap();  
        }
        
        let mut data = Arc::try_unwrap(data).unwrap().into_inner().unwrap();
        log::info!("Collected data from pairs file.");

        // let mut data: HashMap<(u32, u32), Vec<SmallIntVec>> = HashMap::new();
        // rdr.records().enumerate().for_each(|(idx, record_result)| {
        //     match record_result {
        //         Ok(record) => {
        //             if filter_mapq && record.len() >= 8 {
        //                 let mapq = record.get(7).unwrap_or(&"60").parse::<u8>().unwrap_or_default();
        //                 if mapq <= min_quality {
        //                     return; 
        //                 }
        //             }

        //             let pos1 = record[2].parse::<u32>().unwrap_or_default();
        //             let pos2 = record[4].parse::<u32>().unwrap_or_default();

        //             let (contig1, contig2, pos1, pos2) = if record[1] > record[3] {
        //                 (&record[3], &record[1],  pos2, pos1)
        //             } else {
        //                 (&record[1], &record[3], pos1, pos2)
        //             };

        //             if let (Some(&idx1), Some(&idx2)) = (contig_idx.get(contig1), contig_idx.get(contig2)) {
        //                 data.entry((idx1, idx2))
        //                     .and_modify(|e| e.push(smallvec![pos1, pos2]))
        //                     .or_insert_with(|| vec![smallvec![pos1, pos2]]);
        //             }
        //         },
        //         Err(e) => {
        //             eprintln!("{:?}", e);
        //         }
        //     }
        // });

        
        if output_split_contacts{
            let writer = common_writer(format!("{}.split.contacts", output_prefix.to_string()).as_str());
            let mut writer = Arc::new(Mutex::new(writer));
            data.par_iter().for_each(|(cp, vec) | {
                if vec.len() < min_contacts as usize {
                    return;
                }
                let length1 = idx_sizes.get(&cp.0).unwrap();
                let length2 = idx_sizes.get(&cp.1).unwrap();
                let contig1 = idx_contig.get(&cp.0).unwrap();
                let contig2 = idx_contig.get(&cp.1).unwrap();
                let res = vec.iter().map(
                    |x| {
                        let pos1 = x[0];
                        let pos2 = x[1];
                        let split_index1 = (pos1 / (length1 / 2)) as u8;
                        let split_index2 = (pos2 / (length2 / 2)) as u8;
                    
                        (split_index1, split_index2)
                    }
                ).collect::<Vec<_>>();

                let mut contact_hash = HashMap::with_capacity(4);
                res.iter().for_each(|(split_idx1, split_idx2)| {
                    *contact_hash.entry((split_idx1, split_idx2)).or_insert(0) += 1;
                });
                
                let mut buffer = Vec::with_capacity(4);
                contact_hash.iter().for_each(|(cp, count)| {
                    if count >= &min_contacts {
                        buffer.push(format!("{}_{}\t{}_{}\t{}\n", contig1, cp.0, contig2, cp.1, count));
                    }
                    
                });
                let buffer = buffer.join("");
                let mut writer = writer.lock().unwrap();
                writer.write_all(buffer.as_bytes()).unwrap();

            });

            log::info!("Successful output split contacts file `{}`", 
                                &format!("{}.split.contacts", output_prefix.to_string()));

            drop(split_contigsizes);
        }

        log::info!("Calculating the distance between contigs");
        
        
        match low_memory {
            true => {
                let mut wtr = common_writer(output);
                let wtr = Arc::new(Mutex::new(wtr));
                data.par_iter().for_each(|(cp, vec)| {
                    if vec.len() < min_contacts as usize {
                        return;
                    }
                    let length1 = idx_sizes.get(&cp.0).unwrap();
                    let length2 = idx_sizes.get(&cp.1).unwrap();
                
                    let res = vec.iter().map(
                        |x| {
                            let pos1 = x[0];
                            let pos2 = x[1];
                            smallvec![
                                length1.wrapping_sub(pos1).wrapping_add(pos2),                      // ctg1+ ctg2+
                                length1.wrapping_sub(pos1).wrapping_add(*length2).wrapping_sub(pos2), // ctg1+ ctg2-
                                pos1.wrapping_add(pos2),                                             // ctg1- ctg2+
                                pos1.wrapping_add(*length2).wrapping_sub(pos2),                      // ctg1- ctg2-
                            ]
                        }
                    ).collect::<Vec<_>>();
                    
                    let zipped: Vec<Vec<u32>> = res[0].iter().enumerate().map(|(i, _)| {
                        res.iter().map(|x: &SmallIntVec4| x[i]).collect::<Vec<_>>()
                    }).collect::<Vec<_>>();

                    let count = zipped[0].len();
                    if count < min_contacts as usize {
                        return
                    }
                    let contig1 = idx_contig.get(&cp.0).unwrap();
                    let contig2 = idx_contig.get(&cp.1).unwrap();

                    let mut buffer = Vec::with_capacity(4);
                    for (i, res1) in zipped.iter().enumerate() {
                        let res1 = res1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ");
                        let line = match i {
                            0 => format!("{}+ {}+\t{}\t{}\n", contig1, contig2, count, res1),
                            1 => format!("{}+ {}-\t{}\t{}\n", contig1, contig2, count, res1),
                            2 => format!("{}- {}+\t{}\t{}\n", contig1, contig2, count, res1),
                            3 => format!("{}- {}-\t{}\t{}\n", contig1, contig2, count, res1),
                            _ => panic!("Error: Invalidindex"),
                        };
                    
                        buffer.extend_from_slice(line.as_bytes());
                    }
                    {
                        let mut wtr = wtr.lock().unwrap();
                        wtr.write_all(&buffer).unwrap();
                    }

                });
            },
                false => {

                    let result = data.par_iter_mut(
                        ).filter_map(|(cp, vec)| {
                            if vec.len() >= min_contacts as usize {
                                Some((cp, vec))
                            } else {
                                None
                            }
                        }
                    ).map_init(|| Vec::new(), |res, (cp, vec)|{
                        let length1 = idx_sizes.get(&cp.0).unwrap();
                        let length2 = idx_sizes.get(&cp.1).unwrap();
                        
                        res.reserve(vec.len());
                        for x in vec {
                            let pos1 = x[0];
                            let pos2 = x[1];
                            
                            res.push([
                                length1.wrapping_sub(pos1).wrapping_add(pos2),                      // ctg1+ ctg2+
                                length1.wrapping_sub(pos1).wrapping_add(*length2).wrapping_sub(pos2), // ctg1+ ctg2-
                                pos1.wrapping_add(pos2),                                             // ctg1- ctg2+
                                pos1.wrapping_add(*length2).wrapping_sub(pos2),                      // ctg1- ctg2-
                            ]);
                        }
                    
                        let zipped: Vec<Vec<u32>> = (0..res[0].len())
                            .map(|i| res.iter().map(|x| x[i]).collect::<Vec<_>>())
                            .collect();

                        res.clear(); 
                        (*cp, zipped)
                        
                    }).collect::<HashMap<_, Vec<Vec<u32>>>>();
            

                    log::info!("Starting to output clm.");
                    
                    let wtr = Mutex::new(common_writer(output));

                    result.par_iter()
                        .filter(|(_, res)| res[0].len() >= min_contacts as usize)
                        .for_each(|(cp_idx, res)| {
                            let contig1 = idx_contig.get(&cp_idx.0).unwrap();
                            let contig2 = idx_contig.get(&cp_idx.1).unwrap();
                            let cp = ContigPair2 { Contig1: contig1, Contig2: contig2 };
                            let res_len = res[0].len();
                            let mut buffer = String::new();
                            for (i, res1) in res.iter().enumerate() {
                                let res1 = res1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ");

                                use std::fmt::Write;
                                write!(
                                    buffer,
                                    "{}\t{}\t{}\n",
                                        match i {
                                            0 => format!("{}+ {}+", cp.Contig1, cp.Contig2),
                                            1 => format!("{}+ {}-", cp.Contig1, cp.Contig2),
                                            2 => format!("{}- {}+", cp.Contig1, cp.Contig2),
                                            3 => format!("{}- {}-", cp.Contig1, cp.Contig2),
                                            _ => panic!("Error: Invalid index"),
                                        },
                                        res_len,
                                        res1 
                            
                                ).unwrap();
                            }

                            let mut wtr = wtr.lock().unwrap();
                            wtr.write_all(buffer.as_bytes()).unwrap();
                    });
                
            }
        }
        log::info!("Successful output clm file `{}`", output);
    }

    pub fn filter_by_mapq_single_threads(&mut self, min_quality: u8, output: &String) {
        let mut wtr = common_writer(output);
        self.header = PairHeader::new();
        self.header.from_pairs(&self.file);
        wtr.write_all(self.header.to_string().as_bytes()).unwrap();

        let mut reader = common_reader(&self.file);
        let mut buffer = String::new();
        let filter_mapq = min_quality > 0;
    
        let mut record_len = 0;
        
        for line in reader.lines().filter_map(|line| line.ok()) {
            if line.starts_with('#') {
                continue;
            }

            let record = line.trim_end_matches('\n').split('\t').collect::<Vec<&str>>();
            
            if record_len == 0 {
                record_len = record.len();
            }

            if record_len < 8 && filter_mapq {
                buffer.push_str(&line);
                buffer.push('\n');
                if buffer.len() >= 8192 {
                    wtr.write_all(buffer.as_bytes()).unwrap();
                    buffer.clear();  
                }
                continue
            }
            
            if let Ok(mapq) = record[7].parse::<u8>() {
                if mapq >= min_quality {
                    buffer.push_str(&line);
                    buffer.push('\n');
                    if buffer.len() >= 8192 {
                        wtr.write_all(buffer.as_bytes()).unwrap();
                        buffer.clear();  
                    }
                }
            } else {
                eprintln!("Warning: could not parse mapq value: {}", record[7]);
            }
        }
        
        if !buffer.is_empty() {
            wtr.write_all(buffer.as_bytes()).unwrap();
        }

           
        log::info!("Successful output new pairs file `{}`", output);
    }

    pub fn filter_by_mapq(&mut self, min_quality: u8, output: &String) {
        let mut wtr = common_writer(output);
        self.header = PairHeader::new();
        self.header.from_pairs(&self.file);
        wtr.write_all(self.header.to_string().as_bytes()).unwrap();

        let mut reader = common_reader(&self.file);

        let (sender, receiver) = bounded::<Vec<String>>(1000);

        let mut handles = vec![];
        let filter_mapq = min_quality > 0;
        let mut wtr = Arc::new(Mutex::new(wtr));

        for _ in 0..4 {
            let receiver = receiver.clone();
            let wtr = Arc::clone(&wtr);
            handles.push(thread::spawn(move || {
                while let Ok((records)) = receiver.recv() {
                    // let filter_needed = filter_mapq && records.get(0).map_or(false, |record| record.len() >= 8);
                    let no_filter = !filter_mapq || records.get(0).map_or(false, |record| record.len() < 8);
                
                    let mut buffer: Vec<&String> = records.par_iter().filter_map(|line| {
                        let record = line.trim_end_matches('\n').split('\t').collect::<Vec<&str>>();
                        
                        let record_len = record.len();

                        if no_filter  {
                            // buffer.push_str(&line);
                            // buffer.push('\n');
                            return Some(line);
                        }

                        if record_len > 7{
                            if let Ok(mapq) = record[7].parse::<u8>() {
                                if mapq >= min_quality {
                                    return Some(line);
                                    
                                }
                            } 
                        } else {
                            return Some(line);
                        }
                        

                        None
                    }).collect();

                    {
                        let mut wtr = wtr.lock().unwrap();
                        for line in &buffer {
                            if !line.trim().is_empty(){
                                wtr.write_all(line.as_bytes()).unwrap();
                                wtr.write_all(b"\n").unwrap()
                            }
                            
                        }
                    }
                }
            }));

        }

        let mut batch = Vec::with_capacity(1000);
        for line in reader.lines().filter_map(|line| line.ok()) {
            if line.starts_with('#') {
                continue;
            }
            batch.push(line);
            if batch.len() >= 1000 {
                sender.send(std::mem::take(&mut batch)).unwrap();
            }
        }

        if !batch.is_empty() {
            sender.send(batch).unwrap();
        }

        drop(sender);
        for handle in handles {
            handle.join().unwrap();
        }

        log::info!("Successful output new pairs file `{}`", output);
    }

    pub fn intersect(&mut self, hcr_bed: &String, invert: bool, min_quality: u8, output: &String) {
        type IvU8 = Interval<usize, u8>;
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        let mut writer = common_writer(output);
        self.header = PairHeader::new();
        self.header.from_pairs(&self.file);
        writer.write_all(self.header.to_string().as_bytes()).unwrap();
        let mut wtr = csv::WriterBuilder::new()
                        .has_headers(false)
                        .delimiter(b'\t')
                        .from_writer(writer);
        
        let filter_mapq = min_quality > 0;
        for (idx, record) in self.parse().unwrap().records().enumerate() {
            let record = record.unwrap();

            if filter_mapq {
                if record.len() >= 8 {
                    let mapq = record.get(7).unwrap_or(&"60").parse::<u8>().unwrap_or_default();
                    if mapq < min_quality {
                        continue
                    }
                }
            }
            let chrom1 = &record[1];
            let pos1 = record[2].parse::<usize>().unwrap();
            let chrom2 = &record[3];
            let pos2 = record[4].parse::<usize>().unwrap();

            let is_in_regions = interval_hash.get(chrom1).map_or(false, |interval1| {
                interval_hash.get(chrom2).map_or(false, |interval2| {
                    interval1.count(pos1, pos1+1) > 0 && interval2.count(pos2, pos2+1) > 0
                })
            });

            if is_in_regions ^ invert {
                wtr.write_record(&record);
            }
              
        }


        log::info!("Successful output new pairs file into `{}`", output);
    }


    pub fn intersect_multi_threads(&mut self, hcr_bed: &String, invert: bool, 
                                    min_quality: u8, output: &String) {
       
        type IvU8 = Interval<usize, u8>;
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        
        let mut parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };


        let mut wtr = common_writer(output);
        wtr.write_all(self.header.to_string().as_bytes()).unwrap();
        
        log::info!("Parsing pairs file...");
        let (sender, receiver) = bounded::<Vec<(usize, StringRecord)>>(1000);
        let mut handles = vec![];
        let filter_mapq = min_quality > 0;
        let mut wtr = Arc::new(Mutex::new(wtr));
        for _ in 0..4 {
            let receiver = receiver.clone();
            let interval_hash = interval_hash.clone();
            let wtr = Arc::clone(&wtr);
            handles.push(thread::spawn(move || {
                while let Ok((records)) = receiver.recv() {
                    let filter_needed = filter_mapq && records.get(0).map_or(false, |(_, record)| record.len() >= 8);
                    let mut data = vec![];
                    for (idx, record) in records.iter().filter(|(idx, record)| {
                        !filter_needed || record.get(7)
                            .and_then(|mapq_str| mapq_str.parse::<u8>().ok())
                            .map_or(true, |mapq| mapq >= min_quality)
                    }) {
                        let chrom1 = &record[1];
                        let pos1 = record[2].parse::<usize>().unwrap();
                        let chrom2 = &record[3];
                        let pos2 = record[4].parse::<usize>().unwrap();
                        let is_in_regions = interval_hash.get(chrom1).map_or(false, |interval1| {
                            interval_hash.get(chrom2).map_or(false, |interval2| {
                                interval1.count(pos1, pos1+1) > 0 && interval2.count(pos2, pos2+1) > 0
                            })
                        });

                        if is_in_regions ^ invert {
                            let record = record.iter().map(|x| x.to_string()).collect::<Vec<String>>();
                            if !record.is_empty() {
                                data.push(record.join("\t"));
                            }
                        }
                    }
                    
                    if !data.is_empty() {
                        let data = data.join("\n") + "\n";
                        let mut wtr = wtr.lock().unwrap();
                        wtr.write_all(data.as_bytes()).unwrap(); 
                    }
                }
            }));
        }


        let mut batch = Vec::with_capacity(1000);
        for (idx, record) in rdr.records().enumerate() {
            let record = record.unwrap();
            batch.push((idx, record));
            if batch.len() >= 1000 {
                // sender.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).unwrap();
                sender.send(std::mem::take(&mut batch)).unwrap();
            }
        }

        if !batch.is_empty() {
            sender.send(batch).unwrap();
        }

        drop(sender);
        for handle in handles {
            handle.join().unwrap();
        }

        log::info!("Successful output new pairs file into `{}`", output);
        

    }

    pub fn to_depth(&mut self, binsize: u32, min_quality: u8, output: &String) {
        use hashbrown::HashMap as BrownHashMap;
        
        let mut parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };
        let pair_header = &self.header;
        let contigsizes_data = &pair_header.chromsizes;
        let contigsizes_data: BrownHashMap<String, u64> = contigsizes_data.iter().map(|x| (x.chrom.clone(), x.size)).collect();
        
        let filter_mapq = min_quality > 0;

        let mut depth: Arc<Mutex<BrownHashMap<String, Vec<u32>>>> = Arc::new(Mutex::new(BrownHashMap::new()));
        let (sender, receiver) = bounded::<Vec<(usize, StringRecord)>>(1000);
        let mut handles = vec![];
        for _ in 0..4 {
            let receiver = receiver.clone();
            let depth = Arc::clone(&depth);
            let contigsizes_data = contigsizes_data.clone();
            handles.push(thread::spawn(move || {

                while let Ok((records)) = receiver.recv() {

                    let mut updates: BrownHashMap<String, Vec<u32>> = BrownHashMap::new();
                    let filter_needed = filter_mapq &&  records.get(0).map_or(false, |(_, record)| record.len() >= 8);

                    for (idx, record) in records.iter().filter(|(idx, record)| {
                        !filter_needed || record.get(7)
                            .and_then(|mapq_str| mapq_str.parse::<u8>().ok())
                            .map_or(true, |mapq| mapq > min_quality)
                    }) {
                        let contig1 = &record[1];
                        let contig2 = &record[3];
                        let pos1 = record[2].parse::<u32>().unwrap();
                        let pos2 = record[4].parse::<u32>().unwrap();
                        for contig in [&contig1, &contig2] {
                            let size = *contigsizes_data.get(*contig).unwrap_or(&0);
                            if size != 0 {
                                let bin_index = if contig == &contig1 { pos1 } else { pos2 } as usize / binsize as usize;
                                updates.entry(contig.to_string())
                                       .or_insert_with(|| vec![0; (size / binsize as u64 + 1) as usize])
                                       [bin_index] += 1;
                            }
                        }

                    }
                    {
                        let mut depth = depth.lock().unwrap();
                        for (key, local_values) in updates {
                            if let Some(depth_values) = depth.get_mut(&key) {
                                for (i, value) in local_values.iter().enumerate() {
                                    depth_values[i] += value;
                                }
                            } else {
                                depth.insert(key, local_values);
                            }
                        }
                    }
                }

            }));
        }

        let mut batch = Vec::with_capacity(1000);
        for (idx, record) in rdr.records().enumerate() {
            let record = record.unwrap();
            
            batch.push((idx, record));
            if batch.len() >= 1000 {
                // sender.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).unwrap();
                sender.send(std::mem::take(&mut batch)).unwrap();
            }

        }

        if !batch.is_empty() {
            sender.send(batch).unwrap();
        }

        drop(sender);
        for handle in handles {
            handle.join().unwrap();
        }

        let mut depth =  Arc::try_unwrap(depth).unwrap().into_inner().unwrap();

        // let mut depth: BrownHashMap<String, Vec<u32>> = BrownHashMap::new();
        // for (i, record) in rdr.records().enumerate() {
        //     let record = record.unwrap();
        //     if filter_mapq {
        //         if record.len() >= 8 {
        //             let mapq = record.get(7).unwrap_or(&"60").parse::<u8>().unwrap_or_default();
        //             if mapq <= min_quality {
        //                 continue
        //             }
        //         }
        //     }

        //     let contig1 = &record[1];
        //     let contig2 = &record[3];
        //     let pos1 = record[2].parse::<u32>().unwrap();
        //     let pos2 = record[4].parse::<u32>().unwrap();

        //     if !depth.contains_key(&record[1]) {
        //         let size1 = *contigsizes_data.get(contig1).unwrap_or(&0);
        //         if size1 != 0 {
        //             let bins1 = depth.entry(contig1.to_string()).or_insert_with(|| vec![0; (size1 / binsize as u64 + 1) as usize]);
        //             bins1[pos1 as usize / binsize as usize] += 1;
        //         }
        
        //     } 
        //     if !depth.contains_key(&record[3]) {
        //         let size2 = *contigsizes_data.get(contig2).unwrap_or(&0);
        //         if size2 != 0 {
        //             let bins2 = depth.entry(contig2.to_string()).or_insert_with(|| vec![0; (size2 / binsize as u64 + 1) as usize]);
        //             bins2[pos2 as usize / binsize as usize] += 1;
        //         }
               
        //     }

        //     if let Some(bins1) = depth.get_mut(&record[1]) {
        //         bins1[pos1 as usize / binsize as usize] += 1;
        //     }
            
        //     if let Some(bins2) = depth.get_mut(&record[3]) {
        //         bins2[pos2 as usize / binsize as usize] += 1;
        //     }
        
        // }

        log::info!("Collected data from pairs file.");
        let depth: BTreeMap<_, _> = depth.into_iter().collect();
        let mut wtr = common_writer(output);
        let wtr = Arc::new(Mutex::new(wtr));

        depth.par_iter().for_each(|(contig, bins)| {
            let size = contigsizes_data.get(contig).unwrap_or(&0);
            let mut buffer = Vec::with_capacity(bins.len() * 50);
            for (bin, count) in bins.iter().enumerate() {
                        let bin_start = bin * binsize as usize;
                        let mut bin_end = bin_start + binsize as usize;
        
                        if bin_end > (*size).try_into().unwrap() {
                            bin_end = *size as usize;
                        }
                        
                        buffer.extend_from_slice(format!("{}\t{}\t{}\t{}\n", contig, bin_start, bin_end, count).as_bytes());
            }
            {
                let mut wtr = wtr.lock().unwrap();
                wtr.write_all(&buffer).unwrap();
            }
            
        });

        log::info!("Successful output depth file into `{}`", output);
        
    }

    pub fn break_contigs(&mut self, break_bed: &String, output: &String) {
        type IvString = Interval<usize, String>;
        let bed = Bed4::new(break_bed);
        let interval_hash = bed.to_interval_hash();
        
        let mut writer = common_writer(output);

        let mut parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let pair_header = &self.header;
        let contigsizes = &pair_header.chromsizes;
        let contigsizes_data: HashMap<String, u64> = contigsizes.iter().map(|x| (x.chrom.clone(), x.size)).collect();
        let mut new_contigsizes_data = contigsizes_data.clone();
        for contig in contigsizes_data.keys() {
            if interval_hash.contains_key(contig) {
                new_contigsizes_data.remove(contig);
                let interval = interval_hash.get(contig).unwrap();
                let mut res = interval.iter().collect::<Vec<_>>();
                for sub_res in res.iter() {
                    // let new_contig = format!("{}:{}-{}", contig, sub_res.start, sub_res.stop);
                    let new_contigs = &sub_res.val;
                    let size = sub_res.stop - sub_res.start + 1;
                    new_contigsizes_data.insert(new_contigs.clone(), size.try_into().unwrap());
                }
            }
        }
        // HashMap to Vec<ChromSizeRecord>
        let mut new_contigsizes: Vec<ChromSizeRecord> = new_contigsizes_data.iter().map(
                    |(chrom, size)| ChromSizeRecord{chrom: chrom.clone(), size: *size}).collect();

        self.header = PairHeader::new();
        self.header.from_chromsizes(new_contigsizes);
        writer.write_all(self.header.to_string().as_bytes()).unwrap();
        
        let mut wtr = csv::WriterBuilder::new()
                        .has_headers(false)
                        .delimiter(b'\t')
                        .from_writer(writer);

        for (idx, record) in rdr.records().enumerate() {
            let record = record.unwrap();
            let chrom1 = &record[1];
            let chrom2 = &record[3];
        
            let is_break_contig1 = interval_hash.contains_key(chrom1);
            let is_break_contig2 = interval_hash.contains_key(chrom2);

            if is_break_contig1 || is_break_contig2 {
                let pos1 = record[2].parse::<usize>().unwrap();
                let pos2 = record[4].parse::<usize>().unwrap();
                let mut new_record = csv::StringRecord::new();
                let (new_chrom1, new_pos1): (&str, usize) = match is_break_contig1 {
                    true => {
                        let interval = interval_hash.get(chrom1).unwrap();
                        let mut res = interval.find(pos1, pos1+1).collect::<Vec<_>>();
                        if res.len() > 0 {
                            let new_pos = pos1 - res[0].start + 1;
                            let new_chrom1 = &res[0].val;
                            (new_chrom1, new_pos )
                        } else {
                            (chrom1, pos1 )
                        }
                    }
                    false => {
                        (chrom1, pos1 )
                    }
                };

                let (new_chrom2, new_pos2): (&str, usize) = match is_break_contig2 {
                    true => {
                        let interval = interval_hash.get(chrom2).unwrap();
                        let mut res = interval.find(pos2, pos2+1).collect::<Vec<_>>();
                        if res.len() > 0 {
                            let new_pos = pos2 - res[0].start + 1;
                            // let new_chrom2 = format!("{}:{}-{}", chrom2, res[0].start, res[0].stop);
                            let new_chrom2 = &res[0].val;
                            (new_chrom2, new_pos)
                        } else {
                            (chrom2, pos2)
                        }
                    }
                    false => {
                        (chrom2, pos2 )
                    }
                };

                new_record.push_field(&record[0]);
                new_record.push_field(&new_chrom1);
                new_record.push_field(&new_pos1.to_string());
                new_record.push_field(&new_chrom2);
                new_record.push_field(&new_pos2.to_string());
                new_record.push_field(&record[5]);
                new_record.push_field(&record[6]);

                if record.len() >= 8 {
                    new_record.push_field(&record[7]);
                } 

                if record.len() >= 9 {
                    new_record.push_field(&record[8]);
                }

                wtr.write_record(&new_record).unwrap();


            } else {
                wtr.write_record(&record).unwrap();
            }

        }

        log::info!("Successful output new pairs file into `{}`", output);
    }

    pub fn dup(&mut self, collapsed_list: &String, seed: usize, output: &String) {
        // random assign a collapsed contig to duplicated contig by copy number

        let reader = common_reader(collapsed_list);
        let mut collapsed_contigs: HashMap<String, Vec<String>> = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let s: Vec<&str> = record.split("\t").collect();
            if s.len() != 2 {
                log::warn!("Invalid record: {}", record);
                continue;
            }
            let contig1 = s[0].to_string();
            let contig2 = s[1].to_string();
            collapsed_contigs.entry(contig1.clone()).or_insert(vec![contig1]).push(contig2);
        }
        
        let seed_bytes = seed.to_ne_bytes();
        let mut seed_array = [0u8; 32];
        for (i, byte) in seed_bytes.iter().enumerate() {
            seed_array[i] = *byte;
        }

        let reader = common_reader(&self.file);

        let mut wtr = common_writer(output);
       
        let mut output_counts = 0;
        let mut rng = StdRng::from_seed(seed_array);
        for (idx, record) in reader.lines().enumerate() {
            let record = record.unwrap();
            if record.starts_with("#") {
                if record.starts_with("#chromsize") {
                    writeln!(wtr, "{}", record);
                    let s: Vec<&str> = record.trim().split(" ").collect();
                    let size: u64 = s[2].clone().parse::<u64>().unwrap();
                    let chrom: String = s[1].to_string();
                    if collapsed_contigs.contains_key(&chrom) {
                        let collapsed_contigs = collapsed_contigs.get(&chrom).unwrap();
                        for contig in collapsed_contigs {
                            writeln!(wtr, "#chromsize: {} {}", contig, size);
                        }
                    }
                } else {
                    writeln!(wtr, "{}", record);
                }
                    
                continue;
            }

            let s = record.trim().split("\t").collect::<Vec<&str>>();
            let mut contig1 = s[1].to_string();
            let mut contig2 = s[3].to_string();
            if collapsed_contigs.contains_key(&contig1) {
                let idx = rng.gen_range(0..collapsed_contigs.get(&contig1).unwrap().len());
                contig1 = collapsed_contigs.get(&contig1).unwrap()[idx].clone();
            }
            if collapsed_contigs.contains_key(&contig2) {
                let idx = rng.gen_range(0..collapsed_contigs.get(&contig2).unwrap().len());
                contig2 = collapsed_contigs.get(&contig2).unwrap()[idx].clone();
            }
            let record = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", s[0], contig1, s[2], contig2, s[4], s[5], s[6], s[7]);
            writeln!(wtr, "{}", record);
        
        }

    }

    pub fn downsample(&mut self, n: usize, p: f64, seed: usize, output: &String) {

        let seed_bytes = seed.to_ne_bytes();
        let mut seed_array = [0u8; 32];
        for (i, byte) in seed_bytes.iter().enumerate() {
            seed_array[i] = *byte;
        }
        
        let percent = if p > 0.0 {
            p
        } else {
            let reader = common_reader(&self.file);
            let mut total_count = 0;
            let mut skip_count = 0;
            for chunk in reader.bytes() {
                let bytes = chunk.unwrap();
                total_count += bytecount::count(&[bytes], b'\n');
                skip_count += bytecount::count(&[bytes], b'#');
            }
            let record_counts = total_count - skip_count;
            log::info!("Total records: {}", record_counts);
            n as f64 / record_counts as f64

        };
        
        let reader = common_reader(&self.file);

        let mut wtr = common_writer(output);
       
        let mut output_counts = 0;
        let mut rng = StdRng::from_seed(seed_array);
        for (idx, record) in reader.lines().enumerate() {
            let record = record.unwrap();
            if record.starts_with("#") {
                writeln!(wtr, "{}", record);
                continue;
            }
            let random_number: f64 = rng.gen();
            if random_number < percent {
                output_counts += 1;
                if ((output_counts > n) && (p == 0.0)) {
                    break;
                }
                writeln!(wtr, "{}", record);
            }
        }
        
        // let (tx, rx) = mpsc::channel();
        // let num_threads = 10;
        // reader.lines().enumerate().par_bridge()
        //     .for_each_with(tx.clone(), |s, (idx, record)| {
        //         let record = record.unwrap();
        //         if record.starts_with("#") {
        //             tx.send(record).unwrap();

        //         } else {
        //             let mut rng = StdRng::from_seed(seed_array);
        //             let random_number: f64 = rng.gen();
        //             if random_number < percent {
        //                 tx.send(record).unwrap();
        //             }
        //         }
        //     });
        
        // let mut output_counts = 0;
        // for record in rx {
        //     if output_counts > n && p == 0.0 {
        //         break;
        //     }

        //     writeln!(wtr, "{}", record);
        //     output_counts += 1
        // }




        log::info!("Successful output new pairs file into `{}`", output);
        
    }
}


pub fn merge_pairs(input: Vec<&String>, output: &String) {
    // let mut wtr = common_writer(output);

    // for (i, file) in input.iter().enumerate() {
    //     let mut reader = common_reader(file);
    //     for line in reader.lines() {
    //         let line = line.unwrap();
    //         if line.starts_with("#") {
    //             if i == 0  {
    //                 writeln!(wtr, "{}", line);
    //             }  else {
    //                 continue
    //             }
    //         } else {
    //             writeln!(wtr, "{}", line);
    //         }

    //     }
    // }
    // wtr.flush().unwrap()

    let mut wtr = common_writer(output);
    
    if let Some(first_file) = input.first() {
        let reader = common_reader(first_file);
        for line in reader.lines() {
            let line = line.unwrap();
            if line.starts_with("#") {
                writeln!(wtr, "{}", line).unwrap();
            } else {
                break; 
            }
        }
    }

    let mut wtr = Arc::new(Mutex::new(wtr));
    input.par_iter().enumerate().for_each(|(i, file)| {
        let reader = common_reader(file);
        let mut local_buffer = Vec::new();

        let mut j = 0;
        for line in reader.lines() {
            let line = line.unwrap();
            if line.starts_with("#") {
                continue
            }

            local_buffer.push(line);
            j += 1;

            if j >= 20000 {
                let mut wtr = wtr.lock().unwrap();
                for line in local_buffer.drain(..) {
                    writeln!(wtr, "{}", line).unwrap();
                }
                j = 0;
            }
        }
    });


}


