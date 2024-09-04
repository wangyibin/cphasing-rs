
use std::collections::HashMap;
use std::fs::File;
use std::io::{ Write, BufRead };
use std::path::Path;
use serde::{Deserialize, Serialize};
use rust_lapper::{Interval, Lapper};

use crate::core::{common_reader, common_writer};
use crate::methy::ModRecord;


type IvU8 = Interval<usize, u8>;
type IvString = Interval<usize, String>;


pub trait BedIterator {
    fn next(&mut self) -> Option<ModRecord>;

}

#[derive(Debug, Deserialize, Serialize)]
pub struct Bed3Record {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
}

pub struct Bed3 {
    pub file: File,
    pub reader: csv::Reader<File>,
}

impl Bed3 {
    pub fn new(bed: &String) -> Self {
        if Path::new(bed).exists() {
            let file = File::open(bed).unwrap();
            let reader = csv::ReaderBuilder::new()
                .flexible(true)
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(file.try_clone().unwrap());
            Self {
                file: file,
                reader: reader,
            }
        } else {
            Self {
                file: File::create(bed).unwrap(),
                reader: csv::ReaderBuilder::new()
                    .flexible(true)
                    .delimiter(b'\t')
                    .has_headers(false)
                    .from_reader(File::create(bed).unwrap().try_clone().unwrap()),
            }
        }
        
    }
    
    pub fn to_interval_hash(self) -> HashMap<String, Lapper<usize, u8>> {
        let mut hcr: HashMap<String, Vec<IvU8>> = HashMap::new();
        for i in self {
            if hcr.contains_key(&i.chrom) {
                hcr.get_mut(&i.chrom).unwrap().push(IvU8 {
                    start: i.start,
                    stop: i.end,
                    val: 0,
                });
            } else {
                hcr.insert(i.chrom, vec![IvU8 {
                    start: i.start,
                    stop: i.end,
                    val: 0,
                }]);
            }
        }
        let mut hcr_lapper: HashMap<String, Lapper<usize, u8>> = HashMap::new();
        for (k, v) in hcr.iter() {
            hcr_lapper.insert(k.to_string(), Lapper::new(v.to_vec()));
        }
        hcr_lapper
    }
}


impl Iterator for Bed3 {
    type Item = Bed3Record;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.reader.records().next();
        match record {
            Some(Ok(record)) => {
                let chrom = record[0].to_string();
                let start = record[1].parse::<usize>().unwrap();
                let end = record[2].parse::<usize>().unwrap();
                Some(Bed3Record {
                    chrom: chrom,
                    start: start,
                    end: end,
                })
            }
            _ => None,
        }
    }
}


#[derive(Debug, Deserialize, Serialize)]
pub struct Bed4Record {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub gene: String,
}

pub struct Bed4 {
    pub file: File,
    pub reader: csv::Reader<File>,
}

impl Bed4 {
    pub fn new(bed: &String) -> Self {
        if Path::new(bed).exists() {
            let file = File::open(bed).unwrap();
            let reader = csv::ReaderBuilder::new()
                .flexible(true)
                .delimiter(b'\t')
                .has_headers(false)
                .from_reader(file.try_clone().unwrap());
            Self {
                file: file,
                reader: reader,
            }
        } else {
            Self {
                file: File::create(bed).unwrap(),
                reader: csv::ReaderBuilder::new()
                    .flexible(true)
                    .delimiter(b'\t')
                    .has_headers(false)
                    .from_reader(File::create(bed).unwrap().try_clone().unwrap()),
            }
        }
        
    }
    
    pub fn to_interval_hash(self) -> HashMap<String, Lapper<usize, String>> {
        let mut hcr: HashMap<String, Vec<IvString>> = HashMap::new();
        for i in self {
            if hcr.contains_key(&i.chrom) {
                hcr.get_mut(&i.chrom).unwrap().push(IvString {
                    start: i.start,
                    stop: i.end,
                    val: i.gene,
                });
            } else {
                hcr.insert(i.chrom, vec![IvString {
                    start: i.start,
                    stop: i.end,
                    val: i.gene,
                }]);
            }
        }
        let mut hcr_lapper: HashMap<String, Lapper<usize, String>> = HashMap::new();
        for (k, v) in hcr.iter() {
            hcr_lapper.insert(k.to_string(), Lapper::new(v.to_vec()));
        }
        hcr_lapper
    }
}


impl Iterator for Bed4 {
    type Item = Bed4Record;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.reader.records().next();
        match record {
            Some(Ok(record)) => {
                let chrom = record[0].to_string();
                let start = record[1].parse::<usize>().unwrap();
                let end = record[2].parse::<usize>().unwrap();
                let gene = record[3].to_string();
                Some(Bed4Record {
                    chrom: chrom,
                    start: start,
                    end: end,
                    gene: gene
                })
            }
            _ => None,
        }
    }
}


#[derive(Debug, Deserialize, Serialize)]
pub struct BedCpGRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub frac: f32,
    pub score: u32,
}

pub struct BedCpG {
    pub file: File,
    pub reader: csv::Reader<File>,
}

impl BedCpG {
    pub fn new(bed: &String) -> Self {
        let file = File::open(bed).unwrap();
        let reader = csv::ReaderBuilder::new()
            .flexible(true)
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(file.try_clone().unwrap());
        Self {
            file: file,
            reader: reader,
        }
    }
    
}

impl BedIterator for BedCpG {

    fn next(&mut self) -> Option<ModRecord> {
        let record = self.reader.records().next();
        match record {
            Some(Ok(record)) => {
                let chrom = record[0].to_string();
                let start = record[1].parse::<usize>().unwrap();
                let end = record[2].parse::<usize>().unwrap();
                let frac = record[3].parse::<f32>().unwrap();
                let score = record[6].parse::<u32>().unwrap();
                Some(ModRecord {
                    chrom: chrom,
                    start: start,
                    end: end,
                    frac: frac,
                    score: score,
                })
            }
            _ => None,
        }
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct BedMethylSimpleRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub mod_base: String,
    pub score: u32,
    pub strand: char,
    pub frac: f32,
}

pub struct BedMethylSimple {
    pub file: String,
    pub reader: csv::Reader<Box<dyn BufRead + Send>>,
}

impl BedMethylSimple {
    pub fn new(bed: &String) -> Self {
        let _file = common_reader(&bed);
        let reader = csv::ReaderBuilder::new()
            .flexible(true)
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(_file);
        Self {
            file: bed.clone(),
            reader: reader,
        }
    }
    
}

impl BedIterator for BedMethylSimple {

    fn next(&mut self) -> Option<ModRecord> {
        let record = self.reader.records().next();
        match record {
            Some(Ok(record)) => {
                let chrom = record[0].to_string();
                let start = record[1].parse::<usize>().unwrap();
                let end = record[2].parse::<usize>().unwrap();
                let mod_base = record[3].to_string();
                let score = record[4].parse::<u32>().unwrap();
                let strand = record[5].parse::<char>().unwrap();
                let record_9_str = record[9].to_string();
                let frac = match record_9_str.parse::<f32>() {
                    Ok(value) => value,
                    Err(_) => {
                        let other_records: Vec<&str> = record_9_str.split(' ').collect();
                        match other_records.get(1) {
                            Some(&value) => value.parse::<f32>().unwrap_or(0.0),
                            None => 0.0,
                        }
                    }
                };
                

                Some(ModRecord {
                    chrom: chrom,
                    start: start,
                    end: end,
                    frac: frac,
                    score: score
                })
            }
            _ => None,
        }
    }
}


impl Iterator for Box<dyn BedIterator> {
    type Item = ModRecord;

    fn next(&mut self) -> Option<Self::Item> {
        (**self).next()
    }
}