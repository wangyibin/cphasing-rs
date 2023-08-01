
use std::fs::File;
use std::io::Write;
use std::path::Path;
use serde::{Deserialize, Serialize};

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
    pub file: File,
    pub reader: csv::Reader<File>,
}

impl BedMethylSimple {
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

impl Iterator for BedMethylSimple {
    type Item = BedMethylSimpleRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.reader.records().next();
        match record {
            Some(Ok(record)) => {
                let chrom = record[0].to_string();
                let start = record[1].parse::<usize>().unwrap();
                let end = record[2].parse::<usize>().unwrap();
                let mod_base = record[3].to_string();
                let score = record[4].parse::<u32>().unwrap();
                let strand = record[5].parse::<char>().unwrap();
                let frac = record[10].parse::<f32>().unwrap();

                Some(BedMethylSimpleRecord {
                    chrom: chrom,
                    start: start,
                    end: end,
                    mod_base: mod_base,
                    score: score,
                    strand: strand,
                    frac: frac
                })
            }
            _ => None,
        }
    }
}