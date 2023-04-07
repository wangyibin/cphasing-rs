use std::collections::HashMap;
use std::borrow::Cow;
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ChromSizeRecord };
use crate::porec::PoreCRecord;



#[derive(Debug, Clone)]
pub struct PairHeader {
    pub chromsizes: Vec<ChromSizeRecord>,
    pub header: Vec<String>,
}

impl PairHeader {
    pub fn new() -> PairHeader {
        let mut chromsizes: Vec<ChromSizeRecord> = Vec::new();
        let mut header: Vec<String> = Vec::new();

        PairHeader { chromsizes: chromsizes, header: header }
    }

    pub fn from_pairs(&mut self, name: &String) {
        let reader = common_reader(name);
        for line in reader.lines() {
            if let Ok(line) = line {
                if line.starts_with("#") {
                    if &line[0..10] == "#chromsize" {
                        let s: String = line.replace("#chromsize: ", "");
                        let s: Vec<&str> = s.split(" ").collect();
                        
                        let size: u32 = s[1].clone().parse::<u32>().unwrap();
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
        let PAIR_HEADER: Vec<String> = vec!["readID".to_string(), 
                                    "chrom1".to_string(), 
                                    "pos1".to_string(), 
                                    "chrom2".to_string(), 
                                    "pos2".to_string(), 
                                    "strand1".to_string(),
                                    "strand2".to_string()
                                    ];
        self.chromsizes = chromsizes;
        self.header = PAIR_HEADER;
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
    pub pos1: u32,
    pub chrom2: String,
    pub pos2: u32,
    pub strand1: char,
    pub strand2: char,
}

impl PairRecord {
    pub fn from_pore_c_pair(pair: Vec<PoreCRecord>, readID: u64) -> PairRecord {
        let (pair1, pair2) = (pair[0].clone(), pair[1].clone());
        
        let pos1: u32 = (pair1.target_end + pair1.target_start) / 2;
        let pos2: u32 = (pair2.target_end + pair2.target_start) / 2;

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


