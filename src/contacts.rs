use anyhow::Result as AnyResult;
use std::borrow::Cow;
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::io::{ Read, Write };
use serde::{ Serialize, Deserialize };
use rayon::prelude::*;

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ContactRecord {
    pub chrom1: String,
    pub chrom2: String,
    pub count: u32,
}

impl ContactRecord {
    pub fn new() -> Self {
        Self {
            chrom1: String::new(),
            chrom2: String::new(),
            count: 0,
        }
    }

    pub fn is_some(&self) -> bool {
        match self.count {
            0 => false,
            _ => true,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Contacts {
    pub file: String,
    pub records: Vec<ContactRecord>,
}

impl BaseTable for Contacts {
    fn new(name: &String) -> Contacts {
        Contacts {
            file: name.clone(),
            records: Vec::new(),
        }
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


impl Contacts {
    pub fn parse(&mut self) {

       let input = common_reader(&self.file);
         let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(false)
                .from_reader(input);
    
          let mut records: Vec<ContactRecord> = Vec::new();
          for result in rdr.deserialize() {
                let record: ContactRecord = result.unwrap();
                self.records.push(record);
          }
       
    }

    pub fn to_data(&self, re_count: HashMap<String, u32>, symmetric: bool) -> HashMap<ContigPair, f64> {
        
        let mut data: HashMap<ContigPair, f64> = self.records.par_iter(
            ).map(|record| {
                let mut contig_pair = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
                if !re_count.contains_key(&record.chrom1) || !re_count.contains_key(&record.chrom2) {
                    (contig_pair, record.count as f64)

                } else {
                    let re_count1 = re_count.get(&record.chrom1).unwrap();
                    let re_count2 = re_count.get(&record.chrom2).unwrap();
                    let re_count = (*re_count1 + *re_count2) as f64;
                    let count = record.count as f64;
                    let ratio = count / re_count;
                    (contig_pair, ratio)
                }
            }).collect();
        data 
    }


    pub fn write(&self, output: &String) {
        let mut wtr = common_writer(output);
        for record in &self.records {
            wtr.write_all(format!("{}\t{}\t{}\n", record.chrom1, record.chrom2, record.count).as_bytes()).unwrap();
        }
    }
}