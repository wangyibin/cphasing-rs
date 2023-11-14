use anyhow::Result as AnyResult;
use std::borrow::Cow;
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::io::{ Read, Write };
use serde::{ Serialize, Deserialize };
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PixelRecord {
    pub chrom1: String,
    pub start1: u32,
    pub end1: u32,
    pub chrom2: String,
    pub start2: u32,
    pub end2: u32,
    pub count: u32,
}

impl PixelRecord {
    pub fn new() -> Self {
        Self {
            chrom1: String::new(),
            start1: 0,
            end1: 0,
            chrom2: String::new(),
            start2: 0,
            end2: 0,
            count: 0,
        }
    }

}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Pixels {
    pub file: String,
    pub records: Vec<PixelRecord>,
}

impl BaseTable for Pixels {
    fn new(name: &String) -> Pixels {
        Pixels {
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


impl Pixels {
    pub fn parse(&mut self) {

       let input = common_reader(&self.file);
         let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(false)
                .from_reader(input);
    
          let mut records: Vec<PixelRecord> = Vec::new();
          for result in rdr.deserialize() {
                let record: PixelRecord = result.unwrap();
                self.records.push(record);
          }
       
    }

    pub fn to_data(&self, re_count: HashMap<String, u32>, symmetric: bool) -> HashMap<ContigPair, f64> {
        let mut data: HashMap<ContigPair, f64> = HashMap::new();
        for record in &self.records {
            let mut contig_pair = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
            if !re_count.contains_key(&record.chrom1) || !re_count.contains_key(&record.chrom2) {
                data.insert(contig_pair, record.count as f64);

            } else {
                let re_count1 = re_count.get(&record.chrom1).unwrap();
                let re_count2 = re_count.get(&record.chrom2).unwrap();
                let re_count = (*re_count1 + *re_count2) as f64;
                let count = record.count as f64;
                let ratio = count / re_count;
                data.insert(contig_pair, ratio);
            }
            
        }

        data 
    }
}