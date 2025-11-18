#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
use anyhow::Result as AnyResult;
use indexmap::IndexMap;
use std::borrow::Cow;
use std::collections::HashMap;
use std::path::Path;
use std::io::BufRead;
use serde::{ Deserialize, Serialize };

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CountReRecord {
    Contig: String, 
    RECounts: u32,
    Length: u32,
}

#[derive(Debug, Clone)]
pub struct CountRE {
    file: String,
    records: Vec<CountReRecord>,
}

impl BaseTable for CountRE {
    fn new(name: &String) -> CountRE {
        CountRE {
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
        file_prefix.to_string()
    }

}

impl CountRE {
    pub fn parse(&mut self) {
        let mut input = common_reader(&self.file);
        for (line_number, line) in input.lines().enumerate() {
            match line {
                Ok(content) => {
      
                    if content.starts_with('#') {
                        continue;
                    }
    
        
                    let fields: Vec<&str> = content.split('\t').collect();
    
                    if fields.len() != 3 {
                        log::info!(
                            "Skipping line {}: expected 3 fields, found {}",
                            line_number + 1,
                            fields.len()
                        );
                        continue;
                    }
    
                    let contig = fields[0].to_string();
                    let re_counts = match fields[1].parse::<u32>() {
                        Ok(value) => value,
                        Err(_) => {
                            log::info!(
                                "Skipping line {}: invalid RECounts value '{}'",
                                line_number + 1,
                                fields[1]
                            );
                            continue;
                        }
                    };
                    let length = match fields[2].parse::<u32>() {
                        Ok(value) => value,
                        Err(_) => {
                            log::info!(
                                "Skipping line {}: invalid Length value '{}'",
                                line_number + 1,
                                fields[2]
                            );
                            continue;
                        }
                    };
    
                    let record = CountReRecord {
                        Contig: contig,
                        RECounts: re_counts,
                        Length: length,
                    };
                    self.records.push(record);
                }
                Err(e) => {
                    log::info!("Failed to read line {}: {}", line_number + 1, e);
                }
            }
        }
    }

    pub fn from_hashmap(&mut self, counts: IndexMap<String, u64>, 
                        chromsizes: HashMap<String, u64>) {
        for (contig, count) in counts {
            let length = chromsizes.get(&contig).unwrap().clone();
            let record = CountReRecord {
                Contig: contig,
                RECounts: count as u32,
                Length: length as u32,
            };
            self.records.push(record);
        }
    } 

    pub fn to_data(&self) -> HashMap<String, u32> {
        let mut data: HashMap<String, u32> = HashMap::new();
        for record in &self.records {
            data.insert(record.Contig.clone(), record.RECounts);
        } 

        data
    }

    pub fn to_lengths(&self) -> HashMap<String, u32> {
        let mut data: HashMap<String, u32> = HashMap::new();
        for record in &self.records {
            data.insert(record.Contig.clone(), record.Length);
        } 

        data
    }

    pub fn write(&self, output: &String) {
        let mut writer = common_writer(output);
        writer.write_all(b"#Contig\tRECounts\tLength\n").unwrap();
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(writer);
        
        
        for record in &self.records {
            wtr.serialize(record).unwrap();
        }
        log::info!("Successful written count RE file to `{}`", output);
    }
}