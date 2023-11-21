use anyhow::Result as anyResult;
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::process::exit;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use crate::core::{ common_reader, common_writer };

use crate::core::{ BaseTable, ContigPair };



#[derive(Debug, Deserialize, Serialize)]
pub struct AlleleRecord {
    pub idx1: u64,
    pub idx2: u64,
    pub contig1: String,
    pub contig2: String,
    pub length: u64,
    pub mz: u64,
    pub mz_unique: u64,
    pub similarity: f64,
}


#[derive(Debug)]
pub struct AlleleTable {
    file: String,
}

impl BaseTable for AlleleTable {
    fn new(name: &String) -> AlleleTable {
        AlleleTable { file: name.clone() }
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

impl AlleleTable {
    pub fn parse(&self) -> anyResult<csv::Reader<Box<dyn BufRead + Send>>> {
        let input = common_reader(&self.file);
        let rdr = csv::ReaderBuilder::new()
                            .flexible(true)
                            .has_headers(false)
                            .comment(Some(b'#'))
                            .delimiter(b'\t')
                            .from_reader(input);
        
        Ok(rdr)
    }

    pub fn allele_records(&self) -> anyResult<Vec<AlleleRecord>> {
        let parse_result = self.parse();
        let mut records: Vec<AlleleRecord> = Vec::new();

        match parse_result {
            Ok(mut rdr) => {
                for result in rdr.deserialize() {
                    let record: AlleleRecord = result?;
                    records.push(record);
                }
            },
            Err(e) => {
                eprintln!("Error parsing allele table: {}", e);
                exit(1);
            }
        }

        Ok(records)
    }
    
    pub fn get_allelic_contig_pairs(&self) -> HashSet<ContigPair> {
        let mut contig_pairs: HashSet<ContigPair> = HashSet::new();
        let records = self.allele_records().unwrap();
        for record in records {
            let contig1 = record.contig1;
            let contig2 = record.contig2;
            let mut contig_pair = ContigPair::new(contig1, contig2);
            contig_pair.order();
            contig_pairs.insert(contig_pair);
        }

        contig_pairs
    }

    pub fn get_allelic_record_by_contig_pairs(&self) -> HashMap<ContigPair, AlleleRecord> {
        let mut data: HashMap<ContigPair, AlleleRecord> = HashMap::new();
        let records = self.allele_records().unwrap();
        for record in records {
            let contig1 = record.contig1.clone();
            let contig2 = record.contig2.clone();
            let mut contig_pair = ContigPair::new(contig1, contig2);
            contig_pair.order();
            data.insert(contig_pair, record);
        }

        data
    }

    pub fn get_allelic_contigs(&self, method: &str) -> HashMap<String, Vec<String>> {
        let mut data: HashMap<String, Vec<String>> = HashMap::new();
        let records = self.allele_records().unwrap();
        for record in records {
            let contig1 = record.contig1;
            let contig2 = record.contig2;
            if !data.contains_key(&contig1) {
                data.insert(contig1.clone(), Vec::new());
                data.get_mut(&contig1).unwrap().push(contig2.clone());
            } else {
                if method == "fast" {
                    if data.contains_key(&contig1) {
                        continue;
                    }
                }
                data.get_mut(&contig1).unwrap().push(contig2.clone());
            }
     
            
            
        }

        data
    }

    pub fn write(&self, records: &Vec<AlleleRecord>) -> anyResult<()> {
        let writer = common_writer(&self.file);
        let mut wtr = csv::WriterBuilder::new()
                            .delimiter(b'\t')
                            .from_writer(writer);

        for record in records {
            wtr.serialize(record)?;
        }

        wtr.flush()?;
        Ok(())
    }
}

