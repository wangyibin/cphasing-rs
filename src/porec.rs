use anyhow::Result as anyResult;
use itertools::{Itertools, Combinations};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::borrow::Cow;
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ChromSize, ChromSizeRecord };
use crate::pairs::{ PairRecord, PairHeader };
use crate::paf::PAFLine;



#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PoreCRecord {
    pub read_idx: u64,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub query_strand: char,
    pub target: String, 
    pub target_start: u64, 
    pub target_end: u64,
    pub mapq: u8,
    pub identity: f32, 
    pub filter_reason: String,
} 



impl PoreCRecord {
    pub fn from_paf_record(record: PAFLine, read_idx: u64,
                            identity: f32, filter_reason: String) -> PoreCRecord{
        PoreCRecord {
            read_idx: read_idx,
            query_length: record.query_length,
            query_start: record.query_start,
            query_end: record.query_end,
            query_strand: record.query_strand,
            target: record.target,
            target_start: record.target_start,
            target_end: record.target_end,
            mapq: record.mapq,
            identity: identity,
            filter_reason: filter_reason
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "Read index: {}\nQuery length: {}\nQuery start: {}\nQuery end: {}\nQuery strand: {}\nTarget: {}\nTarget start: {}\nTarget end: {}\nMapping quality: {}\nIdentity: {}\nFilter reason: {}\n",
            self.read_idx,
            self.query_length,
            self.query_start,
            self.query_end,
            self.query_strand,
            self.target,
            self.target_start,
            self.target_end,
            self.mapq,
            self.identity,
            self.filter_reason,
        )
    }
}

impl PartialOrd for PoreCRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.query_start.partial_cmp(&other.query_start)
    }
}

impl PartialEq for PoreCRecord {
    fn eq(&self, other: &Self) -> bool {
        self.query_start == other.query_start
    }
}

#[derive(Debug)]
pub struct Concatemer {
    pub records: Vec<PoreCRecord>,
}

impl Concatemer {

    pub fn new() -> Concatemer{
        Concatemer {
            records: Vec::new(),
        }
    }

    pub fn push(&mut self, pcr: PoreCRecord) {
        self.records.push(pcr);
    }

    pub fn clear(&mut self) {
        self.records.clear();
    }

    pub fn sort(&mut self) {
        // self.records.sort_by(| a, b | a.partial_cmp(&b).unwrap());
        self.records.sort_unstable_by_key(|x| (x.target.clone(), x.target_start));
    }

    pub fn decompose(&mut self) -> Combinations<std::vec::IntoIter<PoreCRecord>> {
        let r = self.records.clone();
        r.into_iter().combinations(2)
        
    }

}

#[derive(Debug, Clone)]
pub struct ConcatemerSummary {
    pub summary: HashMap<u32, u64>,
}

impl ConcatemerSummary {
    pub fn new() -> ConcatemerSummary {
        ConcatemerSummary { summary: HashMap::<u32, u64>::new() }
    }
    
    pub fn count(&mut self, concatemer: &Concatemer) {

        let concatemer_count: u32 = concatemer.records.len().try_into().unwrap();
        *self.summary.entry(concatemer_count).or_insert(0) += 1;

    }

    pub fn to_string(&self) -> String {
        let mut vec = self.summary.iter().collect::<Vec<_>>();
        vec.sort_by_key(|(key, _)| *key);

        let string = vec.iter()
                        .map(|(key, value)| format!("{}\t{}", key, value))
                        .collect::<Vec<_>>()
                        .join("\n");

        string
    }

    pub fn save(&self, output: &String) {
        let mut wtr = common_writer(output);    

        let result: String = self.to_string();

        wtr.write_all(result.as_bytes()).unwrap();
        log::info!("Successful output summary of concatemer `{}`", output);
    }
    // pub fn collapse(&self) -> HashMap<&str, u64> {
        
    // }
}

#[derive(Debug)]
pub struct PoreCTable {
    file: String,
}

impl BaseTable for PoreCTable {
    fn new(name: &String) -> PoreCTable {
        PoreCTable { file: name.clone() }
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


impl PoreCTable {
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

    pub fn to_pairs(&self, chromsizes: &String, output: &String) -> Result<(), Box<dyn Error>> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
        };

        let chromsizes: ChromSize = ChromSize::new(chromsizes);
        let chromsizes_data: Vec<ChromSizeRecord> = chromsizes.to_vec().unwrap();

        let mut ph: PairHeader = PairHeader::new();
        ph.from_chromsizes(chromsizes_data);


        let mut writer = common_writer(output);
        writer.write_all(ph.to_string().as_bytes()).unwrap();

        let mut wtr = csv::WriterBuilder::new()
                            .has_headers(false)
                            .delimiter(b'\t')
                            .from_writer(writer);

        let mut concatemer: Concatemer = Concatemer::new();  
        let mut concatemer_summary: ConcatemerSummary = ConcatemerSummary::new();

        let mut old_read_idx: u64 = 0; 
        let mut flag: bool = false; 
        let mut read_id: u64 = 0;
        
        for (i, line) in rdr.deserialize().enumerate() {
            let record: PoreCRecord = match line {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}", i + 1);
                    continue
                },
            };
            if record.read_idx != old_read_idx && flag == true {
                concatemer.sort();
                
                concatemer_summary.count(&concatemer);

                for pair in concatemer.decompose() {
                    wtr.serialize(PairRecord::from_pore_c_pair(pair, read_id));
                    read_id += 1;
                }
                concatemer.clear();
            }
            flag = true;
            old_read_idx = record.read_idx;
            concatemer.push(record);
           
        }

        // process last concatemer
        concatemer.sort();
        for pair in concatemer.decompose() {
            wtr.serialize(PairRecord::from_pore_c_pair(pair, read_id));
            read_id += 1;
        }

        log::info!("Successful output pairs `{}`", output);

        concatemer_summary.save(&format!("{}.concatemer.summary", self.prefix()));
        Ok(())
    }
}