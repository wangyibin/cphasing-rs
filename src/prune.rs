use anyhow::Result as anyResult;
use std::borrow::Cow;
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };

pub struct PruneRecord {
    pub contig1: String,
    pub contig2: String,
    pub mz1: u64,
    pub mz2: u64,
    pub mzShared: u64,
    pub similarity: f64,
    pub allele_type: u8,
}

#[derive(Debug)]
pub struct PruneTable {
    file: String,
}

impl BaseTable for PruneTable {
    fn new(name: &String) -> PruneTable {
        PruneTable { file: name.clone() }
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

impl PruneTable {
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

    pub fn contig_pairs(&self) -> anyResult<Vec<ContigPair>> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
        };

        let mut contig_pairs: Vec<ContigPair> = Vec::new();
        for result in rdr.records() {
            let record = result?;
            let contig1 = record[0].to_string();
            let contig2 = record[1].to_string();
           

            let contig_pair = ContigPair::new(contig1, contig2);
            contig_pairs.push(contig_pair);
        }

        Ok(contig_pairs)
    }

}
