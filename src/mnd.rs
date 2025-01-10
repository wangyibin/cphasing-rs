// juicer mnd file format

use std::borrow::Cow;
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ChromSize, ChromSizeRecord };

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct MndRecord {
    pub strand1: i32,
    pub chrom1: String,
    pub pos1: u64, 
    pub frag2: u32,
    pub strand2: i32,
    pub chrom2: String,
    pub pos2: u64,
    pub frag1: u32,
    pub mapq1: u32,
    pub cigar1: char, 
    pub sequence1: char,
    pub mapq2: u32,
    pub cigar2: char,
    pub sequence2: char,
    pub readname1: char,
    pub readname2: char,
}

impl Default for MndRecord {
    fn default() -> Self {
        Self {
            strand1: 0,
            chrom1: String::from("."),
            pos1: 0,
            frag2: 1,
            strand2: -1,
            chrom2: String::from("."),
            pos2: 0,
            frag1: 0,
            mapq1: 2,
            cigar1: '-',
            sequence1: '-',
            mapq2: 2,
            cigar2: '-',
            sequence2: '-',
            readname1: '-',
            readname2: '-',
        }
    }
}

#[derive(Debug, Clone)]
pub struct MndTable {
    file: String,
}

impl BaseTable for MndTable {
    fn new(name: &String) -> MndTable {
        MndTable { file: name.clone() }
    }

    fn file_name(&self) ->  Cow<'_, str> {
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

