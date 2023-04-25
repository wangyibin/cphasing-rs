use anyhow::Result as AnyResult;
use std::collections::HashMap;
use std::error::Error;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::io;
use std::path::Path;

use seq_io::prelude::*;
use seq_io::fastx::Reader as FastxReader;
use seq_io::parallel::read_process_fastx_records;

use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;

pub struct Fastx {
    file: String,
}

impl BaseTable for Fastx {
    fn new(name: &String) -> Fastx {
        Fastx { file: name.clone() }
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

impl Fastx {
    pub fn get_chrom_size(&self) -> AnyResult<HashMap<String, u64>>  {
        let reader = common_reader(&self.file);
        let reader = FastxReader::new(reader);
        let mut chrom_size: HashMap<String, u64> = HashMap::new();

        read_process_fastx_records(reader, 4, 2,
            |record, length| { // runs in worker
                *length = record.seq().len();
            },
            |record, length| { // runs in main thread
                chrom_size.insert(record.id().unwrap().to_owned(), *length as u64); 
                None::<()>
            }).unwrap();
    
        Ok(chrom_size)
    }
}

