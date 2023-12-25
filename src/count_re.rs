use anyhow::Result as AnyResult;
use indexmap::IndexMap;
use std::borrow::Cow;
use std::collections::HashMap;
use std::path::Path;
use serde::{ Deserialize, Serialize };

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CountReRecord {
    Contig: String, 
    RECounts: u32,
    Length: u32,
}

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
        let input = common_reader(&self.file);
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .has_headers(false)
            .from_reader(input);

        for result in rdr.deserialize() {
            let record: CountReRecord = result.unwrap();
            self.records.push(record);
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