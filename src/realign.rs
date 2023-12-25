#[warn(unused_assignments)]
use anyhow::Result as anyResult;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::path::Path;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;

#[derive(Debug, Clone)]
pub struct PAFLine {
    pub query: String,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub query_strand: char,
    pub target: String, 
    pub target_length: u64,
    pub target_start: u64, 
    pub target_end: u64,
    pub match_n: u32,
    pub alignment_length: u32,
    pub mapq: u8,
    pub cigar: String,
    pub tags: Vec<String>,
}

impl PAFLine {
    fn new(fields: Vec<String>) -> Self {
        PAFLine {
            query: fields[0].clone(),
            query_length: fields[1].parse::<u32>().unwrap(),
            query_start: fields[2].parse::<u32>().unwrap(),
            query_end: fields[3].parse::<u32>().unwrap(),
            query_strand: fields[4].parse::<char>().unwrap(),
            target: fields[5].clone(),
            target_length: fields[6].parse::<u64>().unwrap(),
            target_start: fields[7].parse::<u64>().unwrap(),
            target_end: fields[8].parse::<u64>().unwrap(),
            match_n: fields[9].parse::<u32>().unwrap(),
            alignment_length: fields[10].parse::<u32>().unwrap(),
            mapq: fields[11].parse::<u8>().unwrap(),
            cigar: fields[12].clone(),
            tags: fields[13..].to_vec(),
        }
    }

    fn get_tag(&self, tag: &str) -> String {
        let mut value = String::new();
        for t in &self.tags {
            if t.starts_with(tag) {
                value = t.to_string();
                break;
            }
        }
        value
    }

    fn is_secondary(&self) -> bool {
        let tag = self.get_tag("tp:A:");
        if tag == "tp:A:S" {
            return true;
        } else {
            return false;
        }
    }

    fn is_primary(&self) -> bool {
        let tag = self.get_tag("tp:A:");
        if tag == "tp:A:P" {
            return true;
        } else {
            return false;
        }
    }

    fn from(paf_line: &PAFLine) -> Self {
        PAFLine {
            query: paf_line.query.clone(),
            query_length: paf_line.query_length,
            query_start: paf_line.query_start,
            query_end: paf_line.query_end,
            query_strand: paf_line.query_strand,
            target: paf_line.target.clone(),
            target_length: paf_line.target_length,
            target_start: paf_line.target_start,
            target_end: paf_line.target_end,
            match_n: paf_line.match_n,
            alignment_length: paf_line.alignment_length,
            mapq: paf_line.mapq,
            cigar: paf_line.cigar.clone(),
            tags: paf_line.tags.clone(),
        }
    }
    
    fn to_string(&self) -> String {
        let mut fields: Vec<String> = Vec::new();
        fields.push(self.query.clone());
        fields.push(self.query_length.to_string());
        fields.push(self.query_start.to_string());
        fields.push(self.query_end.to_string());
        fields.push(self.query_strand.to_string());
        fields.push(self.target.clone());
        fields.push(self.target_length.to_string());
        fields.push(self.target_start.to_string());
        fields.push(self.target_end.to_string());
        fields.push(self.match_n.to_string());
        fields.push(self.alignment_length.to_string());
        fields.push(self.mapq.to_string());
        fields.push(self.cigar.clone());
        for t in &self.tags {
            fields.push(t.to_string());
        }
        fields.join("\t")
    }
}


#[derive(Debug)]
pub struct PAFTable {
    file: String,
}


impl BaseTable for PAFTable {
    fn new(name: &String) -> PAFTable {
        PAFTable { file: name.clone() }
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

impl PAFTable {
    pub fn parse(&self) -> anyResult<Box<dyn BufRead + Send + 'static>> {
        let mut reader = common_reader(&self.file);
        Ok(reader)
    }

}

#[derive(Debug, Clone)]
pub struct PAFReadUnit {
    data: Vec<PAFLine>,
}

impl PAFReadUnit {
    fn new() -> Self {
        PAFReadUnit {
            data: Vec::new(),
        }
    }

    fn clear (&mut self) {
        self.data.clear();
    }
}

pub struct PAFAlignmentUnit {
    Primary: Vec<PAFLine>,
    Secondary: Vec<Vec<PAFLine>>,
}

impl PAFAlignmentUnit {
    fn new() -> Self {
        PAFAlignmentUnit {
            Primary: Vec::new(),
            Secondary: Vec::new(),
        }
    }

    fn clear (&mut self) {
        self.Primary.clear();
        self.Secondary.clear();
    }

    fn add_primary(&mut self, record: PAFLine) {
        self.Primary.push(record);
        self.Secondary.push(Vec::new());
    }

    fn add_secondary(&mut self, record: PAFLine, idx: usize) {
        // idx: index of primary alignment
        self.Secondary[idx].push(record);
    }

    fn is_empty(&self) -> bool {
        self.Primary.is_empty()
    }

    fn read_id(&self) -> String {
        self.Primary[0].query.clone()
    }

    fn rescue(&mut self, mapq: u8) {
        // count map quality >= 1 {"target1": 2, "target2": 1}
        
        let mut high_mapq: HashMap<String, u32> = HashMap::new();
        let mut high_high_mapq: HashMap<String, u32> = HashMap::new();
        for r in &self.Primary {
            let target = r.target.clone();
            if r.mapq >= mapq {
                let count = high_mapq.entry(target.clone()).or_insert(0);
                *count += 1;
            }

            if r.mapq > 1 {
                let count = high_high_mapq.entry(target).or_insert(0);
                *count += 1;
            }
        }
        

        if high_mapq.len() == 0 {
            return;
        }

        {
        'outer: for (mut p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
            if p.mapq > 1 {
                continue;
            }

            let mut res: HashSet<String> = HashSet::new();
            let mut res_record_idx: HashMap<String, usize> = HashMap::new();
            let target = p.target.clone();
    
            if high_high_mapq.contains_key(&target) {
                res_record_idx.insert(target.clone(), 0);
                res.insert(target);
            }

            for (j, r) in s.iter().enumerate() {
                let target = r.target.clone();
                if high_high_mapq.contains_key(&target) {
                    res_record_idx.insert(target.clone(), j+1);
                    res.insert(target);
                }
            }

            let mut max_target = String::new();
            if res.len() == 0 {
                continue;
            } else if res.len() == 1 {
                let target = res.iter().next().unwrap();
                max_target = target.clone();
            } else {
                // max in res 
                let mut max = 0;
                for target in res {
                    let count = high_mapq.get(&target).unwrap();
                    if *count > max {
                        max = *count;
                        max_target = target;
                    }
                }
        
            }

            let record_idx = res_record_idx.get(&max_target).unwrap();

            high_mapq.entry(max_target).and_modify(|x| *x += 1);
            if record_idx == &0 {
                p.mapq = 1;
            } else {
                let mut new_record = PAFLine::from(&s[*record_idx - 1]);
                // println!("{:?}", p);
                // println!("{:?}", new_record);
                new_record.mapq = 1;
                new_record.tags = new_record.tags.iter().map(|x| {
                    if x.starts_with("tp:A:") {
                        "tp:A:P".to_string()
                    } else {
                        x.to_string()
                    }
                }).collect();
                *p = new_record;
            
            }
        }

        
        'outer: for (mut p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
            if p.mapq >= mapq {
                continue;
            }

            let mut res: HashSet<String> = HashSet::new();
            let mut res_record_idx: HashMap<String, usize> = HashMap::new();
            let target = p.target.clone();
    
            if high_mapq.contains_key(&target) {
                res_record_idx.insert(target.clone(), 0);
                res.insert(target);
            }

            for (j, r) in s.iter().enumerate() {
                let target = r.target.clone();
                if high_mapq.contains_key(&target) {
                    res_record_idx.insert(target.clone(), j+1);
                    res.insert(target);
                }
            }

            let mut max_target = String::new();
            if res.len() == 0 {
                continue;
            } else if res.len() == 1 {
                let target = res.iter().next().unwrap();
                max_target = target.clone();
            } else {
                // max in res 
                let mut max = 0;
                for target in res {
                    let count = high_mapq.get(&target).unwrap();
                    if *count > max {
                        max = *count;
                        max_target = target;
                    }
                }
            }
            
            let record_idx = res_record_idx.get(&max_target).unwrap();
            if record_idx == &0 {
                p.mapq = 1;
            } else {
                let mut new_record = PAFLine::from(&s[*record_idx - 1]);
                new_record.mapq = 1;
                new_record.tags = new_record.tags.iter().map(|x| {
                    if x.starts_with("tp:A:") {
                        "tp:A:P".to_string()
                    } else {
                        x.to_string()
                    }
                }).collect();
                *p = new_record;
            }

    }
    }
    }

}


pub fn parse_paf_read_unit(read_unit: &PAFReadUnit) -> PAFAlignmentUnit {
    let mut idx: u64 = 0;
    let mut au = PAFAlignmentUnit::new();
    for r in &read_unit.data {
        let read_id = r.query.clone();
        

        if !r.is_secondary() {
            au.add_primary(r.clone());
            idx += 1;
        } else {
            if idx == 0 {
                log::warn!("Secondary alignment `{:?}` could not found primary, \
                            skipped.\
                            The input bam should be sorted by read name.", read_id);
                continue;
            }
            let idx2 = idx - 1;
            au.add_secondary(r.clone(), idx2 as usize);
        }
        
    }
    au
}


pub fn read_paf(input_paf: &String, mapq: u8, output: &String) {
    let mut paf = PAFTable::new(input_paf);
    let parse_result = paf.parse();
    let mut rdr = match parse_result {
        Ok(v) => v,
        Err(error) => panic!("Error: Could not parse input file: {:?}", paf.file_name()),
    };

    let mut total_reads: u64 = 0;
    let mut total_alignments: u64 = 0;
    let mut total_unmapped: u64 = 0;
    let mut old_read_id = String::from("");

    let wtr = common_writer(output);
    let mut writer = BufWriter::new(wtr);
    
    let mut read_unit = PAFReadUnit::new();

    for line in rdr.lines() {
        let fields: Vec<String> = line.unwrap().split('\t').map(|x| x.to_string()).collect();
        
        assert!(fields.len() > 12, "Error: PAF file should have at least 12 columns");
        let paf_line: PAFLine = PAFLine::new(fields);

        if paf_line.target == "*" {
            total_unmapped += 1;
            total_reads += 1;
            continue;
        }

        let read_id = paf_line.query.clone();

        if old_read_id != paf_line.query {
            if old_read_id != "" {
                let mut au = parse_paf_read_unit(&read_unit);
                au.rescue(mapq);
                for r in au.Primary {
                    
                    writer.write(r.to_string().as_bytes()).unwrap();
                    writer.write(b"\n").unwrap();
                }
            }

            total_reads += 1;
            read_unit.clear();
            read_unit.data.push(paf_line);
        } else {
            read_unit.data.push(paf_line);
        }

        old_read_id = read_id;

        

    }
}




#[derive(Clone, Debug)]
pub struct ReadUnit {
    data: Vec<Record>,
}

impl ReadUnit {
    fn new() -> Self {
        ReadUnit {
            data: Vec::new(),
        }
    }

    fn clear (&mut self) {
        self.data.clear();
    }
}


#[derive(Debug)]
pub struct AlignmentUnit {
    Primary: Vec<Record>,
    Secondary: Vec<Vec<Record>>,
}


impl AlignmentUnit {
    fn new() -> Self {
        AlignmentUnit {
            Primary: Vec::new(),
            Secondary: Vec::new(),
        }
    }

    fn clear (&mut self) {
        self.Primary.clear();
        self.Secondary.clear();
    }

    fn add_primary(&mut self, record: Record) {
        self.Primary.push(record);
        self.Secondary.push(Vec::new());
    }

    fn add_secondary(&mut self, record: Record, idx: usize) {
        // idx: index of primary alignment
        self.Secondary[idx].push(record);
    }

    fn is_empty(&self) -> bool {
        self.Primary.is_empty()
    }

    fn read_id(&self) -> String {
        String::from_utf8(self.Primary[0].qname().to_vec())
                                .unwrap().to_string()
    }

    fn rescue(&mut self, mapq: u8) {
        
        let mut high_mapq: HashMap<u64, u32> = HashMap::new();
        let mut high_high_mapq: HashMap<u64, u32> = HashMap::new();
        for r in &self.Primary {
            let target: u64 = r.tid().try_into().unwrap();
            if r.mapq() >= mapq {
                let count = high_mapq.entry(target).or_insert(0);
                *count += 1;
            } 

            if r.mapq() > 1 {
                let count = high_high_mapq.entry(target).or_insert(0);
                *count += 1;
            }
        }
        

        if high_mapq.len() == 0 {
            return;
        }
      
        'outer: for (mut p, s) in self.Primary.iter_mut().zip(self.Secondary.iter()) {
            if p.mapq() >= mapq {
                continue;
            }

            let mut res: HashSet<u64> = HashSet::new();
            let mut res_record_idx: HashMap<u64, usize> = HashMap::new();
            let target: u64 = p.tid().try_into().unwrap();
            if high_mapq.contains_key(&target) {
                res_record_idx.insert(target, 0);
                res.insert(target);


            }
            
            'inger: for (j, r) in s.iter().enumerate() {
                let target: u64 = r.tid().try_into().unwrap();
                if high_mapq.contains_key(&target) {
                    res_record_idx.insert(target, j);
                    res.insert(target);
                }
            }
            
            let mut max_target = 0;
            if res.len() == 0 {
                continue;
            } else if res.len() == 1 {
                let target = res.iter().next().unwrap();
                max_target = *target;
            } else {
                // max in res 
                let mut max = 0;
                for target in res {
                    let count = high_mapq.get(&target).unwrap();
                    if *count > max {
                        max = *count;
                        max_target = target;
                    }
                }
            
            }
            let record_idx = res_record_idx.get(&max_target).unwrap();
            
            if record_idx == &0 {
                p.set_tid(max_target.try_into().unwrap());
                p.set_mapq(1);
            } else {
                let flag = p.flags();
                let r = &s[*record_idx - 1];
                let mut new_record = Record::from(r.clone());
                p = &mut new_record;
                p.set_mapq(1);
                p.set_flags(flag);
            } 
            
        }
    
    
    }
 

}

fn parse_read_unit(read_unit: &ReadUnit) -> AlignmentUnit {
    let mut idx: u64 = 0;
    let mut au = AlignmentUnit::new();
    for r in &read_unit.data {
        let read_id = String::from_utf8(r.qname().to_vec())
                                .unwrap().to_string();
        

        if !r.is_secondary() {
            au.add_primary(r.clone());
            idx += 1;
        } else {
            if idx == 0 {
                log::warn!("Secondary alignment `{:?}` could not found primary, \
                            skipped.\
                            The input bam should be sorted by read name.", read_id);
                continue;
            }
            let idx2 = idx - 1;
            au.add_secondary(r.clone(), idx2 as usize);
        }
        
    }
    au
}


pub fn read_bam(input_bam: &String, mapq: u8, output: &String) {
    let mut bam = Reader::from_path(input_bam).unwrap();
    let bam_header = Header::from_template(bam.header());
    let bam_header = HeaderView::from_header(&bam_header);

    let mut total_reads: u64 = 0;
    let mut total_alignments: u64 = 0;
    let mut total_unmapped: u64 = 0;
    let mut old_read_id = String::from("");

    let mut read_unit = ReadUnit::new();

    let header = Header::from_template(&bam_header);
    let mut writer = Writer::from_path(output, &header, bam::Format::Bam).unwrap();

    for r in bam.records() {
        let record = r.unwrap();
        let read_id = String::from_utf8(record.qname().to_vec())
                                .unwrap().to_string();

        if record.is_unmapped() {
            total_unmapped += 1;
            total_reads += 1;
            continue;
        }

        if old_read_id != read_id {
            if old_read_id != "" {
                let mut au = parse_read_unit(&read_unit);
                au.rescue(mapq);
                for r in au.Primary {
                    writer.write(&r).unwrap();
                }
            }


            total_reads += 1;
            read_unit.clear();
            read_unit.data.push(record);
        } else {
            read_unit.data.push(record);
        }
      
        old_read_id = read_id;
    

    }

}




