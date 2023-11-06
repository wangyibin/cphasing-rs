use anyhow::Result as anyResult;
use std::borrow::Cow;
use std::collections::HashSet;
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ChromSizeRecord, ContigPair};
use crate::mnd::MndRecord;
use crate::porec::PoreCRecord;



#[derive(Debug, Clone)]
pub struct PairHeader {
    pub chromsizes: Vec<ChromSizeRecord>,
    pub header: Vec<String>,
}

impl PairHeader {
    pub fn new() -> PairHeader {
        let chromsizes: Vec<ChromSizeRecord> = Vec::new();
        let header: Vec<String> = Vec::new();

        PairHeader { chromsizes: chromsizes, header: header }
    }

    pub fn from_pairs(&mut self, name: &String) {
        let reader = common_reader(name);
        for line in reader.lines() {
            if let Ok(line) = line {
                if line.starts_with("#") {
                    if &line[0..10] == "#chromsize" {
                        let s: String = line.replace("#chromsize: ", "");
                        let s: Vec<&str> = s.split(" ").collect();
                        
                        let size: u64 = s[1].clone().parse::<u64>().unwrap();
                        let chrom: String = s[0].to_string();
                        let c: ChromSizeRecord = ChromSizeRecord { chrom: chrom, size: size };
                        self.chromsizes.push(c);

                    } else if &line[0..8] == "#columns" {
                        
                        let s: String = line.replace("#columns: ", "");
                        let s: Vec<&str> = s.split(" ").collect();
                        
                        let s: Vec<String> = s.iter().map(|&x| x.to_owned()).collect();
                        self.header = s;

                    } else {
                        continue
                    }
                } else {
                    break
                }
            }
        }
    }

    pub fn from_chromsizes(&mut self, chromsizes: Vec<ChromSizeRecord>) {
        let pair_header: Vec<String> = vec!["readID".to_string(), 
                                    "chrom1".to_string(), 
                                    "pos1".to_string(), 
                                    "chrom2".to_string(), 
                                    "pos2".to_string(), 
                                    "strand1".to_string(),
                                    "strand2".to_string()
                                    ];
        self.chromsizes = chromsizes;
        self.header = pair_header;
    }

    pub fn to_string(&self) -> String {
        let mut result = String::new();
        
        result.push_str("## pairs format 1.0\n");
        result.push_str("#shape: upper triangle\n");

        for chromsize in &self.chromsizes {
            result.push_str(&format!("#chromsize: {} {}\n", 
                                chromsize.chrom, chromsize.size));
        }
        result.push_str(&format!("#columns: {}\n", &self.header.join(" ")));

        result
    }

}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PairRecord {
    pub readID: u64,
    pub chrom1: String,
    pub pos1: u64,
    pub chrom2: String,
    pub pos2: u64,
    pub strand1: char,
    pub strand2: char,
}

impl PairRecord {
    pub fn from_pore_c_pair(pair: Vec<PoreCRecord>, readID: u64) -> PairRecord {
        let (pair1, pair2) = (pair[0].clone(), pair[1].clone());
        
        let pos1: u64 = (pair1.target_end + pair1.target_start) / 2;
        let pos2: u64 = (pair2.target_end + pair2.target_start) / 2;

        PairRecord {
            readID: readID,
            chrom1: pair1.target,
            pos1: pos1, 
            chrom2: pair2.target,
            pos2: pos2,
            strand1: pair1.query_strand,
            strand2: pair2.query_strand,
        }
    }

}

#[derive(Debug, Clone)]
pub struct Pairs {
    file: String,
    header: PairHeader,
}

impl BaseTable for Pairs {
    fn new(name: &String) -> Pairs {
        Pairs { file: name.clone(), header: PairHeader::new() }
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

impl Pairs {
    pub fn parse(&mut self) -> anyResult<csv::Reader<Box<dyn BufRead + Send>>> {
        let input = common_reader(&self.file);
        self.header = PairHeader::new();
        self.header.from_pairs(&self.file);

        let rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .from_reader(input);

        Ok(rdr)
    }

    pub fn remove_by_contig_pairs(&mut self, contigs: HashSet<ContigPair>, output: &String) -> anyResult<()> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let ph = self.header.clone();
        let mut writer = common_writer(output);
        writer.write_all(ph.to_string().as_bytes());
        let mut wtr = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(writer);

        for result in rdr.deserialize() {
            if result.is_err() {
                println!("{:?}", result);
                continue
            }
            let record: PairRecord = result?;
            let cp = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
            if !contigs.contains(&cp) {
                wtr.serialize(record)?;
            }
        }

        Ok(())
    }

    pub fn to_mnd(&mut self, output: &String) -> anyResult<()>{
        let parse_result = self.parse();
        let mnd_defualt = MndRecord::default();

        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let wtr = common_writer(output);
        let mut wtr = csv::WriterBuilder::new()
                .has_headers(false)
                .delimiter(b' ')
                .from_writer(wtr);

        // for result in rdr.deserialize() {
        //     if result.is_err() {
        //         println!("{:?}", result);
        //         continue
        //     }
        //     let record: PairRecord = result?;
        //     let strand1 = if record.strand1 == '+' { 0 } else { -1 };
        //     let strand2 = if record.strand2 == '+' { 0 } else { -1 };
        //     let mnd_record = MndRecord {
        //         strand1: strand1,
        //         chrom1: record.chrom1,
        //         pos1: record.pos1,
        //         frag1: mnd_defualt.frag1,
        //         strand2: strand2,
        //         chrom2: record.chrom2,
        //         pos2: record.pos2,
        //         frag2: mnd_defualt.frag2,
        //         mapq1: mnd_defualt.mapq1,
        //         cigar1: mnd_defualt.cigar1,
        //         sequence1: mnd_defualt.sequence1,
        //         mapq2: mnd_defualt.mapq2,
        //         cigar2: mnd_defualt.cigar2,
        //         sequence2: mnd_defualt.sequence2,
        //         readname1: mnd_defualt.readname1,
        //         readname2: mnd_defualt.readname2,
        //     };
        
        for (idx, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    
                    let record = PairRecord {
                        readID: idx as u64,
                        chrom1: record[1].to_string(),
                        pos1: record[2].parse::<u64>().unwrap(),
                        chrom2: record[3].to_string(),
                        pos2: record[4].parse::<u64>().unwrap(),
                        strand1: record[5].parse::<char>().unwrap(),
                        strand2: record[6].parse::<char>().unwrap(),
                    };

                    let strand1 = if record.strand1 == '+' { 0 } else { -1 };
                    let strand2 = if record.strand2 == '+' { 0 } else { -1 };
                    let mnd_record = MndRecord {
                        strand1: strand1,
                        chrom1: record.chrom1,
                        pos1: record.pos1,
                        frag1: mnd_defualt.frag1,
                        strand2: strand2,
                        chrom2: record.chrom2,
                        pos2: record.pos2,
                        frag2: mnd_defualt.frag2,
                        mapq1: mnd_defualt.mapq1,
                        cigar1: mnd_defualt.cigar1,
                        sequence1: mnd_defualt.sequence1,
                        mapq2: mnd_defualt.mapq2,
                        cigar2: mnd_defualt.cigar2,
                        sequence2: mnd_defualt.sequence2,
                        readname1: mnd_defualt.readname1,
                        readname2: mnd_defualt.readname2,
                    };
                    wtr.serialize(mnd_record)?;
                }
                Err(e) => {
                    println!("{:?}", e);
                    continue
            }
        }

            // wtr.serialize(mnd_record)?;
        }

        wtr.flush()?;

        log::info!("Successful output mnd file `{}`", output);
        
        Ok(())
    }

    pub fn to_bam(&self, output: &String) {
        
        let mut wtr = bam::Writer::from_path(output, &bam::Header::new()).unwrap();
    }
}


