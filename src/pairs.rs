use anyhow::Result as anyResult;
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::sync::{ Arc, Mutex };
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ChromSizeRecord, ContigPair};
use crate::mnd::MndRecord;
use crate::porec::PoreCRecord;
use crate::contacts::{ Contacts, ContactRecord };



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
        writer.write_all(ph.to_string().as_bytes()).unwrap();
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
                    eprintln!("{:?}", e);
                    continue
            }
        }

            // wtr.serialize(mnd_record)?;
        }

        wtr.flush()?;

        log::info!("Successful output mnd file `{}`", output);
        
        Ok(())
    }

    // convert pairs to pseudo bam file 
    pub fn to_bam(&mut self, output: &String) {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };
        let pair_header = self.header.clone();
        let mut bam_header = Header::new();
        let mut record = HeaderRecord::new(b"HD");
        record.push_tag(b"VN", &"1.6") ;
        record.push_tag(b"SO", &"coordinate");    
        bam_header.push_record(&record);

        for chromsize in &pair_header.chromsizes {
            let mut record = HeaderRecord::new(b"SQ");
            record.push_tag(b"SN", &chromsize.chrom);
            record.push_tag(b"LN", &chromsize.size);
            bam_header.push_record(&record);
        }
        let mut record = HeaderRecord::new(b"PG");
        record.push_tag(b"ID", &"pairs2bam");
        record.push_tag(b"PN", &"cphasing-rs");
        let version = env!("CARGO_PKG_VERSION");
        record.push_tag(b"VN", &version);
        record.push_tag(b"CL", &"cphasing-rs pairs2bam");
        bam_header.push_record(&record);
        
        let bam_header_view = HeaderView::from_header(&bam_header);

        let mut wtr = Writer::from_path(output, &bam_header, bam::Format::Bam).unwrap();

        for (id, record) in rdr.records().enumerate() {
            match record {
                Ok(record) => {
                    let flag1 = match record[5].parse::<char>().unwrap() {
                        '+' => 0,
                        '-' => 16,
                        _ => 0 
                    };

                    let flag2 = match record[6].parse::<char>().unwrap() {
                        '+' => 0,
                        '-' => 16,
                        _ => 0 
                    };
                    {
                        let mut bam_record1 = Record::new();
                        let cigar = CigarString(vec![Cigar::Match(150)]);
                        let cigar: Option<&CigarString> = Some(&cigar);
                        bam_record1.set(record[0].as_bytes(), cigar, &[], &[]);
                        bam_record1.set_tid(bam_header_view.tid(record[1].as_bytes()).unwrap().try_into().unwrap());
                        bam_record1.set_pos(record[2].parse::<i64>().unwrap());
                        bam_record1.set_mapq(60);
                        bam_record1.set_mtid(bam_header_view.tid(record[3].as_bytes()).unwrap().try_into().unwrap());
                        bam_record1.set_mpos(record[4].parse::<i64>().unwrap());
                        bam_record1.set_insert_size(0);
                        // bam_record1.set_qname(record[0].as_bytes());
                        bam_record1.set_flags(flag1);
                        bam_record1.push_aux(b"RG", Aux::String("1")).unwrap();
                        bam_record1.push_aux(b"AS", Aux::U8(0)).unwrap();
                        bam_record1.push_aux(b"NM", Aux::U8(0)).unwrap();
                        bam_record1.push_aux(b"MQ", Aux::U8(60)).unwrap();
                    
                        wtr.write(&bam_record1).unwrap();
                    }
                    
                    {
                        let mut bam_record2 = Record::new();
                        let cigar = CigarString(vec![Cigar::Match(150)]);
                        let cigar: Option<&CigarString> = Some(&cigar);
                        bam_record2.set(record[0].as_bytes(), cigar, &[], &[]);
                        bam_record2.set_tid(bam_header_view.tid(record[3].as_bytes()).unwrap().try_into().unwrap());
                        bam_record2.set_pos(record[4].parse::<i64>().unwrap());
                        bam_record2.set_mapq(60);
                        bam_record2.set_mtid(bam_header_view.tid(record[1].as_bytes()).unwrap().try_into().unwrap());
                        bam_record2.set_mpos(record[2].parse::<i64>().unwrap());
                        bam_record2.set_insert_size(0);
                        // bam_record2.set_qname(record[0].as_bytes());
                        bam_record2.set_flags(flag2);
                        bam_record2.push_aux(b"RG", Aux::String("1")).unwrap();
                        bam_record2.push_aux(b"AS", Aux::U8(0)).unwrap();
                        bam_record2.push_aux(b"NM", Aux::U8(0)).unwrap();
                        bam_record2.push_aux(b"MQ", Aux::U8(60)).unwrap();

                        wtr.write(&bam_record2).unwrap();
                    }

                }
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue
                }
            }
                
            }
            log::info!("Successful output bam file `{}`", output);

    }  

    pub fn to_contacts(&mut self, max_contacts: u32) -> anyResult<Contacts> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let mut contacts = Contacts::new(&format!("{}.pixels", self.prefix()).to_string());
        
        let mut contact_hash: HashMap<ContigPair,u32>  = HashMap::new();
        for result in rdr.deserialize() {
            let record: PairRecord = result?;
            let mut cp = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
            cp.order();
            if !contact_hash.contains_key(&cp) {
                contact_hash.insert(cp, 1);
            } else {
                let count = contact_hash.get(&cp).unwrap() + 1;
                contact_hash.insert(cp, count);
            }
            
        }
        
        
        // let contact_hash: Arc<Mutex<HashMap<ContigPair, u32>>> = Arc::new(Mutex::new(HashMap::new()));

        // rdr.deserialize::<PairRecord>()
        //     .par_bridge()
        //     .try_for_each(|result: Result<PairRecord, csv::Error>| -> Result<(), csv::Error> {
        //         let record = result?;
        //         let mut cp = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
        //         cp.order();

        //         let mut hash = contact_hash.lock().unwrap();
        //         if let Some(count) = hash.get_mut(&cp) {
        //             *count += 1;
        //         } else {
        //             hash.insert(cp, 1);
        //         }

        //         Ok(())
        // })?;
        
        let mut contact_records: Vec<ContactRecord> = contact_hash.par_iter(
            ).map(|(cp, count)| {
                if count >= &max_contacts{
                    let record = ContactRecord {
                        chrom1: cp.Contig1.to_owned(),
                        chrom2: cp.Contig2.to_owned(),
                        count: *count,
                    };
                    record
                } else {
                    let record = ContactRecord {
                        chrom1: "".to_string(),
                        chrom2: "".to_string(),
                        count: 0,
                    };
                    record
                }
                
            }).collect();
        
        contact_records.retain(|x| x.is_some());
        
        // let contact_records: Vec<ContactRecord> = contact_hash
        //     .lock()
        //     .unwrap()
        //     .iter()
        //     .map(|(cp, count)| {
        //         ContactRecord {
        //             chrom1: cp.Contig1.to_owned(),
        //             chrom2: cp.Contig2.to_owned(),
        //             count: *count,
        //         }
        //     })
        //     .collect();
        contacts.records = contact_records;
        Ok(contacts)

    }   

    pub fn to_clm(&mut self, max_contacts: u32, output: &String) {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(r) => r,
            Err(e) => panic!("Error: Could not parse input file: {:?}", self.file_name()),
        };

        let pair_header = self.header.clone();
        let chromsizes: Vec<ChromSizeRecord> = pair_header.chromsizes.clone();

        let mut wtr = common_writer(output);

        let mut data: HashMap<ContigPair, Vec<Vec<u64>>> = HashMap::new();
        for result in rdr.deserialize() {
            let record: PairRecord = result.unwrap();
            let mut contig_pair = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
            contig_pair.order();
            if let Some(vec) = data.get_mut(&contig_pair) {
                if record.chrom1 > record.chrom2 {
                    vec.push(vec![record.pos2, record.pos1]);
                } else {
                    vec.push(vec![record.pos1, record.pos2]);
                }
            } else {
                let mut vec = Vec::new();
                if record.chrom1 > record.chrom2 {
                    vec.push(vec![record.pos2, record.pos1]);
                } else {
                    vec.push(vec![record.pos1, record.pos2]);
                }
                data.insert(contig_pair, vec);
            }
        }

        let result = data.par_iter(
            ).filter_map(
                |(cp, vec)| {
                    if vec.len() >= max_contacts as usize {
                        Some((cp, vec))
                    } else {
                        None
                    }
                }
            ).map(|(cp, vec)|{
                
                let length1 = chromsizes.iter().find(|&x| x.chrom == cp.Contig1).unwrap().size;
                let length2 = chromsizes.iter().find(|&x| x.chrom == cp.Contig2).unwrap().size;
                
                // ctg1+ ctg2+
                let res1 = vec.par_iter(
                    ).map(|x| {
                        let pos1 = x[0];
                        let pos2 = x[1];
                        let res1 = length1 - pos1 + pos2;
                        res1
                    }).collect::<Vec<_>>();
                
                // ctg1+ ctg2-
                let res2 = vec.par_iter(
                ).map(|x| {
                    let pos1 = x[0];
                    let pos2 = x[1];
                    let res2 = length1 - pos1 + length2 - pos2;
                    res2
                }).collect::<Vec<_>>();

                // ctg1- ctg2+
                let res3 = vec.par_iter(
                ).map(|x| {
                    let pos1 = x[0];
                    let pos2 = x[1];
                    let res3 = pos1 + pos2;
                    res3
                }).collect::<Vec<_>>();
                
                // ctg1- ctg2-
                let res4 = vec.par_iter(
                ).map(|x| {
                    let pos1 = x[0];
                    let pos2 = x[1];
                    let res4 = pos1 + length2 - pos2;
                    res4
                }).collect::<Vec<_>>();

                (cp, vec![res1, res2, res3, res4])
                
            }).collect::<HashMap<_, Vec<Vec<u64>>>>();

        
        for (cp, res) in result {
            for (i, res1) in res.into_iter().enumerate() {
                let count = res1.len();
                let res1 = res1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ");
                let contig_pairs = match i {
                    0 => format!("{}+ {}+", cp.Contig1, cp.Contig2),
                    1 => format!("{}+ {}-", cp.Contig1, cp.Contig2),
                    2 => format!("{}- {}+", cp.Contig1, cp.Contig2),
                    3 => format!("{}- {}-", cp.Contig1, cp.Contig2),
                    _ => panic!("Error: Invalid index"),
                };
                let line = format!("{}\t{}\t{}\n", contig_pairs, count, res1);
                wtr.write_all(line.as_bytes()).unwrap();
            }
        }
        
        log::info!("Successful output clm file `{}`", output);
    }
}


