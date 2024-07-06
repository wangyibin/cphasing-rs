use anyhow::Result as anyResult;
use bio::io::fastq;
use flate2::read;
use flate2::write::GzEncoder;
use gzp::deflate::Gzip;
use gzp::{ZBuilder, Compression};
use gzp::{par::compress::{ParCompress, ParCompressBuilder}};
use rust_htslib::bam::{
    self, 
    record::Aux, record::Cigar, 
    record::CigarStringView,
    Header, HeaderView,
    Read, Reader, Record,
    Writer, 
};
use std::io::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::error::Error;
use std::result::Result;
use std::path::{Path, PathBuf};
use std::process::{Command, exit};
use serde::{ Deserialize, Serialize };
use rayon::prelude::*;


const BUFFER_SIZE: usize = 128 * 1024;
type DynResult<T> = anyResult<T, Box<dyn Error + 'static>>;

pub trait BaseTable {
    fn new(name: &String) -> Self;

    fn file_name(&self) -> Cow<'_, str>;

    fn prefix(&self) -> String;
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ChromSizeRecord {
    pub chrom: String,
    pub size: u64,
}

#[derive(Debug, Clone)]
pub struct ChromSize {
    pub file: String,
}

impl BaseTable for ChromSize {
    fn new(name: &String) -> ChromSize {
        ChromSize { file: name.clone() }
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

impl ChromSize {
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

    pub fn data(&self) -> Result<HashMap<String, u64>, Box<dyn Error>> {
        let parse_result = self.parse();

        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
        };

        let mut db: HashMap<String, u64> = HashMap::new();

        for line in rdr.deserialize() {
            let record: ChromSizeRecord = line?;

            db.insert(record.chrom, record.size);
        }

        Ok(db) 
    }

    pub fn to_vec(&self) -> Result<Vec<ChromSizeRecord>, Box<dyn Error>> {
        let chromsizes = self.data().unwrap();
        let mut vec: Vec<ChromSizeRecord> = chromsizes
                                            .iter()
                                            .map(|(k, v)| ChromSizeRecord { 
                                                                chrom: k.to_string(), 
                                                                size: *v })
                                            .collect();
        vec.sort_unstable_by_key(|x| x.chrom.clone());
        Ok(vec)
    }

   
}


#[derive(PartialEq, Eq, Hash, Debug, Clone)]
pub struct ContigPair {
    pub Contig1: String,
    pub Contig2: String,
}

impl ContigPair {
    pub fn new(contig1: String, contig2: String) -> ContigPair {

        ContigPair { Contig1: contig1, Contig2: contig2 }
    }

    pub fn from_vec(vec: Vec<&String>) -> ContigPair {
        ContigPair { Contig1: (*vec[0].clone()).to_string(), 
                    Contig2: (*vec[1].clone()).to_string()}
    }

    pub fn swap(&mut self) {
        std::mem::swap(&mut self.Contig1, &mut self.Contig2);
    }

    pub fn order(&mut self) {
        if self.Contig1 > self.Contig2 {
            self.swap();
        }
    }


}

#[derive(PartialEq, Eq, Hash, Debug, Clone)]
pub struct ContigPair2<'a> {
    pub Contig1: &'a String,
    pub Contig2: &'a String,
}

impl ContigPair2<'_> {
    pub fn new<'a>(contig1: &'a String, contig2: &'a String) -> ContigPair2<'a> {

        ContigPair2 { Contig1: contig1, Contig2: contig2 }
    }

    // pub fn from_vec(vec: Vec<&String>) -> ContigPair2 {
    //     ContigPair2 { Contig1: &(*vec[0]).to_string(), 
    //                 Contig2: &(*vec[1]).to_string()}
    // }

    pub fn swap(&mut self) {
        std::mem::swap(&mut self.Contig1, &mut self.Contig2);
    }

    pub fn order(&mut self) {
        if *self.Contig1 > *self.Contig2 {
            self.swap();
        }
    }


}

#[derive(PartialEq, Eq, Hash, Debug, Clone)]
pub struct ContigPair3<'a> {
    pub Contig1: &'a str,
    pub Contig2: &'a str,
}

impl ContigPair3<'_> {
    pub fn new<'a>(contig1: &'a str, contig2: &'a str) -> ContigPair3<'a> {

        ContigPair3 { Contig1: contig1, Contig2: contig2 }
    }

    pub fn swap(&mut self) {
        std::mem::swap(&mut self.Contig1, &mut self.Contig2);
    }

    pub fn order(&mut self) {
        if *self.Contig1 > *self.Contig2 {
            self.swap();
        }
    }
}



// {parse_input, parse_output, common_reader, common_writer} learn from https://github.com/mrvollger/rustybam/blob/main/src/myio.rs
pub fn parse_input(path: Option<PathBuf>) -> DynResult<Box<dyn BufRead + Send + 'static>> {
    let fp: Box<dyn BufRead + Send + 'static> = match path {
        Some(path) => {
            if path.as_os_str() == "-" {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, io::stdin()))
            } else {
                Box::new(BufReader::with_capacity(BUFFER_SIZE, File::open(path)?))
            }
        }
        None => Box::new(BufReader::with_capacity(BUFFER_SIZE, io::stdin())),
    };
    Ok(fp)
}

pub fn common_reader(file: &str) -> Box<dyn BufRead + Send + 'static> {
    log::info!("Load `{}`", &file);
    let suffix = Path::new(file).extension();
    let file_path = PathBuf::from(file);

    if suffix == Some(OsStr::new("gz")) {
        let fp = match File::open(&file_path) {
            Err(error) => panic!("No such of file `{}`: {}", file_path.display(), error),
            Ok(fp) => fp,
        };
        Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            read::MultiGzDecoder::new(fp),
        ))
    } else {

        parse_input(Some(file_path)).expect("No such of file")
    }
}


pub fn parse_output(path: Option<PathBuf>) -> anyResult<Box<dyn Write + Send + 'static>> {
    let op: Box<dyn Write + Send + 'static> = match path {
        Some(path) => {
            if path.as_os_str() == "-" {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, io::stdout()))
            } else {
                Box::new(BufWriter::with_capacity(BUFFER_SIZE, File::create(path)?))
            }
        }
        None => Box::new(BufWriter::with_capacity(BUFFER_SIZE, io::stdout()))
    };
    Ok(op)
}

pub fn common_writer(file: &str) -> Box<dyn Write + Send + 'static> {
    let suffix = Path::new(file).extension();
    let file_path = PathBuf::from(file);

    let buffered = parse_output(Some(file_path)).expect("Error: failed to create output file ");
    
    if suffix == Some(OsStr::new("gz")) {
        let writer = ParCompressBuilder::<Gzip>::new()
            .num_threads(8).expect("REASON")
            // .compression_level(Compression::new(6))
            .from_writer(buffered);

        Box::new(writer)
    } else {
        buffered
    }

}

pub fn which(program: &str) -> Option<PathBuf> {
    let paths = std::env::var("PATH").unwrap();
    for path in paths.split(":") {
        let path = Path::new(path).join(program);
        if path.exists() {
            return Some(path);
        }
    }
    None
}

pub fn check_program(program: &str) {
    if which(program).is_none() {
        eprintln!("Error: {} is not installed", program);
        exit(1);
    }
}


// split contig size by binsize
pub fn binify(contigsizes: &HashMap<String, u64>, binsize: u32) -> anyResult<HashMap<String, Vec<u32>>> {
    
    let bins_db: HashMap<String, Vec<u32>> = contigsizes.par_iter().map(|(contig, size)| {
        let n_bins: u32 = (size / binsize as u64).try_into().unwrap();
        let mut bins = Vec::new();
        for i in 0..(n_bins + 1) {
            bins.push(i * binsize);
        }


        if let Some(last) = bins.last_mut() {
            *last = *size as u32;
        }

        (contig.to_string(), bins)
    }).collect();

    Ok(bins_db)

}