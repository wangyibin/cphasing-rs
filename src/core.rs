#![allow(non_snake_case)]
use anyhow::Result as anyResult;
// use async_compression::tokio::bufread::GzipDecoder;
// use tokio::io::TBufReader;

use bio::io::fastq;
use flate2::read;
use flate2::write::GzEncoder;
use gzp::deflate::{Gzip, Mgzip};
use gzp::{ZBuilder, Compression};
use gzp::{par::compress::{ParCompress, ParCompressBuilder}};
use gzp::{par::decompress::ParDecompressBuilder};
use rust_htslib::bam::{
    self, 
    record::Aux, record::Cigar, 
    record::CigarStringView,
    Header, HeaderView,
    Read, Reader, Record,
    Writer, 
};
use std::env;
use std::io::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;
use std::ffi::OsStr;
use std::io::Cursor;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read as StdRead, Write};
use std::error::Error;
use std::result::Result;
use std::path::{Path, PathBuf};
use std::process::{Command, exit};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;
use serde::{ Deserialize, Serialize };
use rayon::prelude::*;



const BUFFER_SIZE: usize = 256 * 1024;
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

impl std::fmt::Display for ChromSizeRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}\t{}", self.chrom, self.size)
    }
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

struct PrependThenReader<R: StdRead> {
    head: Cursor<Vec<u8>>,
    tail: R,
}
impl<R: StdRead> StdRead for PrependThenReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.head.position() < self.head.get_ref().len() as u64 {
            self.head.read(buf)
        } else {
            self.tail.read(buf)
        }
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


pub fn is_gzip_file(file_path: &str) -> io::Result<bool> {
    let mut file = File::open(file_path)?;
    let mut magic_number = [0u8; 2];
    file.read_exact(&mut magic_number)?;

    Ok(magic_number == [0x1F, 0x8B])
}

fn valid_gzip_header(buf: &[u8], off: usize) -> bool {
    if off + 10 > buf.len() { return false; }
    if buf[off] != 0x1F || buf[off + 1] != 0x8B || buf[off + 2] != 0x08 { return false; }
    let flg = buf[off + 3];

    if (flg & 0xE0) != 0 { return false; }

    let mut p = off + 10;


    if (flg & 0x04) != 0 {
        if p + 2 > buf.len() { return false; }
        let xlen = u16::from_le_bytes([buf[p], buf[p + 1]]) as usize;
        p += 2;
        if p + xlen > buf.len() { return false; }
        p += xlen;
    }

    if (flg & 0x08) != 0 {
        while p < buf.len() && buf[p] != 0 { p += 1; }
        if p >= buf.len() { return false; }
        p += 1;
    }

    if (flg & 0x10) != 0 {
        while p < buf.len() && buf[p] != 0 { p += 1; }
        if p >= buf.len() { return false; }
        p += 1;
    }

    if (flg & 0x02) != 0 {
        if p + 2 > buf.len() { return false; }
        p += 2;
    }
    true
}

pub fn is_mgzip_file(file_path: &str) -> io::Result<bool> {
    if file_path == "-" { return Ok(false); }
    let path = std::path::Path::new(file_path);
    let mut f = File::open(path)?;

    let mut buf = vec![0u8; 8 * 1024 * 1024];
    let n = f.read(&mut buf)?;
    if n < 20 { return Ok(false); }
    buf.truncate(n);


    if !valid_gzip_header(&buf, 0) { return Ok(false); }


    let mut i = 1;
    while i + 10 <= buf.len() {
        if buf[i] == 0x1F && buf[i + 1] == 0x8B && buf[i + 2] == 0x08 {
           
            if i >= 8 {
                let trailer = &buf[i - 8..i];
               
                let isize = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
                if isize != 0 && valid_gzip_header(&buf, i) {
                    return Ok(true);
                }
            }
        }
        i += 1;
    }
    // while i + 10 <= buf.len() {
    //     if buf[i] == 0x1F && buf[i + 1] == 0x8B && buf[i + 2] == 0x08 {
    //         if valid_gzip_header(&buf, i) {
    //             return Ok(true);
    //         }
    //     }
    //     i += 1;
    // }

    Ok(false)
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

    let threads = env::var("CPHASING_IO_THREADS")
        .unwrap_or_else(|_| "8".to_string())
        .parse::<usize>()
        .unwrap_or(8);

    if suffix == Some(OsStr::new("gz")) && 
     ( is_gzip_file(&file_path.to_string_lossy()).unwrap_or(false) | 
        is_mgzip_file(&file_path.to_string_lossy()).unwrap_or(false) ) {
        let fp = match File::open(&file_path) {
            Err(error) => panic!("No such of file `{}`: {}", file_path.display(), error),
            Ok(fp) => fp,
        };

        let try_mgzip = is_mgzip_file(&file_path.to_string_lossy()).unwrap_or(false);
        if try_mgzip {
            log::info!("`{}` detected multi-member (try Mgzip)…", file_path.display());
            let mgz = ParDecompressBuilder::<Mgzip>::new()
                .num_threads(threads).expect("set num_threads failed")
                .from_reader(BufReader::with_capacity(BUFFER_SIZE, File::open(&file_path).expect("No such of file")));
            const PROBE: usize = 64 * 1024;
            let mut mgz_probe = mgz;
            let mut head = Vec::with_capacity(PROBE);
            match (&mut mgz_probe).take(PROBE as u64).read_to_end(&mut head) {
                Ok(_) => {
                   
                    let reader = PrependThenReader { head: Cursor::new(head), tail: mgz_probe };
                    return Box::new(BufReader::with_capacity(BUFFER_SIZE, reader));
                }
                Err(e) => {
                    log::warn!("Mgzip probe failed on `{}`: {}. Fallback to plain gzip.", file_path.display(), e);
                }
            }
        }

        log::info!("`{}` treated as gzip (flate2 MultiGzDecoder)", file_path.display());
        return Box::new(BufReader::with_capacity(BUFFER_SIZE, read::MultiGzDecoder::new(fp)));
        
    } else if suffix == Some(OsStr::new("mgz")) && is_mgzip_file(&file_path.to_string_lossy()).unwrap_or(false)  {
        log::info!("`{}` is detected as mgzip file", file_path.display());
        Box::new(BufReader::with_capacity(
            BUFFER_SIZE,
            ParDecompressBuilder::<Mgzip>::new()
                .num_threads(threads).expect("set num_threads failed")
                .from_reader(BufReader::with_capacity(BUFFER_SIZE, File::open(&file_path).expect("No such of file"))),
        ))

    } else {

        parse_input(Some(file_path.clone())).expect(format!("No such of file, {}", file_path.display()).as_str())
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
    let threads = env::var("CPHASING_IO_THREADS")
        .unwrap_or_else(|_| "8".to_string())
        .parse::<usize>()
        .unwrap_or(8);
    if suffix == Some(OsStr::new("gz")) {
        let writer = ParCompressBuilder::<Mgzip>::new()
            .num_threads(threads).expect("REASON")
            // .compression_level(Compression::new(6))
            .from_writer(buffered);

        Box::new(writer)
    } else if suffix == Some(OsStr::new("mgz")) {
        let writer = ParCompressBuilder::<Mgzip>::new()
            .num_threads(threads).expect("REASON")
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

pub struct Prof {
    pub on: bool,
    pub read_ns: AtomicU64,
    pub worker_proc_ns: AtomicU64,
    pub worker_fmt_ns: AtomicU64,
    pub writer_ns: AtomicU64,
    pub out_bytes: AtomicU64,
}
impl Prof {
    pub fn new() -> Self {
        let on = std::env::var_os("CPHASING_PROF").is_some();
        Self {
            on,
            read_ns: AtomicU64::new(0),
            worker_proc_ns: AtomicU64::new(0),
            worker_fmt_ns: AtomicU64::new(0),
            writer_ns: AtomicU64::new(0),
            out_bytes: AtomicU64::new(0),
        }
    }
    pub fn add(ns: &AtomicU64, dur: u128) { ns.fetch_add(dur as u64, Ordering::Relaxed); }
}
