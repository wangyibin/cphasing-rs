#[allow(dead_code)]
use anyhow::Result as AnyResult;
use bio::io::fastq::{Reader, Record, Writer};
// use bio::io::fasta::{Reader as FastaReader,
//                     Writer as FastaReader}
use bio::pattern_matching::horspool::Horspool;
use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use indexmap::IndexMap;
use hashbrown::HashMap as hashHashMap;
use rayon::prelude::*;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::io;
use std::path::Path;
use std::ops::AddAssign;
use std::sync::{Arc, Mutex};
use std::thread;
use seq_io::prelude::*;
use seq_io::fastq;
use seq_io::fastx::Reader as FastxReader;
use seq_io::parallel::{ read_process_fastx_records, read_process_recordsets };

use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;
use crate::sketch::hash;


#[derive(Debug, Clone)]
enum FileType {
    Fasta,
    Fastq,
    Unknown,
}

fn get_file_type(filename: &str) -> io::Result<FileType>  {
    let reader = common_reader(&filename);
    let mut lines = reader.lines();
    let first_line = match lines.next() {
        Some(line) => line?,
        None => return Ok(FileType::Unknown),
    };

    if first_line.starts_with(">") {
        Ok(FileType::Fasta)
    } else if first_line.starts_with("@") {
        Ok(FileType::Fastq)
    } else {
        Ok(FileType::Unknown)
    }
}


pub struct Fastx {
    pub file: String,
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
                *length = record.seq_lines()
                                .fold(0, |l, seq| l + seq.len());
            },
            |record, length| { // runs in main thread
                chrom_size.insert(record.id().unwrap().to_owned(), *length as u64); 
                None::<()>
            }).unwrap();
    
        Ok(chrom_size)
    }

    pub fn get_chrom_seqs (&self) -> AnyResult<HashMap<String, String>> {
        let reader = common_reader(&self.file);
        let reader = FastxReader::new(reader);
        let mut chrom_seqs: HashMap<String, String> = HashMap::new();

        read_process_fastx_records(reader, 8, 4,
            |record, seq| { // runs in worker
                *seq = record.seq_lines()
                                .fold(String::new(), |mut s, seq| {
                                    s.push_str(&String::from_utf8(seq.to_vec()).unwrap());
                                    s
                });
            },
            |record, seq| { // runs in main thread
                chrom_seqs.insert(record.id().unwrap().to_owned(), seq.to_owned()); 
                None::<()>
            }).unwrap();
    
        Ok(chrom_seqs)
    }

    pub fn count_re(&self, patterns: &String) -> AnyResult<(hashHashMap<String, u64>, hashHashMap<String, u64>)> {
        use aho_corasick::AhoCorasick;
        let chrom_count: hashHashMap<String, u64> = hashHashMap::new();
        let pattern_vec = patterns.split(",").collect::<Vec<&str>>();
        let pattern_vec_uppercase: Vec<String> = pattern_vec.iter().map(|x| x.to_uppercase()).collect();
        // let pattern_vec_lowercase: Vec<String> = pattern_vec.iter().map(|x| x.to_lowercase()).collect();
        // let all_patterns: Vec<&str> = pattern_vec_uppercase
        //                                 .iter()
        //                                 .chain(pattern_vec_lowercase.iter())
        //                                 .map(|x| x.as_str())
        //                                 .collect();
        
        let chrom_count_mutex = Mutex::new(hashHashMap::new());
        let chrom_length_mutex = Mutex::new(hashHashMap::new());
        
        
        let ac = AhoCorasick::new(&pattern_vec_uppercase).unwrap();

        log::set_max_level(log::LevelFilter::Off);
        let reader = common_reader(&self.file);
        log::set_max_level(log::LevelFilter::Info);
        let reader = FastxReader::new(reader);

        
        read_process_fastx_records(
            reader,
            4, // Number of worker threads
            2, // Number of partitions
            |record, (count, length)| {
                // Worker thread: Count matches for all patterns
                let sequence = record.seq().to_ascii_uppercase();
                let mut local_count = 0u64;

                for mat in ac.find_iter(&sequence) {
                    local_count += 1;
                }

                *length = record.seq_lines().fold(0, |l, seq| l + seq.len());
                
                *count = local_count;
            },
            |record, (count, length)| {
                let mut chrom_count = chrom_count_mutex.lock().unwrap();
                chrom_count
                    .entry(record.id().unwrap().to_owned())
                    .or_insert(0)
                    .add_assign(*count as u64);
                let mut chrom_length = chrom_length_mutex.lock().unwrap();
                chrom_length
                    .entry(record.id().unwrap().to_owned())
                    .or_insert(0)
                    .add_assign(*length as u64);
                None::<()>
            },
        ).unwrap();
        
        let chrom_count = chrom_count_mutex.into_inner().unwrap();
        let chrom_length = chrom_length_mutex.into_inner().unwrap();
        let pattern_count = pattern_vec_uppercase.len() as u64;

        Ok((chrom_count, chrom_length))
    }

    pub fn count_re2(&self, patterns: &String) -> AnyResult<IndexMap<String, u64>> {
        
        let mut chrom_count: IndexMap<String, u64> = IndexMap::new();
        let pattern_vec = patterns.split(",").collect::<Vec<&str>>();
    
        // convert to uppercase
        let pattern_vec_uppercase = pattern_vec.iter().map(|x| x.to_uppercase()).collect::<Vec<String>>();
        for pattern in &pattern_vec_uppercase {
            log::info!("Counting restriction site `{}` in `{}`", pattern, self.file);
        }
        // convert to lowercase
        let pattern_vec_lowercase = pattern_vec.iter().map(|x| x.to_lowercase()).collect::<Vec<String>>();
        // merge two vectors
        let pattern_vec = pattern_vec_uppercase.iter().chain(pattern_vec_lowercase.iter()).collect::<Vec<&String>>();


        for pattern in pattern_vec {
            log::set_max_level(log::LevelFilter::Off);
            let reader = common_reader(&self.file);
            log::set_max_level(log::LevelFilter::Info);
            let reader = FastxReader::new(reader);
            read_process_fastx_records(reader, 4, 2,
                |record, count| { // runs in worker
                    *count = record.seq().windows(pattern.len()).filter(|&x| x == pattern.as_bytes()).count();
                                    
                },
                |record, count| { // runs in main thread
                    chrom_count.entry(record.id().unwrap().to_owned()).or_insert(0).add_assign(*count as u64);
                    None::<()>
                }).unwrap();
        }
        
    
        Ok(chrom_count)
    }

    pub fn digest(&self, patterns: &String, slope: i64) -> AnyResult<HashMap<String, Vec<Vec<i64>>>> {
        
        let pattern_vec = patterns.split(",").collect::<Vec<&str>>();
        
        
        let mut positions: HashMap<String, Vec<Vec<i64>>> = HashMap::new();

        for pattern in pattern_vec {
            let mut pattern_uppercase = pattern.to_owned();
            pattern_uppercase.make_ascii_uppercase();
            log::info!("Identifing restriction site positions `{}` in `{}`", pattern_uppercase, self.file);
            let horspool = Horspool::new(pattern_uppercase.as_bytes());
            let reader = common_reader(&self.file);
            let reader = bio::io::fasta::Reader::new(reader);
            for result in reader.records() {
                let record = result?;
                let name = record.id().to_owned();
                let seq = record.seq();
                // seq uppercase

                let mut seq = seq.to_owned();
                seq.make_ascii_uppercase();
                let seq_length = seq.len();
                let seq_str = std::str::from_utf8(&seq).unwrap();
                let occ: Vec<usize> = horspool.find_all(&seq).collect();
                let len_vec = vec![seq_length];

                if occ.len() >= 1 {
                    let end_pos = seq_length;
                    for (i, pos) in occ.into_iter().enumerate() {
                        
                        let pos_mid = if pos >= seq_length {
                            pos as i64
                        } else {
                            (pos + pattern.len() / 2) as i64
                        };
                        let mut pos_start = pos_mid - slope;
                        let mut pos_end = pos_mid + slope;
                        if pos_start < 0 {
                            pos_start = 0
                        }
                        if pos_end > end_pos as i64 {
                            pos_end = end_pos as i64
                        }
                        
                        positions.entry(name.clone()).or_insert(Vec::new()).push(vec![pos_start, pos_end]);

                    }
                }
            }
        }
        
        Ok(positions)
    }
    pub fn slidefasta(&self, output: &String, window: u64, step: u64)  -> AnyResult<()>{
        use std::fmt::Write;
        let window = window as usize;
        let step = step as usize;
        let step = if step == 0 { window } else { step };
        let reader = common_reader(&self.file);

        let reader = bio::io::fasta::Reader::new(reader);
        
        
        let writer = common_writer(output);
        let wtr = bio::io::fasta::Writer::new(writer);
        let wtr = Arc::new(Mutex::new(wtr));
        let (sender, receiver) = bounded::<(String, Vec<u8>)>(100);
        let consumer_count = 8;
        let mut handles = Vec::new();
        for _ in 0..consumer_count {
            let receiver = receiver.clone();
            let writer = Arc::clone(&wtr);

            handles.push(thread::spawn(move || {
                
                while let Ok((id, seq_bytes)) = receiver.recv() {
                    let seq_length = seq_bytes.len();
                    let seq_str = match std::str::from_utf8(&seq_bytes) {
                        Ok(s) => s,
                        Err(e) => {
                            log::error!("Invalid UTF-8: {}", e);
                            continue;
                        }
                    };
                    let mut start = 0;
                    let mut end = window.min(seq_length);
                    let mut vec = Vec::new();
                    let mut name_buf = String::with_capacity(64);
                    while start < seq_length {
                        name_buf.clear();
                        let _ = write!(name_buf, "{}:{}-{}", id, start, end);
                        let record = bio::io::fasta::Record::with_attrs(
                            &name_buf,
                            None,
                            &seq_str[start..end].as_bytes()
                        );
                        
                        vec.push(record);

                        start += step;
                        end += step;
                        if end > seq_length {
                            end = seq_length;
                        }
                    }
                    {
                        let mut wtr = writer.lock().unwrap();
                        for record in vec {
                            match wtr.write_record(&record) {
                                Ok(_) => {},
                                Err(e) => {
                                    log::error!("Error writing record: {}", e);
                                }
                            }
                        }
                    }
                }
            }));
        }

        for result in reader.records() {
            match result {
                Ok(rec) => {
                    
                    if sender.send((rec.id().to_owned(), rec.seq().to_vec())).is_err() {
                        break;
                    }
                }
                Err(e) => {
                    log::error!("Error reading record: {}", e);
                    continue;
                }
            }
        }

        drop(sender);
    
        for handle in handles {
            handle.join().unwrap();
        }

        
        // for result in reader.records() {
        //     let record = match result {
        //         Ok(record) => record,
        //         Err(e) => {
        //             log::error!("Error reading record: {}", e);
        //             continue
        //         }
        //     };
        //     let seq = record.seq();
        //     let seq_length = seq.len();
        //     let seq_str = std::str::from_utf8(&seq).unwrap();
        //     let mut start = 0;
        //     let mut end = window;
        //     let mut i = 0;
        //     if end >= seq_length {
        //         end = seq_length;
        //     }
            
        //     while start < seq_length {
        //         let seq_window = &seq_str[start..end as usize];
        //         let record = bio::io::fasta::Record::with_attrs(format!("{}:{}-{}", record.id(), start, end).as_str(),
        //                                                     None, seq_window.as_bytes());
        //         match wtr.write_record(&record) {
        //             Ok(_) => {},
        //             Err(e) => {
        //                 log::error!("Error writing record: {}", e);
        //             }
                
        //         };
        //         start += step;
        //         end += step;
        //         if end >= seq_length {
        //             end = seq_length;
        //         }
                
        //         i += 1;
        //     } 
        // }

        Ok(())


    }
    pub fn slide(&self, output: &String, 
                    window: u64, step: u64, 
                    min_length: u64, filetype: &str,
                    coordinate_suffix: bool) -> AnyResult<()>{
        let window = window as usize;
        let step = step as usize;
        let min_length = min_length as usize;
        let step = if step == 0 { window } else { step };
        let total_reads = 0;
        let file_type = match filetype {
            "fasta" => FileType::Fasta,
            "fastq" => FileType::Fastq,
            _ => {
                get_file_type(&self.file).unwrap()
            }
        };

        let reader = common_reader(&self.file);
        match file_type {
            FileType::Fastq => {
                let reader = Reader::with_capacity(100000, reader);
                let writer = common_writer(output);
                let mut wtr = Writer::new(writer);
                let mut output_counts = 0;
                let mut slided_counts = 0;
                let mut filter_counts = 0;
                for result in reader.records() {
                    let record = match result {
                        Ok(record) => record,
                        Err(e) => {
                            log::error!("Error reading record: {}", e);
                            continue
                        }
                    };
                    let seq = record.seq();
                    let seq_length = seq.len();
                    
                    if seq_length < min_length {
                        filter_counts += 1;
                        continue
                    }
                    output_counts += 1;
                    
                    let qual = record.qual();

                    let mut start = 0;
                    let mut end = window;
                    let mut i = 0;
                    if end >= seq_length {
                        end = seq_length;
                    }
                    
                    while start < seq_length {
                        let seq_window = &seq[start..end as usize];
                        let qual_window = &qual[start..end as usize];
                        let record = match coordinate_suffix {
                            false => Record::with_attrs(format!("{}_{}", record.id(), i).as_str(),
                                                        None, seq_window, qual_window),
                            true => Record::with_attrs(format!("{}:{}-{}", record.id(), start, end).as_str(),
                                        None, seq_window, qual_window)
                        };
                        match wtr.write_record(&record) {
                            Ok(_) => {},
                            Err(e) => {
                                log::error!("Error writing record: {}", e);
                            }
                        
                        };
                        start += step;
                        end += step;
                        if end >= seq_length {
                            end = seq_length;
                        }
                        
                        i += 1;
                    } 

                    slided_counts += i;

                    
                }
                log::info!("Filtered {} reads.", filter_counts);
                log::info!("Slide {} reads into {} reads", output_counts, slided_counts);
            },
            FileType::Fasta => {
                let reader = bio::io::fasta::Reader::with_capacity(1024, reader);
                let writer = common_writer(output);
                let mut wtr = bio::io::fasta::Writer::new(writer);
                
                let mut output_counts = 0;
                let mut slided_counts = 0;
                let mut filter_counts = 0;
                for result in reader.records() {
                    let record = match result {
                        Ok(record) => record,
                        Err(e) => {
                            log::error!("Error reading record: {}", e);
                            continue
                        }
                    };
                    let seq = record.seq();
                    let seq_length = seq.len();
                    
                    if seq_length < min_length {
                        filter_counts += 1;
                        continue
                    }
                    output_counts += 1;
            
                    let mut start = 0;
                    let mut end = window;
                    let mut i = 0;
                    if end >= seq_length {
                        end = seq_length;
                    }
                    
                    while start < seq_length {
                        let seq_window = &seq[start..end as usize];
                       
                        let record = match coordinate_suffix {
                            false => bio::io::fasta::Record::with_attrs(format!("{}_{}", record.id(), i).as_str(),
                                                                        None, seq_window),
                            true => bio::io::fasta::Record::with_attrs(format!("{}:{}-{}", record.id(), start, end).as_str(),
                                                                        None, seq_window),
                        };
                        match wtr.write_record(&record) {
                            Ok(_) => {},
                            Err(e) => {
                                log::error!("Error writing record: {}", e);
                            }
                        
                        };
                        start += step;
                        end += step;
                        if end >= seq_length {
                            end = seq_length;
                        }
                        
                        i += 1;
                    } 

                    slided_counts += i;

                    
                }
                log::info!("Filtered {} reads.", filter_counts);
                log::info!("Slide {} reads into {} reads", output_counts, slided_counts);
            },
            FileType::Unknown => {
                log::error!("Unknown type of input file. must a fastq or a fasta.");
            }
        }
        
        Ok(())
    }

    pub fn kmer_count(&self, k: usize) -> AnyResult<HashMap<String, u64>> {
        let reader = common_reader(&self.file);
        let reader = bio::io::fasta::Reader::new(reader);
        log::info!("Counting kmers of length {}", k);
        let mut kmer_count: HashMap<String, u64> = HashMap::new();
        for result in reader.records() {
            let record = match result {
                Ok(record) => record,
                Err(e) => {
                    log::error!("Error reading record: {}", e);
                    continue
                }
            };
            let seq = record.seq();
            let seq_length = seq.len();
            let seq_str = std::str::from_utf8(&seq).unwrap();
            for i in 0..(seq_length - k + 1) {
                let kmer = &seq_str[i..(i+k)];
                *kmer_count.entry(kmer.to_string()).or_insert(0) += 1;
                
                let rev_kmer = bio::alphabets::dna::revcomp(kmer.bytes());
                *kmer_count.entry(String::from_utf8(rev_kmer).unwrap()).or_insert(0) += 1;
            }
        }

        log::info!("{} kmers counted", kmer_count.len());
        Ok(kmer_count)
    }

    pub fn kmer_positions(&self, k: usize, kmer_list: &String, output: &String) -> AnyResult<()> {
        let reader = common_reader(kmer_list);
        let reader = BufReader::new(reader);
        let mut kmer_set = HashSet::new();
        for line in reader.lines() {
            let line = line?;
            let kmer_list = line.split("\t").collect::<Vec<&str>>();
            kmer_set.insert(kmer_list[0].to_string());
        }
     
        let reader = common_reader(&self.file);
        let reader = bio::io::fasta::Reader::new(reader);
        let mut writer = common_writer(output);

        for result in reader.records() {
            let record = match result {
                Ok(record) => record,
                Err(e) => {
                    log::error!("Error reading record: {}", e);
                    continue
                }
            };
            let seq = record.seq();
            let seq_length = seq.len();
            let seq_str = std::str::from_utf8(&seq).unwrap();
            for i in 0..(seq_length - k + 1) {
                let kmer = &seq_str[i..(i+k)];
                if kmer_set.contains(kmer) {
                    let _ = writer.write_all(format!("{}\t{}\t{}\n", record.id(), i, i+k).as_bytes());
                } else {
                    let rev_kmer = String::from_utf8(bio::alphabets::dna::revcomp(kmer.bytes())).unwrap();
                    if kmer_set.contains(&rev_kmer) {
                        let _ = writer.write_all(format!("{}\t{}\t{}\n", record.id(), i, i+k).as_bytes());
                    }
                }

            }
        }
      
        Ok(())
    }

    pub fn mask_high_frequency_kmer(&self, k: usize, threshold: u64, output: &String) -> AnyResult<()> {
        let kmer_count = self.kmer_count(k)?;
        let mut high_freq_kmer: Vec<&String> = Vec::new();
        for (kmer, count) in kmer_count.iter() {
            if *count > threshold {
                high_freq_kmer.push(kmer);
            }
        }

        // k number of "N"
        let mut n = String::new();
        for _ in 0..k {
            n.push('N');
        }

        let reader = common_reader(&self.file);
        let reader = bio::io::fasta::Reader::new(reader);
        let mut writer = common_writer(output);
        for result in reader.records() {
            let record = match result {
                Ok(record) => record,
                Err(e) => {
                    log::error!("Error reading record: {}", e);
                    continue
                }
            };
            
            let seq = record.seq();
            let seq_length = seq.len();

            let seq_str = std::str::from_utf8(&seq).unwrap();
            let mut masked_seq = seq_str.to_string();
            for kmer in &high_freq_kmer {
                let re = regex::Regex::new(kmer).unwrap();
                masked_seq = re.replace_all(&masked_seq, n.clone()).to_string();
            }
            let _ = writer.write_all(format!(">{}\n{}\n", record.id(), masked_seq).as_bytes());
           
        }
        Ok(())
    }


    pub fn split_by_cluster(&self, cluster_file: &String, trim_length: usize) -> AnyResult<()> {
        let reader = common_reader(cluster_file);
        let mut cluster_map: HashMap<String, String> = HashMap::new();
        let mut groups = Vec::new();
        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();
            let line = line.split("\t").collect::<Vec<&str>>();
            let group = line[0].to_string();
            groups.push(group.clone());

            let contigs = line[2].split(" ").collect::<Vec<&str>>();
            for contig in contigs {

                cluster_map.insert(contig.to_string(), group.clone());
            }
        }

        let reader = common_reader(&self.file);
        let reader = bio::io::fasta::Reader::new(reader);
        let mut writer_map: HashMap<String, Box<dyn Write + Send>> = HashMap::new();
        for group in groups {
            let writer = common_writer(&format!("{}.contigs.fasta", group));
            writer_map.insert(group, writer);
        }

        let trim_threshold = trim_length * 3;

        for result in reader.records() {
            let record = match result {
                Ok(record) => record,
                Err(e) => {
                    log::error!("Error reading record: {}", e);
                    continue
                }
            };
            let id = record.id();
            
            if let Some(cluster) = cluster_map.get(id) {
                if let Some(writer) = writer_map.get_mut(cluster) {
                    let mut seq = record.seq();
                    let length = seq.len();
                    if length > trim_threshold {
                        seq = &seq[trim_length..length - trim_length];
                    }
                    let seq_str = std::str::from_utf8(&seq).unwrap();
                    writer.write_all(format!(">{}\n{}\n", id, seq_str).as_bytes()).unwrap();
                }
            }
        }

        Ok(())
    }
}

// split fastq into several files by record number
pub fn split_fastq(input_fastq: &String, output_prefix: &String, 
             record_num: usize) -> Result<(), Box<dyn std::error::Error>> {
    let buf = common_reader(input_fastq);
    let fastq = Reader::from_bufread(buf);
    let mut i = 0;
    let mut j = 0;
    log::info!("write {} records to {}", &record_num, format!("{}_{}.fastq.gz", output_prefix, j));
    let writer = common_writer(&format!("{}_{}.fastq.gz", output_prefix, j));
    let mut wtr = Writer::new(writer);
    for r in fastq.records() {
        let record = r?;
        wtr.write_record(&record)?;
        i += 1;
        if i == record_num {
            j += 1;
            let writer = common_writer(&format!("{}_{}.fastq.gz", output_prefix, j));
            wtr = Writer::new(writer);
            log::info!("write {} records to {}", i, format!("{}_{}.fastq.gz", output_prefix, j));
            i = 0;
        }
    }

    Ok(())
}

