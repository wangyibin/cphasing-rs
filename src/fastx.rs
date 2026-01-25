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
use needletail::{parse_fastx_file, Sequence};

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
    // pub fn get_chrom_size(&self) -> AnyResult<HashMap<String, u64>>  {
    //     use rayon::iter::{ParallelBridge, ParallelIterator};
    //     let reader = common_reader(&self.file);
    //     let fasta_reader = bio::io::fasta::Reader::new(reader);
    //     let chrom_size = Arc::new(Mutex::new(HashMap::new()));

    //     fasta_reader.records().par_bridge().for_each(|result| {
    //         if let Ok(record) = result {
    //             let id = record.id().to_owned();
    //             let len = record.seq().len() as u64;
    //             chrom_size.lock().unwrap().insert(id, len);
    //         }
    //     });
    
    //     Ok(Arc::try_unwrap(chrom_size).unwrap().into_inner().unwrap())
    // }

    // pub fn get_chrom_seqs (&self) -> AnyResult<HashMap<String, String>> {
    //     use rayon::iter::{ParallelBridge, ParallelIterator};
    //     let reader = common_reader(&self.file);
    //     let fasta_reader = bio::io::fasta::Reader::new(reader);
    //     let chrom_seqs = Arc::new(Mutex::new(HashMap::new()));

    //     fasta_reader.records().par_bridge().for_each(|result| {
    //         if let Ok(record) = result {
    //             let id = record.id().to_owned();
    //             let seq = String::from_utf8(record.seq().to_vec())
    //                 .expect("Invalid UTF-8 sequence");
    //             chrom_seqs.lock().unwrap().insert(id, seq);
    //         }
    //     });
    
    //     Ok(Arc::try_unwrap(chrom_seqs).unwrap().into_inner().unwrap())
    // }

    pub fn get_chrom_size(&self) -> AnyResult<HashMap<String, u64>> {
        let mut reader = needletail::parse_fastx_file(&self.file)
            .map_err(|e| anyhow::anyhow!(e.to_string()))?;
        
        let mut chrom_size = HashMap::new();
        while let Some(record) = reader.next() {
            let seq_rec = record.map_err(|e| anyhow::anyhow!(e.to_string()))?;
            let id = std::str::from_utf8(seq_rec.id())?.to_owned();
            
            chrom_size.insert(id, seq_rec.num_bases() as u64);
        }
        Ok(chrom_size)
    }

    pub fn get_chrom_seqs(&self) -> AnyResult<HashMap<String, String>> {
        let mut reader = needletail::parse_fastx_file(&self.file)
            .map_err(|e| anyhow::anyhow!(e.to_string()))?;
        let mut chrom_seqs = HashMap::new();

        while let Some(record) = reader.next() {
            let seq_rec = record.map_err(|e| anyhow::anyhow!(e.to_string()))?;
            let id = std::str::from_utf8(seq_rec.id())?.to_owned();
            let seq = String::from_utf8(seq_rec.seq().to_vec())?;
            chrom_seqs.insert(id, seq);
        }
        Ok(chrom_seqs)
    }

    // pub fn get_chrom_size(&self) -> AnyResult<HashMap<String, u64>>  {
    //     let reader = common_reader(&self.file);
    //     let reader = FastxReader::new(reader);
    //     let mut chrom_size: HashMap<String, u64> = HashMap::new();

    //     read_process_fastx_records(reader, 4, 2,
    //         |record, length| { // runs in worker
    //             *length = record.seq_lines()
    //                             .fold(0, |l, seq| l + seq.len());
    //         },
    //         |record, length| { // runs in main thread
    //             chrom_size.insert(record.id().unwrap().to_owned(), *length as u64); 
    //             None::<()>
    //         }).unwrap();
    
    //     Ok(chrom_size)
    // }

    // pub fn get_chrom_seqs (&self) -> AnyResult<HashMap<String, String>> {
    //     let reader = common_reader(&self.file);
    //     let reader = FastxReader::new(reader);
    //     let mut chrom_seqs: HashMap<String, String> = HashMap::new();

    //     read_process_fastx_records(reader, 8, 4,
    //         |record, seq| { // runs in worker
    //             *seq = record.seq_lines()
    //                             .fold(String::new(), |mut s, seq| {
    //                                 s.push_str(&String::from_utf8(seq.to_vec()).unwrap());
    //                                 s
    //             });
    //         },
    //         |record, seq| { // runs in main thread
    //             chrom_seqs.insert(record.id().unwrap().to_owned(), seq.to_owned()); 
    //             None::<()>
    //         }).unwrap();
    
    //     Ok(chrom_seqs)
    // }

    // pub fn count_re(&self, patterns: &String) -> AnyResult<(hashHashMap<String, u64>, hashHashMap<String, u64>)> {
    //     use aho_corasick::AhoCorasick;
    //     let chrom_count: hashHashMap<String, u64> = hashHashMap::new();
    //     let pattern_vec = patterns.split(",").collect::<Vec<&str>>();
    //     let pattern_vec_uppercase: Vec<String> = pattern_vec.iter().map(|x| x.to_uppercase()).collect();
    //     // let pattern_vec_lowercase: Vec<String> = pattern_vec.iter().map(|x| x.to_lowercase()).collect();
    //     // let all_patterns: Vec<&str> = pattern_vec_uppercase
    //     //                                 .iter()
    //     //                                 .chain(pattern_vec_lowercase.iter())
    //     //                                 .map(|x| x.as_str())
    //     //                                 .collect();
        
    //     let chrom_count_mutex = Mutex::new(hashHashMap::new());
    //     let chrom_length_mutex = Mutex::new(hashHashMap::new());
        
        
    //     let ac = AhoCorasick::new(&pattern_vec_uppercase).unwrap();

    //     log::set_max_level(log::LevelFilter::Off);
    //     let reader = common_reader(&self.file);
    //     log::set_max_level(log::LevelFilter::Info);
    //     let reader = FastxReader::new(reader);

        
    //     read_process_fastx_records(
    //         reader,
    //         4, // Number of worker threads
    //         2, // Number of partitions
    //         |record, (count, length)| {
    //             // Worker thread: Count matches for all patterns
    //             let sequence = record.seq().to_ascii_uppercase();
    //             let mut local_count = 0u64;

    //             for mat in ac.find_iter(&sequence) {
    //                 local_count += 1;
    //             }

    //             *length = record.seq_lines().fold(0, |l, seq| l + seq.len());
                
    //             *count = local_count;
    //         },
    //         |record, (count, length)| {
    //             let mut chrom_count = chrom_count_mutex.lock().unwrap();
    //             chrom_count
    //                 .entry(record.id().unwrap().to_owned())
    //                 .or_insert(0)
    //                 .add_assign(*count as u64);
    //             let mut chrom_length = chrom_length_mutex.lock().unwrap();
    //             chrom_length
    //                 .entry(record.id().unwrap().to_owned())
    //                 .or_insert(0)
    //                 .add_assign(*length as u64);
    //             None::<()>
    //         },
    //     ).unwrap();
        
    //     let chrom_count = chrom_count_mutex.into_inner().unwrap();
    //     let chrom_length = chrom_length_mutex.into_inner().unwrap();
    //     let pattern_count = pattern_vec_uppercase.len() as u64;

    //     Ok((chrom_count, chrom_length))
    // }

    pub fn count_re(&self, patterns: &String) -> AnyResult<(hashHashMap<String, u64>, hashHashMap<String, u64>)> {
        use aho_corasick::AhoCorasick;
        use rayon::iter::{ParallelBridge, ParallelIterator};

        let pattern_vec: Vec<String> = patterns.split(',').map(|s| s.to_uppercase()).collect();
        let ac = AhoCorasick::new(&pattern_vec).unwrap();

        let chrom_count = Arc::new(Mutex::new(hashHashMap::new()));
        let chrom_length = Arc::new(Mutex::new(hashHashMap::new()));

        let reader = common_reader(&self.file);
        let fasta_reader = bio::io::fasta::Reader::new(reader);

        fasta_reader.records().par_bridge().for_each(|result| {
            if let Ok(record) = result {
                let id = record.id().to_owned();
                let mut seq = record.seq().to_vec();
                seq.make_ascii_uppercase();

                let mut local_count = 0u64;
                for _ in ac.find_iter(&seq) {
                    local_count += 1;
                }

                let seq_len = seq.len() as u64;

                {
                    let mut c_map = chrom_count.lock().unwrap();
                    c_map.entry(id.clone()).or_insert(0).add_assign(local_count);
                    
                    let mut l_map = chrom_length.lock().unwrap();
                    l_map.entry(id).or_insert(0).add_assign(seq_len);
                }
            }
        });

        let final_count = Arc::try_unwrap(chrom_count).unwrap().into_inner().unwrap();
        let final_len = Arc::try_unwrap(chrom_length).unwrap().into_inner().unwrap();

        Ok((final_count, final_len))
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


    pub fn slidefasta(&self, output: &String, window: u64, step: u64) -> AnyResult<()> {
        use itoa;
        use std::io::Write as _;

        let window = window as usize;
        let step = if step == 0 { window } else { step as usize };

        let num_workers = rayon::current_num_threads().max(4).min(32);
        let writer_capacity = 64; 
        let pool_size = writer_capacity + num_workers + 8; 
    
        let (tx_writer, rx_writer) = bounded::<Vec<u8>>(writer_capacity); 
        let (tx_free, rx_free) = unbounded::<Vec<u8>>();
        for _ in 0..pool_size { 
            tx_free.send(Vec::with_capacity(16 * 1024 * 1024)).unwrap(); 
        }


        let output_clone = output.clone();
        let tx_free_clone = tx_free.clone();
        let writer_handle = thread::spawn(move || {
            let mut writer = common_writer(&output_clone);
            while let Ok(mut data) = rx_writer.recv() {
                let _ = writer.write_all(&data);
                data.clear();
                let _ = tx_free_clone.send(data); 
            }
        });

        let (tx_worker, rx_worker) = bounded::<(String, Arc<Vec<u8>>, usize, usize)>(200);
        let mut worker_handles = Vec::new();
        let num_workers = rayon::current_num_threads().max(4); 

        for _ in 0..num_workers {
            let rx = rx_worker.clone();
            let tx_w = tx_writer.clone();
            let rx_f = rx_free.clone();
            worker_handles.push(thread::spawn(move || {
                let mut itoa_b = itoa::Buffer::new();
                let mut batch_buffer = rx_f.recv().unwrap(); 

                while let Ok((id, seq, start_bound, end_bound)) = rx.recv() {
                    let mut curr = start_bound;
                    while curr < end_bound {
                        let end = (curr + window).min(seq.len());

                        batch_buffer.push(b'>');
                        batch_buffer.extend_from_slice(id.as_bytes());
                        batch_buffer.push(b':');
                        batch_buffer.extend_from_slice(itoa_b.format(curr).as_bytes());
                        batch_buffer.push(b'-');
                        batch_buffer.extend_from_slice(itoa_b.format(end).as_bytes());
                        batch_buffer.push(b'\n');
                        batch_buffer.extend_from_slice(&seq[curr..end]);
                        batch_buffer.push(b'\n');

                        if batch_buffer.len() > 14 * 1024 * 1024 {
                            if tx_w.send(batch_buffer).is_err() { return; }
                            batch_buffer = rx_f.recv().unwrap(); 
                        }

                        if end == seq.len() || curr + step >= end_bound { break; }
                        curr += step;
                    }
                }
                if !batch_buffer.is_empty() { let _ = tx_w.send(batch_buffer); }
            }));
        }

        let mut reader = needletail::parse_fastx_file(&self.file)
            .map_err(|e| anyhow::anyhow!(e.to_string()))?;
        
        while let Some(record) = reader.next() {
            let seq_rec = record.map_err(|e| anyhow::anyhow!(e.to_string()))?;
            let id = std::str::from_utf8(seq_rec.id())?.to_owned();
            let seq = Arc::new(seq_rec.normalize(false).into_owned());
            
            let seq_len = seq.len();
            let chunk_size = 5000 * step; 
            let mut chunk_start = 0;
            while chunk_start < seq_len {
                let chunk_end = (chunk_start + chunk_size).min(seq_len);
                tx_worker.send((id.clone(), Arc::clone(&seq), chunk_start, chunk_end))?;
                chunk_start += chunk_size;
            }
        }

        drop(tx_worker);
        drop(tx_writer);
        for h in worker_handles { h.join().unwrap(); }
    
        writer_handle.join().unwrap();
        Ok(())
    }


    // pub fn slidefasta(&self, output: &String, window: u64, step: u64) -> AnyResult<()> {
    //     use itoa;
    //     use std::io::Write as _;

    //     let window = window as usize;
    //     let step = step as usize;
    //     let step = if step == 0 { window } else { step };

    //     let (tx_writer, rx_writer) = bounded::<Vec<u8>>(40); 
    //     let (tx_free, rx_free) = unbounded::<Vec<u8>>();
    //     for _ in 0..80 { 
    //         tx_free.send(Vec::with_capacity(16 * 1024 * 1024)).unwrap(); 
    //     }

    //     let output_clone = output.clone();
    //     let tx_free_clone = tx_free.clone();
    //     let writer_handle = thread::spawn(move || {
    //         let mut writer = common_writer(&output_clone);
    //         while let Ok(mut data) = rx_writer.recv() {
    //             let _ = writer.write_all(&data);
    //             data.clear();
    //             let _ = tx_free_clone.send(data); 
    //         }
    //     });

    //     let (tx_worker, rx_worker) = bounded::<(String, Vec<u8>)>(100);
    //     let mut worker_handles = Vec::new();
    //     let num_workers = rayon::current_num_threads().max(4).min(16); 
    //     for _ in 0..num_workers {
    //         let rx = rx_worker.clone();
    //         let tx_w = tx_writer.clone();
    //         let rx_f = rx_free.clone();
    //         worker_handles.push(thread::spawn(move || {
    //             let mut itoa_b = itoa::Buffer::new();
    //             let mut batch_buffer = rx_f.recv().unwrap(); 
                
    //             while let Ok((id, seq)) = rx.recv() {
    //                 let seq_len = seq.len();
    //                 let id_bytes = id.as_bytes();
    //                 let mut start = 0;

    //                 while start < seq_len {
    //                     let end = (start + window).min(seq_len);

    //                     batch_buffer.push(b'>');
    //                     batch_buffer.extend_from_slice(id_bytes);
    //                     batch_buffer.push(b':');
    //                     batch_buffer.extend_from_slice(itoa_b.format(start).as_bytes());
    //                     batch_buffer.push(b'-');
    //                     batch_buffer.extend_from_slice(itoa_b.format(end).as_bytes());
    //                     batch_buffer.push(b'\n');
                        
    //                     batch_buffer.extend_from_slice(&seq[start..end]);
    //                     batch_buffer.push(b'\n');

    //                     if batch_buffer.len() > 14 * 1024 * 1024 {
    //                         if tx_w.send(batch_buffer).is_err() { return; }
    //                         batch_buffer = rx_f.recv().unwrap(); 
    //                     }

    //                     if end == seq_len { break; }
    //                     start += step;
    //                 }
    //             }
    //             if !batch_buffer.is_empty() {
    //                 let _ = tx_w.send(batch_buffer);
    //             }
    //         }));
    //     }


    //     let mut reader = needletail::parse_fastx_file(&self.file)
    //         .map_err(|e| anyhow::anyhow!(e.to_string()))?;
        
    //     while let Some(record) = reader.next() {
    //         let seq_rec = record.map_err(|e| anyhow::anyhow!(e.to_string()))?;
    //         let id = std::str::from_utf8(seq_rec.id())?.to_owned();
    //         let seq = seq_rec.normalize(false).into_owned(); 
    //         if tx_worker.send((id, seq)).is_err() { break; }
    //     }

    //     drop(tx_worker);
    //     for h in worker_handles { h.join().unwrap(); }
    //     drop(tx_writer);
    //     writer_handle.join().unwrap();
    //     Ok(())
    // }

    // pub fn slidefasta(&self, output: &String, window: u64, step: u64) -> AnyResult<()> {
    //     use std::fmt::Write as _; 
    //     use std::io::Write as _;  

    //     let window = window as usize;
    //     let step = step as usize;
    //     let step = if step == 0 { window } else { step };

    //     let (tx_writer, rx_writer) = bounded::<Vec<u8>>(1000);
    //     let output_clone = output.clone();
    //     let writer_handle = thread::spawn(move || {
    //         let mut writer = common_writer(&output_clone);
    //         while let Ok(data) = rx_writer.recv() {
    //             let _ = writer.write_all(&data);
    //         }
    //     });

    //     let (tx_worker, rx_worker) = bounded::<(String, Vec<u8>)>(100);
    //     let mut worker_handles = Vec::new();
    //     let num_workers = 8;

    //     for _ in 0..num_workers {
    //         let rx = rx_worker.clone();
    //         let tx_w = tx_writer.clone();
    //         worker_handles.push(thread::spawn(move || {
    //             let mut name_buf = String::with_capacity(128);
                
    //             while let Ok((id, seq)) = rx.recv() {
    //                 let seq_len = seq.len();
    //                 let mut start = 0;
                    
    //                 let mut batch_buffer = Vec::with_capacity(1024 * 1024); 

    //                 while start < seq_len {
    //                     let mut end = start + window;
    //                     if end > seq_len { end = seq_len; }

    //                     name_buf.clear();
    //                     let _ = write!(name_buf, ">{}:{}-{}\n", id, start, end);
                        
    //                     batch_buffer.extend_from_slice(name_buf.as_bytes());
    //                     batch_buffer.extend_from_slice(&seq[start..end]);
    //                     batch_buffer.push(b'\n');

    //                     if batch_buffer.len() > 8 * 1024 * 1024 { 
    //                         if tx_w.send(batch_buffer.split_off(0)).is_err() { break; }
    //                     }

    //                     if end == seq_len { break; }
    //                     start += step;
    //                 }
    //                 if !batch_buffer.is_empty() {
    //                     let _ = tx_w.send(batch_buffer);
    //                 }
    //             }
    //         }));
    //     }


    //     let mut reader = needletail::parse_fastx_file(&self.file)
    //         .map_err(|e| anyhow::anyhow!(e.to_string()))?;
        
    //     while let Some(record) = reader.next() {
    //         let seq_rec = record.map_err(|e| anyhow::anyhow!(e.to_string()))?;
    //         let id = std::str::from_utf8(seq_rec.id())?.to_owned();
    //         let seq = seq_rec.normalize(false).into_owned(); 
    //         if tx_worker.send((id, seq)).is_err() { break; }
    //     }


    //     drop(tx_worker);
    //     for h in worker_handles { h.join().unwrap(); }
    //     drop(tx_writer);
    //     writer_handle.join().unwrap();

    //     Ok(())
    // }

    // pub fn slidefasta(&self, output: &String, window: u64, step: u64)  -> AnyResult<()>{
    //     use std::fmt::Write;
    //     let window = window as usize;
    //     let step = step as usize;
    //     let step = if step == 0 { window } else { step };
    //     let reader = common_reader(&self.file);

    //     let reader = bio::io::fasta::Reader::new(reader);
        
    //     let writer = common_writer(output);
    //     let wtr = bio::io::fasta::Writer::new(writer);
    //     let wtr = Arc::new(Mutex::new(wtr));
    //     let (sender, receiver) = bounded::<(String, Vec<u8>)>(100);
    //     let consumer_count = 8;
    //     let mut handles = Vec::new();
    //     for _ in 0..consumer_count {
    //         let receiver = receiver.clone();
    //         let writer = Arc::clone(&wtr);

    //         handles.push(thread::spawn(move || {
                
    //             while let Ok((id, seq_bytes)) = receiver.recv() {
    //                 let seq_length = seq_bytes.len();
    //                 let seq_str = match std::str::from_utf8(&seq_bytes) {
    //                     Ok(s) => s,
    //                     Err(e) => {
    //                         log::error!("Invalid UTF-8: {}", e);
    //                         continue;
    //                     }
    //                 };
    //                 let mut start = 0;
    //                 let mut end = window.min(seq_length);
    //                 let mut vec = Vec::new();
    //                 let mut name_buf = String::with_capacity(64);
    //                 while start < seq_length {
    //                     name_buf.clear();
    //                     let _ = write!(name_buf, "{}:{}-{}", id, start, end);
    //                     let record = bio::io::fasta::Record::with_attrs(
    //                         &name_buf,
    //                         None,
    //                         &seq_str[start..end].as_bytes()
    //                     );
                        
    //                     vec.push(record);

    //                     start += step;
    //                     end += step;
    //                     if end > seq_length {
    //                         end = seq_length;
    //                     }
    //                 }
    //                 {
    //                     let mut wtr = writer.lock().unwrap();
    //                     for record in vec {
    //                         match wtr.write_record(&record) {
    //                             Ok(_) => {},
    //                             Err(e) => {
    //                                 log::error!("Error writing record: {}", e);
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }));
    //     }

    //     for result in reader.records() {
    //         match result {
    //             Ok(rec) => {
                    
    //                 if sender.send((rec.id().to_owned(), rec.seq().to_vec())).is_err() {
    //                     break;
    //                 }
    //             }
    //             Err(e) => {
    //                 log::error!("Error reading record: {}", e);
    //                 continue;
    //             }
    //         }
    //     }

    //     drop(sender);
    
    //     for handle in handles {
    //         handle.join().unwrap();
    //     }

    //     Ok(())

    // }

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


    // pub fn split_by_cluster(&self, cluster_file: &String, trim_length: usize) -> AnyResult<()> {
    //     // let reader = common_reader(cluster_file);
    //     // let mut cluster_map: HashMap<String, String> = HashMap::new();
    //     // let mut groups = Vec::new();
    //     // for line in reader.lines() {
    //     //     let line = line.unwrap();
    //     //     let line = line.trim();
    //     //     let line = line.split("\t").collect::<Vec<&str>>();
    //     //     let group = line[0].to_string();
    //     //     groups.push(group.clone());

    //     //     let contigs = line[2].split(" ").collect::<Vec<&str>>();
    //     //     for contig in contigs {

    //     //         cluster_map.insert(contig.to_string(), group.clone());
    //     //     }
    //     // }

    //     // let reader = common_reader(&self.file);
    //     // let reader = bio::io::fasta::Reader::new(reader);
    //     // let mut writer_map: HashMap<String, Box<dyn Write + Send>> = HashMap::new();
    //     // for group in groups {
    //     //     let writer = common_writer(&format!("{}.contigs.fasta", group));
    //     //     writer_map.insert(group, writer);
    //     // }

    //     // let trim_threshold = trim_length * 3;

    //     // for result in reader.records() {
    //     //     let record = match result {
    //     //         Ok(record) => record,
    //     //         Err(e) => {
    //     //             log::error!("Error reading record: {}", e);
    //     //             continue
    //     //         }
    //     //     };
    //     //     let id = record.id();
            
    //     //     if let Some(cluster) = cluster_map.get(id) {
    //     //         if let Some(writer) = writer_map.get_mut(cluster) {
    //     //             let mut seq = record.seq();
    //     //             let length = seq.len();
    //     //             if length <= 0 {
    //     //                 continue;
    //     //             }
    //     //             if length > trim_threshold {
    //     //                 seq = &seq[trim_length..length - trim_length];
    //     //             }
    //     //             if seq.len() <= 0 {
    //     //                 continue;
    //     //             }
    //     //             let seq_str = std::str::from_utf8(&seq).unwrap();
    //     //             writer.write_all(format!(">{}\n{}\n", id, seq_str).as_bytes()).unwrap();
    //     //         }
    //     //     }
    //     // }

    //     let reader = common_reader(cluster_file);
      
    //     let mut cluster_map: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();
    //     let mut groups = HashSet::new(); // Use HashSet to avoid duplicates

    //     for line in reader.lines() {
    //         let line = line.unwrap();
    //         let line = line.trim();
    //         if line.is_empty() { continue; }
    //         let parts = line.split("\t").collect::<Vec<&str>>();
    //         if parts.len() < 3 { continue; }
            
    //         let group = parts[0].to_string();
    //         groups.insert(group.clone());

    //         let contigs = parts[2].split(" ").collect::<Vec<&str>>();
    //         for contig in contigs {
              
    //             if let Some(idx) = contig.rfind('|') {
    //                 let raw_id = &contig[..idx];
    //                 let range_part = &contig[idx+1..];
    //                 if let Some(sep_idx) = range_part.find('_') {
    //                     let start_str = &range_part[..sep_idx];
    //                     let end_str = &range_part[sep_idx+1..];
    //                     if let (Ok(start), Ok(end)) = (start_str.parse::<usize>(), end_str.parse::<usize>()) {
    //                         cluster_map.entry(raw_id.to_string())
    //                             .or_insert_with(Vec::new)
    //                             .push((group.clone(), start, end));
    //                     }
    //                 }
    //             } else {
                  
    //                 cluster_map.entry(contig.to_string())
    //                     .or_insert_with(Vec::new)
    //                     .push((group.clone(), 0, usize::MAX));
    //             }
    //         }
    //     }

    //     let reader = common_reader(&self.file);
    //     let reader = bio::io::fasta::Reader::new(reader);
    //     let mut writer_map: HashMap<String, Box<dyn Write + Send>> = HashMap::new();
    //     for group in groups {
    //         let writer = common_writer(&format!("{}.contigs.fasta", group));
    //         writer_map.insert(group, writer);
    //     }

    //     let trim_threshold = trim_length * 3;

    //     for result in reader.records() {
    //         let record = match result {
    //             Ok(record) => record,
    //             Err(e) => {
    //                 log::error!("Error reading record: {}", e);
    //                 continue
    //             }
    //         };
    //         let id = record.id();
            
    //         if let Some(targets) = cluster_map.get(id) {
    //             let full_seq = record.seq();
    //             let full_len = full_seq.len();

    //             for (group, req_start, req_end) in targets {
    //                 if let Some(writer) = writer_map.get_mut(group) {
                        
    //                     let mut start = *req_start;
    //                     let mut end = *req_end;

                       
    //                     if end == usize::MAX {
    //                         end = full_len;
    //                     }

                        
    //                     if start >= full_len { continue; }
    //                     if end > full_len { end = full_len; }
    //                     if start >= end { continue; }

    //                     let slice_len = end - start;
    //                     if slice_len > trim_threshold {
    //                         if start < trim_length {
    //                             start += trim_length;
    //                         }
                           
    //                         if slice_len > trim_threshold {
        
    //                             let mut effective_end = end;
    //                             if (effective_end - start) > trim_threshold {
    //                                  effective_end -= trim_length;
    //                             }
    //                             end = effective_end;
    //                         }
    //                     }
    //                     if start >= end { continue; }

    //                     let seq_slice = &full_seq[start..end];
    //                     let seq_str = std::str::from_utf8(seq_slice).unwrap();
                        
                    
    //                     let header = if *req_end == usize::MAX {
    //                         id.to_string()
    //                     } else {
    //                         format!("{}|{}_{}", id, start, end)
    //                     };

    //                     writer.write_all(format!(">{}\n{}\n", header, seq_str).as_bytes()).unwrap();
    //                 }
    //             }
    //         }
    //     }

    //     Ok(())
    // }

    // pub fn split_by_cluster(&self, cluster_file: &String, trim_length: usize) -> AnyResult<()> {
    //     let reader = common_reader(cluster_file);
        
    //     let mut cluster_map: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();
    //     let mut groups = HashSet::new(); // Use HashSet to avoid duplicates

    //     for line in reader.lines() {
    //         let line = line.unwrap();
    //         let line = line.trim();
    //         if line.is_empty() { continue; }
    //         let parts = line.split("\t").collect::<Vec<&str>>();
    //         if parts.len() < 3 { continue; }
            
    //         let group = parts[0].to_string();
    //         groups.insert(group.clone());

    //         let contigs = parts[2].split(" ").collect::<Vec<&str>>();
    //         for contig in contigs {
                
    //             if let Some(idx) = contig.rfind('|') {
    //                 let raw_id = &contig[..idx];
    //                 let range_part = &contig[idx+1..];
    //                 if let Some(sep_idx) = range_part.find('_') {
    //                     let start_str = &range_part[..sep_idx];
    //                     let end_str = &range_part[sep_idx+1..];
    //                     if let (Ok(start), Ok(end)) = (start_str.parse::<usize>(), end_str.parse::<usize>()) {
    //                         cluster_map.entry(raw_id.to_string())
    //                             .or_insert_with(Vec::new)
    //                             .push((group.clone(), start, end));
    //                     }
    //                 }
    //             } else {
                    
    //                 cluster_map.entry(contig.to_string())
    //                     .or_insert_with(Vec::new)
    //                     .push((group.clone(), 0, usize::MAX));
    //             }
    //         }
    //     }

    //     // Wrap cluster_map in Arc for parallel access
    //     let cluster_map = Arc::new(cluster_map);

    //     let mut writer_map: HashMap<String, Box<dyn Write + Send>> = HashMap::new();
    //     for group in groups {
    //         let writer = common_writer(&format!("{}.contigs.fasta", group));
    //         writer_map.insert(group, writer);
    //     }

    //     let trim_threshold = trim_length * 3;

    //     // Use seq_io for parallel processing
    //     let reader = common_reader(&self.file);
    //     let reader = FastxReader::new(reader);

    //     read_process_fastx_records(
    //         reader,
    //         4, // threads
    //         2, // queue size
    //         |record, buffer: &mut Vec<(String, String)>| { // Worker thread: Process and Format
    //             buffer.clear();
    //             let id = match record.id() {
    //                 Ok(s) => s,
    //                 Err(_) => return,
    //             };
                
    //             if let Some(targets) = cluster_map.get(id) {
    //                 let full_seq = record.seq();
    //                 let full_len = full_seq.len();

    //                 for (group, req_start, req_end) in targets {
    //                     let mut start = *req_start;
    //                     let mut end = *req_end;

    //                     if end == usize::MAX {
    //                         end = full_len;
    //                     }

    //                     if start >= full_len { continue; }
    //                     if end > full_len { end = full_len; }
    //                     if start >= end { continue; }

    //                     let slice_len = end - start;
    //                     if slice_len > trim_threshold {
    //                         if start < trim_length {
    //                             start += trim_length;
    //                         }
                            
    //                         if slice_len > trim_threshold {
    //                             let mut effective_end = end;
    //                             if (effective_end - start) > trim_threshold {
    //                                     effective_end -= trim_length;
    //                             }
    //                             end = effective_end;
    //                         }
    //                     }
    //                     if start >= end { continue; }

    //                     let seq_slice = &full_seq[start..end];
    //                     let seq_str = String::from_utf8_lossy(seq_slice);
                        
    //                     let header = if *req_end == usize::MAX {
    //                         id.to_string()
    //                     } else {
    //                         format!("{}|{}_{}", id, start, end)
    //                     };

    //                     // Format the output string here in parallel
    //                     let output = format!(">{}\n{}\n", header, seq_str);
    //                     buffer.push((group.clone(), output));
    //                 }
    //             }
    //         },
    //         |_, buffer: &mut Vec<(String, String)>| { // Main thread: Write to files
    //             for (group, content) in buffer.iter() {
    //                 if let Some(writer) = writer_map.get_mut(group) {
    //                     writer.write_all(content.as_bytes()).unwrap();
    //                 }
    //             }
    //             None::<()>
    //         }
    //     ).unwrap();

    //     Ok(())
    // }
    
    pub fn split_by_cluster(&self, cluster_file: &String, trim_length: usize) -> AnyResult<()> {
        let reader = common_reader(cluster_file);
        
        let mut cluster_map: HashMap<String, Vec<(String, usize, usize)>> = HashMap::new();
        let mut groups = HashSet::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() { continue; }
            let parts = line.split('\t').collect::<Vec<&str>>();
            if parts.len() < 3 { continue; }
            
            let group = parts[0].to_string();
            groups.insert(group.clone());

            let contigs = parts[2].split(' ').collect::<Vec<&str>>();
            for contig in contigs {
                if let Some(idx) = contig.rfind('|') {
                    let raw_id = &contig[..idx];
                    let range_part = &contig[idx+1..];
                    if let Some(sep_idx) = range_part.find('_') {
                        let start_str = &range_part[..sep_idx];
                        let end_str = &range_part[sep_idx+1..];
                        if let (Ok(start), Ok(end)) = (start_str.parse::<usize>(), end_str.parse::<usize>()) {
                            cluster_map.entry(raw_id.to_string())
                                .or_insert_with(Vec::new)
                                .push((group.clone(), start, end));
                        }
                    }
                } else {
                    cluster_map.entry(contig.to_string())
                        .or_insert_with(Vec::new)
                        .push((group.clone(), 0, usize::MAX));
                }
            }
        }

        let cluster_map = Arc::new(cluster_map);
        let mut writer_map: HashMap<String, Box<dyn Write + Send>> = HashMap::new();
        for group in groups {
            let writer = common_writer(&format!("{}.contigs.fasta", group));
            writer_map.insert(group, writer);
        }

        let trim_threshold = trim_length * 3;
        
        let (tx, rx) = bounded::<(String, String)>(2000);

        let writer_handle = thread::spawn(move || {
            let mut writer_map = writer_map;
            while let Ok((group, content)) = rx.recv() {
                if let Some(writer) = writer_map.get_mut(&group) {
                    let _ = writer.write_all(content.as_bytes());
                }
            }
        });

        let reader = common_reader(&self.file);
        let fasta_reader = bio::io::fasta::Reader::new(reader);

        fasta_reader.records().par_bridge().for_each_with(tx, |tx, result| {
            if let Ok(record) = result {
                let id = record.id();
                
                if let Some(targets) = cluster_map.get(id) {
                    let full_seq = record.seq();
                    let full_len = full_seq.len();

                    for (group, req_start, req_end) in targets {
                        let mut start = *req_start;
                        let mut end = *req_end;

                        if end == usize::MAX { end = full_len; }
                        if start >= full_len { continue; }
                        if end > full_len { end = full_len; }
                        if start >= end { continue; }

                        let slice_len = end - start;
                        if slice_len > trim_threshold {
                            if start < trim_length { start += trim_length; }
                            
                            if (end - start) > trim_threshold {
                                end -= trim_length;
                            }
                        }
                        if start >= end { continue; }

                        let seq_slice = &full_seq[start..end];
                        let seq_str = String::from_utf8_lossy(seq_slice);
                        
                        let header = if *req_end == usize::MAX {
                            id.to_string()
                        } else {
                            format!("{}|{}_{}", id, start, end)
                        };

                        let output = format!(">{}\n{}\n", header, seq_str);
                        let _ = tx.send((group.clone(), output));
                    }
                }
            }
        });

        writer_handle.join().unwrap();

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

