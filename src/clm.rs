use anyhow::Result as anyResult;
// use indicatif::ProgressBar;
use std::collections::{HashMap, HashSet};
use std::io::{ Write, BufReader, BufRead };
use std::borrow::Cow;
use std::path::Path;
use std::io::BufWriter;
use std::fs::File;

use std::sync::{ Arc, Mutex };
use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use std::thread;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use memmap2::Mmap;
use indexmap::IndexMap;
use std::collections::BTreeMap;

use crate::core::BaseTable;
use crate::core::{ ContigPair2, ContigPair3 };
use crate::core::{ common_reader, common_writer };

#[derive(Debug, Clone)]
pub struct Clm {
    file: String
}

impl BaseTable for Clm {
    fn new(name: &String) -> Self {
        Clm {file: name.clone()}
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

impl Clm {
    // pub fn split_clm(&self, cluster_file: &String, output_dir: &String) -> anyResult<()> {
    //     let mut cluster_map: HashMap<String, Vec<String>> = HashMap::new();
        
    //     std::fs::create_dir_all(output_dir)?;

    //     let mut cluster_file = BufReader::new(common_reader(cluster_file));
    //     let mut line = String::new();
    //     while cluster_file.read_line(&mut line)? > 0 {
    //         let mut iter = line.split_whitespace();
    //         let cluster = iter.next().unwrap();
    //         let cluster_set = cluster_map.entry(cluster.to_string()).or_insert(Vec::new());
    //         for item in iter {
    //             cluster_set.push(item.to_string());
    //         }
    //         line.clear();
    //     }

    //     let mut cluster_paired_map: HashMap<ContigPair3, &String> = HashMap::new();
    //     for (cluster, contigs) in &cluster_map {
    //         for i in 0..contigs.len() {
    //             for j in i+1..contigs.len() {
    //                 let (a, b) = if contigs[i] <= contigs[j] {
    //                     (&contigs[i], &contigs[j])
    //                 } else {
    //                     (&contigs[j], &contigs[i])
    //                 };
    //                 let contig_pair = ContigPair3::new(a, b);
                    
    //                 cluster_paired_map.insert(contig_pair, cluster);
    //             }
    //         }
           
    //     }

    //     let reader = common_reader(&self.file);
        
    //     let mut writer_map: HashMap<&String, Box<dyn Write + Send>> = HashMap::with_capacity(cluster_map.len());
    //     for (cluster, _) in &cluster_map {
    //         let file_name = format!("{}/{}.clm", output_dir, cluster);
    //         let writer = common_writer(&file_name);
    //         writer_map.insert(cluster, writer);
    //     }
    //     fn strip_orient(s: &str) -> &str {
    //         if let Some(last) = s.as_bytes().last() {
    //             if *last == b'+' || *last == b'-' {
    //                 return &s[..s.len() - 1];
    //             }
    //         }
    //         s
    //     }
    //     for line in reader.lines() {
    //         let line = line?;
    //         let mut iter = line.split_whitespace();
    //         let c1_raw = match iter.next() {
    //             Some(x) => x,
    //             None => continue,
    //         };
    //         let c2_raw = match iter.next() {
    //             Some(x) => x,
    //             None => continue,
    //         };

    //         let contig1 = strip_orient(c1_raw);
    //         let contig2 = strip_orient(c2_raw);
    //         let (a, b) = if contig1 <= contig2 { (contig1, contig2) } else { (contig2, contig1) };
    //         let contig_pair = ContigPair3::new(a, b);

    //         let cluster = match cluster_paired_map.get(&contig_pair) {
    //             Some(cluster) => cluster,
    //             None => continue
    //         };
           

    //         let writer = writer_map.get_mut(cluster).unwrap();
    //         writeln!(writer, "{}", line)?;
        
    //     }
        

    //     Ok(())
    // }
    pub fn split_clm(&self, cluster_file: &String, output_dir: &String) -> anyResult<()> {

        let mut cluster_map: HashMap<String, Vec<String>> = HashMap::new();
        std::fs::create_dir_all(output_dir)?;

        let mut cluster_reader = BufReader::new(common_reader(cluster_file));
        let mut line = String::new();
        while cluster_reader.read_line(&mut line)? > 0 {
            let mut iter = line.split_whitespace();
            if let Some(cluster) = iter.next() {
                let cluster_set = cluster_map.entry(cluster.to_string()).or_default();
                for item in iter {
                    cluster_set.push(item.to_string());
                }
            }
            line.clear();
        }

        let mut contig_to_cids: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut writers = Vec::new();

        for (cluster, contigs) in cluster_map {
            let cid = writers.len();
            for c in contigs {
                contig_to_cids.entry(c).or_default().push(cid);
            }
            let file_path = format!("{}/{}.clm", output_dir, cluster);
            let writer = Arc::new(Mutex::new(BufWriter::new(common_writer(&file_path))));
            writers.push(writer);
        }

        let contig_to_cids = Arc::new(contig_to_cids);
        let writers = Arc::new(writers);
        let (tx, rx) = bounded::<Vec<String>>(200);

        let num_workers = 8; 
        let mut worker_handles = Vec::new();
        for _ in 0..num_workers {
            let rx = rx.clone();
            let contig_to_cids = Arc::clone(&contig_to_cids);
            let writers = Arc::clone(&writers);
            
            worker_handles.push(thread::spawn(move || {
                #[inline]
                fn strip_orient(s: &str) -> &str {
                    let b = s.as_bytes();
                    if !b.is_empty() && (b[b.len()-1] == b'+' || b[b.len()-1] == b'-') {
                        &s[..s.len() - 1]
                    } else {
                        s
                    }
                }

                #[inline]
                fn write_to_common_clusters(
                    cids1: &[usize],
                    cids2: &[usize],
                    writers: &[Arc<Mutex<BufWriter<Box<dyn Write + Send>>>>],
                    line: &str,
                ) {
                    for &cid in cids1 {
                        if cids2.iter().any(|&x| x == cid) {
                            let mut w = writers[cid].lock().unwrap();
                            writeln!(w, "{}", line).unwrap();
                        }
                    }
                }

                for batch in rx {
                    for line in batch {
                        let mut parts = line.splitn(3, |c: char| c == '\t' || c == ' ');
                        let c1_raw = match parts.next() {
                            Some(x) => x,
                            None => continue,
                        };
                        let c2_raw = match parts.next() {
                            Some(x) => x,
                            None => continue,
                        };

                        let contig1 = strip_orient(c1_raw);
                        let contig2 = strip_orient(c2_raw);

                        let cids1 = contig_to_cids.get(contig1);
                        let cids2 = contig_to_cids.get(contig2);
                        if let (Some(c1s), Some(c2s)) = (cids1, cids2) {
                            write_to_common_clusters(c1s, c2s, &writers, &line);
                        }
                    }
                }
            }));
        }

        let reader = common_reader(&self.file);
        let mut batch = Vec::with_capacity(2000);
        for line in reader.lines() {
            batch.push(line?);
            if batch.len() >= 2000 {
                tx.send(batch).unwrap();
                batch = Vec::with_capacity(2000);
            }
        }
        if !batch.is_empty() { tx.send(batch).unwrap(); }
        drop(tx); 
        for h in worker_handles {
            h.join().unwrap();
        }
        
        for w in writers.iter() {
            w.lock().unwrap().flush()?;
        }

        Ok(())
    }

    pub fn parse_clm(&self) -> anyResult<IndexMap<String, Vec<u32>>> {
        let reader = common_reader(&self.file);
        let mut data: IndexMap<String, Vec<u32>> = IndexMap::new();

        for line in reader.lines() {
            let line = line?;
            let mut iter = line.split('\t');
            let contigs = iter.next().unwrap();
            let _ = iter.next().unwrap();
            
            let values = iter.next().unwrap();
            let values = values.split(' ').map(|x| x.parse::<u32>().unwrap()).collect::<Vec<u32>>();
            data.entry(contigs.to_string()).or_insert_with(Vec::new).extend(values);
        }

        Ok(data)
    }

}


pub fn merge_clm(clm_files: Vec<String>, output: &String) -> anyResult<()> {
    let mut wtr = BufWriter::new(common_writer(output));

    let (sender, receiver) = bounded(100);
    let mut data: IndexMap<String, Vec<u32>> = IndexMap::new();

    let producer_handle = thread::spawn(move || {
        for clm_file in clm_files {
            let reader = common_reader(&clm_file);
            for line in reader.lines() {
                let line = line.unwrap();
                sender.send(line).unwrap();
            }
        }
    });

    let consumer_handles: Vec<_> = (0..8).map(|_| {
        let receiver = receiver.clone();
        thread::spawn(move || {
            let mut local_data: IndexMap<String, Vec<u32>> = IndexMap::new();
            for line in receiver.iter() {
                let mut iter = line.split('\t');
                let contigs = iter.next().unwrap();
                let _ = iter.next().unwrap();
                
                let values = iter.next().unwrap();
                let values = values.split(' ').map(|x| x.parse::<u32>().unwrap()).collect::<Vec<u32>>();
                local_data.entry(contigs.to_string()).or_insert_with(Vec::new).extend(values);
            }
            local_data
        })
    }).collect();

    producer_handle.join().unwrap();

    for handle in consumer_handles {
        let local_data = handle.join().unwrap();
        for (contigs, values) in local_data {
            data.entry(contigs).or_insert_with(Vec::new).extend(values);
        }
    }


    let (sender, receiver) = bounded(100);

    let producer_handle = thread::spawn(move || {
        for (contigs, values) in data {
            sender.send((contigs, values)).unwrap();
        }
    });

    let consumer_handles: Vec<_> = (0..8).map(|_| {
        let receiver = receiver.clone();
        thread::spawn(move || {
            let mut local_buffer = Vec::new();
            for (contigs, values) in receiver.iter() {
                let count = values.len();
                let values_str = values.iter().map(ToString::to_string).collect::<Vec<String>>().join(" ");
                local_buffer.push(format!("{}\t{}\t{}", contigs, count, values_str));
            }
            local_buffer
        })
    }).collect();

    producer_handle.join().unwrap();

    for handle in consumer_handles {
        let local_buffer = handle.join().unwrap();
        for line in local_buffer {
            writeln!(wtr, "{}", line).unwrap();
        }
    }

    wtr.flush().unwrap();

    Ok(())
}
