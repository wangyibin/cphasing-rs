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
    pub fn split_clm(&self, cluster_file: &String, output_dir: &String) -> anyResult<()> {
        let mut cluster_map: HashMap<String, Vec<String>> = HashMap::new();

        let mut cluster_file = BufReader::new(common_reader(cluster_file));
        let mut line = String::new();
        while cluster_file.read_line(&mut line)? > 0 {
            let mut iter = line.split_whitespace();
            let cluster = iter.next().unwrap();
            let mut cluster_set = cluster_map.entry(cluster.to_string()).or_insert(Vec::new());
            for item in iter {
                cluster_set.push(item.to_string());
            }
            line.clear();
        }

        let mut cluster_paired_map: HashMap<ContigPair3, &String> = HashMap::new();
        for (cluster, contigs) in &cluster_map {
            for i in 0..contigs.len() {
                for j in i+1..contigs.len() {
                    let contig_pair = if contigs[i] > contigs[j]{
                        ContigPair3::new(&contigs[j], &contigs[i])
                    } else {
                        ContigPair3::new(&contigs[i], &contigs[j])
                    };
                 
                    cluster_paired_map.insert(contig_pair, cluster);
                }
            }
           
        }

        let reader = common_reader(&self.file);
        
        let mut writer_map: HashMap<&String, Box<dyn Write + Send>> = HashMap::with_capacity(cluster_map.len());
        for (cluster, _) in &cluster_map {
            let file_name = format!("{}/{}.clm", output_dir, cluster);
            let writer = common_writer(&file_name);
            writer_map.insert(&cluster, writer);
        }

        for line in reader.lines() {
            let line = line?;
            let mut iter = line.split_whitespace();
            let contig1 = iter.next().unwrap();
            let contig2 = iter.next().unwrap();
            let contig1 = &contig1[..contig1.len() - 1];
            let contig2 = &contig2[..contig2.len() - 1];
            let contig_pair = ContigPair3::new(contig1, contig2);
            let cluster = match cluster_paired_map.get(&contig_pair) {
                Some(cluster) => cluster,
                None => continue
            };
           

            let writer = writer_map.get_mut(cluster).unwrap();
            writeln!(writer, "{}", line)?;
        
        }
        

        Ok(())
    }
}


pub fn merge_clm(clm_files: Vec<String>, output: &String) -> anyResult<()> {
    let mut wtr = BufWriter::new(common_writer(output));
    

    // let mut data: IndexMap<String, Vec<u32>> = IndexMap::new();
  
    // for clm_file in clm_files {
    //     let reader = common_reader(&clm_file);
    //     for line in reader.lines() {
    //         let line = line.unwrap();
    //         let mut iter = line.split("\t");
    //         let contigs = iter.next().unwrap();
    //         let _ = iter.next().unwrap();
            
    //         let values = iter.next().unwrap();
    //         let values = values.split(" ").map(|x| x.parse::<u32>().unwrap()).collect::<Vec<u32>>();

    //         data.entry(contigs.to_string()).or_insert(Vec::new()).extend(values);
    //     }
    // }

    // for (contigs, values) in data {
    //     let count: usize = values.iter().count();
    //     let values_str = values.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
    //     writeln!(wtr, "{}\t{}\t{}", contigs, count, values_str).unwrap();
    // }



    // let datas = clm_files.par_iter().map(|clm_file| {
    //     let reader = common_reader(&clm_file);
    //     let mut data = IndexMap::new();
    //     for line in reader.lines() {
    //         let line = line.unwrap();
    //         let mut iter = line.split("\t");
    //         let contigs = iter.next().unwrap();
    //         let count = iter.next().unwrap();
    //         let count = count.parse::<usize>().unwrap();
            
    //         let values = iter.next().unwrap();
    //         let values = values.split(" ").map(|x| x.parse::<usize>().unwrap()).collect::<Vec<usize>>();
    //         data.entry(contigs.to_string()).or_insert(Vec::new()).extend(values);
    //     }

    //     data
    // }).collect::<Vec<IndexMap<String, Vec<usize>>>>();

    // for d in datas {
    //     for (contigs, values) in d {
    //         data.entry(contigs.to_string()).or_insert(Vec::new()).extend(values);
    //     }
    // }

    // for d in datas {
    //     for (contigs, values) in d {
    //         data.entry(contigs).or_insert_with(Vec::new).extend(values);
    //     }
    // }


    // let mut data: Arc<Mutex<IndexMap<String, Vec<usize>>>> = Arc::new(Mutex::new(IndexMap::new()));
    // let (sender, receiver) = bounded::<Vec<String>>(1000);

    // let mut handles = vec![];
    // for _ in 0..8 {
    //     let receiver = receiver.clone();
    //     let data = Arc::clone(&data);
    //     handles.push(thread::spawn(move || {
    //         let mut updates: IndexMap<String, Vec<usize>> = IndexMap::new();
    //         while let Ok(records) = receiver.recv() {
    //             for line in records {
    //                 let mut iter = line.split("\t");
    //                 let contigs = iter.next().unwrap();
    //                 let _ = iter.next().unwrap();
                 
    //                 let values = iter.next().unwrap();
    //                 let values = values.split(" ").map(|x| x.parse::<usize>().unwrap()).collect::<Vec<usize>>();
    //                 updates.entry(contigs.to_string()).or_insert(Vec::new()).extend(values);
    //             }

    //             {
    //                 let mut data = data.lock().unwrap();
    //                 for (contigs, values) in &updates {
    //                     data.entry(contigs.to_string()).or_insert_with(Vec::new).extend(values);
    //                 }
    //             }
    //         }
    //     })
            
    //     )
        
    // }

    // let batch_size = 4000;
    // let mut batch = Vec::with_capacity(batch_size);
    // for clm_file in clm_files {
    //     let reader = common_reader(&clm_file);
    //     for line in reader.lines() {
    //         let line = line.unwrap();
            
    //         batch.push(line);
    //         if batch.len() >= batch_size {
    //             sender.send(std::mem::take(&mut batch)).unwrap();
    //         }


    //     }
    // }

    // if !batch.is_empty() {
    //     sender.send(batch).unwrap();
    // }

    // drop(sender);

    // for handle in handles {
    //     handle.join().unwrap();
    // }


    // let data = data.lock().unwrap();
    // for (contigs, values) in data.iter() {
    //     let count: usize = values.iter().count();
    //     let values_str = values.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
    //     writeln!(wtr, "{}\t{}\t{}", contigs, count, values_str).unwrap();
    // }



    // let (sender, receiver) = bounded(100);
    // let mut data: IndexMap<String, Vec<u32>> = IndexMap::new();
    // let handle = thread::spawn(move || {
    //     for clm_file in clm_files {
    //         let reader = common_reader(&clm_file);
    //         for line in reader.lines() {
    //             let line = line.unwrap();
    //             sender.send(line).unwrap();
    //         }
    //     }
    // });

    // for line in receiver.iter() {
    //     let mut iter = line.split('\t');
    //     let contigs = iter.next().unwrap();
    //     let _ = iter.next().unwrap();
        
    //     let values = iter.next().unwrap();
    //     let values = values.split(' ').map(|x| x.parse::<u32>().unwrap()).collect::<Vec<u32>>();
    //     data.entry(contigs.to_string()).or_insert_with(Vec::new).extend(values);
    // }

    // handle.join().unwrap();

    // for (contigs, values) in data {
    //     let count: usize = values.iter().count();
    //     let values_str = values.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
    //     writeln!(wtr, "{}\t{}\t{}", contigs, count, values_str).unwrap();
    // }



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

    // for (contigs, values) in data {
    //     let count: usize = values.iter().count();
    //     let values_str = values.iter().map(|x| x.to_string()).collect::<Vec<String>>().join(" ");
    //     writeln!(wtr, "{}\t{}\t{}", contigs, count, values_str).unwrap();
    // }

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
