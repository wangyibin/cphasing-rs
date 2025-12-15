#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as AnyResult;
use std::borrow::Cow;
use std::collections::{ HashSet };
use hashbrown::{ HashMap };
use indexmap::IndexMap;
use std::error::Error;
use std::path::Path;
use std::io::{ Read, Write, BufReader, BufRead, BufWriter };
use rayon::prelude::*;

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair, ContigPair2, ContigPair3 };


#[derive(Debug, Clone)]
pub struct SplitContacts {
    pub file: String, 
    pub contigs: HashSet<String>,
    pub data: HashMap<ContigPair, Vec<f64>>,
}

impl BaseTable for SplitContacts {
    fn new(name: &String) -> SplitContacts {
        SplitContacts {
            file: name.clone(),
            contigs: HashSet::new(),
            data: HashMap::new(),
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

        (*file_prefix).to_string()
    }
    
}


impl SplitContacts {
    // pub fn read_from_file(file: &String) -> AnyResult<SplitContacts> {
    //     let mut split_contacts = SplitContacts::new(file);
    //     let mut buf = String::with_capacity(1024);
    //     let mut reader = common_reader(file);
       
    //     while reader.read_line(&mut buf)? > 0 {
    //         let line = buf.trim_end();
    //         if line.is_empty() {
    //             buf.clear();
    //             continue;
    //         }

           
    //         let mut fields = line.split('\t');
                
    //         let raw_contig1 = fields.next();
    //         let raw_contig2 = fields.next();
    //         let raw_count = fields.next();
        

    //         if let (Some(s1), Some(s2), Some(s_count)) = (raw_contig1, raw_contig2, raw_count) {
    //             let count: f64 = s_count.parse().unwrap_or(0.0);


    //             let (c1_name, idx1) = s1.rsplit_once('_').unwrap_or((s1, "0"));
    //             let (c2_name, idx2) = s2.rsplit_once('_').unwrap_or((s2, "0"));

    //             let pair = if c1_name < c2_name {
    //                 ContigPair {
    //                 Contig1: c1_name.to_string(),
    //                 Contig2: c2_name.to_string(),
    //                 } 
    //             } else {
    //                 ContigPair {
    //                 Contig1: c2_name.to_string(),
    //                 Contig2: c1_name.to_string(),
    //                 } 
    //             };

            

    //             split_contacts.contigs.insert(pair.Contig1.clone());
    //             split_contacts.contigs.insert(pair.Contig2.clone());


    //             let counts = split_contacts.data.entry(pair).or_insert_with(|| vec![0.0; 4]);

    //             match (idx1, idx2) {
    //                 ("0", "0") => counts[0] += count,
    //                 ("0", "1") => counts[1] += count,
    //                 ("1", "0") => counts[2] += count,
    //                 ("1", "1") => counts[3] += count,
    //                 _ => {},
    //             }
    //         }


    //         buf.clear();
    //     }

    //     Ok(split_contacts)
    // }
    pub fn read_from_file(file: &String, whitelist: Option<&HashSet<String>>) -> AnyResult<SplitContacts> {
        let mut split_contacts = SplitContacts::new(file);
        let mut buf = String::with_capacity(1024);
        let mut reader = common_reader(file);
       
        while reader.read_line(&mut buf)? > 0 {
            let line = buf.trim_end();
            if line.is_empty() {
                buf.clear();
                continue;
            }

           
            let mut fields = line.split('\t');
                
            let raw_contig1 = fields.next();
            let raw_contig2 = fields.next();
            let raw_count = fields.next();
        

            if let (Some(s1), Some(s2), Some(s_count)) = (raw_contig1, raw_contig2, raw_count) {
                let (c1_name, idx1) = s1.rsplit_once('_').unwrap_or((s1, "0"));
                let (c2_name, idx2) = s2.rsplit_once('_').unwrap_or((s2, "0"));

                if let Some(valid_contigs) = whitelist {
                    if !valid_contigs.contains(c1_name) || !valid_contigs.contains(c2_name) {
                        buf.clear();
                        continue;
                    }
                }

                let count: f64 = s_count.parse().unwrap_or(0.0);

                let pair = if c1_name < c2_name {
                    ContigPair {
                    Contig1: c1_name.to_string(),
                    Contig2: c2_name.to_string(),
                    } 
                } else {
                    ContigPair {
                    Contig1: c2_name.to_string(),
                    Contig2: c1_name.to_string(),
                    } 
                };

            

                split_contacts.contigs.insert(pair.Contig1.clone());
                split_contacts.contigs.insert(pair.Contig2.clone());


                let counts = split_contacts.data.entry(pair).or_insert_with(|| vec![0.0; 4]);

                match (idx1, idx2) {
                    ("0", "0") => counts[0] += count,
                    ("0", "1") => counts[1] += count,
                    ("1", "0") => counts[2] += count,
                    ("1", "1") => counts[3] += count,
                    _ => {},
                }
            }


            buf.clear();
        }

        Ok(split_contacts)
    }

    // pub fn read_from_file(file: &String, whitelist: Option<&HashSet<String>>) -> AnyResult<SplitContacts> {
    //     let mut split_contacts = SplitContacts::new(file);
        
    //     let file_handle = common_reader(file);
    //     let reader = BufReader::new(file_handle);

    //     let (final_data, final_contigs) = reader.lines()
    //         .par_bridge() 
    //         .map(|line| line.unwrap_or_default())
    //         .filter(|line| !line.is_empty())
    //         .fold(
    //             || (HashMap::new(), HashSet::new()), 
    //             |(mut local_data, mut local_contigs), line| {
    //                 let mut parts = line.split('\t');

    //                 if let (Some(s1), Some(s2), Some(s_count_str)) = (parts.next(), parts.next(), parts.next()) {
                        
    //                     let (c1_name, idx1) = s1.rsplit_once('_').unwrap_or((s1, "0"));
    //                     let (c2_name, idx2) = s2.rsplit_once('_').unwrap_or((s2, "0"));

    //                     if let Some(valid_contigs) = whitelist {
    //                         if !valid_contigs.contains(c1_name) || !valid_contigs.contains(c2_name) {
    //                             return (local_data, local_contigs);
    //                         }
    //                     }


    //                     let count: f64 = s_count_str.parse().unwrap_or(0.0);

    //                     if count > 0.0 {
    //                         let pair = if c1_name < c2_name {
    //                             ContigPair {
    //                                 Contig1: c1_name.to_string(),
    //                                 Contig2: c2_name.to_string(),
    //                             } 
    //                         } else {
    //                             ContigPair {
    //                                 Contig1: c2_name.to_string(),
    //                                 Contig2: c1_name.to_string(),
    //                             } 
    //                         };

    //                         local_contigs.insert(pair.Contig1.clone());
    //                         local_contigs.insert(pair.Contig2.clone());

    //                         let counts = local_data.entry(pair).or_insert_with(|| vec![0.0; 4]);
    //                         match (idx1, idx2) {
    //                             ("0", "0") => counts[0] += count,
    //                             ("0", "1") => counts[1] += count,
    //                             ("1", "0") => counts[2] += count,
    //                             ("1", "1") => counts[3] += count,
    //                             _ => {},
    //                         }
    //                     }
    //                 }
    //                 (local_data, local_contigs)
    //             }
    //         )
    //         .reduce(
    //             || (HashMap::new(), HashSet::new()), 
    //             |(mut data1, mut contigs1), (data2, contigs2)| {
    //                 for (k, v) in data2 {
    //                     let entry = data1.entry(k).or_insert(vec![0.0; 4]);
    //                     for i in 0..4 {
    //                         entry[i] += v[i];
    //                     }
    //                 }
    //                 contigs1.extend(contigs2);
    //                 (data1, contigs1)
    //             }
    //         );

    //     split_contacts.data = final_data;
    //     split_contacts.contigs = final_contigs;

    //     Ok(split_contacts)
    // }

    pub fn to_contacts(&self, contig2idx: &HashMap<String, usize>) -> HashMap<(usize, usize), f64> {
        let mut contacts: HashMap<(usize, usize), f64> = HashMap::new();

        for (contigpair, counts) in self.data.iter() {
            let contig1 = &contigpair.Contig1;
            let contig2 = &contigpair.Contig2;
            if !contig2idx.contains_key(contig1) || !contig2idx.contains_key(contig2) {
                continue;
            }
            let pair = (
                *contig2idx.get(&contigpair.Contig1).unwrap(),
                *contig2idx.get(&contigpair.Contig2).unwrap(),
            );
            

            let count = counts.iter().sum();
            contacts.insert(pair, count);
        }

        contacts

    }

    pub fn normalize_by_contig_sizes(&mut self, contig_sizes: &IndexMap<String, usize>) {
        for (contig_pair, counts) in self.data.iter_mut() {
            let size1 = *contig_sizes.get(&contig_pair.Contig1).unwrap_or(&1) as f64;
            let size2 = *contig_sizes.get(&contig_pair.Contig2).unwrap_or(&1) as f64;

            let norm_factor = (size1 * size2).sqrt();

            for count in counts.iter_mut() {
                *count /= norm_factor;
            }
        }
    }

    pub fn normalize_by_cis(&mut self) {
        let cis_data: HashMap<String, (f64, f64)> = self.data.iter()
                    .filter(|(pair, _counts)| pair.Contig1 == pair.Contig2)
                    .map(|(pair, counts)| (pair.Contig1.clone(), (counts[0], counts[3])))
                    .collect();

        let smoothing = 1.0; 

        for (contig_pair, counts) in self.data.iter_mut() {
            let contig1 = &contig_pair.Contig1;
            let contig2 = &contig_pair.Contig2;

            let (c1_0, c1_1) = cis_data.get(contig1).unwrap_or(&(0.0, 0.0));
            let (c2_0, c2_1) = cis_data.get(contig2).unwrap_or(&(0.0, 0.0));

            let norm_00 = ((c1_0 + smoothing) * (c2_0 + smoothing)).sqrt();
            let norm_01 = ((c1_0 + smoothing) * (c2_1 + smoothing)).sqrt();
            let norm_10 = ((c1_1 + smoothing) * (c2_0 + smoothing)).sqrt();
            let norm_11 = ((c1_1 + smoothing) * (c2_1 + smoothing)).sqrt();

            counts[0] /= norm_00;
            counts[1] /= norm_01;
            counts[2] /= norm_10;
            counts[3] /= norm_11;
        }

    }

    pub fn to_detailed_contact_matrix(&self, contig2idx: &HashMap<String, usize>) -> HashMap<(usize, usize), (f64, f64, f64, f64)> {
        let mut detailed_contacts: HashMap<(usize, usize), (f64, f64, f64, f64)> = HashMap::new();

        for (contig_pair, counts) in self.data.iter() {
            let pair = contig_pair;
            if !contig2idx.contains_key(&pair.Contig1) || !contig2idx.contains_key(&pair.Contig2) {
                continue;
            }
            let (idx1, idx2) = (
                *contig2idx.get(&pair.Contig1).unwrap(),
                *contig2idx.get(&pair.Contig2).unwrap(),
            );

            let head_head = (counts[0] + 1.0).ln();
            let head_tail = (counts[1] + 1.0).ln();
            let tail_head = (counts[2] + 1.0).ln();
            let tail_tail = (counts[3] + 1.0).ln();

    

            // let head_head = counts[0];
            // let head_tail = counts[1];
            // let tail_head = counts[2];
            // let tail_tail = counts[3];
            
            detailed_contacts.insert((idx1, idx2), (head_head, head_tail, tail_head, tail_tail));
            detailed_contacts.insert((idx2, idx1), (head_head, tail_head, head_tail, tail_tail));
        }

        detailed_contacts
    }

    pub fn extract_by_contigs(&mut self, contig_indices: &Vec<String>) {
        let contig_indices = contig_indices.iter().cloned().collect::<HashSet<String>>();
        self.data.retain(|contig_pair, _count| {
            contig_indices.contains(&contig_pair.Contig1) &&
            contig_indices.contains(&contig_pair.Contig2)
        });

        self.contigs.retain(|contig| {
            contig_indices.contains(contig)
        });
    }

}


pub fn split_contacts_by_clusters(
    split_contacts: &String,
    cluster_file: &String,
    output_dir: &String
) -> AnyResult<()> {
    let mut cluster_map: HashMap<String, Vec<String>> = HashMap::new();

    let mut cluster_file = BufReader::new(common_reader(cluster_file));
    let mut line = String::new();
    while cluster_file.read_line(&mut line)? > 0 {
        let mut iter = line.split_whitespace();
        let cluster = iter.next().unwrap();
        let cluster_set = cluster_map.entry(cluster.to_string()).or_insert(Vec::new());
        for item in iter {
            cluster_set.push(item.to_string());
        }
        line.clear();
    }

    let mut cluster_paired_map: HashMap<ContigPair3, &String> = HashMap::new();
    for (cluster, contigs) in &cluster_map {
        for i in 0..contigs.len() {
            for j in i+1..contigs.len() {
                let (a, b) = if contigs[i] <= contigs[j] {
                    (&contigs[i], &contigs[j])
                } else {
                    (&contigs[j], &contigs[i])
                };
                let contig_pair = ContigPair3::new(a, b);
                
                cluster_paired_map.insert(contig_pair, cluster);
            }
        }
        
    }


    let mut writer_map: HashMap<&String, Box<dyn Write + Send>> = HashMap::with_capacity(cluster_map.len());
    for (cluster, _) in &cluster_map {
        let file_name = format!("{}/{}.split.contacts.gz", output_dir, cluster);
        let writer = common_writer(&file_name);
        writer_map.insert(cluster, writer);
    }

    let reader = common_reader(&split_contacts);

    for line in reader.lines() {
        let line = line?;
        let mut fields = line.split('\t');
        let raw_contig1 = fields.next();
        let raw_contig2 = fields.next();
        let raw_count = fields.next();

        let (contig1, contig2) = if let (Some(s1), Some(s2)) = (raw_contig1, raw_contig2) {
            let (c1_name, _idx1) = s1.rsplit_once('_').unwrap_or((s1, "0"));
            let (c2_name, _idx2) = s2.rsplit_once('_').unwrap_or((s2, "0"));
            (c1_name, c2_name)
        } else {
            continue;
        };

        let (a, b) = if contig1 <= contig2 {
            (contig1, contig2)
        } else {
            (contig2, contig1)
        };

        let contig_pair = ContigPair3::new(a, b);

        if let Some(cluster) = cluster_paired_map.get(&contig_pair) {
            if let Some(writer) = writer_map.get_mut(cluster) {
                writeln!(writer, "{}", line)?;
            }
        }

    }


    Ok(())
    
}