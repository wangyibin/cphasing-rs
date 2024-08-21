use anyhow::Result as anyResult;
// use indicatif::ProgressBar;
use std::collections::{HashMap, HashSet};
use std::io::{ Write, BufReader, BufRead };
use std::borrow::Cow;
use std::path::Path;
use std::io::BufWriter;
use std::fs::File;

use std::sync::{ Arc, Mutex };
use rayon::prelude::*;

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