use anyhow::Result as anyResult;
// use rdst::RadixSort;
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::process::exit;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use petgraph::prelude::*;
use petgraph::visit::NodeIndexable;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use itertools::Itertools;
use seq_io::prelude::*;
use seq_io::fastq;
use seq_io::fastx::Reader as FastxReader;
use seq_io::parallel::{ read_process_fastx_records, read_process_recordsets };

use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair, ContigPair2 };
use crate::sketch::{ sketch, MinimizerInfo, MinimizerData };

// maximal cliques 
pub fn bron_kerbosch(
    graph: &Graph<(), u64, Undirected>,
    mut r: HashSet<NodeIndex>,
    mut p: HashSet<NodeIndex>,
    mut x: HashSet<NodeIndex>,
    cliques: &mut Vec<HashSet<NodeIndex>>,
) {
    if p.is_empty() && x.is_empty() {
        cliques.push(r.clone());
        return;
    }

    let mut p_candidates = p.clone();
    p_candidates.retain(|&v| {
        let mut neighbors = graph.neighbors(v);
        neighbors.all(|n| p.contains(&n))
    });

    let p_cloned = p.clone();
    for v in p_cloned {
        let mut neighbors = graph.neighbors(v);
        let mut r_new = r.clone();
        r_new.insert(v);

        let p_new = p.intersection(&neighbors.clone().collect::<HashSet<_>>()).cloned().collect();
        let x_new = x.intersection(&neighbors.collect::<HashSet<_>>()).cloned().collect();

        bron_kerbosch(graph, r_new, p_new, x_new, cliques);

        p.remove(&v);
        x.insert(v);
    }
}

pub fn find_cliques(graph: &Graph<(), u64, Undirected>) -> Vec<HashSet<NodeIndex>> {
    let mut cliques = Vec::new();
    let n = graph.node_count();
    let mut all_nodes = (0..n).map(NodeIndex::new).collect::<HashSet<NodeIndex>>();
    let mut r = HashSet::new();

    let mut x = HashSet::new();

    bron_kerbosch(graph, r, all_nodes, x, &mut cliques);
    cliques
}

#[derive(Debug, Clone)]
pub struct AlleleHeader {
    pub header: Vec<String>,
    pub contigsizes: HashMap<String, u64>,
    pub contigs: HashSet<String>,
    pub minimizer: HashMap<String, u64>,
    pub unique_minimizer: HashMap<String, u64>,
}

impl AlleleHeader {
    pub fn new() -> Self {
        Self {
            header: Vec::new(),
            contigsizes: HashMap::new(),
            contigs: HashSet::new(),
            minimizer: HashMap::new(),
            unique_minimizer: HashMap::new(),
        }
    }

    pub fn from_file(&mut self, file: &str) -> anyResult<()> {
        let input = common_reader(file);
        let reader = BufReader::new(input);
        let mut header: Vec<String> = Vec::new();
        let mut contigsizes: HashMap<String, u64> = HashMap::new();
        let mut contigs: HashSet<String> = HashSet::new();
        let mut minimizer: HashMap<String, u64> = HashMap::new();
        let mut unique_minimizer: HashMap<String, u64> = HashMap::new();

        for result in reader.lines() {
            let line = result?;
            if !line.starts_with("#") {
                break;
            }
            let line_vec = line.strip_prefix("#").unwrap().split(" ").collect::<Vec<&str>>();
            let contig = line_vec[0].to_string();
            let size = line_vec[1].parse::<u64>().unwrap();
            let min = line_vec[2].parse::<u64>().unwrap();
            let unique_min = line_vec[3].parse::<u64>().unwrap();
            header.push(contig.clone());
            contigsizes.insert(contig.clone(), size);
            contigs.insert(contig.clone());
            minimizer.insert(contig.clone(), min);
            unique_minimizer.insert(contig.clone(), unique_min);
        }

        self.header = header;
        self.contigsizes = contigsizes;
        self.contigs = contigs;
        self.minimizer = minimizer;
        self.unique_minimizer = unique_minimizer;

        Ok(())
    }

    pub fn to_unique_minimizer_density(&mut self) -> HashMap<String, f64> {
        let mut density_hash: HashMap<String, f64> = HashMap::new();
        for (contig, min) in self.unique_minimizer.iter() {
            let minimizer = self.minimizer.get(contig).unwrap();
            let density = *min as f64 / *minimizer as f64;
            density_hash.insert(contig.clone(), density);
        }

        density_hash
    }

}

#[derive(Debug, Deserialize, Serialize)]
pub struct AlleleRecord2 {
    pub idx1: u32,
    pub idx2: u32,
    pub contig1: String,
    pub contig2: String,
    pub mz1: u32,
    pub mz2: u32,
    pub mz_shared: u32,
    pub similarity: f64,
    pub strand: i8,
}


#[derive(Debug)]
pub struct AlleleTable2 {
    pub file: String,
    pub header: AlleleHeader,
    pub allele_records: Vec<AlleleRecord2>,
}

impl BaseTable for AlleleTable2 {
    fn new(name: &String) -> AlleleTable2 {
        let mut header = AlleleHeader::new();
        let _ = header.from_file(name);
        
        AlleleTable2 { file: name.clone(),
                        header: header,
                        allele_records: Vec::new() }
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

impl AlleleTable2 {
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

    pub fn allele_records(&self) -> anyResult<Vec<AlleleRecord2>> {
        let parse_result = self.parse();
        let mut records: Vec<AlleleRecord2> = Vec::new();

        match parse_result {
            Ok(mut rdr) => {
                for result in rdr.deserialize() {
                    let record: AlleleRecord2 = result?;
                    records.push(record);
                }
            },
            Err(e) => {
                eprintln!("Error parsing allele table: {}", e);
                exit(1);
            }
        }

        Ok(records)
    }
    
    pub fn get_allelic_contig_pairs(&self) -> HashSet<ContigPair2> {
        let mut contig_pairs: HashSet<ContigPair2> = HashSet::new();
        let mut records = &self.allele_records;
        // sort records by mz_shared in ascending order
        // records.sort_by(|a, b| a.mz_shared.partial_cmp(&b.mz_shared).unwrap());

        for record in records {
            let contig1 = &record.contig1;
            let contig2 = &record.contig2;
            let mut contig_pair = ContigPair2::new(&contig1, &contig2);
            contig_pair.order();
            contig_pairs.insert(contig_pair);
        }

        contig_pairs
    }

    // pub fn get_allelic_contig_groups_by_cliques(&self, whitehash: &HashSet<&String>) -> HashMap<&String, Vec<Vec<&String>>> {
    //     let mut records = self.allele_records().unwrap();
    //     records.sort_by(|a, b| a.mz_shared.partial_cmp(&b.mz_shared).unwrap());
        
    //     let mut graph = Graph::<(), u64, Undirected>::new_undirected();
    //     let mut nodes: HashMap<String, NodeIndex> = HashMap::new();
    //     let mut idx_nodes: HashMap<NodeIndex, String> = HashMap::new();
    //     for record in records {
    //         let contig1 = record.contig1;
    //         let contig2 = record.contig2;
    //         if whitehash.len() > 0 {
    //             if !whitehash.contains(&contig1) || !whitehash.contains(&contig2) {
    //                 continue;
    //             }
    //         }
    //         let mz_shared = record.mz_shared;
            
    //         let node1 = *nodes.entry(contig1.clone()).or_insert_with(|| graph.add_node(()));
    //         let node2 = *nodes.entry(contig2.clone()).or_insert_with(|| graph.add_node(()));
    //         idx_nodes.insert(node1, contig1.clone());
    //         idx_nodes.insert(node2, contig2.clone());

    //         graph.add_edge(node1, node2, 1);
    //     }

    //     let cliques = find_cliques(&graph);

    //     let new_cliques = cliques.par_iter().map(|clique| {
    //         let mut new_clique = Vec::new();
    //         for node in clique {
    //             let contig = idx_nodes.get(node).unwrap();
    //             new_clique.push(contig.clone());
    //         }
    //         // sort clique by contig name
    //         new_clique.sort();
    //         new_clique
    //     }).collect::<Vec<Vec<String>>>();

    //     // remove duplicate 
    //     let mut unique_cliques: Vec<Vec<String>> = Vec::new();
    //     for clique in new_cliques {
    //         if !unique_cliques.contains(&clique) {
    //             unique_cliques.push(clique);
    //         }
    //     }
        
    //     let mut res = HashMap::new();
    //     for clique in unique_cliques.iter() {
    //         for contig in clique.iter() {
    //             res.entry(contig.clone()).or_insert_with(Vec::new);
    //             res.get_mut(contig).unwrap().push(clique.clone());
    //         }
           
    //     }
        
    //     res 

    // }

    pub fn get_allelic_record_by_contig_pairs(&self) -> HashMap<ContigPair2, &AlleleRecord2> {
        let mut data: HashMap<ContigPair2, &AlleleRecord2> = HashMap::new();
        let records = &self.allele_records;
        for record in records {
            
            let mut contig_pair = ContigPair2::new(&record.contig1, &record.contig2);
            contig_pair.order();
            data.insert(contig_pair, &record);
        }

        data
    }
    
    pub fn get_allelic_contigs_precise(&self, whitehash: &HashSet<&String>) -> HashMap<&String, Vec<Vec<&String>>> {
        let mut data: HashMap<&String, Vec<Vec<&String>>> = HashMap::new();
        let records = &self.allele_records;
        let check_whitehash = !whitehash.is_empty();

        for record in records {
            
            if check_whitehash && (!whitehash.contains(&record.contig1) || !whitehash.contains(&record.contig2)) {
                continue;
            }

            if !data.contains_key(&record.contig1) {
                data.insert(&record.contig1, Vec::new());
                data.get_mut(&record.contig1).unwrap().push(vec![&record.contig1, &record.contig2]);
            } else {
                data.get_mut(&record.contig1).unwrap()[0].push(&record.contig2);
            }
        }

        data
    }

    pub fn get_allelic_contigs(&self, method: &str, whitehash: &HashSet<&String>) -> HashMap<&String, Vec<Vec<&String>>> {
        let mut data: HashMap<&String, Vec<Vec<&String>>> = HashMap::new();
        let mut records = &self.allele_records;
        // sort records by mz_shared in ascending order
        // records.sort_by(|a, b| a.mz_shared.partial_cmp(&b.mz_shared).unwrap());
        for record in records {
            
            if whitehash.len() > 0 {
                if !whitehash.contains(&record.contig1) || !whitehash.contains(&record.contig2) {
                    continue;
                }
            }
            if !data.contains_key(&record.contig1) {
                data.insert(&record.contig1, vec![vec![&record.contig1, &record.contig2]]);
                // data.get_mut(&contig1).unwrap().push(contig2.clone());
            } else {
                if method == "fast" {
                    if data.contains_key(&record.contig1) {
                        continue;
                    }
                }
                data.get_mut(&record.contig1).unwrap().push(vec![&record.contig1, &record.contig2]);
            }
            
        }

        // data to HashMap<String, Vec<Vec<String>>>
        // let mut res = HashMap::new();
        // for (contig, contigs) in data.iter() {
        //     for contig in contigs.iter() {
        //         res.entry(contig.clone()).or_insert_with(Vec::new);
        //         res.get_mut(contig).unwrap().push(contigs.clone());
        //     }
        // }
        
        // res
        data
    }

    pub fn write(&self, records: &Vec<AlleleRecord2>) -> anyResult<()> {
        let writer = common_writer(&self.file);
        let mut wtr = csv::WriterBuilder::new()
                            .delimiter(b'\t')
                            .from_writer(writer);

        for record in records {
            wtr.serialize(record)?;
        }

        wtr.flush()?;
        Ok(())
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AlleleStrandRecord {
    pub idx1: u32,
    pub idx2: u32,
    pub contig1: String,
    pub contig2: String,
    pub length1: u32, 
    pub length2: u32,
    pub matches: u32,
    pub identity: f32, 
    pub strand: i8, 
}

#[derive(Debug, Clone,)]
pub struct AlleleStrandTable {
    pub file: String, 
    pub allele_strand_records: Vec<AlleleStrandRecord>,
}

impl BaseTable for AlleleStrandTable {
    fn new(name: &String) -> AlleleStrandTable {
        AlleleStrandTable {
            file: name.clone(), 
            allele_strand_records: Vec::new(),
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

impl AlleleStrandTable {
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
    
    pub fn get_allele_strand_records(&mut self) {
        let parse_result = self.parse();
      

        match parse_result {
            Ok(mut rdr) => {
                for result in rdr.deserialize() {
                    let record: AlleleStrandRecord = result.unwrap();
                    self.allele_strand_records.push(record);
                }
            },
            Err(e) => {
                eprintln!("Error parsing allele table: {}", e);
                exit(1);
            }
        }

    }

    pub fn get_info(&self) -> HashMap<ContigPair2, &AlleleStrandRecord> {

        let mut data: HashMap<ContigPair2, &AlleleStrandRecord> = HashMap::new();
        
        for record in &self.allele_strand_records {
            
            if record.contig1 >= record.contig2 {
                continue;
            }

            let mut contig_pair = ContigPair2::new(&record.contig1, &record.contig2);
    
            data.insert(contig_pair, &record);
        }

        data 

    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AlleleRecord {
    pub data: Vec<String>,
}

#[derive(Debug)]
pub struct AlleleTable {
    pub file: String,
    pub allele_records: Vec<AlleleRecord>,
}

impl AlleleTable {
    pub fn new(name: &String) -> AlleleTable {
        AlleleTable { file: name.clone(), allele_records: Vec::new() }
    }

    pub fn parse(&self) -> anyResult<Box<dyn BufRead + Send>> {
        let input = common_reader(&self.file);
        
        Ok(input)
    }

    pub fn allele_records(&self) -> anyResult<Vec<AlleleRecord>> {
        let parse_result = self.parse();
        let mut records: Vec<AlleleRecord> = Vec::new();

        let reader = BufReader::new(parse_result.unwrap());
        for result in reader.lines() {
            let line = result?;
            if line.starts_with("#") {
                continue;
            }
            //remove whitespace suffix 
            let line = line.trim_end();
            let line_vec = line.split("\t").collect::<Vec<&str>>();
            
            // convert Vec<&str> to Vec<String>
            let record = line_vec[2..].iter().map(|x| x.to_string()).collect::<Vec<String>>();
            if record.len() < 2 {
                continue;
            }
            records.push(AlleleRecord { data: record });

        }

        Ok(records)
    }

    pub fn to_contig_db(&self, method: &str,  whitehash: &HashSet<&String>) -> HashMap<&String, Vec<&AlleleRecord>> {
        let mut data: HashMap<&String, Vec<&AlleleRecord>> = HashMap::new();
        for record in &self.allele_records {
           
            'inner: for contig in &record.data {
                if whitehash.len() > 0 {
                    if !whitehash.contains(contig) {
                        continue 'inner;
                    }
                }
                
                match method {
                    "fast" => {
                       
                        if data.contains_key(contig) {
                            continue
                        }
                        data.entry(contig).or_insert_with(Vec::new);
                        data.get_mut(contig).unwrap().push(record);
                    },
                    "precise" => {
                        data.entry(contig).or_insert_with(Vec::new);
                        data.get_mut(contig).unwrap().push(record);
                    },
                    _ => todo!()
                }
               
            }
        }
        
        data
    }

    pub fn get_allelic_contig_pairs(&self, whitehash: &HashSet<&String>) -> HashSet<ContigPair2> {
        let mut contig_pairs: HashSet<ContigPair2> = HashSet::new();
        
        for record in &self.allele_records {
            for i in 0..record.data.len() - 1 {
                'inner: for j in (i+1)..record.data.len() {
                    let contig1 = &record.data[i];
                    let contig2 = &record.data[j];
                    if whitehash.len() == 0 {
                        if !whitehash.contains(contig1) || !whitehash.contains(contig2) {
                            continue 'inner;
                        }
                    }
                    let contig_pair = ContigPair2::new(&contig1, &contig2);
                    
                    contig_pairs.insert(contig_pair);
                }
                
            }
        }

        contig_pairs
    }

   


}



// https://github.com/lh3/partig

pub struct Uinfo {
    pub cnt1: u32,
    pub cnt2: u32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Anchor {
    rid1: u32,
    pos1: u32,
    rid2: u32,
    pos2: u32,
    rev: u8,
}

impl Anchor {
    fn swap_mut(&mut self) {
        std::mem::swap(&mut self.rid1, &mut self.rid2);
        std::mem::swap(&mut self.pos1, &mut self.pos2);
    }

    fn swap(&self) -> Self {
        Anchor {
            rid1: self.rid2,
            pos1: self.pos2,
            rid2: self.rid1,
            pos2: self.pos1,
            rev: self.rev
        }
    }
    
}

impl Ord for Anchor {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rid1.cmp(&other.rid1).then(self.rid2.cmp(&other.rid2)).then(self.rev.cmp(&other.rev))    
    }
}

impl PartialOrd for Anchor {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Copy, Clone, Debug)]
pub struct MatchRecord {
    pub rid1: u32,
    pub rid2: u32,
    pub rev: u8,
    pub mz1: u32,
    pub mz2: u32,
    pub mz_shared: u32,
    pub similarity: f64,
}

impl PartialEq for MatchRecord {
    fn eq(&self, other: &Self) -> bool {
        self.rid1 == other.rid1 && self.rid2 == other.rid2 && self.rev == other.rev
    }
}

impl Eq for MatchRecord {}

impl Ord for MatchRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rid1.cmp(&other.rid1).then(self.rid2.cmp(&other.rid2)).then(self.rev.cmp(&other.rev))
    
    }
}

impl PartialOrd for MatchRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


pub struct AllelesFasta {
    pub file: String,
    pub contigs: Vec<String>,
    pub contig_lengths: HashMap<String, u32>,
}

impl BaseTable for AllelesFasta {
    fn new(name: &String) -> AllelesFasta {
        AllelesFasta { file: name.clone(), contigs: Vec::new(),
                        contig_lengths: HashMap::new()}
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

impl AllelesFasta {
    pub fn seqs(&mut self) -> anyResult<Vec<Vec<u8>>> {
        let reader = common_reader(&self.file);
        let reader = FastxReader::new(reader);
        let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();

        read_process_fastx_records(reader, 4, 2,
            |record, seq| { // runs in worker
                *seq = record.seq_lines()
                            .fold(Vec::new(), |mut s, seq| {s.extend_from_slice(seq); s});
                                // .fold(Vec::new(), |s, seq| s +  &String::from_utf8(seq.to_vec()).unwrap());
            },
            |record, seq| { // runs in main thread
                seqs.insert(record.id().unwrap().to_owned(), seq.to_owned()); 
                None::<()>
            }).unwrap();
        
        // sort

        self.contigs = seqs.keys().map(|x| x.to_string()).sorted().collect();
        self.contig_lengths = seqs.iter().map(|(k, v)| (k.clone(), v.len() as u32)).collect();
        let seqs = self.contigs.par_iter().map(|x| seqs.get(x).unwrap().clone()).collect::<Vec<Vec<u8>>>();
        Ok(seqs)
    }

    fn collect_minimizers(&mut self, k: usize, w: usize) -> anyResult<Vec<MinimizerData>> {
        let seqs = self.seqs().unwrap();
        log::info!("Load `{}` sequences.", seqs.len());
        let mut minimizers: Vec<MinimizerData> = seqs.par_iter().enumerate().flat_map(|(i, seq)| {
            let rid: u32 = i as u32;
            sketch(seq, rid, k, w)
        }).collect();

        log::info!("Collected {} minimizers.", minimizers.len());

        minimizers.par_sort();
        Ok(minimizers)
    }

    fn collect_anchors<'a>(&'a self, k: usize, w: usize, mut minimizers: Vec<MinimizerData>
                            ) -> anyResult<(Vec<Anchor>, HashMap<u32, usize>, HashMap<u32, usize>)> {

        let mut anchor: Vec<Vec<u64>> = Vec::new();
        let mz_num = minimizers.len();
        let max_occ = 100;
        let max_minimizer_count = 1000;
        let mut unique_minimizer_count_db = HashMap::new();
        let mut contig_minimizer_count_db = HashMap::new();

        let groups = minimizers.iter_mut().group_by(move |x| x.minimizer)
                    .into_iter()
                    .filter_map(|(minimizer, group)| {
                        let group: Vec<&mut MinimizerData> = group.collect();
                        let group_length = group.len();
                        
                        let contig_counts = group.iter().fold(HashMap::new(), |mut acc, x| {
                            *acc.entry(x.info.rid).or_insert(0) += 1;
                            acc
                        });

                        for (rid, count) in contig_counts {
                            *contig_minimizer_count_db.entry(rid).or_insert(0) += count;
                            if count == 1 {
                                *unique_minimizer_count_db.entry(rid).or_insert(0) += 1;
                            }
                        }

                        if group_length == 1 || group_length > max_minimizer_count {
                            return None;
                        }

                      
                        Some( if group_length > max_occ {
                            group.into_iter().take(max_occ).collect()
                        } else {
                            group
                        })

                    }).collect::<Vec<Vec<&mut MinimizerData>>>();
        log::info!("Collected {} groups.", groups.len());
        let mut anchor = groups.par_iter()
                            .flat_map(|group| {
                                let group_length = group.len();
                                let mut pairs = Vec::with_capacity(group_length * (group_length - 1));
                            

                                for i in 0..group_length {
                                    for j in (i+1)..group_length {
                                        if group[i].info.rid == group[j].info.rid {
                                            continue;
                                        }

                                        let rev_xor = group[i].info.rev ^ group[j].info.rev;
                                        let anchor1 = Anchor { rid1: group[i].info.rid, 
                                                                    pos1: group[i].info.pos,
                                                                    rid2: group[j].info.rid,
                                                                    pos2: group[j].info.pos, 
                                                                    rev: rev_xor };
                                        pairs.push(anchor1);

                                    }
                                }
                                pairs.into_par_iter()
                                
                            }).collect::<Vec<_>>();
        drop(minimizers);                 

        let anchor2 = anchor.par_iter().map(|group| {
            group.swap()
        }).collect::<Vec<_>>();
        anchor.extend(anchor2);

        anchor.par_sort();
        log::info!("Collected {} anchors.", anchor.len());
        Ok((anchor, unique_minimizer_count_db, contig_minimizer_count_db))
                        
    }


    // longest increasing sequences 
    pub fn lis<'a>(&'a self, anchor: &'a Vec<&Anchor>) -> usize {
        let n = anchor.len();
        let mut dp: Vec<usize> = vec![1; n];
        for i in 0..n {
            for j in 0..i {
                if anchor[i].pos2 > anchor[j].pos2 {
                    dp[i] = dp[i].max(dp[j] + 1);
                }
            }
        }

        dp.into_iter().max().unwrap()

    }

    // https://rosettacode.org/wiki/Longest_increasing_subsequence
    pub fn lis_optimized(&self, anchors: &Vec<&Anchor>) -> usize {
       
        let x: Vec<_> = anchors.iter().map(|x| x.pos2).collect();
        let n = x.len();
        let mut m = vec![usize::MAX; n + 1];
        let mut p = vec![usize::MAX; n];
        let mut l = 0;

        for i in 0..n {
            let mut lo = 1;
            let mut hi = l;

            while lo <= hi {
                let mid = (lo + hi) / 2;

                let mid_index = m[mid];
                if mid_index != usize::MAX {
                    let mid_value = x[mid_index];
                    if mid_value < x[i] {
                        lo = mid + 1;
                    } else {
                        hi = mid - 1;
                    }
                } else {
                    break; 
                }
            }

            let new_l = lo;
            p[i] = m.get(new_l - 1).copied().unwrap_or(usize::MAX);
            m[new_l] = i;

            if new_l > l {
                l = new_l;
            }
        }

        l
    }

    fn calculate_simularity<'a>(&'a self, anchor: &'a Vec<Anchor>, 
                                contig_minimizer_count_db: &HashMap<u32, usize>,
                                k: usize, min_sim: f64) -> Vec<MatchRecord> {
        let min_cnt = 5;
        
        let mut anchors = anchor.iter().group_by(|x| (x.rid1, x.rid2, x.rev))
                    .into_iter()
                    .filter_map(|(rid, group)|{
                        let mut group: Vec<_> = group.collect();
                        if group.len() < min_cnt {
                            return None;
                        }

                        match group[0].rev {
                            1 => group.sort_unstable_by_key(|x| std::cmp::Reverse(x.pos1)),
                            _ => group.sort_unstable_by_key(|x| x.pos1),
                        }
                        
                        Some(group)
                    } ).collect::<Vec<Vec<&Anchor>>>();
       
        let res = anchors.par_iter().filter_map(|group| {
            let m = self.lis_optimized(group);
            if m < min_cnt {
                return None;
            }

            let rid1 = group[0].rid1;
            let rid2 = group[0].rid2;

            let n1 = contig_minimizer_count_db.get(&rid1).unwrap();
            let n2 = contig_minimizer_count_db.get(&rid2).unwrap();

            let similarity = (2.0 * (m as f64 / (n1 + n2) as f64)).powf(1.0 / k as f64);
            if similarity >= min_sim {
                
                Some(MatchRecord { rid1: rid1, rid2: rid2, rev: group[0].rev, 
                            mz1: *n1 as u32, mz2: *n2 as u32, 
                            mz_shared: m as u32, similarity: similarity })
            } else {
                return None;
            }
            
            }).collect::<Vec<MatchRecord>>();
        
        log::info!("Collected {} matches.", res.len());

        res 
    }


    pub fn symmetric_and_filter(&self, matches: &mut Vec<MatchRecord>) -> Vec<MatchRecord>  {
        
        let mut filtered_matches: HashMap::<(u32, u32), &MatchRecord> = HashMap::new();

        for record in matches.iter() {
            let key = (record.rid1, record.rid2);
            filtered_matches.entry(key).and_modify(|e| {
                if e.similarity < record.similarity {
                    *e = record;
                }
            }).or_insert(record);
        }

        let mut res = Vec::with_capacity(filtered_matches.len() * 2);
        for record in filtered_matches.values() {
            res.push(**record);
            let record2 = MatchRecord { rid1: record.rid2, rid2: record.rid1, rev: record.rev, 
                                        mz1: record.mz2, mz2: record.mz1, mz_shared: record.mz_shared, 
                                        similarity: record.similarity };
            res.push(record2);

        }

        res 
    }

    pub fn run(&mut self, k: usize, w: usize,
                m: f64, output: &String) {
        // get time 
        let start = std::time::Instant::now();
        let minimizers = self.collect_minimizers(k, w).unwrap();
        
        let (anchor, unique_minimizer_count_db, contig_minimizer_count_db) = self.collect_anchors(k, w, minimizers).unwrap();
    
        let mut matches = self.calculate_simularity(&anchor, &contig_minimizer_count_db, k, m);
        drop(anchor);
        let matches = self.symmetric_and_filter(&mut matches);

        let mut writer = common_writer(output);
        
        for (i, record) in matches.iter().enumerate() {
            let contig1 = self.contigs.get(record.rid1 as usize).unwrap();
            let contig2 = self.contigs.get(record.rid2 as usize).unwrap();
            let strand = if record.rev == 1 { "-1" } else { "1" };
            writer.write_all(format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                                i, i, 
                                contig1, contig2,  
                                record.mz1, record.mz2, 
                                record.mz_shared,
                                record.similarity, strand).as_bytes()).unwrap();
        }
        writer.flush().unwrap();
        log::info!("Successful output allele table to `{}`.", output);
        
    }
}


