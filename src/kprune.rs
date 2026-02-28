#[warn(unused_assignments)]
use anyhow::Result as anyResult;
use no_panic::no_panic;
use ordered_float::OrderedFloat;
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use rustc_hash::FxHashMap; 
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;
use pathfinding::prelude::{kuhn_munkres, Matrix, Weights};


use crate::alleles::AlleleTable2;
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair, ContigPair2 };
use crate::contacts::{ Contacts, Contacts2 };
use crate::count_re::CountRE;



// kuhn munkres algorithm
pub fn maximum_bipartite_matching(matrix: Matrix<OrderedFloat<f64>>) -> Vec<usize>  {

    let (cash_flow, assignments) = kuhn_munkres(&matrix);

    assignments
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PruneRecord {
    pub contig1: String,
    pub contig2: String,
    pub mz1: u64,
    pub mz2: u64,
    pub mz_shared: u64,
    pub similarity: f64,
    pub allele_type: u8,
}

#[derive(Debug, Clone)]
pub struct PruneTable {
    pub file: String,
    pub records: Vec<PruneRecord>,
}

impl BaseTable for PruneTable {
    fn new(name: &String) -> PruneTable {
        PruneTable { file: name.clone(),
                     records: Vec::new(), }
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

impl PruneTable {
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

    pub fn contig_pairs(&self) -> anyResult<Vec<ContigPair>> {
        let parse_result = self.parse();
        let mut rdr = match parse_result {
            Ok(v) => v,
            Err(error) => panic!("Could not parse input file: {:?}", self.file_name()),
        };

        let mut contig_pairs: Vec<ContigPair> = Vec::new();
        for result in rdr.records() {
            let record = result?;
            let contig1 = record[0].to_string();
            let contig2 = record[1].to_string();
           
            let contig_pair = ContigPair::new(contig1, contig2);
            contig_pairs.push(contig_pair);
        }

        Ok(contig_pairs)
    }

    pub fn write(&self, wtr: &mut csv::Writer<Box<dyn std::io::Write>>) {

        for record in &self.records {
            wtr.serialize(record).expect("Could not write record")
        } 
        
    }

}



pub struct KPruner {
    pub alleletable: AlleleTable2,
    pub contacts: Contacts2,
    // pub contacts_data: HashMap<ContigPair2<'a>, f64>,
    pub prunetable: PruneTable,
    pub countre: Option<CountRE>,
    // pub contig_pairs: Vec<ContigPair2<'a>>,
    pub normalization_method: String,
    pub allelic_counts: u32,
    pub potential_cross_allelic_counts: u32,
    pub cross_allelic_counts: u32,

}

impl KPruner {
    pub fn new(alleletable: &String, contacts: &String, prunetable: &String,
                    count_re: &Option<String>,
                    normalization_method: &String) -> KPruner {
        
        let count_re = match count_re {
            Some(v) => Some(CountRE::new(v)),
            None => None,
        };
        // let mut count_re = CountRE::new(count_re);
        // count_re.parse();
        let mut contacts = Contacts2::new(contacts);
        contacts.parse();
        
        log::set_max_level(log::LevelFilter::Off);
        let mut alleletable = AlleleTable2::new(alleletable);
        alleletable.allele_records = alleletable.allele_records().unwrap();
        alleletable.header = alleletable.parse_header().unwrap();
        log::set_max_level(log::LevelFilter::Info);
         
        // let contacts_data = contacts.to_data(&unique_min, normalization_method);
        // let contig_pairs: Vec<ContigPair2> = contacts_data.keys().cloned().collect();
    
        KPruner {
            alleletable: alleletable,
            contacts: contacts,
            // contacts_data: contacts_data,
            prunetable: PruneTable::new(prunetable),
            countre: count_re,
            // contig_pairs: contig_pairs,
            normalization_method: normalization_method.clone(),
            allelic_counts: 0,
            potential_cross_allelic_counts: 0,
            cross_allelic_counts: 0,
        }
    }

    pub fn prune(&mut self, method: &str, 
            whitehash: &HashSet<&String>,
            writer: &mut Box<dyn Write + Send>,
            partial_whitelist: bool
        ) {

        log::info!("Starting allelic identification ...");
        
        log::set_max_level(log::LevelFilter::Off);
        let unique_min = self.alleletable.header.to_unique_minimizer_density();
        log::set_max_level(log::LevelFilter::Info);

        let length_hash = &self.alleletable.header.contigsizes; 
        let mut contacts_data = self.contacts.to_data(&unique_min, &self.normalization_method, 
                                                    &self.countre, Some(length_hash));

        if !whitehash.is_empty() {
            if !partial_whitelist {
                contacts_data.retain(|x, _| whitehash.contains(x.Contig1) && whitehash.contains(x.Contig2));
            }
        }
     
        let mut contig_pairs: Vec<&ContigPair2> = contacts_data.keys().collect();
        let mut allelic_contig_pairs = self.alleletable.get_allelic_contig_pairs();

        contig_pairs.retain(|x| !allelic_contig_pairs.contains(x));
        if !whitehash.is_empty() {
            if !partial_whitelist {
                allelic_contig_pairs.retain(|x| whitehash.contains(x.Contig1) && whitehash.contains(x.Contig2));
            }
        }


        allelic_contig_pairs.retain(|x| contacts_data.contains_key(x));
        self.allelic_counts += allelic_contig_pairs.len() as u32;
        log::info!("Allelic contig pairs: {}", self.allelic_counts);
    
        log::info!("Starting cross-allelic identification ...");
        
        let allelic_contigs = match method {
            "precise" => 
                if !partial_whitelist {
                    self.alleletable.get_allelic_contigs_precise(whitehash)
                } else{
                    let _whitehash: HashSet<&String> = HashSet::new();
                    self.alleletable.get_allelic_contigs_precise(&_whitehash)
                },
            _ => 
                if !partial_whitelist {
                    self.alleletable.get_allelic_contigs(method, whitehash)
                } else {
                    let _whitehash: HashSet<&String> = HashSet::new();
                    self.alleletable.get_allelic_contigs(method, &_whitehash)
                },

        }; 

        let filter_func = |x: &&ContigPair2| -> bool {
            if !whitehash.is_empty() {
                if partial_whitelist {
                    if !whitehash.contains(&x.Contig1) && !whitehash.contains(&x.Contig2) { return false; }
                } else {
                    if !whitehash.contains(&x.Contig1) || !whitehash.contains(&x.Contig2) { return false; }
                }
            }
            allelic_contigs.contains_key(&x.Contig1) && 
            allelic_contigs.contains_key(&x.Contig2) && 
            x.Contig1 != x.Contig2
        };
        contig_pairs.retain(filter_func);

        let mut all_names = HashSet::new();
        for pair in contacts_data.keys() {
            all_names.insert(&pair.Contig1);
            all_names.insert(&pair.Contig2);
        }
        for (name, groups) in &allelic_contigs {
            all_names.insert(name);
            for group in groups {
                for n in group { all_names.insert(n); }
            }
        }
        
        let name_to_id: HashMap<&String, u32> = all_names.into_iter()
            .enumerate()
            .map(|(i, name)| (*name, i as u32))
            .collect();
            
        let get_id = |name: &String| *name_to_id.get(name).expect("ID mapping failed");


        let mut contacts_id_map: FxHashMap<(u32, u32), f64> = FxHashMap::with_capacity_and_hasher(
                                                                    contacts_data.len(), Default::default());
        for (pair, &val) in &contacts_data {
            let u = get_id(&pair.Contig1);
            let v = get_id(&pair.Contig2);
            contacts_id_map.insert(if u < v { (u, v) } else { (v, u) }, val);
        }

    
        let max_id = name_to_id.len();
        let mut allelic_id_vec: Vec<Vec<Vec<u32>>> = vec![vec![]; max_id];
        for (name, groups) in &allelic_contigs {
            let id = get_id(name);
            let groups_id: Vec<Vec<u32>> = groups.iter()
                .map(|g| g.iter().map(|n| get_id(n)).collect())
                .collect();
            allelic_id_vec[id as usize] = groups_id;
        }

        let cross_allelic: Vec<&&ContigPair2> = contig_pairs
            .par_iter()
            .filter_map(|contig_pair| {
                let id1 = get_id(&contig_pair.Contig1);
                let id2 = get_id(&contig_pair.Contig2);

                let alleles1 = &allelic_id_vec[id1 as usize];
                let alleles2 = &allelic_id_vec[id2 as usize];

                if alleles1.is_empty() || alleles2.is_empty() { return None; }

                let mut is_weak = false;

                'outer: for g1 in alleles1 {
                    for g2 in alleles2 {
                        let r = g1.len();
                        let c = g2.len();
                        
                        if r <= 1 && c <= 1 { continue; }

                        let swapped = r > c;
                        let (row_grp, col_grp) = if swapped { (g2, g1) } else { (g1, g2) };
                        let (row_target, col_target) = if swapped { (id2, id1) } else { (id1, id2) };

                        let row_idx_opt = row_grp.iter().position(|&x| x == row_target);
                        let col_idx_opt = col_grp.iter().position(|&x| x == col_target);

                        let row_idx = match row_idx_opt { Some(i) => i, None => continue };
                        let col_idx = match col_idx_opt { Some(i) => i, None => continue };

                        let rows = row_grp.len();
                        let cols = col_grp.len();
                        
                        let mut matrix = Matrix::new(rows, cols, OrderedFloat(0.0));
                        
                        for i in 0..rows {
                            let u = row_grp[i];
                            for j in 0..cols {
                                let v = col_grp[j];
                                let key = if u < v { (u, v) } else { (v, u) };
                                if let Some(&val) = contacts_id_map.get(&key) {
                                    matrix[(i, j)] = OrderedFloat(val);
                                }
                            }
                        }
                       

                        let assignments = maximum_bipartite_matching(matrix);
                        
                        if assignments[row_idx] != col_idx {
                            is_weak = true;
                            break 'outer;
                        }
                    }
                }

                if is_weak { Some(contig_pair) } else { None }
            })
            .collect();


        self.cross_allelic_counts += cross_allelic.len() as u32;
        log::info!("Cross allelic contig pairs: {}", self.cross_allelic_counts);
        
        let allelic_record_hashmap = self.alleletable.get_allelic_record_by_contig_pairs();
        let mut buffer = Vec::new();
        
        for contig_pair in allelic_contig_pairs.iter() {
            if let Some(record) = allelic_record_hashmap.get(contig_pair) {
                 let (m1, m2) = if record.contig1 > record.contig2 {
                     (record.mz2, record.mz1)
                 } else {
                     (record.mz1, record.mz2)
                 };
                 writeln!(buffer, "{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                            contig_pair.Contig1, contig_pair.Contig2,
                            m1, m2, record.mz_shared, record.similarity, 0).unwrap();
            }
        }

        for contig_pair in cross_allelic.iter() {
            writeln!(buffer, "{}\t{}\t0\t0\t0\t0\t1", 
            contig_pair.Contig1, contig_pair.Contig2).unwrap();
        }

        writer.write_all(&buffer).unwrap();
    }

    pub fn prune_bak(&mut self, method: &str, 
            whitehash: &HashSet<&String>,
            writer: &mut Box<dyn Write + Send> ) {

        log::info!("Starting allelic identification ...");
        
        log::set_max_level(log::LevelFilter::Off);
        let unique_min = self.alleletable.header.to_unique_minimizer_density();
        log::set_max_level(log::LevelFilter::Info);

        let length_hash = &self.alleletable.header.contigsizes; 
        let mut contacts_data = self.contacts.to_data(&unique_min, &self.normalization_method, 
                                                    &self.countre, Some(length_hash));

        if !whitehash.is_empty() {
            contacts_data.retain(|x, _| whitehash.contains(x.Contig1) &&  whitehash.contains(x.Contig2));
        }
     
        let mut contig_pairs: Vec<&ContigPair2> = contacts_data.keys().collect();
        let mut allelic_contig_pairs = self.alleletable.get_allelic_contig_pairs();

        contig_pairs.retain(|x| !allelic_contig_pairs.contains(x));
        if !whitehash.is_empty() {
            allelic_contig_pairs.retain(|x| whitehash.contains(x.Contig1) && whitehash.contains(x.Contig2));
        }


        allelic_contig_pairs.retain(|x| contacts_data.contains_key(x));
        self.allelic_counts += allelic_contig_pairs.len() as u32;
        log::info!("Allelic contig pairs: {}", self.allelic_counts);
    
        log::info!("Starting cross-allelic identification ...");
        
        let allelic_contigs = match method {
            "precise" => self.alleletable.get_allelic_contigs_precise(whitehash),
            _ => self.alleletable.get_allelic_contigs(method, whitehash),
        }; 

        let filter_func = |x: &&ContigPair2| -> bool {
            if !whitehash.is_empty() {
                if !whitehash.contains(&x.Contig1) || !whitehash.contains(&x.Contig2) { return false; }
            }
            allelic_contigs.contains_key(&x.Contig1) && 
            allelic_contigs.contains_key(&x.Contig2) && 
            x.Contig1 != x.Contig2
        };
        contig_pairs.retain(filter_func);

        let mut all_names = HashSet::new();
        for pair in contacts_data.keys() {
            all_names.insert(&pair.Contig1);
            all_names.insert(&pair.Contig2);
        }
        for (name, groups) in &allelic_contigs {
            all_names.insert(name);
            for group in groups {
                for n in group { all_names.insert(n); }
            }
        }
        
        let name_to_id: HashMap<&String, u32> = all_names.into_iter()
            .enumerate()
            .map(|(i, name)| (*name, i as u32))
            .collect();
            
        let get_id = |name: &String| *name_to_id.get(name).expect("ID mapping failed");

        let mut contacts_id_map: FxHashMap<(u32, u32), f64> = FxHashMap::with_capacity_and_hasher(
            contacts_data.len(), Default::default());
        for (pair, &val) in &contacts_data {
            let u = get_id(&pair.Contig1);
            let v = get_id(&pair.Contig2);
            contacts_id_map.insert(if u < v { (u, v) } else { (v, u) }, val);
        }

        let max_id = name_to_id.len();
        let mut allelic_id_vec: Vec<Vec<Vec<u32>>> = vec![vec![]; max_id];
        for (name, groups) in &allelic_contigs {
            let id = get_id(name);
            let groups_id: Vec<Vec<u32>> = groups.iter()
                .map(|g| g.iter().map(|n| get_id(n)).collect())
                .collect();
            allelic_id_vec[id as usize] = groups_id;
        }

        let mut group_pairs: FxHashMap<(usize, usize), (Vec<u32>, Vec<u32>)> = FxHashMap::default();
        
        for contig_pair in &contig_pairs {
            let id1 = get_id(&contig_pair.Contig1);
            let id2 = get_id(&contig_pair.Contig2);
            let alleles1 = &allelic_id_vec[id1 as usize];
            let alleles2 = &allelic_id_vec[id2 as usize];

            for g1 in alleles1 {
                for g2 in alleles2 {
                    if g1.len() <= 1 && g2.len() <= 1 { continue; }
                    let g1_ptr = g1.as_ptr() as usize;
                    let g2_ptr = g2.as_ptr() as usize;
                    let key = if g1_ptr < g2_ptr { (g1_ptr, g2_ptr) } else { (g2_ptr, g1_ptr) };
                    group_pairs.entry(key).or_insert_with(|| (g1.clone(), g2.clone()));
                }
            }
        }

    
        let group_matching_results: FxHashMap<(usize, usize), Vec<usize>> = group_pairs
            .par_iter()
            .map(|(&key, (g1, g2))| {
                let swapped = g1.len() > g2.len();
                let (r_grp, c_grp) = if !swapped { (g1, g2) } else { (g2, g1) };

                let rows = r_grp.len();
                let cols = c_grp.len();
                let mut matrix = Matrix::new(rows, cols, OrderedFloat(0.0));
                for i in 0..rows {
                    let u = r_grp[i];
                    for j in 0..cols {
                        let v = c_grp[j];
                        let k = if u < v { (u, v) } else { (v, u) };
                        if let Some(&val) = contacts_id_map.get(&k) {
                            matrix[(i, j)] = OrderedFloat(val);
                        }
                    }
                }
                (key, maximum_bipartite_matching(matrix))
            })
            .collect();

        let cross_allelic: Vec<&&ContigPair2> = contig_pairs
            .par_iter()
            .filter(|contig_pair| {
                let id1 = get_id(&contig_pair.Contig1);
                let id2 = get_id(&contig_pair.Contig2);
                let alleles1 = &allelic_id_vec[id1 as usize];
                let alleles2 = &allelic_id_vec[id2 as usize];

                for g1 in alleles1 {
                    for g2 in alleles2 {
                        if g1.len() <= 1 && g2.len() <= 1 { continue; }
                        
                        let g1_ptr = g1.as_ptr() as usize;
                        let g2_ptr = g2.as_ptr() as usize;
                        let is_swapped = g1_ptr > g2_ptr;
                        let key = if !is_swapped { (g1_ptr, g2_ptr) } else { (g2_ptr, g1_ptr) };
                        
                        if let Some(assignments) = group_matching_results.get(&key) {
                            let swapped_in_matrix = g1.len() > g2.len();
                            let (row_grp, col_grp) = if !swapped_in_matrix { (g1, g2) } else { (g2, g1) };
                           
                            let (row_target, col_target) = if !swapped_in_matrix { (id1, id2) } else { (id2, id1) };

                            let row_idx = row_grp.iter().position(|&x| x == row_target).unwrap();
                            let col_idx = col_grp.iter().position(|&x| x == col_target).unwrap();

                            if assignments[row_idx] != col_idx {
                                return true; 
                            }
                        }
                    }
                }
                false
            })
            .collect();


        self.cross_allelic_counts += cross_allelic.len() as u32;
        log::info!("Cross allelic contig pairs: {}", self.cross_allelic_counts);
        
        let allelic_record_hashmap = self.alleletable.get_allelic_record_by_contig_pairs();
        let mut buffer = Vec::new();
        
        for contig_pair in allelic_contig_pairs.iter() {
            if let Some(record) = allelic_record_hashmap.get(contig_pair) {
                 let (m1, m2) = if record.contig1 > record.contig2 {
                     (record.mz2, record.mz1)
                 } else {
                     (record.mz1, record.mz2)
                 };
                 writeln!(buffer, "{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                            contig_pair.Contig1, contig_pair.Contig2,
                            m1, m2, record.mz_shared, record.similarity, 0).unwrap();
            }
        }

        for contig_pair in cross_allelic.iter() {
            writeln!(buffer, "{}\t{}\t0\t0\t0\t0\t1", 
            contig_pair.Contig1, contig_pair.Contig2).unwrap();
        }

        writer.write_all(&buffer).unwrap();
    }

    
}

