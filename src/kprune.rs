#[warn(unused_assignments)]
use anyhow::Result as anyResult;
use no_panic::no_panic;
use ordered_float::OrderedFloat;
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;
use rayon::{ ThreadPoolBuilder, ThreadPool };
use pathfinding::prelude::{kuhn_munkres, Matrix, Weights};

use crate::alleles::AlleleTable;
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };
use crate::contacts::Contacts;
use crate::count_re::CountRE;

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

#[derive(Debug)]
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
    pub alleletable: AlleleTable,
    pub contacts: HashMap<ContigPair, f64>,
    pub prunetable: PruneTable,
    pub contig_pairs: Vec<ContigPair>,
    pub allelic_counts: u32,
    pub potential_cross_allelic_counts: u32,
    pub cross_allelic_counts: u32,
}

impl KPruner {
    pub fn new(alleletable: &String, contacts: &String, prunetable: &String,
                    nomalization_method: &String) -> KPruner {
        
        // let mut count_re = CountRE::new(count_re);
        // count_re.parse();
        let mut contacts = Contacts::new(contacts);
        contacts.parse();
      
        let mut alleletable = AlleleTable::new(alleletable);
        let unique_min = alleletable.header.to_unique_minimizer_density();
        let contact_data = contacts.to_data(&unique_min, nomalization_method);
        let contig_pairs: Vec<ContigPair> = contact_data.keys().cloned().collect();
    
        KPruner {
            alleletable: alleletable,
            contacts: contact_data,
            prunetable: PruneTable::new(prunetable),
            contig_pairs: contig_pairs,
            allelic_counts: 0,
            potential_cross_allelic_counts: 0,
            cross_allelic_counts: 0,
        }
    }

    pub fn allelic(&mut self, whitehash: &HashSet<String>) -> HashSet<ContigPair> {
        log::info!("Starting allelic identification");
        let mut allelic_contig_pairs = self.alleletable.get_allelic_contig_pairs(); 
        // let prune_allelic_contig_pairs = self.contig_pairs.iter()
        //                                                 .filter(
        //                                                     |x| allelic_contig_pairs.contains(x))
        //                                                 .cloned()
        //                                                 .collect::<Vec<ContigPair>>();

        
        self.contig_pairs.retain(|x| !allelic_contig_pairs.contains(x));
        if whitehash.len() > 0 {
            allelic_contig_pairs.retain(|x| whitehash.contains(&x.Contig1) && whitehash.contains(&x.Contig2));
        }

        // remove allelic_contig_pairs from self.contacts, which not contain it
        allelic_contig_pairs.retain(|x| self.contacts.contains_key(x));

        
        for contig_pair in allelic_contig_pairs.iter() {
            self.contacts.remove(contig_pair);
        }
        self.allelic_counts += allelic_contig_pairs.len() as u32;
        log::info!("Allelic contig pairs: {}", allelic_contig_pairs.len());
        allelic_contig_pairs

    }

    // deprecated
    pub fn potential_cross_allelic(&mut self, whitehash: &HashSet<String>) -> HashSet<ContigPair> {
        let cis_data = self.contacts.par_iter().filter(|(contig_pair, _)| {
            contig_pair.Contig1 == contig_pair.Contig2
        }).map(|(contig_pair, count)| {
            (contig_pair.Contig1.clone(), *count)
        }).collect::<HashMap<String, f64>>();
        // filter contig pairs that contig1 == contig2 
        self.contig_pairs.retain(|x| x.Contig1 != x.Contig2);
        if whitehash.len() > 0 {
            self.contig_pairs.retain(|x| whitehash.contains(&x.Contig1) && whitehash.contains(&x.Contig2));
        }

        let potential_cross_allelic: HashSet<ContigPair> = self.contig_pairs
                                                    .par_iter()
                                                    .filter_map(|contig_pair| {
            let count1 = cis_data.get(&contig_pair.Contig1).unwrap_or(&0.0);
            let count2 = cis_data.get(&contig_pair.Contig2).unwrap_or(&0.0);
            if count1 == &0.0 && count2 == &0.0 {
                Some(contig_pair.clone())
            } else {
                None
            }
          
        }).collect();

        self.contig_pairs.retain(|x| !potential_cross_allelic.contains(x));
        // remove potential_cross_allelic from self.contacts
        // for contig_pair in potential_cross_allelic.iter() {
        //     self.contacts.remove(contig_pair);
        // }
        log::info!("Potential cross allelic contig pairs: {}", potential_cross_allelic.len());
        self.potential_cross_allelic_counts += potential_cross_allelic.len() as u32; 
        potential_cross_allelic
    }

    pub fn cross_allelic(&mut self, method: &str, whitehash: &HashSet<String>) -> Vec<ContigPair> {
        
        let mut allelic_contigs = self.alleletable.get_allelic_contigs(method, whitehash);
    

        // remove both contig1 and contig2 not in whitehash
        if whitehash.len() > 0 {
            self.contig_pairs.retain(|x| whitehash.contains(&x.Contig1) && whitehash.contains(&x.Contig2));
            self.contig_pairs.retain(|x| allelic_contigs.contains_key(&x.Contig1) && allelic_contigs.contains_key(&x.Contig2));
        }
        
      
        // filter contig pairs that contig1 == contig2 
        self.contig_pairs.retain(|x| x.Contig1 != x.Contig2);

        let cross_allelic: Vec<&ContigPair> =  self.contig_pairs
                                                    .par_iter()
                                                    .filter_map(|contig_pair| {

            let alleles1 = allelic_contigs.get(&contig_pair.Contig1)?;
            let alleles2 = allelic_contigs.get(&contig_pair.Contig2)?;
                                     
                            
            let mut group1 = std::iter::once(&contig_pair.Contig1).chain(alleles1.iter()).collect::<Vec<&String>>();
            let mut group2 = std::iter::once(&contig_pair.Contig2).chain(alleles2.iter()).collect::<Vec<&String>>();
            if group1.len() > group2.len() {
                std::mem::swap(&mut group1, &mut group2);
            }
            let mut matrix = Matrix::new(group1.len(), group2.len(), OrderedFloat(0.0));
            matrix.par_iter_mut().enumerate().for_each(|(index, element) | {
                let i = index / group2.len();
                let j = index % group2.len();
    
                let mut tmp_contig_pair = ContigPair::new(group1[i].to_string(), group2[j].to_string());
                tmp_contig_pair.order();
                let score = self.contacts.get(&tmp_contig_pair).unwrap_or(&0.0);

                *element = OrderedFloat(*score);
            });
                    
           

            // if &contig_pair.Contig1 == &"5G.ctg44".to_string() && &contig_pair.Contig2 == &"5H.ctg72".to_string() {
            //     println!("{:?}", &contig_pair);
            //     dbg!("Same contig: {:?}", &matrix);

            // }
            
            
            // match matrix[(0, 0)] == OrderedFloat(0.0) {
            //     true => {
            //         Some(contig_pair)
            //     },
            //     false => {
            //         let assignments = maximum_bipartite_matching(matrix);
            //         if assignments[0] != 0 {
            //             Some(contig_pair)
            //         } else {
            //             None
            //         }
            //     }
            // }
            // if &contig_pair.Contig1 == &"1A.ctg151".to_string() && &contig_pair.Contig2 == &"1B.ctg189".to_string() {
            //     dbg!("{:?}", group1);
            //     dbg!("{:?}", group2);
            //     dbg!("Same contig: {:?}", &matrix);
            // }
            let assignments = maximum_bipartite_matching(matrix);
            // if &contig_pair.Contig1 == &"1A.ctg151".to_string() && &contig_pair.Contig2 == &"1B.ctg189".to_string() {
            //     dbg!(&assignments);
              
            // }     
            if assignments[0] != 0 {
                Some(contig_pair)
            } else {
                None
            }

            }).collect();


        log::info!("Cross allelic contig pairs: {}", cross_allelic.len());
        let cross_allelic = cross_allelic.into_iter().cloned().collect::<Vec<ContigPair>>();
        self.cross_allelic_counts += cross_allelic.len() as u32;
        cross_allelic
    }

    pub fn prune(&mut self, method: &str, whitehash: &HashSet<String>) {
        log::info!("Using {} method to prune", method);

        let allelic = self.allelic(whitehash);
        let potential_cross_allelic = self.potential_cross_allelic(whitehash);
        let cross_allelic = self.cross_allelic(method, whitehash);
        
        let allelic_record_hashmap = self.alleletable.get_allelic_record_by_contig_pairs();

        for contig_pair in allelic.iter() {
            let record = allelic_record_hashmap.get(contig_pair).unwrap();
            let prune_record = PruneRecord {
                contig1: contig_pair.Contig1.clone(),
                contig2: contig_pair.Contig2.clone(),
                mz1: 0,
                mz2: 0,
                mz_shared: record.mz_unique,
                similarity: record.similarity,
                allele_type: 0,
            };
            self.prunetable.records.push(prune_record);
        }

        for contig_pair in potential_cross_allelic.into_iter() {
            let prune_record = PruneRecord {
                contig1: contig_pair.Contig1.clone(),
                contig2: contig_pair.Contig2.clone(),
                mz1: 0,
                mz2: 0,
                mz_shared: 0,
                similarity: 0.0,
                allele_type: 1,
            };
            self.prunetable.records.push(prune_record);
        }

        for contig_pair in cross_allelic.into_iter() {
            let prune_record = PruneRecord {
                contig1: contig_pair.Contig1.clone(),
                contig2: contig_pair.Contig2.clone(),
                mz1: 0,
                mz2: 0,
                mz_shared: 0,
                similarity: 0.0,
                allele_type: 1,
            };
            self.prunetable.records.push(prune_record);

        }
    }
}

// kuhn munkres algorithm
pub fn maximum_bipartite_matching(matrix: Matrix<OrderedFloat<f64>>) -> Vec<usize>  {

    let (cash_flow, assignments) = kuhn_munkres(&matrix);

    assignments
}