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
    // pub contig_pairs: Vec<ContigPair2<'a>>,
    pub normalization_method: String,
    pub allelic_counts: u32,
    pub potential_cross_allelic_counts: u32,
    pub cross_allelic_counts: u32,

}

impl KPruner {
    pub fn new(alleletable: &String, contacts: &String, prunetable: &String,
                    normalization_method: &String) -> KPruner {
        
        // let mut count_re = CountRE::new(count_re);
        // count_re.parse();
        let mut contacts = Contacts2::new(contacts);
        contacts.parse();
      
        let mut alleletable = AlleleTable2::new(alleletable);
        alleletable.allele_records = alleletable.allele_records().unwrap();
        
        
        // let contacts_data = contacts.to_data(&unique_min, normalization_method);
        // let contig_pairs: Vec<ContigPair2> = contacts_data.keys().cloned().collect();
    
        KPruner {
            alleletable: alleletable,
            contacts: contacts,
            // contacts_data: contacts_data,
            prunetable: PruneTable::new(prunetable),
            // contig_pairs: contig_pairs,
            normalization_method: normalization_method.clone(),
            allelic_counts: 0,
            potential_cross_allelic_counts: 0,
            cross_allelic_counts: 0,
        }
    }

    // pub fn allelic(&mut self, whitehash: &HashSet<&String>) -> HashSet<ContigPair2> {
    //     log::info!("Starting allelic identification ...");
    //     log::set_max_level(log::LevelFilter::Off);
    //     let mut allelic_contig_pairs = self.alleletable.get_allelic_contig_pairs(); 
    //     log::set_max_level(log::LevelFilter::Info);
        // let prune_allelic_contig_pairs = self.contig_pairs.iter()
        //                                                 .filter(
        //                                                     |x| allelic_contig_pairs.contains(x))
        //                                                 .cloned()
        //                                                 .collect::<Vec<ContigPair>>();

        
        // self.contig_pairs.retain(|x| !allelic_contig_pairs.contains(x));
        // if whitehash.len() > 0 {
        //     allelic_contig_pairs.retain(|x| whitehash.contains(x.Contig1) && whitehash.contains(x.Contig2));
        // }

        // remove allelic_contig_pairs from self.contacts, which not contain it
        // allelic_contig_pairs.retain(|x| self.contacts_data.contains_key(x));

        
        // for contig_pair in allelic_contig_pairs.iter() {
        //     self.contacts_data.remove(contig_pair);
        // }
    //     self.allelic_counts += allelic_contig_pairs.len() as u32;
    //     log::info!("Allelic contig pairs: {}", self.allelic_counts);
    //     allelic_contig_pairs

    // }

    // deprecated
    // pub fn potential_cross_allelic(&mut self, whitehash: &HashSet<&String>) -> HashSet<ContigPair2> {
    //     let cis_data = self.contacts_data.par_iter().filter(|(contig_pair, _)| {
    //         contig_pair.Contig1 == contig_pair.Contig2
    //     }).map(|(contig_pair, count)| {
    //         (contig_pair.Contig1, *count)
    //     }).collect::<HashMap<&String, f64>>();
    //     // filter contig pairs that contig1 == contig2 
    //     self.contig_pairs.retain(|x| x.Contig1 != x.Contig2);
    //     if whitehash.len() > 0 {
    //         self.contig_pairs.retain(|x| whitehash.contains(&x.Contig1) && whitehash.contains(&x.Contig2));
    //     }

    //     let potential_cross_allelic: HashSet<ContigPair2> = self.contig_pairs
    //                                                 .par_iter()
    //                                                 .filter_map(|contig_pair| {
    //         let count1 = cis_data.get(&contig_pair.Contig1).unwrap_or(&0.0);
    //         let count2 = cis_data.get(&contig_pair.Contig2).unwrap_or(&0.0);
    //         if count1 == &0.0 && count2 == &0.0 {
    //             Some(contig_pair.clone())
    //         } else {
    //             None
    //         }
          
    //     }).collect();

    //     self.contig_pairs.retain(|x| !potential_cross_allelic.contains(x));
    //     // remove potential_cross_allelic from self.contacts
    //     // for contig_pair in potential_cross_allelic.iter() {
    //     //     self.contacts.remove(contig_pair);
    //     // }
    //     log::info!("Potential cross allelic contig pairs: {}", potential_cross_allelic.len());
    //     self.potential_cross_allelic_counts += potential_cross_allelic.len() as u32; 
    //     potential_cross_allelic
    // }

    pub fn prune(&mut self, method: &str, 
                whitehash: &HashSet<&String>,
                mut writer: &mut Box<dyn Write + Send> ) {

        log::info!("Starting allelic identification ...");
        log::set_max_level(log::LevelFilter::Off);
        let unique_min = self.alleletable.header.to_unique_minimizer_density();
        log::set_max_level(log::LevelFilter::Info);
        let mut contacts_data = self.contacts.to_data(&unique_min, &self.normalization_method);

        // filter contact data
        if !whitehash.is_empty() {
            contacts_data.retain(|x, y| whitehash.contains(x.Contig1) &&  whitehash.contains(x.Contig2));
        }
        
        let mut contig_pairs: Vec<&ContigPair2> = contacts_data.keys().collect();
        let mut allelic_contig_pairs = self.alleletable.get_allelic_contig_pairs();

        contig_pairs.retain(|x| !allelic_contig_pairs.contains(x));

        if !whitehash.is_empty() {
            allelic_contig_pairs.retain(|x| whitehash.contains(x.Contig1) && whitehash.contains(x.Contig2));
        }

        

        // remove allelic_contig_pairs from self.contacts, which not contain it
        allelic_contig_pairs.retain(|x| contacts_data.contains_key(x));
        // for contig_pair in allelic_contig_pairs.iter() {
        //     contacts_data.remove(contig_pair);
        // }

        self.allelic_counts += allelic_contig_pairs.len() as u32;
        log::info!("Allelic contig pairs: {}", self.allelic_counts);
    
        log::info!("Starting cross-allelic identification ...");
        
        let mut allelic_contigs = match method {
            "precise" => self.alleletable.get_allelic_contigs_precise(whitehash),
            // "multiple" => self.alleletable.get_allelic_contig_groups_by_cliques(whitehash),
            _ => self.alleletable.get_allelic_contigs(method, whitehash),
        }; 



        // // remove both contig1 and contig2 not in whitehash
        // if whitehash.len() > 0 {
        //     contig_pairs.retain(|x| whitehash.contains(&x.Contig1) 
        //                                 && whitehash.contains(&x.Contig2));
        //     contig_pairs.retain(|x| allelic_contigs.contains_key(x.Contig1) 
        //                                 && allelic_contigs.contains_key(x.Contig2));
        // }
        
        // // filter contig pairs that contig1 == contig2 
        // contig_pairs.retain(|x| x.Contig1 != x.Contig2);

        if !whitehash.is_empty() {
            contig_pairs.retain(|x| 
                whitehash.contains(&x.Contig1) && 
                whitehash.contains(&x.Contig2) &&
                allelic_contigs.contains_key(&x.Contig1) &&
                allelic_contigs.contains_key(&x.Contig2) &&
                x.Contig1 != x.Contig2
            );
        } else {
            contig_pairs.retain(|x| 
                allelic_contigs.contains_key(&x.Contig1) &&
                allelic_contigs.contains_key(&x.Contig2) &&
                x.Contig1 != x.Contig2
            );
        }
        let cross_allelic: Vec<&&ContigPair2> =  contig_pairs
                                                    .par_iter()
                                                    .filter_map(|contig_pair| {

            let alleles1 = allelic_contigs.get(contig_pair.Contig1)?;
            let alleles2 = allelic_contigs.get(contig_pair.Contig2)?;
         
            let mut is_weak = false;
                                                        
            'outer: for init_group1 in alleles1.iter() {
                 
                'inner: for init_group2 in alleles2.iter() {
                    let mut group1_length = init_group1.len();
                    let mut group2_length = init_group2.len();
                    
                    if group1_length <= 1 && group2_length <= 1 {
                        break 'outer;
                    }

                    let longer_group1 = group1_length > group2_length;
                    let (contig1, contig2) = if longer_group1 {
                        (&contig_pair.Contig2, &contig_pair.Contig1)
                    } else {
                        (&contig_pair.Contig1, &contig_pair.Contig2)
                    };
    
                    let (group1, group2) = if longer_group1 {
                        (&init_group2, &init_group1)
                    } else {
                        (&init_group1, &init_group2)
                    };
                   
                    // index of contig1 in group1 
                    let idx1 = match group1.iter().position(|x| x == contig1) {
                        Some(v) => v,
                        None => continue 'inner,
                    };
                        
                    // index of contig2 in group2
                    let idx2 = match group2.iter().position(|x| x == contig2) {
                        Some(v) => v,
                        None => continue 'inner,
                    };

                    if group1_length > group2_length {
                        std::mem::swap(&mut group1_length, &mut group2_length)
                    }
                  

                    let mut matrix = Matrix::new(group1_length, group2_length, OrderedFloat(0.0));
                    matrix.iter_mut().enumerate().for_each(|(index, element) | {
                        let i = index / group2_length;
                        let j = index % group2_length;
            
                        let mut tmp_contig_pair = ContigPair2::new(&group1[i], &group2[j]);
                        tmp_contig_pair.order();

                        if let Some(value) = contacts_data.get(&tmp_contig_pair) {
                            *element = OrderedFloat(*value);
                        }
                    });

                    let assignments = maximum_bipartite_matching(matrix);
                
                    if assignments[idx1] != idx2 {
                        is_weak = true;
                        break 'outer;
                    } 
                
                }
            } 
             
           
            if is_weak {
                Some(contig_pair)
            } else {
                None
            }
                                                                               
        }).collect();

            
        
        // let cross_allelic = cross_allelic.into_iter().cloned().collect::<Vec<ContigPair>>();
        self.cross_allelic_counts += cross_allelic.len() as u32;
        log::info!("Cross allelic contig pairs: {}", self.cross_allelic_counts);
        
        let mut allelic_record_hashmap = self.alleletable.get_allelic_record_by_contig_pairs();
        
        let mut buffer = Vec::new();
        for contig_pair in allelic_contig_pairs.iter() {
            let record = allelic_record_hashmap.get(contig_pair).unwrap();
            write!(buffer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
                            contig_pair.Contig1, contig_pair.Contig2,
                            record.mz1, record.mz2, record.mz_shared,
                            record.similarity, 0).unwrap();
        }

        for contig_pair in cross_allelic.iter() {
            write!(buffer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\n", 
            contig_pair.Contig1, contig_pair.Contig2,
            0, 0, 0, 0, 1).unwrap();
        }

        writer.write_all(&buffer).unwrap();
        
    }

    // pub fn prune(&mut self, method: &str, whitehash: &HashSet<&String>) {
    //     log::info!("Using `{}` method to prune", method);

    //     let allelic = self.allelic(whitehash);
    //     // let potential_cross_allelic = self.potential_cross_allelic(whitehash);
    //     // let potential_cross_allelic: HashSet<ContigPair> = HashSet::new();
    //     let mut cross_allelic = self.cross_allelic(method, whitehash);
        
    //     let mut allelic_record_hashmap = self.alleletable.get_allelic_record_by_contig_pairs();

    //     for contig_pair in allelic.iter() {
    //         let record = allelic_record_hashmap.get(contig_pair).unwrap();
    //         let prune_record = PruneRecord {
    //             contig1: contig_pair.Contig1.clone(),
    //             contig2: contig_pair.Contig2.clone(),
    //             mz1: record.mz1,
    //             mz2: record.mz2,
    //             mz_shared: record.mz_shared,
    //             similarity: record.similarity,
    //             allele_type: 0,
    //         };
    //         self.prunetable.records.push(prune_record);
    //     }

        

    //     // for contig_pair in potential_cross_allelic.into_iter() {
    //     //     let prune_record = PruneRecord {
    //     //         contig1: contig_pair.Contig1.clone(),
    //     //         contig2: contig_pair.Contig2.clone(),
    //     //         mz1: 0,
    //     //         mz2: 0,
    //     //         mz_shared: 0,
    //     //         similarity: 0.0,
    //     //         allele_type: 1,
    //     //     };
    //     //     self.prunetable.records.push(prune_record);
    //     // }

    //     for contig_pair in cross_allelic.into_iter() {
    //         let prune_record = PruneRecord {
    //             contig1: contig_pair.Contig1.clone(),
    //             contig2: contig_pair.Contig2.clone(),
    //             mz1: 0,
    //             mz2: 0,
    //             mz_shared: 0,
    //             similarity: 0.0,
    //             allele_type: 1,
    //         };
    //         self.prunetable.records.push(prune_record);

    //     }
    // }
}

