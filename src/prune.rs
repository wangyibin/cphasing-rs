#[warn(unused_assignments)]
use anyhow::Result as anyResult;
use ordered_float::OrderedFloat;
use pathfinding::prelude::{ kuhn_munkres, Matrix };
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::sync::Mutex;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;

use crate::alleles::{ AlleleTable, AlleleRecord, AlleleStrandRecord, AlleleStrandTable };
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair, ContigPair2 };
use crate::contacts::{ Contacts, Contacts2 };

use crate::kprune::{ PruneRecord, PruneTable, maximum_bipartite_matching };


pub struct PruneLineRecord<'a> {
    pub allelic_data: HashMap<&'a String, HashSet<&'a String>>,
    pub cross_allelic_data: HashMap<&'a String, HashSet<&'a String>>
}

#[derive(Debug)]
pub struct Pruner {
    pub alleletable: AlleleTable,
    pub allele_strand_table: AlleleStrandTable,
    // pub prunetable: PruneTable,
    pub contacts: Contacts2,
    pub allcontigs: HashSet<String>,
    // pub contig1_to_contig2: HashMap<String, HashSet<String>>,
    // pub contacts_data: HashMap<ContigPair2<'a>, f64>,
    pub normalization_method: String,
    pub allelic_counts: u32,
    pub cross_allelic_counts: u32, 
}


impl Pruner {
    pub fn new(alleletable: &String, 
                    allele_strand_table: &String, 
                    contacts: &String, 
                    // prunetable: &String, 
                    normalization_method: &String) -> Self {
        let mut contacts = Contacts2::new(contacts);
        contacts.parse(); 
        let mut allcontigs = contacts.contigs();
        
        // let mut contig1_to_contig2 = contacts.contig1_to_contig2();
        let mut alleletable = AlleleTable::new(alleletable);
        alleletable.allele_records = alleletable.allele_records().unwrap();
        let mut allele_strand_table = AlleleStrandTable::new(allele_strand_table);
       
        // let mut prunetable = PruneTable::new(prunetable);
        // let unique_min: HashMap<&String, f64> = HashMap::new();
        // let contacts_data = contacts.to_data(&unique_min, normalization_method);
       
        Pruner {
            alleletable: alleletable,
            allele_strand_table: allele_strand_table,
            // prunetable: prunetable,
            contacts: contacts,
            allcontigs: allcontigs,
            // contig1_to_contig2: contig1_to_contig2,
            // contacts_data: contacts_data,
            normalization_method: normalization_method.clone(),
            allelic_counts: 0,
            cross_allelic_counts: 0,
        }
    }

    pub fn prune(&mut self, whitehash: &HashSet<&String>, mut writer: &mut Box<dyn Write> ) {
        let unique_min: HashMap<String, f64> = HashMap::new();
        let contact_data = self.contacts.to_data(&unique_min, &self.normalization_method);

        // self.allcontigs to sorted Vec 
        let mut allcontigs = self.allcontigs.iter().collect::<Vec<&String>>();
        allcontigs.sort();
        
        self.allele_strand_table.get_allele_strand_records();
        let allele_strand_records = &self.allele_strand_table.allele_strand_records;
        let allele_strand_info = &self.allele_strand_table.get_info();


        let all_prune_line_records = self.alleletable.allele_records.par_iter().map(|record| {
            let mut remove_db: HashMap<&String, HashSet<&String>> = HashMap::new();
            let mut retain_db: HashMap<&String, &String> = HashMap::new();
            let mut count_db: HashMap<&String, f64> = HashMap::new();
            let mut prune_line_record = PruneLineRecord {
                allelic_data: HashMap::new(),
                cross_allelic_data: HashMap::new(),
            };
            for i in 0..(record.data.len() - 1) {
                'inner: for j in (i+1)..record.data.len() {
                    
                    let mut contig1 = &record.data[i]; 
                    let mut contig2 = &record.data[j];
                    if contig1 == contig2 {
                        continue 'inner;
                    }

                    if whitehash.len() != 0 {
                        if !whitehash.contains(contig1) || !whitehash.contains(contig2) {
                            continue 'inner;
                        }
                    }

                    if *contig1 > *contig2 {
                        std::mem::swap(&mut contig1, &mut contig2);
                    }
                    
                    prune_line_record.allelic_data.entry(contig1).or_insert_with(HashSet::new).insert(contig2);
                    remove_db.entry(contig1).or_insert_with(HashSet::new).insert(contig2);
                    
                }
            }
           
            if remove_db.len() == 0 {
                return prune_line_record;
            }

            for contig in record.data.iter() {
                
                'inner: for mut query_contig in allcontigs.iter() {
                    if &contig == query_contig {
                        continue 'inner;
                    }

                    if whitehash.len() != 0 {
                        if !whitehash.contains(contig) || !whitehash.contains(query_contig) {
                            continue 'inner;
                        }
                    }

                    let mut contig1 = &contig;
                    let mut contig2 = query_contig;
                    let mut flag = false;
                     
                    if **contig1 > **contig2 {
                        // println!("{:?} {:?}", contig1, contig2);
                        std::mem::swap(&mut contig1, &mut contig2);
                        // println!("{:?} {:?}", contig1, contig2);
                        flag = true;
                    }
                    
                    let mut contig_pair = ContigPair2::new(contig1, contig2);
                    // println!("{:?}", contig_pair);
                    // contig_pair.order();
                    // println!("{:?}", contig_pair);
                    if !contact_data.contains_key(&contig_pair) {
                        continue 'inner;
                    }

                    if remove_db.get(contig_pair.Contig1).map_or(false, |v| v.contains(contig_pair.Contig2)) {
                        continue 'inner;
                    }
                    let count = contact_data[&contig_pair];
                   
                    if !retain_db.contains_key(query_contig) {
                        retain_db.entry(query_contig).or_insert(contig);
                        count_db.entry(query_contig).or_insert(count);

                    } else {
                        let old_contig1 = retain_db[query_contig];
                        let old_count = count_db[query_contig];
                        if count > old_count {
                            
                            if query_contig > &old_contig1 {
                                prune_line_record.cross_allelic_data.entry(
                                    query_contig).or_insert_with(HashSet::new).insert(old_contig1);
                            } else {
                                prune_line_record.cross_allelic_data.entry(
                                    old_contig1).or_insert_with(HashSet::new).insert(query_contig);
                            }
                           
                            // if (old_contig1 == &String::from("1A.ctg10") || contig == &String::from("1A.ctg10")) && query_contig == &&String::from("1A.ctg3") {

                            //     println!("old_contig1: {:?}", old_contig1);
                            //     println!("old_count: {:?}", old_count);
                            //     println!("count: {:?}", count);
                            //     println!("new_contig1: {:?}", contig);
                            //     println!("{:?}", prune_line_record.cross_allelic_data.get(&String::from("1A.ctg3")));
                            
                            // }
                                
                            retain_db.insert(query_contig, contig);
                            count_db.insert(query_contig, count);
                        } else {
                            // if (old_contig1 == &String::from("1A.ctg10") || contig == &String::from("1A.ctg10")) && query_contig == &&String::from("1A.ctg3") {
                            //     println!("old_contig1: {:?}", old_contig1);
                            //     println!("contig1: {:?}", contig);
                            //     println!("old_count: {:?}", old_count);
                            //     println!("count: {:?}", count);
                            //     println!("{:?}", prune_line_record.cross_allelic_data.get(&String::from("1A.ctg3")));
                                
                            
                            // }
                            if &contig < query_contig {
                                prune_line_record.cross_allelic_data.entry(
                                    query_contig).or_insert_with(HashSet::new).insert(contig);
                            } else {
                                prune_line_record.cross_allelic_data.entry(
                                    contig).or_insert_with(HashSet::new).insert(query_contig);
                            }
                          
                        }
                    }
                    
                }

            }
          
            prune_line_record

        }).collect::<Vec<PruneLineRecord>>();


        let mut allelic_contig_pairs: HashSet<ContigPair2> = HashSet::new();
        let mut cross_allelic_contig_pairs: HashSet<ContigPair2> = HashSet::new();
        for prune_line_record in all_prune_line_records.iter() {
            for (contig1, contig2s) in prune_line_record.allelic_data.iter() {
                for contig2 in contig2s.iter() {
                    let mut contig_pair = ContigPair2::new(contig1, contig2);
                    contig_pair.order();
                    allelic_contig_pairs.insert(contig_pair);
            
                }
            }
            for (contig1, contig2s) in prune_line_record.cross_allelic_data.iter() {
                for contig2 in contig2s.iter() {
                    let mut contig_pair = ContigPair2::new(contig1, contig2);
                    contig_pair.order();
                    cross_allelic_contig_pairs.insert(contig_pair);
                   
                }
            }
        }
        
        cross_allelic_contig_pairs.retain(|contig_pair| !allelic_contig_pairs.contains(contig_pair));
        
        for contig_pair in allelic_contig_pairs.into_iter() {
            self.allelic_counts += 1;
            match allele_strand_info.get(&contig_pair) {
                Some(allele_strand_record) => {
                    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            contig_pair.Contig1, contig_pair.Contig2, 
                            allele_strand_record.length1, allele_strand_record.length2,
                            allele_strand_record.matches, allele_strand_record.identity, 
                            0);
                },
                None => {
                    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            contig_pair.Contig1, contig_pair.Contig2, 
                            0, 0, 0, 0.98, 0);
                }
            }
            
        }

        for contig_pair in cross_allelic_contig_pairs.into_iter() {
            self.cross_allelic_counts += 1;
            writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    contig_pair.Contig1, contig_pair.Contig2, 
                    0, 0, 0, 0, 1);
        }
    }


    pub fn kprune(&mut self, whitehash: &HashSet<&String>, method: &str, writer: &mut Box<dyn Write> ) {
        let unique_min: HashMap<String, f64> = HashMap::new();
        let contact_data = self.contacts.to_data(&unique_min, &self.normalization_method);
        let mut allelic_contig_pairs = self.alleletable.get_allelic_contig_pairs(whitehash);
        if whitehash.len() > 0 {
            allelic_contig_pairs.retain(|x| whitehash.contains(x.Contig1) && whitehash.contains(x.Contig2))
        }

        let allelic_contig_db = self.alleletable.to_contig_db(method, whitehash);

        let mut allelic_contig_db2: HashMap<&String, Vec<AlleleRecord>> = HashMap::new();
        if whitehash.len() > 0{
            for (contig, record) in allelic_contig_db.iter() {
                let mut tmp_record: Vec<AlleleRecord> = Vec::new();
                
                for alleles in record.iter() {
                    let mut tmp_alleles = alleles.clone().clone();
                    tmp_alleles.data.retain(|x| whitehash.contains(x));
                    if tmp_alleles.data.len() > 1 {
                        tmp_record.push(tmp_alleles);
                    }
                }

                if tmp_record.len() > 0 {
                    if let Some(old_value) = allelic_contig_db2.insert(contig, tmp_record) {
                        
                    }
                }
            }

        }

        let mut allelic_contig_db:HashMap<&String, Vec<&AlleleRecord>> = HashMap::new();
        for (contig, record) in allelic_contig_db2.iter() {
            let mut tmp_record: Vec<&AlleleRecord> = Vec::new();
            for alleles in record.iter() {
                tmp_record.push(alleles);
            }
            allelic_contig_db.insert(contig, tmp_record);
        }
        
        self.allele_strand_table.get_allele_strand_records();
        let allele_strand_info = self.allele_strand_table.get_info();
        
        let mut contig_pairs = contact_data.keys().collect::<HashSet<&ContigPair2>>();
        
        if whitehash.len() != 0 {
            contig_pairs.retain(|contig_pair| whitehash.contains(contig_pair.Contig1) && whitehash.contains(contig_pair.Contig2));
        }
        
    
        
        // filtered the contig pairs that are not in the allelic contig pairs

        // allelic_contig_pairs.retain(|contig_pair| contig_pairs.contains(contig_pair));

        // allelic_contig_pairs.retain(|x| allele_strand_info.contains_key(x));
        contig_pairs.retain(|x| !allelic_contig_pairs.contains(x));
        contig_pairs.retain(|x| x.Contig1 != x.Contig2);

        let cross_allelic: Vec<&&ContigPair2> = contig_pairs.par_iter().filter_map(|contig_pair| {
      
            let alleles1 = allelic_contig_db.get(contig_pair.Contig1)?;
            let alleles2 = allelic_contig_db.get(contig_pair.Contig2)?;
            
            let mut is_weak = false; 

            'outer: for group1 in alleles1.iter() {
           
                let init_group1 = &group1.data;
                
                'inner: for group2 in alleles2.iter() {
                    let init_group2 = &group2.data;
                    let mut group1_length = init_group1.len();
                    let mut group2_length = init_group2.len();


                    if group1_length <= 1 && group2_length <= 1 {
                        // is_weak = true;
                        break 'outer;
                    }
                    let contig1 = match group1_length > group2_length {
                        true => contig_pair.Contig2, 
                        false => contig_pair.Contig1,
                    };
                    
                    let contig2 = match group1_length > group2_length {
                        true => contig_pair.Contig1, 
                        false => contig_pair.Contig2,
                    };

                    let group1 = match group1_length > group2_length {
                        true => &init_group2,
                        false => &init_group1,
                    };
                    let group2 = match group1_length > group2_length {
                        true => &init_group1,
                        false => &init_group2,
                    };

                    let idx1 = match group1.iter().position(|x| x == contig1) {
                        Some(idx) => idx,
                        None => continue 'inner,
                    };


                    let idx2 = match group2.iter().position(|x| x == contig2) {
                        Some(idx) => idx,
                        None => continue 'inner,
                    };
           
                    // if group1_length > group2_length {
                    //     // std::mem::swap(&mut group1, &mut group2);
                    //     // std::mem::swap(&mut idx1, &mut idx2);
                    //     std::mem::swap(&mut group1_length, &mut group2_length)
                    // }
                    let group1_length = group1.len();
                    let group2_length = group2.len();
    

                    let mut matrix = Matrix::new(group1_length, group2_length, OrderedFloat(0.0));
                    matrix.iter_mut().enumerate().for_each(|(index, element)| {
                        let i = index / group2_length;
                        let j = index % group2_length;
                        let mut contig_pair1 = ContigPair2::new(&group2[j], &group1[i]);
                        // let mut contig_pair2 = ContigPair2::new(&group2[j], &group1[i]);
                        contig_pair1.order();
                        // contig_pair2.order();
                        if let Some(value) = contact_data.get(&contig_pair1) {
                            *element = OrderedFloat(*value);
                        } 
                        // else if let Some(value) = contact_data.get(&contig_pair2) {
                        //     *element = OrderedFloat(*value);
                        // }
                    });

                    let assignments = maximum_bipartite_matching(matrix);
                    if assignments[idx1] != idx2 {
                        is_weak = true;
                    } else {
                        is_weak = false; 
                    }
                    if is_weak {
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

        
        for contig_pair in allelic_contig_pairs.into_iter() {
            self.allelic_counts += 1;
            match allele_strand_info.get(&contig_pair) {
                Some(allele_strand_record) => {
                    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            contig_pair.Contig1, contig_pair.Contig2, 
                            allele_strand_record.length1, allele_strand_record.length2,
                            allele_strand_record.matches, allele_strand_record.identity, 
                            0).unwrap();
                },
                None => {
                    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            contig_pair.Contig1, contig_pair.Contig2, 
                            0, 0, 0, 0.75, 0).unwrap();
                }
            }
            
        }

        for contig_pair in cross_allelic.into_iter() {
            self.cross_allelic_counts += 1;
            writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    contig_pair.Contig1, contig_pair.Contig2, 
                    0, 0, 0, 0, 1).unwrap();
        }

    }

}


