#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::Result as AnyResult;
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::io::{ Read, Write, BufReader, BufRead, BufWriter };
use serde::{ Serialize, Deserialize };
use rayon::prelude::*;

use crate::alleles::{ AlleleTable2, AlleleHeader };
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair, ContigPair2 };
use crate::count_re::CountRE;


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ContactRecord {
    pub chrom1: String,
    pub chrom2: String,
    pub count: f64,
}

impl ContactRecord {
    pub fn new() -> Self {
        Self {
            chrom1: String::new(),
            chrom2: String::new(),
            count: 0.0,
        }
    }

    pub fn is_some(&self) -> bool {
        match self.count {
            0.0 => false,
            _ => true,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Contacts {
    pub file: String,
    pub records: Vec<ContactRecord>,
}

impl BaseTable for Contacts {
    fn new(name: &String) -> Contacts {
        Contacts {
            file: name.clone(),
            records: Vec::new(),
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


impl Contacts {
    pub fn parse(&mut self) {

       let input = common_reader(&self.file);
         let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(false)
                .from_reader(input);
    
          for result in rdr.deserialize() {
                let record: ContactRecord = result.unwrap();
                self.records.push(record);
          }
       
    }

    pub fn from_clm(clm: &String) -> Self {
        let buffered = common_reader(clm);
        let reader = BufReader::new(buffered);
        let mut records: Vec<ContactRecord> = Vec::new();
        for (i, record) in reader.lines().enumerate() {
            if i % 4 != 0 {
                continue
            }
            let record = record.unwrap();
            let mut record_iter = record.split("\t");
            let contig_pair = record_iter.next().unwrap().to_string();
            let count = record_iter.next().unwrap().parse::<f64>().unwrap();
            let mut contig_pair_iter = contig_pair.split(" ");
            let chrom1 = contig_pair_iter.next().unwrap().to_string();
            let chrom2 = contig_pair_iter.next().unwrap().to_string();

            let mut contact_record = ContactRecord::new();
            contact_record.chrom1 = chrom1;
            contact_record.chrom1.pop();
            contact_record.chrom2 = chrom2;
            contact_record.chrom2.pop();
            contact_record.count = count;

            records.push(contact_record);
        }

        let mut contacts = Contacts {
            file: clm.clone(),
            records: records,
        };
        contacts.file = contacts.prefix() + ".contacts";
        contacts
    }

    pub fn to_data(&self, unique_min: &HashMap<String, f64>, 
                    normalization_method: &String) -> HashMap<ContigPair, f64> {//, re_count: HashMap<String, u32>, lengths: HashMap<String, u32>) -> HashMap<ContigPair, f64> {
        // let total_re_count = re_count.values().sum::<u32>();
        // let total_length = lengths.values().sum::<u32>();
        // let re_density = total_re_count as f64 / total_length as f64;
        // let longest_re = re_count.values().max().unwrap();
        // let longest_re_square = (longest_re * longest_re) as f64;
        let mut data: HashMap<ContigPair, f64> = self.records.par_iter(
            ).map(|record| {
                let contig_pair = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
                let count = record.count;

                (contig_pair, count)
            }).collect();
        
        
        // get contig1 == contig2 data
        let cis_data = data.par_iter().filter(|(contig_pair, _)| {
            contig_pair.Contig1 == contig_pair.Contig2
        }).map(|(contig_pair, count)| {
            (contig_pair.Contig1.clone(), *count)
        }).collect::<HashMap<String, f64>>();

        let normalization_method = normalization_method.as_str();
        data.par_iter_mut().for_each(|(contig_pair, count)| {
            let mut ratio = 0.0;
            let count1 = cis_data.get(&contig_pair.Contig1).unwrap_or(&0.0);
            let count2 = cis_data.get(&contig_pair.Contig2).unwrap_or(&0.0);

            if contig_pair.Contig1 == contig_pair.Contig2 {
                ratio = match normalization_method {
                    "none" => *count as f64,
                    "cis" => {
                        
                        match count1 * count2 {
                            0.0 => 0.0,
                            _ => *count / ((count1 * count2).sqrt())
                        }
                    },
                    _ => {
                        let m1 = unique_min.get(&contig_pair.Contig1).unwrap_or(&0.0);
                        -(m1 + 1.0).log2()
                    }
                };
            } else {
                
                let m1_log = match normalization_method   {
                    "none" => 0.0,
                    "cis" => 0.0,
                    _ => {
                        let m1 = unique_min.get(&contig_pair.Contig1).unwrap_or(&0.0);
                       
                        -(m1 + 1.0).log2()
                    }
                };
                let m2_log = match normalization_method {
                    "none" => 0.0,
                    "cis" => 0.0,
                    _ => {
                        let m2 = unique_min.get(&contig_pair.Contig2).unwrap_or(&0.0);
                        -(m2 + 1.0).log2()
                    }
                };

  
                ratio = match count1 * count2 {
                    0.0 => 0.0,
                    _ => {
                        if normalization_method == "none" {
                            *count  
                        } 
                        else if normalization_method == "cis" {
                            *count / ((count1 * count2).sqrt())
                        }
                        else if normalization_method == "cis_unique" {
                            match m1_log * m2_log {
                                0.0 => 0.0,
                                _ => *count / ((count1 * count2).sqrt()) * (m1_log * m2_log)
                            }
                        }
                        else {
                            *count
                        }
                    }
                    // _ => *count / ((count1) * (count2)).sqrt(),
                    // _ => *count,
                    //     // _ => *count / ((count1 / (m1_log.powf(2.0))) * (count2 / (m2_log.powf(2.0)))).sqrt(),
                    // _ => match m1_log * m2_log {
                    //     0.0 => 0.0,
                    //     _ => *count / ((count1 * count2).sqrt()) * (m1_log * m2_log),

                    // }
                
                };
            }
            
            // replace NaN with 0.0
            if ratio.is_nan() {
                ratio = 0.0;
            }

            if ratio < 0.0 {
                ratio = 0.0;
            }
            
            *count = ratio;
        
        });

        // let mut data: HashMap<ContigPair, f64> = self.records.par_iter(
        //     ).map(|record| {
               
        //         let contig_pair = ContigPair::new(record.chrom1.clone(), record.chrom2.clone());
        //         if !re_count.contains_key(&record.chrom1) || !re_count.contains_key(&record.chrom2) {
        //             (contig_pair, record.count as f64)

        //         } else {
        //             let contig1_length = lengths.get(&record.chrom1).unwrap();
        //             let contig2_length = lengths.get(&record.chrom2).unwrap();
        //             let re_count1 = re_count.get(&record.chrom1).unwrap_or(&0);
        //             let re_count2 = re_count.get(&record.chrom2).unwrap_or(&0);

        //             // let r1 = *re_count1 as f64 / ((*contig1_length as f64) * re_density);
        //             // let r2 = *re_count2 as f64 / ((*contig2_length as f64) * re_density);
        //             let count = record.count as f64;
        //             // let ratio = count * 10000 / ((*contig1_length as f64 / 10000.0) * (*contig2_length as f64 / 10000.0));
        //             // best in 20231118
        //             let ratio = match re_count1 * re_count2 {
        //                 0 => 0.0,
        //                 // _ => count / ((contig1_length * contig2_length) as f64).log(10.0) 
        //                 _ =>  count / (re_count1 * re_count2) as f64,
        //             };

        //             (contig_pair, ratio)
        //         }
        //     }).collect();
        
        // min-max normalization
        // let max = data.values().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        // let min = data.values().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
        // let range = max - min;
        
        // let new_data: HashMap<ContigPair, f64> = data.par_iter().map(|(contig_pair, value)| {
        //     let new_value = (*value - min) / range;
        //     (contig_pair.clone(), new_value)
        // }).collect();

        // new_data  
        data 

    }


    pub fn write(&self, output: &String) {
        let mut wtr = common_writer(output);
        for record in &self.records {
            wtr.write_all(format!("{}\t{}\t{}\n", record.chrom1, record.chrom2, record.count).as_bytes()).unwrap();
        }
    }
}


#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ContactMatrix {
    pub file: String,
    pub records: Vec<Vec<f64>>,
}

impl BaseTable for ContactMatrix {
    fn new(name: &String) -> ContactMatrix {
        ContactMatrix {
            file: name.clone(),
            records: Vec::new(),
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





#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Contacts2 {
    pub file: String,
    pub records: Vec<ContactRecord>,
}

impl BaseTable for Contacts2 {
    fn new(name: &String) -> Contacts2 {
        Contacts2 {
            file: name.clone(),
            records: Vec::new(),
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


impl Contacts2 {
    pub fn parse(&mut self) {

       let input = common_reader(&self.file);
         let mut rdr = csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .comment(Some(b'#'))
                .has_headers(false)
                .from_reader(input);
    
          for result in rdr.deserialize() {
                let record: ContactRecord = result.unwrap();
                self.records.push(record);
          }
       
    }

    pub fn from_clm(clm: &String) -> Self {
        let buffered = common_reader(clm);
        let reader = BufReader::new(buffered);
        let mut records: Vec<ContactRecord> = Vec::new();
        for (i, record) in reader.lines().enumerate() {
            if i % 4 != 0 {
                continue
            }
            let record = record.unwrap();
            let mut record_iter = record.split("\t");
            let contig_pair = record_iter.next().unwrap().to_string();
            let count = record_iter.next().unwrap().parse::<f64>().unwrap();
            let mut contig_pair_iter = contig_pair.split(" ");
            let chrom1 = contig_pair_iter.next().unwrap().to_string();
            let chrom2 = contig_pair_iter.next().unwrap().to_string();

            let mut contact_record = ContactRecord::new();
            contact_record.chrom1 = chrom1;
            contact_record.chrom1.pop();
            contact_record.chrom2 = chrom2;
            contact_record.chrom2.pop();
            contact_record.count = count;

            records.push(contact_record);
        }

        let mut contacts = Contacts2 {
            file: clm.clone(),
            records: records,
        };
        contacts.file = contacts.prefix() + ".contacts";
        contacts
    }

    pub fn to_data(&self, unique_min: &HashMap<String, f64>, 
                    normalization_method: &String,
                    re_count: &Option<CountRE>
                ) -> HashMap<ContigPair2, f64> {//, re_count: HashMap<String, u32>, lengths: HashMap<String, u32>) -> HashMap<ContigPair, f64> {

        let _re_count = if let Some(v) = re_count {
          
            let mut re_count = v.clone();
            re_count.parse();
            let _re_count = re_count.to_data();
         
            _re_count.into_iter()
                .map(|(key, value)| (key, value as f64))
                .collect::<HashMap<_, _>>()
        } else {
           HashMap::new()
        };


        let re_count: HashMap<&String, f64> = _re_count.iter()
            .map(|(key, value)| (key, *value))
            .collect();
        // let total_re_count = re_count.values().sum::<u32>();
        // let total_length = lengths.values().sum::<u32>();
        // let re_density = total_re_count as f64 / total_length as f64;
        // let longest_re = re_count.values().max().unwrap();
        // let longest_re_square = (longest_re * longest_re) as f64;
        
        let mut data: HashMap<ContigPair2, f64> = self.records.par_iter(
            ).map(|record| {
                let mut contig_pair = match record.chrom1 > record.chrom2 {
                    true => ContigPair2::new(&record.chrom2, &record.chrom1),
                    false => ContigPair2::new(&record.chrom1, &record.chrom2)
                };
                contig_pair.order();
                let count = record.count;

                (contig_pair, count)
            }).collect();
    
    
        let normalization_method = match (re_count.len() == 0) {
            true => {
                log::info!("No RE counts provided, switching normalization method to '{}'", normalization_method);
                normalization_method.clone()

            },
            false => {
                log::info!("RE counts provided, using normalization method: {}", "re");
                "re".to_string()
            }

        };
        log::info!("Using normalization method: {}", normalization_method);
        

        // get contig1 == contig2 data
        let cis_data = data.par_iter().filter(|(contig_pair, _)| {
            *contig_pair.Contig1 == *contig_pair.Contig2
        }).map(|(contig_pair, count)| {
            (contig_pair.Contig1, *count)
        }).collect::<HashMap<&String, f64>>();

        let mut total_contacts: HashMap<&String, f64> = HashMap::new();
        if normalization_method == "vc" || normalization_method == "tweight" || normalization_method == "hybrid" {
            for (pair, count) in data.iter() {
                if *pair.Contig1 != *pair.Contig2 {
                    *total_contacts.entry(pair.Contig1).or_insert(0.0) += count;
                    *total_contacts.entry(pair.Contig2).or_insert(0.0) += count;
                }
            }
        }

        let mut t_row_sums: HashMap<&String, f64> = HashMap::new();
        if normalization_method == "tweight" {
            // Create a temporary map for T[i,j] values to calculate row sums
            let mut t_matrix_entries: HashMap<ContigPair2, f64> = HashMap::new();
            for (pair, count) in data.iter() {
                if *pair.Contig1 != *pair.Contig2 {
                    let total1 = total_contacts.get(pair.Contig1).unwrap_or(&1.0);
                    let total2 = total_contacts.get(pair.Contig2).unwrap_or(&1.0);
                    if *total1 > 0.0 && *total2 > 0.0 {
                        let t_value = count / (total1 * total2);
                        t_matrix_entries.insert(pair.clone(), t_value);
                    }
                }
            }
            // Calculate T_row_sum for each contig
            for (pair, t_value) in t_matrix_entries.iter() {
                *t_row_sums.entry(pair.Contig1).or_insert(0.0) += t_value;
                *t_row_sums.entry(pair.Contig2).or_insert(0.0) += t_value;
            }
        }

        let mut scaled_cis: HashMap<&String, f64> = HashMap::new();
        let mut scaled_total: HashMap<&String, f64> = HashMap::new();

        if normalization_method == "hybrid" {
            // 1. Calculate sum for scaling
            let cis_sum: f64 = cis_data.values().sum();
            let total_sum: f64 = total_contacts.values().sum();

            // 2. Create scaled bias vectors
            if cis_sum > 0.0 {
                for (contig, count) in cis_data.iter() {
                    scaled_cis.insert(contig, count / cis_sum);
                }
            }
            if total_sum > 0.0 {
                for (contig, count) in total_contacts.iter() {
                    scaled_total.insert(contig, count / total_sum);
                }
            }
        }

        let normalization_method = normalization_method.as_str();

        data.par_iter_mut().for_each(|(contig_pair, count)| {
            let mut ratio = 0.0;
            let count1 = cis_data.get(&contig_pair.Contig1).unwrap_or(&0.0);
            let count2 = cis_data.get(&contig_pair.Contig2).unwrap_or(&0.0);
            
            if *contig_pair.Contig1 == *contig_pair.Contig2 {
                ratio = match normalization_method {
                    "none" => *count as f64,
                    "cis" => {
                        
                        match count1 * count2 {
                            0.0 => 0.0,
                            _ => *count / ((count1 * count2).sqrt())
                        }
                    },
                    "re" => {
                        let re1 = re_count.get(&contig_pair.Contig1).unwrap_or(&1.0);
                        let re2 = re_count.get(&contig_pair.Contig2).unwrap_or(&1.0);
                        match re1 * re2 {
                            0.0 => 0.0,
                            _ => *count / ((re1 * re2).sqrt())
                        }
                    }
                    _ => { 0.0
                        // let m1 = unique_min.get(&contig_pair.Contig1.to_string()).unwrap_or(&0.0);
                        // -(m1 + 1.0).log2()
                    }
                };
            } else {
                
                let m1_log = match normalization_method   {
                    // "none" => 0.0,
                    // "cis" => 0.0,
                    _ => { 0.0
                        // let m1 = unique_min.get(&contig_pair.Contig1.to_string()).unwrap_or(&0.0);
                       
                        // -(m1 + 1.0).log2()
                    }
                };
                let m2_log = match normalization_method {
                    // "none" => 0.0,
                    // "cis" => 0.0,
                    _ => { 0.0
                        // let m2 = unique_min.get(&contig_pair.Contig2.to_string()).unwrap_or(&0.0);
                        // -(m2 + 1.0).log2()
                    }
                };

                
                ratio = match count1 * count2 {
                    -1.0 => 0.0,
                    _ => {
                        if normalization_method == "none" {
                            *count  
                        } 
                        else if normalization_method == "cis" {
                            *count / (((count1 + 1.0) * (count2 + 1.0)).sqrt())
                        }
                        else if normalization_method == "re" {
                            let re1 = re_count.get(&contig_pair.Contig1).unwrap_or(&1.0);
                            
                            let re2 = re_count.get(&contig_pair.Contig2).unwrap_or(&1.0);
                            match re1 * re2 {
                                0.0 => 0.0,
                                _ => *count / ((re1 * re2).sqrt())
                            }
                        }
                        else if normalization_method == "vc" {
                            let total1 = total_contacts.get(&contig_pair.Contig1).unwrap_or(&0.0); // Use 1.0 to avoid division by zero
                            let total2 = total_contacts.get(&contig_pair.Contig2).unwrap_or(&0.0);
                            
                            *count / ((total1 + 1.0) * (total2 + 1.0))
                            
                        }
                        else if normalization_method == "hybrid" {
                            let w_cis = 0.5;
                            let w_total = 0.5;
                            let pseudo_count = 1e-9; // Use a very small pseudo_count for scaled data
    
                            // Get the SCALED bias values
                            let s_cis1 = scaled_cis.get(&contig_pair.Contig1).unwrap_or(&0.0);
                            let s_cis2 = scaled_cis.get(&contig_pair.Contig2).unwrap_or(&0.0);
    
                            let s_total1 = scaled_total.get(&contig_pair.Contig1).unwrap_or(&0.0);
                            let s_total2 = scaled_total.get(&contig_pair.Contig2).unwrap_or(&0.0);
    
                            // Calculate the hybrid bias score using SCALED values
                            let bias1 = w_cis * s_cis1 + w_total * s_total1 + pseudo_count;
                            let bias2 = w_cis * s_cis2 + w_total * s_total2 + pseudo_count;
                            
                            // The denominator is now a product of small, scaled numbers.
                            // The raw count needs to be divided by this.
                            if bias1 > 0.0 && bias2 > 0.0 {
                                *count / (bias1 * bias2)
                            } else {
                                0.0
                            }
                        }
                        else if normalization_method == "tweight" {
                            let total1 = total_contacts.get(&contig_pair.Contig1).unwrap_or(&1.0);
                            let total2 = total_contacts.get(&contig_pair.Contig2).unwrap_or(&1.0);
                            let t_sum1 = t_row_sums.get(&contig_pair.Contig1).unwrap_or(&1.0);
                            let t_sum2 = t_row_sums.get(&contig_pair.Contig2).unwrap_or(&1.0);
        
                            let bias1 = total1 * t_sum1;
                            let bias2 = total2 * t_sum2;
        
                            if bias1 > 0.0 && bias2 > 0.0 {
                                *count / (bias1 * bias2)
                            } else {
                                0.0
                            }
                        }
                        else if normalization_method == "cis_unique" {
                            match m1_log * m2_log {
                                0.0 => 0.0,
                                _ => *count / ((count1 * count2).sqrt()) * (m1_log * m2_log)
                            }
                        }
                        else {
                            *count
                        }
                    }
                    // _ => *count / ((count1) * (count2)).sqrt(),
                    // _ => *count,
                    //     // _ => *count / ((count1 / (m1_log.powf(2.0))) * (count2 / (m2_log.powf(2.0)))).sqrt(),
                    // _ => match m1_log * m2_log {
                    //     0.0 => 0.0,
                    //     _ => *count / ((count1 * count2).sqrt()) * (m1_log * m2_log),

                    // }


                
                };
            }
            
            if ratio.is_nan() {
                ratio = 0.0;
            }

            if ratio < 0.0 {
                ratio = 0.0;
            }
            
            *count = ratio;
        
        });


        data 

    }

    pub fn contigs(&self) -> HashSet<String> {
        let mut contigs: HashSet<String> = HashSet::new();
        for record in &self.records {
            if record.chrom1 == record.chrom2 {
                continue
            }
            contigs.insert(record.chrom1.clone());
            contigs.insert(record.chrom2.clone());
        }

        contigs
    }

    pub fn contig1_to_contig2(&self) -> HashMap<String, HashSet<String>> {
        let mut contig1_to_contig2: HashMap<String, HashSet<String>> = HashMap::new();
        for record in &self.records {
            if record.chrom1 == record.chrom2 {
                continue;
            }
            contig1_to_contig2.entry(record.chrom1.clone()).or_insert_with(HashSet::new).insert(record.chrom2.clone());
            contig1_to_contig2.entry(record.chrom2.clone()).or_insert_with(HashSet::new).insert(record.chrom1.clone());
        }

        contig1_to_contig2
    }

    pub fn write(&self, output: &String) {
        let mut wtr = common_writer(output);
        for record in &self.records {
            wtr.write_all(format!("{}\t{}\t{}\n", record.chrom1, record.chrom2, record.count).as_bytes()).unwrap();
        }
    }
}



