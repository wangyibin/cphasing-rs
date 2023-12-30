use anyhow::Result as AnyResult;
use std::borrow::Cow;
use std::collections::HashMap;
use std::error::Error;
use std::path::Path;
use std::io::{ Read, Write, BufReader, BufRead, BufWriter };
use serde::{ Serialize, Deserialize };
use rayon::prelude::*;

use crate::alleles::{ AlleleTable, AlleleHeader };
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };


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
            let record: Vec<&str> = record.split("\t").collect();
            let contig_pair = record[0].to_string();
            let count = record[1].parse::<f64>().unwrap();
            let mut contact_record = ContactRecord::new();
            contact_record.chrom1 = contig_pair.split(" ").nth(0).unwrap().to_string();
            contact_record.chrom1.pop();
            contact_record.chrom2 = contig_pair.split(" ").nth(1).unwrap().to_string();
            contact_record.chrom2.pop();
            contact_record.count = count as f64; 

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
            if contig_pair.Contig1 == contig_pair.Contig2 {
                ratio = *count as f64; 
            } else {
                
                let count1 = cis_data.get(&contig_pair.Contig1).unwrap_or(&0.0);
                let count2 = cis_data.get(&contig_pair.Contig2).unwrap_or(&0.0);
                
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

                // if m1_log == 0.0 || m2_log == 0.0{
                //     println!("{} {}: {} {} | {}: {} {}", count, contig_pair.Contig1, 
                //                 m1_log, count1, contig_pair.Contig2, m2_log, count2);
                // }
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




