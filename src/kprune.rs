use anyhow::Result as anyResult;
use ordered_float::OrderedFloat;
use std::borrow::Cow;
use std::collections::{ HashMap, HashSet };
use std::error::Error;
use std::path::Path;
use std::io::{ Write, BufReader, BufRead };
use serde::{ Deserialize, Serialize};
use rayon::prelude::*;
use pathfinding::prelude::{kuhn_munkres, Matrix, Weights};

use crate::alleles::AlleleTable;
use crate::core::{ common_reader, common_writer };
use crate::core::{ BaseTable, ContigPair };
use crate::pixels::Pixels;
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

    pub fn write(&self, output: &String) {
        let mut writer = common_writer(output);
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(writer);

        for record in &self.records {
            wtr.serialize(record).expect("Could not write record")
        } 
        log::info!("Allelic and cross-allelic information written to {}", output);
    }

}



pub struct KPruner {
    pub alleletable: AlleleTable,
    pub contacts: HashMap<ContigPair, f64>,
    pub prunetable: PruneTable,
    pub contig_pairs: Vec<ContigPair>,
}

impl KPruner {
    pub fn new(alleletable: &String, pixels: &String, count_re: &String, prunetable: &String) -> KPruner {
        
        let mut count_re = CountRE::new(count_re);
        count_re.parse();
        let mut pixels = Pixels::new(pixels);
        pixels.parse();
        let contact_data = pixels.to_data(count_re.to_data(), true);
        let mut contig_pairs: Vec<ContigPair> = contact_data.keys().cloned().collect();

        let mut alleletable = AlleleTable::new(alleletable);
    
        KPruner {
            alleletable: alleletable,
            contacts: contact_data,
            prunetable: PruneTable::new(prunetable),
            contig_pairs: contig_pairs,
        }
    }

    pub fn allelic(&mut self) -> HashSet<ContigPair> {
        log::info!("Starting allelic identification");
        let allelic_contig_pairs = self.alleletable.get_allelic_contig_pairs(); 
        // let prune_allelic_contig_pairs = self.contig_pairs.iter()
        //                                                 .filter(
        //                                                     |x| allelic_contig_pairs.contains(x))
        //                                                 .cloned()
        //                                                 .collect::<Vec<ContigPair>>();
                                                    
        self.contig_pairs.retain(|x| !allelic_contig_pairs.contains(x));
        
        log::info!("Allelic contig pairs: {}", allelic_contig_pairs.len());
        allelic_contig_pairs

    }

    pub fn cross_allelic(&mut self) -> Vec<ContigPair> {
        let mut cross_allelic: Vec<ContigPair> = Vec::new();
        let allelic_contigs = self.alleletable.get_allelic_contigs();
        for contig_pair in &self.contig_pairs {
            let alleles1 = allelic_contigs.get(&contig_pair.Contig1).unwrap();
            let alleles2 = allelic_contigs.get(&contig_pair.Contig2).unwrap();
            

            let mut group1 = std::iter::once(&contig_pair.Contig1).chain(alleles1.iter()).collect::<Vec<&String>>();
            let mut group2 = std::iter::once(&contig_pair.Contig2).chain(alleles2.iter()).collect::<Vec<&String>>();

            if group1.len() > group2.len() {
                std::mem::swap(&mut group1, &mut group2);
            }
            

            let mut matrix = Matrix::new(group1.len(), group2.len(), OrderedFloat(0.0));
            for (i, ctg1) in group1.iter().enumerate() {
                for (j, ctg2) in group2.iter().enumerate() {
                    let mut tmp_contig_pair = ContigPair::new(ctg1.to_string(), ctg2.to_string());
                    tmp_contig_pair.order();
                    let score = self.contacts.get(&tmp_contig_pair).unwrap_or(&0.0);

                    matrix[(i, j)] = OrderedFloat(*score);
                }
            }

            let assignments =  maximum_bipartite_matching(matrix);
            if assignments[0] != 0 {
                cross_allelic.push(contig_pair.clone());
            }
            
        }
        log::info!("Cross allelic contig pairs: {}", cross_allelic.len());
        cross_allelic 
    }

    pub fn prune(&mut self) {
        let allelic = self.allelic();
        let mut cross_allelic = self.cross_allelic();

        for contig_pair in allelic.iter() {
            let prune_record = PruneRecord {
                contig1: contig_pair.Contig1.clone(),
                contig2: contig_pair.Contig2.clone(),
                mz1: 0,
                mz2: 0,
                mz_shared: 0,
                similarity: 0.0,
                allele_type: 0,
            };
            self.prunetable.records.push(prune_record);
        }

        for contig_pair in cross_allelic.iter() {
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