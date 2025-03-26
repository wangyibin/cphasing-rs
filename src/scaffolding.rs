#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use rand::{thread_rng, Rng};
use rand::seq::SliceRandom;

use std::collections::{HashMap, HashSet};
use std::collections::VecDeque;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use pathfinding::undirected::kruskal::kruskal;

use crate::core::{common_writer, common_reader};

use crate::core::{ 
    ContigPair,
    ContigPair2,
};


#[derive(Debug, Clone)]
pub struct ContigUnit {
    contig: String,
    index: u32,
    length: u32,
}


#[derive(Debug, Clone, Eq, Hash, PartialEq)]
pub struct SplitContig {
    contig: String,
    index: u8,
}

impl SplitContig {
    pub fn new(contig: &str) -> Self {
        let record = contig.split("_").collect::<Vec<&str>>();
        let contig = record[0].to_string();
        let index = record[1].parse::<u8>().unwrap();

        Self {contig: contig, index: index}

    }   
}


// #[derive(Debug, Clone, Eq, Hash, PartialEq)]
// pub struct SplitContigPair<'a> {
//     pub splitcontig1: SplitContig,
//     pub splitcontig2: SplitContig,
// }

// impl SplitContigPair<'_> {
//     pub fn new<'a>(contig1: &'a SplitContig, contig2: &'a SplitContig) -> Self {
//         Self { splitcontig1: contig1, splitcontig2: contig2 }
//     }
// }

#[derive(Debug, Clone)]
pub struct SplitContacts {
    contigs: HashMap<SplitContig, usize>,
    data: HashMap<(u32, u32), f64>,
}



impl SplitContacts {
    pub fn new() -> Self {
        SplitContacts { contigs: HashMap::new(), data: HashMap::new()}
    }

    pub fn from_path(&mut self, path: &String) {
        let reader = common_reader(path); 

        let mut idx = 0;
        for line in reader.lines().filter_map(|line| line.ok()) {
            let record = line.trim_end_matches('\n').split('\t').collect::<Vec<&str>>();

            let contig1 = record[0];
            let contig2 = record[1];
            let contact = record[2].parse::<f64>().unwrap();
            

            let mut idx1 = *self.contigs.entry(SplitContig::new(contig1)).or_insert_with(|| {
                let current_idx = idx;
                idx += 1;
                current_idx
            });

            let mut idx2 = *self.contigs.entry(SplitContig::new(contig2)).or_insert_with(|| {
                let current_idx = idx;
                idx += 1;
                current_idx
            });

        
            self.data.insert((idx1.try_into().unwrap(), idx2.try_into().unwrap()), contact);

            println!("{:?}", self.contigs);
        }
    }
}



#[derive(Debug, Clone)]
pub struct GeneticAlgorithm {

}

impl GeneticAlgorithm {
    pub fn new() -> Self {
        GeneticAlgorithm {}
    }
}

