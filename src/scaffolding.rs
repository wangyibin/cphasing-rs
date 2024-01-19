use rand::{thread_rng, Rng};
use rand::seq::SliceRandom;

use std::collections::{HashMap, HashSet};
use std::collections::VecDeque;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use crate::core::ContigPair;

#[derive(Debug, Clone)]
pub struct SplitContigUnit {
    pub contig: String,
    pub cis1: f64,
    pub cis2: f64,
    pub trans: f64,
    pub orientation: u8,
    pub order: usize,
}

impl SplitContigUnit {
    pub fn new() -> Self {
        Self {
            contig: String::new(),
            cis1: 0.0,
            cis2: 0.0,
            trans: 0.0,
            orientation: 0,
            order: 0,
        }
    }

}

#[derive(Debug, Clone)]
pub struct SplitTour {
    pub contigs: Vec<SplitContigUnit>,
    pub cost: f64,
    pub split_contact_hash: HashMap<ContigPair, f64>,
}

impl SplitTour {
    pub fn new() -> Self {

        Self {
            contigs: Vec::new(),
            cost: 0.0,
            split_contact_hash: HashMap::new(),
        }
    }

    pub fn from_contacts(&mut self, contact_hash: &HashMap<ContigPair, f64>, contig_list: &Vec<String>) {
        let mut split_contig_hash: HashMap<String, SplitContigUnit> = HashMap::new();
        let mut i = 0;
     
        for (contig_pair, contact) in contact_hash.iter(){
            
            let part1: Vec<&str>  = contig_pair.Contig1.split("_").collect();
            let part2: Vec<&str>  = contig_pair.Contig2.split("_").collect();
            let contig1 = part1[0];
            let contig2 = part2[0];
            let suffix1 = part1[1];
            let suffix2 = part2[1];
            
            let suffix1 = suffix1.parse::<u32>().unwrap();
            let suffix2 = suffix2.parse::<u32>().unwrap();

            
            // contig suffix split by "_" if suffix equal 0 assign contact to cis1 if
            if contig1 != contig2 {
                continue
            }
            
            // init order value by contig_list
            let order = contig_list.iter().position(|x| *x == contig1).unwrap();

            let mut scu = if let Some(split_contig) = split_contig_hash.get_mut(contig1) {
                split_contig
            } else {
                let split_contig = split_contig_hash.insert(contig1.to_string(), SplitContigUnit::new());
                split_contig_hash.get_mut(contig1).unwrap()
            };
            
            match (suffix1, suffix2) {
                (0, 0) => {
                    scu.contig = contig1.to_string();
                    scu.cis1 = *contact;
                    scu.orientation = 0;
                    scu.order = order;
                    
                }, 
                (1, 1) => {
                    scu.contig = contig2.to_string();
                    scu.cis2 = *contact;
                    scu.orientation = 0;
                    scu.order = order;
                },
                (0, 1) => {
                    scu.contig = contig1.to_string();
                    scu.trans = *contact;
                    scu.orientation = 0;
                    scu.orientation = 0;
                    scu.order = order;

                }
                _ => panic!("Only support for contig suffix 0 or 1"),
            }

            
        }

        // sort contigs by order and push to self.contigs 
        let mut contigs: Vec<SplitContigUnit> = Vec::new();
        for (_, scu) in split_contig_hash.iter() {
            contigs.push(scu.clone());
        }
        contigs.sort_by(|a, b| a.order.cmp(&b.order));
        self.contigs = contigs;
    }

    pub fn get_split_contact_hash(&mut self,
        contact_hash: &HashMap<ContigPair, f64>) {

        let mut split_contact_hash: HashMap<ContigPair, f64> = HashMap::new();
        
        // range contigs by 2 step
       


        // combination of contigs
        for i in 0..(self.contigs.len() - 1) {
            let j = i + 1;
            let split_contig1 = &self.contigs[i];
            let split_contig2 = &self.contigs[j];
           
            let contig1 = match split_contig1.orientation {
                0 => format!("{}_0", split_contig1.contig),
                1 => format!("{}_1", split_contig1.contig),
                _ => panic!("Only support for contig suffix 0 or 1"),
            };
            let contig2 = match split_contig2.orientation {
                0 => format!("{}_0", split_contig2.contig),
                1 => format!("{}_1", split_contig2.contig),
                _ => panic!("Only support for contig suffix 0 or 1"),
            };
                

            let mut split_contig_pair = ContigPair::new(contig1, contig2);
            split_contig_pair.order();
            let contig_pair = ContigPair::new(split_contig1.contig.clone(), 
                                                split_contig2.contig.clone());
            
            
            let contact = contact_hash.get(&split_contig_pair).unwrap_or(&0.0);
            match contact {
                0.0 => continue,
                _ => (),
            }

            split_contact_hash.insert(contig_pair.clone(), *contact);
        
        }

        self.split_contact_hash = split_contact_hash;

    }

    pub fn get_cost(&mut self) -> f64 {
        let mut cost = 0.0;
        for (_, contact) in self.split_contact_hash.iter() {
            cost += contact;
        } 

        self.cost = cost;
        cost
    }
    

}


#[derive(Debug, Clone)]
pub struct SplitTourPopulation {
    pub tours: Vec<SplitTour>,
    pub contact_hash: HashMap<ContigPair, f64>,
    pub population_size: usize,
    pub mutation_rate: f64,
    pub crossover_rate: f64,
    pub max_generation: usize,
    pub current_generation: usize,
    pub best_tour: SplitTour,
    pub best_cost: f64,
}


impl SplitTourPopulation {
    pub fn new() -> Self {
        Self {
            tours: Vec::new(),
            contact_hash: HashMap::new(),
            population_size: 100,
            mutation_rate: 0.2,
            crossover_rate: 0.8,
            max_generation: 5000,
            current_generation: 0,
            best_tour: SplitTour::new(),
            best_cost: 0.0,
        }
    
    }

    pub fn init_population(&mut self, split_tour: SplitTour) {
        let mut rng = thread_rng();
        for i in 0..self.population_size {
            let mut tmp_split_tour = split_tour.clone();
            
            tmp_split_tour.contigs.shuffle(&mut rng);
            self.tours.push(tmp_split_tour);
        }
    }

    pub fn best(&mut self) {
        let mut best_cost = 0.0;
        let mut best_tour = SplitTour::new();
        for tour in self.tours.iter_mut() {
            tour.get_split_contact_hash(&self.contact_hash);
            let cost = tour.get_cost();
            if cost > best_cost {
                best_cost = cost;
                best_tour = tour.clone();
            }
        }

        self.best_cost = best_cost;
        self.best_tour = best_tour;
    }

    pub fn selection(&mut self) {
        let mut rng = thread_rng();
        let mut selected_tours: Vec<SplitTour> = Vec::new();
        for _ in 0..self.population_size {
            let tour1 = self.tours.choose(&mut rng).unwrap();
            let tour2 = self.tours.choose(&mut rng).unwrap();
            if tour1.cost > tour2.cost {
                selected_tours.push(tour1.clone());
            } else {
                selected_tours.push(tour2.clone());
            }
        }

        self.tours = selected_tours;
    }

}   


pub fn crossover(mut split_tour1: SplitTour, mut split_tour2: SplitTour, 
                    contact_hash: &HashMap<ContigPair, f64>) -> (SplitTour, SplitTour){
    let mut rng = thread_rng();
    // pick random index
    let mut index = rng.gen_range(0..split_tour1.contigs.len());
    // 
    split_tour1.get_split_contact_hash(contact_hash);
    split_tour2.get_split_contact_hash(contact_hash);

    (split_tour1, split_tour2)
}



fn permute() {

}

fn splice() {

}

fn insertion() {

}

fn inversion() {

}

fn mutation() {
    
}

fn selection() {
    
}

fn fitness() {

}

fn genetic_algorithm() {
    
}