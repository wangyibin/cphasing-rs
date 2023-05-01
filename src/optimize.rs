use lazy_static::lazy_static;
use itertools::Itertools;
use rand::Rng;
use rand::prelude::SliceRandom;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::core::ChromSize;


use crate::core::common_reader;
use crate::core::BaseTable;

lazy_static! {
    #[derive(Debug)]
    pub static ref ORIENTATION_TO_IDX: HashMap<&'static str, u8> = {
        let mut m = HashMap::new();
        m.insert("++", 0);
        m.insert("+-", 1);
        m.insert("-+", 2);
        m.insert("--", 3);
        m
    };
}



lazy_static! {
    pub static ref ORIENTATION_REVERSE_DB: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::new();
        m.insert("++", "--");
        m.insert("+-", "-+");
        m.insert("-+", "+-");
        m.insert("--", "++");
        m
    };
}


#[derive(Debug, Clone)]
pub struct ContigUnit {
    pub contig: String,
    pub orderidx: usize,
    pub orientation: u8,
    pub scores: HashMap<String, Vec<f64>>,
}

impl ContigUnit {
    pub fn new(contig: String, orderidx: usize, orientation: u8, scores: HashMap<String, Vec<f64>>) -> ContigUnit {
        ContigUnit { contig, orderidx, orientation, scores }
    }

    pub fn get_score(&self, contig: &ContigUnit) -> f64 {

        let idx = (self.orientation * 2).pow(1) + (contig.orientation * 2).pow(0);

        let score = self.scores.get(&contig.contig)
            .map(|v| v[idx as usize])
            .unwrap_or(0.0);
        // println!("{} {} {} {} {} {}", self.contig, contig.contig, self.orientation, contig.orientation, idx, score);
        score
    }

    pub fn get_orientation(&self) -> u8 {
        self.orientation.clone()
    }

    pub fn get_orderidx(&self) -> usize {
        self.orderidx
    }

    pub fn get_contig(&self) -> String {
        self.contig.clone()
    }

    pub fn insert_score(&mut self, contig: String, score: f64, orientation: &str) {
        let idx = ORIENTATION_TO_IDX.get(orientation).unwrap().clone();
       
        match self.scores.get_mut(&contig) {
            Some(v) => {
                v[idx as usize] = score;
            },
            None => {
                let mut v = vec![0.0, 0.0, 0.0, 0.0];
                v[idx as usize] = score;
                self.scores.insert(contig, v);
            }
        }
     
    }

    pub fn shuffle(&mut self) {
        let mut rng = rand::thread_rng();
        self.orderidx = rng.gen_range(0..self.orderidx);
    }
}   

#[derive(Debug, Clone)]
pub struct ContigOrder {
    pub contig_units: Vec<ContigUnit>,
}


impl ContigOrder {
    pub fn new(contig_units: Vec<ContigUnit>) -> ContigOrder {
        ContigOrder { contig_units }
    }

    pub fn cost(&self) -> f64 {
        let mut cost = 0.0;

        for (i, contig) in self.contig_units.iter().sorted_by_key(|c| c.orderidx).enumerate() {
            for j in (i + 1)..self.contig_units.len() {
                let score = contig.get_score(&self.contig_units[j]);
                let distance = (j - i) as f64;
                cost += score / distance;

                // println!("{} {} {} {} {}", contig.contig, self.contig_units[j].contig, score, distance, score/distance);
            }
            
        }

        cost
    }
    pub fn shuffle(&mut self) {
        let mut rng = rand::thread_rng();
        self.contig_units.shuffle(&mut rng);
        for (i, contig) in self.contig_units.iter_mut().enumerate() {
            contig.orderidx = i;
        }
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        self.contig_units[i].orderidx = j;
        self.contig_units[j].orderidx = i;
    }

    // pub fn crossover(&mut self, i: usize, j: usize) {
    //     let mut contig_units = self.contig_units.clone();
    //     let mut rng = rand::thread_rng();
    //     let mut i = i;
    //     let mut j = j;

    //     if i == j {
    //         return;
    //     }
        
    //     if i > j {
    //         (i, j) = (j, i);
    //     }

    //     let mut k = rng.gen_range(0..i);
    //     let mut l = rng.gen_range((j + 1)..self.contig_units.len());
        
    //     let mut k_iter = contig_units.splice(k..=i, std::iter::empty());
    //     let mut l_iter = contig_units.splice(l..=j, std::iter::empty());

    //     contig_units.splice(k..=k, l_iter);
    //     contig_units.splice(l..=l, k_iter);

    //     for (i, contig) in contig_units.iter_mut().enumerate() {
    //         contig.orderidx = i;
    //     }

    //     self.contig_units = contig_units;
    // }

    pub fn rotate(&mut self, i: usize) {
        self.contig_units[i].orientation = self.contig_units[i].orientation ^ 1;
    }

    pub fn reverse(&mut self, i: usize, j: usize) {
        for k in i..j + 1 {
            self.contig_units[k].orientation = self.contig_units[k].orientation ^ 1;
            self.contig_units[k].orderidx = j - self.contig_units[k].orderidx + i;
        }
    }


    pub fn contigs(&self) -> Vec<String> {
        self.contig_units.iter().map(|c| c.contig.clone()).collect()
    }
    
    
}

// struct for ContigScoreTable which is a file that contain three columns
#[derive(Debug, Clone)]
pub struct ContigScoreTable {
    file: String,
}

impl BaseTable for ContigScoreTable {
    fn new(name: &String) -> ContigScoreTable {
        ContigScoreTable { file: name.clone() }
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

impl ContigScoreTable {
    pub fn read(&self) -> ContigOrder {
        let mut contig_score_data = HashMap::<String, ContigUnit>::new();

        let reader = common_reader(&self.file);
        let mut reader = BufReader::new(reader);
        let mut line = String::new();

        let mut int_order_idx = 0;
        

        while reader.read_line(&mut line).unwrap() > 0 {
            let mut items = line.split_whitespace();
            let contig1 = items.next().unwrap().to_string();
            let contig2 = items.next().unwrap().to_string();
            let orientation1 = contig1.chars().last().unwrap();
            let orientation2 = contig2.chars().last().unwrap();
            let contig1 = contig1.chars().take(contig1.len()-1).collect::<String>();
            let contig2 = contig2.chars().take(contig2.len()-1).collect::<String>();
            let score = items.next().unwrap().parse::<f64>().unwrap();
            
            let orientation = format!("{}{}", orientation1, orientation2).clone();
            let orientation_reverse = ORIENTATION_REVERSE_DB.get(&orientation.as_str()).unwrap();
            let orientation_reverse = orientation_reverse.to_string();
            
            match contig_score_data.get_mut(&contig1) {
                Some(v) => {

                    v.insert_score(contig2.clone(), score, orientation.as_str());
                  
                },
                None => {
                    let mut scores = HashMap::<String, Vec<f64>>::new();
                    scores.insert(contig2.clone(), vec![0.0, 0.0, 0.0, 0.0]);
                    scores.get_mut(&contig2)
                            .unwrap()[ORIENTATION_TO_IDX.get(
                                                orientation.as_str())
                                                            .unwrap()
                                                            .clone() as usize] = score;
                    let contig_unit = ContigUnit::new(contig1.clone(), int_order_idx, 0 as u8, scores);
                    contig_score_data.insert(contig1.clone(), contig_unit);
                    int_order_idx += 1;
                }
            }

            match contig_score_data.get_mut(&contig2) {
                Some(v) => {
                    v.insert_score(contig1.clone(), score, orientation_reverse.as_str());
                },
                None => {
                    let mut scores = HashMap::<String, Vec<f64>>::new();
                    scores.insert(contig1.clone(), vec![0.0, 0.0, 0.0, 0.0]);
                    scores.get_mut(&contig1)
                            .unwrap()[ORIENTATION_TO_IDX.get(
                                                orientation_reverse.as_str())
                                                        .unwrap().clone() as usize] = score;
                    let contig_unit = ContigUnit::new(contig2.clone(), int_order_idx, 0 as u8, scores);
                    contig_score_data.insert(contig2.clone(), contig_unit);
                    int_order_idx += 1;
                }
            }
            
            line.clear();
            
        }
        let values: Vec<ContigUnit> = contig_score_data.values().cloned().collect();

        // let values: Vec<ContigUnit> = contig_score_data.into_iter()
        //                                 .map(|(k, v)| (k,v))
        //                                 .sorted_by_key(|(_, v)| v.orderidx)
        //                                 .map(|(_, v)| v)
        //                                 .collect();

        ContigOrder::new(values)
    }

}

pub struct SimulatedAnnealing {
    pub contig_order: ContigOrder,
    pub temperature: f64,
    pub cooling_rate: f64,
    pub mimimum_temperature: f64,
    pub max_iteration: usize,
}

impl SimulatedAnnealing {
    pub fn new(contig_order: ContigOrder, 
            temperature: f64, 
            cooling_rate: f64, 
            mimimum_temperature: f64, 
            max_iteration: usize) -> SimulatedAnnealing {
        SimulatedAnnealing {
            contig_order: contig_order,
            temperature: temperature,
            cooling_rate: cooling_rate,
            mimimum_temperature: mimimum_temperature,
            max_iteration: max_iteration,
        }
    }
    

    pub fn run(&mut self) -> ContigOrder {
        let mut best = self.contig_order.clone();
        let mut current = self.contig_order.clone();
        let mut t = self.temperature;
        let mut i = 0;

        fn acceptance_probability(current_cost: f64, new_cost: f64, temperature: f64) -> f64 {
            if new_cost > current_cost {
                1.0
            } else {
                (1.0 * (new_cost - current_cost) / temperature).exp()
            }
        }

        while t > self.mimimum_temperature && i < self.max_iteration {
            let mut new = current.clone();
            let mut i1 = rand::thread_rng().gen_range(0..new.contig_units.len());
            let mut i2 = rand::thread_rng().gen_range(0..new.contig_units.len());
            if i1 == i2 {
                continue;
            }

            if i1 > i2 {
                (i1, i2) = (i2, i1);
            }

            let operate_idx = rand::thread_rng().gen_range(0..100);
            match operate_idx {
                0..=49 => {
                    new.swap(i1, i2);
                }

                // 50..=99 => {
                //     new.crossover(i1, i2);
                // }
                50..=99 => {
                    new.reverse(i1, i2);
                }

                // 90..=99 => {
                //     new.rotate(i1);
                // }
                _ => todo!()
            }
    

            let current_cost = current.cost();
            let new_cost = new.cost();
            let ap = acceptance_probability(current_cost, new_cost, t);
            if ap > rand::thread_rng().gen_range(0.0..1.0) {
                current = new;
            }
            if current_cost > best.cost() {
                println!("Cost: {}, {}", current_cost, best.cost());
                best = current.clone();
            }
            // println!("Cost: {}", current_cost);
            t *= self.cooling_rate;
            i += 1;
        }
        println!("t: {}, i: {}", t, i);
        best
    }
}


pub struct GeneticAlgorithm {
    pub contig_order: ContigOrder,
    pub population_size: usize,
    pub max_iteration: usize,
    pub mutation_rate: f64,
    pub crossover_rate: f64,
}

impl GeneticAlgorithm {
    pub fn new(contig_order: ContigOrder, 
            population_size: usize, 
            max_iteration: usize, 
            mutation_rate: f64, 
            crossover_rate: f64) -> GeneticAlgorithm {
        GeneticAlgorithm {
            contig_order: contig_order,
            population_size: population_size,
            max_iteration: max_iteration,
            mutation_rate: mutation_rate,
            crossover_rate: crossover_rate,
        }
    }

    pub fn select(&self, population: &Vec<ContigOrder>) -> ContigOrder {
        let mut idx = rand::thread_rng().gen_range(0..population.len());
        let mut best = population[idx].clone();
        for _ in 0..10 {
            idx = rand::thread_rng().gen_range(0..population.len());
            if population[idx].cost() < best.cost() {
                best = population[idx].clone();
            }
        }
        best
    }

    fn crossover(&self, parent1: &ContigOrder, parent2: &ContigOrder) -> ContigOrder {
        let mut child = ContigOrder::new(Vec::new());
        let mut i1 = rand::thread_rng().gen_range(0..parent1.contig_units.len());
        let mut i2 = rand::thread_rng().gen_range(0..parent1.contig_units.len());
        if i1 == i2 {
            return parent1.clone();
        }
        if i1 > i2 {
            (i1, i2) = (i2, i1);
        }
        let mut i = 0;
        while i < i1 {
            child.contig_units.push(parent1.contig_units[i].clone());
            i += 1;
        }
        while i < i2 {
            child.contig_units.push(parent2.contig_units[i].clone());
            i += 1;
        }
        while i < parent1.contig_units.len() {
            child.contig_units.push(parent1.contig_units[i].clone());
            i += 1;
        }
        child
    }

    fn mutate(&self, child: &mut ContigOrder) {
        let mut i1 = rand::thread_rng().gen_range(0..child.contig_units.len());
        let mut i2 = rand::thread_rng().gen_range(0..child.contig_units.len());
        if i1 == i2 {
            return;
        }
        if i1 > i2 {
            (i1, i2) = (i2, i1);
        }
        child.reverse(i1, i2);
    }

    fn mutate2(&self, child: &mut ContigOrder) {
        let mut i1 = rand::thread_rng().gen_range(0..child.contig_units.len());
        let mut i2 = rand::thread_rng().gen_range(0..child.contig_units.len());
        if i1 == i2 {
            return;
        }
        if i1 > i2 {
            (i1, i2) = (i2, i1);
        }
        child.swap(i1, i2);
    }

    fn mutate3(&self, child: &mut ContigOrder) {
        let mut i1 = rand::thread_rng().gen_range(0..child.contig_units.len());
        let mut i2 = rand::thread_rng().gen_range(0..child.contig_units.len());
        if i1 == i2 {
            return;
        }
        if i1 > i2 {
            (i1, i2) = (i2, i1);
        }
        child.rotate(i1);
    }



    pub fn run(&mut self) -> ContigOrder {
        let mut population = Vec::<ContigOrder>::new();
        for _ in 0..self.population_size {
            let mut contig_order = self.contig_order.clone();
            contig_order.shuffle();
            population.push(contig_order);
        }

        let mut best = population[0].clone();
        let mut i = 0;
        while i < self.max_iteration {
            let mut new_population = Vec::<ContigOrder>::new();
            for _ in 0..self.population_size {
                let parent1 = self.select(&population);
                let parent2 = self.select(&population);
                let mut child = self.crossover(&parent1, &parent2);
                self.mutate(&mut child);
                new_population.push(child);
            }
            population = new_population;
            for contig_order in &population {
                if contig_order.cost() > best.cost() {
                    best = contig_order.clone();
                    println!("{}", best.cost());
                }
            }
            i += 1;
        }
        best
    }
}