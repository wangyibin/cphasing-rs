#![allow(dead_code)]
use lazy_static::lazy_static;
use ndarray::array;
use ndarray::prelude::*;
use itertools::Itertools;
use rand::Rng;
use rand::prelude::SliceRandom;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::core::ChromSize;


use crate::core::{ common_reader, common_writer };
use crate::core::BaseTable;

lazy_static! {
    #[derive(Debug)]
    pub static ref ORIENTATION_TO_IDX: HashMap<&'static str, u8> = {
        let mut m = HashMap::new();
        m.insert("++", 2);
        m.insert("+-", 3);
        m.insert("-+", 0);
        m.insert("--", 1);
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

lazy_static! {
    pub static ref BASE_MATRIX: Array2<f64> = 
             Array2::from_shape_vec((2, 2), vec![2.0, 3.0, 1.0, 2.0]).unwrap();
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

    pub fn get_score(&self, contig: &ContigUnit, distance: f64) -> Array2<f64> {

        let idx = (self.orientation * 2).pow(1) + (contig.orientation * 2).pow(0);
        // let mut matrix = BASE_MATRIX.clone();
        // matrix.map_inplace(|x| *x = *x as f64 + (distance - 1.0) * 2.0);
        let v1: f64 = ( distance * 2.0 - 1.0 ) * ( distance * 2.0 - 1.0 );
        let v2: f64 = ( distance * 2.0) * ( distance * 2.0 - 1.0);
        let v3: f64 = ( distance * 2.0) * ( distance * 2.0 - 1.0);
        let v4: f64 = ( distance * 2.0 ) * ( distance * 2.0 );

        let v: Vec<f64> = if self.orientation == 0 && contig.orientation == 0 {
            vec![v2, v4, v1, v3]
        } else if self.orientation == 1 && contig.orientation == 0 {
            vec![v1, v3, v2, v4]
        } else if self.orientation == 0 && contig.orientation == 1 {
            vec![v4, v2, v3, v1]
        } else {
            vec![v3, v1, v4, v2]
        }; 
        
        let mut matrix: Array2<f64> = Array2::from_shape_vec((2, 2), v).unwrap();
        let score: Array2<f64> = 
            Array2::from_shape_vec((2, 2), 
                                    self.scores.get(&contig.contig)
                                                    .unwrap_or(&vec![0.0, 0.0, 0.0, 0.0])
                                                    .clone()).unwrap();
        
        let score = score / &matrix;
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
        let contigs: Vec<ContigUnit> = self.contig_units.iter().sorted_by_key(|c| c.orderidx).cloned().collect();
        for (i, contig1) in contigs.iter().enumerate() {
            
            for j in (i + 1)..contigs.len() {
                
                let contig2 = &contigs[j];
                // eprintln!("{} {} {} {}", contig1.contig, contig1.orientation, contig2.contig, contig2.orientation);
                let distance = (j - i) as f64;
                let score = contig1.get_score(&contig2, distance).sum();
                
                cost += score;
            }
        }

        cost
    }

    pub fn shuffle(&mut self) {
        let mut rng = rand::thread_rng();
        // self.contig_units.shuffle(&mut rng);
        for (i, contig) in self.contig_units.iter_mut().enumerate() {
            contig.orderidx = i;
            contig.orientation = rng.gen_range(0..2);

        }
        
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        for contig in self.contig_units.iter_mut() {
            if contig.orderidx == i {
                contig.orderidx = j;
            } else if contig.orderidx == j {
                contig.orderidx = i;
            }
        }
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

    pub fn inversion(&mut self, i: usize, j: usize) {

        for contig in self.contig_units.iter_mut() {
            if contig.orderidx >= i && contig.orderidx <= j {
                contig.orderidx = j - contig.orderidx + i;
                contig.orientation = contig.orientation ^ 1;
            }
        }
    }

    pub fn insertion(&mut self, i: usize, j: usize) {

        if rand::thread_rng().gen_range(0.0..1.0) < 0.5 {
            // insertion i to j
            for contig in self.contig_units.iter_mut() {
                if contig.orderidx == i {
                    contig.orderidx = j;
                } else if contig.orderidx > i {
                    if contig.orderidx <= j {
                        contig.orderidx -= 1;
                }
                }
            }
        } else {
            // insert j to i
            for contig in self.contig_units.iter_mut() {
                if contig.orderidx == j {
                    contig.orderidx = i;
                } else if contig.orderidx >= i {
                    if contig.orderidx < j {
                        contig.orderidx += 1;
                    }
                }
            }
        }
        
        
    }

    pub fn splice(&mut self, i: usize) {
        let total_length: usize = self.contig_units.len();
        let before_length: usize = i + 1;
        let after_length: usize = total_length - before_length; 
        
        for k in 0..(before_length + 1) {
            self.contig_units[k].orderidx = k + after_length;
            
        }
        for k in before_length..total_length {
            self.contig_units[k].orderidx = k - before_length;
            
        }
    }

    pub fn contigs(&self) -> Vec<String> {
        let contigs: Vec<ContigUnit> = self.contig_units.iter().sorted_by_key(|c| c.orderidx).cloned().collect();
        let contigs: Vec<String> = contigs.iter().map(|c| c.contig.clone()).collect();
        contigs
    }

    pub fn orderidxes(&self) -> Vec<usize> {
        let orderidxes: Vec<ContigUnit> = self.contig_units.iter().sorted_by_key(|c| c.orderidx).cloned().collect();
        let orderidxes: Vec<usize> = orderidxes.iter().map(|c| c.orderidx.clone()).collect();
        orderidxes
    }

    pub fn orientations(&self) -> Vec<u8> {
        let orientations: Vec<ContigUnit> = self.contig_units.iter().sorted_by_key(|c| c.orderidx).cloned().collect();
        let orientations: Vec<u8> = orientations.iter().map(|c| c.orientation.clone()).collect();
        orientations
    }

    pub fn save(&self, output: &String) {
        let mut wtr = common_writer(output);    

        let contigs = self.contigs();
    
        wtr.write_all(contigs.join("\n").as_bytes()).unwrap();

        log::info!("Successful output optimize result in `{}`", output);

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
        let tmp = self.contig_order.clone();
        let mut current = self.contig_order.clone();
        current.shuffle();
        let mut best = current.clone();

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
                0..=9 => {
                    new.swap(i1, i2);
                }

                10..=49 => {
                    new.insertion(i1, i2);
                }
                
                50..=69 => {
                    new.inversion(i1, i2);
                }

                70..=79 => {
                    new.splice(i1);
                }

                80..=99 => {
                    new.rotate(i1);
                }
                _ => todo!()
            }
    
            
            let current_cost = current.cost();
            let new_cost = new.cost();
            let ap = acceptance_probability(current_cost, new_cost, t);
            if current_cost > best.cost() {
                eprintln!("Cost: {}, {}", current_cost, best.cost());
                best = current.clone();
            } else {
                if ap > rand::thread_rng().gen_range(0.0..1.0) {
                current = new;
                }
            }
            
            // eprintln!("Cost: {}", current_cost);
            t *= self.cooling_rate;
            i += 1;
        }
        // eprintln!("t: {}, i: {}", t, i);
        eprintln!("{}", tmp.cost());
        eprintln!("{:.?}", tmp.contigs());
        eprintln!("{:.?}", tmp.orientations());
        eprintln!("{:.?}", best.contigs());
        eprintln!("{:.?}", best.orientations());
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
        child.inversion(i1, i2);
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
                    eprintln!("{}", best.cost());
                }
            }
            i += 1;
        }
        best
    }
}