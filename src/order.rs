#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use rand::Rng;
use rand::rngs::SmallRng;
use rand::distributions::{Bernoulli, Distribution};
use rand::prelude::*;
use std::borrow::Cow;
use hashbrown::HashMap;
use indexmap::IndexMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::{Arc, RwLock, Mutex};
use std::fmt::{self, Display, Arguments};
use std::cell::RefCell;
use std::marker::PhantomData;
use std::time::Instant;
use itertools::Itertools;
use genetic_algorithm::strategy::evolve::prelude::*;
use genetic_algorithm::strategy::hill_climb::prelude::*;
use genetic_algorithm::fitness::prelude::*;
use elkai_rs::DistanceMatrix;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

#[derive(Debug, Clone)]
pub struct ContactMatrix<'a> {
    num_contigs: usize,
    length_data: &'a IndexMap<usize, usize>,
    matrix: Vec<HashMap<usize, f64>>, 
    adjacency: Vec<Vec<(usize, f32)>>, 
    flat_edges: Vec<(usize, usize, f32)>,
}

impl<'a> ContactMatrix<'a> {
    pub fn new(length_data: &'a IndexMap<usize, usize>, contacts: Vec<HashMap<usize, u32>>) -> Self {
        let num_contigs = length_data.len();

        let mut adjacency = Vec::with_capacity(num_contigs);
        for contact_map in &contacts {
            let mut row: Vec<(usize, f32)> = contact_map.into_iter()
                .filter(|&(_, w)| *w > 0)
                .map(|(target, weight)| (target.clone(), weight.clone() as f32))
                .collect();
            
            row.sort_unstable_by_key(|k| k.0); 
            
            adjacency.push(row);
        }

        let mut matrix = Vec::with_capacity(num_contigs);
        for contact_map in contacts.clone() {
            let mut row = HashMap::with_capacity(contact_map.len());
            for (target, &weight) in contact_map.iter() {
                if weight > 0 {
                    row.insert(*target, weight as f64);
                }
            }
            matrix.push(row);
        }
        let mut flat_edges = Vec::with_capacity(num_contigs * 10); 
        for (u, contact_map) in contacts.iter().enumerate() {
            let mut edges: Vec<_> = contact_map.iter().collect();
            edges.sort_by_key(|(v, _)| **v); 
            for (v, weight) in edges {
                if *weight > 0 && u < *v { 
                    flat_edges.push((u, *v, *weight as f32));
                }
            }
        }
        ContactMatrix {
            num_contigs,
            length_data,
            matrix,
            adjacency,
            flat_edges
        }
    }


    #[inline(always)]
    fn get(&self, i: usize, j: usize) -> f64 {
        if let Some(row) = self.matrix.get(i) {
            *row.get(&j).unwrap_or(&0.0)
        } else {
            0.0
        }
    }
    #[inline(always)]
    fn get_adjacency(&self, i: usize, j: usize) -> f64 {
        if let Some(row) = self.matrix.get(i) {
            *row.get(&j).unwrap_or(&0.0)
        } else {
            0.0
        }
    }
}

#[derive(Debug, Clone)]
pub struct Tour<T = usize> {
    pub contigs: Vec<T>, 
    pub signs: Vec<bool>, 
}

impl<T: std::fmt::Display> Display for Tour<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, contig) in self.contigs.iter().enumerate() {
            let orientation = if self.signs[i] { "+" } else { "-" };
            write!(f, "{}{} ", contig, orientation)?;
        }
        Ok(())
    }
}

impl<T> Tour<T> {
    pub fn into_iter(self) -> std::vec::IntoIter<T> {
        self.contigs.into_iter()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.contigs.iter()
    }

    pub fn map<U, F>(self, mut f: F) -> Tour<U> 
    where F: FnMut(T) -> U {
        Tour {
            contigs: self.contigs.into_iter().map(f).collect(),
            signs: self.signs,
        }
    }

    pub fn reverse(&mut self) {
        self.contigs.reverse();
        // signs reverse and + to - - to + 
        self.signs.reverse();
        for sign in self.signs.iter_mut() {
            *sign = !*sign;
        }
    }
}


thread_local! {
    static MIDS_BUF: RefCell<Vec<f64>> = RefCell::new(Vec::new());
    static POS_BUF: RefCell<Vec<usize>> = RefCell::new(Vec::new());
}


#[derive(Clone, Debug)]
struct MyFitness<'a> {
    matrix: &'a ContactMatrix<'a>,
}

impl Fitness for MyFitness<'_> {
    type Genotype = UniqueGenotype; 
 
    fn calculate_for_chromosome(&mut self, chromosome: &FitnessChromosome<Self>,
                                _genotype: &FitnessGenotype<Self>,) -> Option<FitnessValue> {
        let genes = &chromosome.genes;
        let matrix = &self.matrix;
        let n = genes.len();
        let num_all_contigs = matrix.num_contigs;
        let mut score = 0.0;
      
        MIDS_BUF.with(|mids_cell| {
            POS_BUF.with(|pos_cell| {
                let mut mids = mids_cell.borrow_mut();
                let mut pos_map = pos_cell.borrow_mut();

                mids.clear();
                if mids.capacity() < n { mids.reserve(n); }
                
                let mut cum_sum = 0.0;

                for &idx in genes {
                    let size = unsafe { *matrix.length_data.get(&idx).unwrap_or(&0) } as f64;
                    mids.push(cum_sum + size * 0.5);
                    cum_sum += size;
                }

                if pos_map.len() < num_all_contigs {
                    pos_map.resize(num_all_contigs, usize::MAX);
                }
                
                for (idx, &id) in genes.iter().enumerate() {
                    pos_map[id] = idx;
                }

                // for (i, &u) in genes.iter().enumerate() {
                //     if let Some(row) = matrix.matrix.get(u) {
                //         for (&v, &contact) in row {
                //             if v >= pos_map.len() { continue; }
                            
                //             let j = pos_map[v];

                //             if j != usize::MAX && j > i {
                //                 let dist = (mids[j] - mids[i]).abs() as f32;
                    
                //                 let log_dist = (dist + 1.0).ln();
                //                 score += contact as f32 * log_dist;
                       
                //             }
                //         }
                //     }
                // }
                // for (i, &u) in genes.iter().enumerate() {
                //     if let Some(row) = matrix.adjacency.get(u) { // 改为 adjacency
                //         for &(v, contact) in row {
                //             if v >= pos_map.len() { continue; }
                //             let j = pos_map[v];
                //             if j != usize::MAX && j > i {
                //                 let dist = (mids[j] - mids[i]).abs();
                //                 let log_dist = (dist + 1.0).ln() as f32;
                //                 score += contact * log_dist;
                //             }
                //         }
                //     }
                // }
                for &(u, v, weight) in &matrix.flat_edges {
                    let idx_u = pos_map[u];
                    let idx_v = pos_map[v];

                    if idx_u != usize::MAX && idx_v != usize::MAX {
                        let dist = (mids[idx_u] - mids[idx_v]).abs();
                        let log_dist = (dist + 1.0).ln() as f32;
                        score += weight * log_dist;
                    }
                }

                for &id in genes {
                    pos_map[id] = usize::MAX;
                }
            })
        });

        
        Some((score * 1000000.0) as isize)     


    }

}


pub fn calculate_fitness(genes: &[usize], matrix: &ContactMatrix) -> isize {
    let n = genes.len();
    let num_all_contigs = matrix.num_contigs;
    let mut score = 0.0;
    
    MIDS_BUF.with(|mids_cell| {
        POS_BUF.with(|pos_cell| {
            let mut mids = mids_cell.borrow_mut();
            let mut pos_map = pos_cell.borrow_mut();

            mids.clear();
            if mids.capacity() < n { mids.reserve(n); }
            
            let mut cum_sum = 0.0;

            for &idx in genes {
                let size = *matrix.length_data.get(&idx).unwrap_or(&0) as f64;
                mids.push(cum_sum + size * 0.5);
                cum_sum += size;
            }

            if pos_map.len() < num_all_contigs {
                pos_map.resize(num_all_contigs, usize::MAX);
            }
            
            for (idx, &id) in genes.iter().enumerate() {
                pos_map[id] = idx;
            }

            // for (i, &u) in genes.iter().enumerate() {
            //     if let Some(row) = matrix.matrix.get(u) {
            //         for (&v, &contact) in row {
            //             if v >= pos_map.len() { continue; }
                        
            //             let j = pos_map[v];

            //             if j != usize::MAX && j > i {
            //                 let dist = (mids[j] - mids[i]).abs() as f64;
            //                 let log_dist = (dist + 1.0).ln();

                            
            //                 score += contact * log_dist;
                          
            //             }
            //         }
            //     }
            // }
            // for (i, &u) in genes.iter().enumerate() {
            //     if let Some(row) = matrix.adjacency.get(u) {
            //         for &(v, contact) in row {
            //             if v >= pos_map.len() { continue; }
                        
            //             let j = pos_map[v];

            //             if j != usize::MAX && j > i {
            //                 let dist = (mids[j] - mids[i]).abs(); // f64
                            
            //                 let log_dist = (dist + 1.0).ln() as f32;

            //                 score += contact * log_dist;
            //             }
            //         }
            //     }
            // }
            for &(u, v, weight) in &matrix.flat_edges {
              
                if let (Some(&idx_u), Some(&idx_v)) = (pos_map.get(u), pos_map.get(v)) {
                    if idx_u != usize::MAX && idx_v != usize::MAX {
                        let dist = (mids[idx_u] - mids[idx_v]).abs();
                        
                        let log_dist = (dist + 1.0).ln() as f32;
                        score += weight * log_dist;
                    }
                }
            }
            for &id in genes {
                pos_map[id] = usize::MAX;
            }
        })
    });

    // let final_score = if score.is_finite() { score } else { 0.0 };
    
    (score * 1000000.0) as isize
}


pub fn calculate_cost(i: usize, j: usize, matrix: &ContactMatrix) -> f64 {
    let len_i = matrix.length_data.get(&i).unwrap_or(&1);
    let len_j = matrix.length_data.get(&j).unwrap_or(&1);
    let dist = (*len_i + *len_j) as f64 / 2.0;
    let log_dist = (dist + 1.0).ln();
    
    let contact = matrix.get(i, j);
    
    if log_dist > 1e-6 {
        contact / log_dist
    } else {
        contact * 1e4
    }
}
    

pub fn two_opt_optimization(genes: &mut Vec<usize>, matrix: &ContactMatrix, max_iterations: usize) {
    let n = genes.len();
    if n < 3 {
        return;
    }
    log::info!("Starting 2-Opt optimization with {} contigs...", n);

   
    let lengths: Vec<f64> = (0..matrix.num_contigs)
        .map(|i| *matrix.length_data.get(&i).unwrap_or(&1) as f64)
        .collect();

    let mut improved = true;
    let mut iterations = 0;

    while improved && iterations < max_iterations {
        improved = false;
        iterations += 1;

        let best_move = (0..n - 2).into_par_iter().map(|i| {
            let mut local_best_delta = 0.0;
            let mut local_best_move = None;

            let n_i = genes[i];
            let n_i1 = genes[i + 1];


            let w_i_i1 = matrix.get(n_i, n_i1);
            let cost_i_i1 = if w_i_i1 > 0.0 {
                let dist = (lengths[n_i] + lengths[n_i1]) * 0.5;
                let log_dist = (dist + 1.0).ln();
                if log_dist > 1e-6 { w_i_i1 / log_dist } else { w_i_i1 }
            } else {
                0.0
            };

            for j in (i + 2)..n - 1 {
                let n_j = genes[j];
                let n_j1 = genes[j + 1];

                let w_i_j = matrix.get(n_i, n_j);
                let w_i1_j1 = matrix.get(n_i1, n_j1);

                if w_i_j == 0.0 && w_i1_j1 == 0.0 {
                    continue;
                }

                let cost_i_j = if w_i_j > 0.0 {
                    let dist = (lengths[n_i] + lengths[n_j]) * 0.5;
                    let log_dist = (dist + 1.0).ln();
                    if log_dist > 1e-6 { w_i_j / log_dist } else { w_i_j }
                } else { 0.0 };

                let cost_i1_j1 = if w_i1_j1 > 0.0 {
                    let dist = (lengths[n_i1] + lengths[n_j1]) * 0.5;
                    let log_dist = (dist + 1.0).ln();
                    if log_dist > 1e-6 { w_i1_j1 / log_dist } else { w_i1_j1 }
                } else { 0.0 };

                let new_cost = cost_i_j + cost_i1_j1;

    
                let w_j_j1 = matrix.get(n_j, n_j1);
                let cost_j_j1 = if w_j_j1 > 0.0 {
                    let dist = (lengths[n_j] + lengths[n_j1]) * 0.5;
                    let log_dist = (dist + 1.0).ln();
                    if log_dist > 1e-6 { w_j_j1 / log_dist } else { w_j_j1 }
                } else { 0.0 };

                let old_cost = cost_i_i1 + cost_j_j1;

                let delta = new_cost - old_cost;

                if delta > 0.0 && delta > local_best_delta {
                    local_best_delta = delta;
                    local_best_move = Some((delta, i, j));
                }
            }
            local_best_move
        })
        .reduce(
            || None,
            |a, b| {
                match (a, b) {
                    (None, None) => None,
                    (Some(v), None) => Some(v),
                    (None, Some(v)) => Some(v),
                    (Some(va), Some(vb)) => {
                        if va.0 > vb.0 { Some(va) } else { Some(vb) }
                    }
                }
            }
        );

        if let Some((_, best_i, best_j)) = best_move {
            genes[best_i + 1..=best_j].reverse();
            improved = true;
        }
    }
}


pub fn block_move_optimization(genes: &mut Vec<usize>, matrix: &ContactMatrix, max_iterations: usize) {
    let n = genes.len();
    if n < 3 { return; }
    log::info!("Starting Block-Move optimization with {} contigs...", n);

    let lengths: Vec<f64> = (0..matrix.num_contigs)
        .map(|i| *matrix.length_data.get(&i).unwrap_or(&1) as f64)
        .collect();

    let mut improved = true;
    let mut iterations = 0;
    
    let get_cost = |i: usize, j: usize| -> f64 {
        let w = matrix.get(i, j);
        if w <= 0.0 { return 0.0; }
        let dist = (lengths[i] + lengths[j]) * 0.5;
        let log_dist = (dist + 1.0).ln();
        if log_dist > 1e-6 { w / log_dist } else { w }
    };

    while improved && iterations < max_iterations {
        improved = false;
        iterations += 1;

        let best_move = (0..n).into_par_iter().map(|i| {
            let mut local_best = None;
            let mut local_max_delta = 0.0;
            
            for j in i..n {

                if i == 0 && j == n - 1 { continue; }
                
                let n_i = genes[i];
                let n_j = genes[j];
                
                let n_prev = if i > 0 { Some(genes[i-1]) } else { None };
                let n_next = if j < n - 1 { Some(genes[j+1]) } else { None };
     
                let mut current_edges_cost = 0.0;
                if let Some(p) = n_prev { current_edges_cost += get_cost(p, n_i); }
                if let Some(nx) = n_next { current_edges_cost += get_cost(n_j, nx); }
                
               
                for k_idx in 0..=n {
       
                    if k_idx >= i && k_idx <= j + 1 { continue; }
                    
                    let target_prev = if k_idx > 0 { Some(genes[k_idx-1]) } else { None };
                    let target_next = if k_idx < n { Some(genes[k_idx]) } else { None };
                    
                    let mut target_edges_cost = 0.0;
                    if let (Some(tp), Some(tn)) = (target_prev, target_next) {
                        target_edges_cost += get_cost(tp, tn);
                    }

                    let mut new_edges_cost = 0.0;

                    if let (Some(p), Some(nx)) = (n_prev, n_next) {
                        new_edges_cost += get_cost(p, nx);
                    }
                    
                    if let Some(tp) = target_prev {
                        new_edges_cost += get_cost(tp, n_i);
                    }
                    if let Some(tn) = target_next {
                        new_edges_cost += get_cost(n_j, tn);
                    }
                    
                    let delta = new_edges_cost - (current_edges_cost + target_edges_cost);
                    
                    if delta > 1e-9 && delta > local_max_delta {
                        local_max_delta = delta;
                        local_best = Some((delta, i, j, k_idx));
                    }
                }
            }
            local_best
        }).reduce(|| None, |a, b| {
            match (a, b) {
                (None, None) => None,
                (Some(v), None) => Some(v),
                (None, Some(v)) => Some(v),
                (Some(va), Some(vb)) => if va.0 > vb.0 { Some(va) } else { Some(vb) },
            }
        });
        
        if let Some((_, best_i, best_j, best_k)) = best_move {
            let block: Vec<_> = genes.drain(best_i..=best_j).collect();
            
            let insert_idx = if best_k > best_i {
                best_k - (best_j - best_i + 1)
            } else {
                best_k
            };
            
            genes.splice(insert_idx..insert_idx, block);
            improved = true;
        }
    }
}


pub fn two_opt_optimization_fitness(genes: &mut Vec<usize>, matrix: &ContactMatrix, max_iterations: usize) {
    let n = genes.len();
    let mut improved = true;
    let mut iterations = 0;

    let mut current_score = calculate_fitness(genes, &matrix);

    while improved && iterations < max_iterations {
        improved = false;
        iterations += 1;

        let best_move = (0..n - 1).into_par_iter().map(|i| {
            let mut local_best_diff = 0;
            let mut local_best_j = 0;
            
            
            let mut best_move_for_i = None;
            let mut max_score_for_i = -1.0; 

            for j in i + 1..n {
                
                let mut new_genes = genes.clone();
                new_genes[i..=j].reverse();
                let new_score = calculate_fitness(&new_genes, &matrix);
                
                if new_score > current_score {
                    if new_score > max_score_for_i as isize {
                        max_score_for_i = new_score as f64;
                        best_move_for_i = Some((i, j, new_score));
                    }
                }
            }
            best_move_for_i
        })
        .reduce(
            || None,
            |a, b| {
                match (a, b) {
                    (Some((_, _, score_a)), Some((_, _, score_b))) => {
                        if score_a > score_b { a } else { b }
                    }
                    (Some(_), None) => a,
                    (None, Some(_)) => b,
                    (None, None) => None,
                }
            }
        );

        if let Some((best_i, best_j, best_new_score)) = best_move {
            if best_new_score > current_score {
                genes[best_i..=best_j].reverse();
                current_score = best_new_score;
                improved = true;
            }
        }
    }
}


pub fn run_lkh_optimizer(
    tour: &Tour,
    contigsizes: IndexMap<usize, usize>,
    matrix: &ContactMatrix,
    iterations: usize,
    seed: u64,
) -> Tour {
    let contacts = &matrix.matrix;
    let num_contigs = contigsizes.len();
    log::info!("Starting LKH optimization with {} contigs...", num_contigs);

    // filter contigsizes by length 
    // let min_length = 20000;
    // let contigsizes: IndexMap<usize, usize> = contigsizes.iter()
    //     .filter(|&(_, &size)| size >= min_length)
    //     .map(|(&id, &size)| (id, size))
    //     .collect();

    let mut id_to_idx = HashMap::new();
    let mut idx_to_id = Vec::with_capacity(num_contigs);
    let mut contigsizes_idx = IndexMap::new();
    for (i, id) in tour.contigs.clone().into_iter().enumerate() {
        id_to_idx.insert(id, i);
        contigsizes_idx.insert(i, *contigsizes.get(&id).unwrap());
        idx_to_id.push(id);
    }

    // Construct distance matrix for LKH
    // Dummy node at the end
    // Dimension: N + 1
    // Distances scaled to [0, 100000]
    let dim = num_contigs + 1;
    let mut dist_matrix = vec![vec![0; dim]; dim];
    
    // Scale factor for distance normalization
    const SCALE_FACTOR: f64 = 100_000.0;
    // Maximum distance for no-contact pairs
    const MAX_DIST: u32 = 10_000_000; 

    // for (u, neighbors) in contacts.iter().enumerate() {
    //     if let Some(&u_idx) = id_to_idx.get(&u) {
    //         let length_u = *contigsizes_idx.get(&u_idx).unwrap_or(&1) as f64;
    //         for (&v, &weight) in neighbors {
    //             if let Some(&v_idx) = id_to_idx.get(&v) {
    //                 if u_idx == v_idx { continue; }
    //                 let length_v = *contigsizes_idx.get(&v_idx).unwrap_or(&1) as f64;

    //                 let dist = if weight > 0.0 {
    //                     let w = weight as f64;
    //                     let min_len = length_u.min(length_v);
    //                     let boost_factor = if min_len < 50_000.0 {
    //                         (50_000.0 / min_len.max(100.0)).sqrt()
       
    //                     } else {
    //                         1.0
    //                     };
    //                     let w_norm = w * boost_factor;

    //                     let log_w = (w_norm + 1.0).ln();
    //                     let max_log = 20.0; 
    //                     let norm = (log_w / max_log).min(1.0);
                        
    //                     ((1.0 - norm) * SCALE_FACTOR) as u32
                  
    //                 } else {
    //                     MAX_DIST
    //                 };
                   
    //                 dist_matrix[u_idx][v_idx] = dist;
    //                 dist_matrix[v_idx][u_idx] = dist;
    //             }
    //         }
    //     }
    // }

    let mut max_weights = vec![1.0; num_contigs];
    for (u, neighbors) in contacts.iter().enumerate() {
        if let Some(&u_idx) = id_to_idx.get(&u) {
            let max_w = neighbors.values().cloned().fold(0u32, |a, b| a.max(b as u32));
            max_weights[u_idx] = max_w as f64;
        }
    }

    for (u, neighbors) in contacts.iter().enumerate() {
        if let Some(&u_idx) = id_to_idx.get(&u) {
            
            for (&v, &weight) in neighbors {
                if let Some(&v_idx) = id_to_idx.get(&v) {
                    if u_idx == v_idx { continue; }
                    
                    
                    let w = weight as f64;
                    let score_u = w / max_weights[u_idx];
                    let score_v = w / max_weights[v_idx];

                    let combined_score = (score_u * score_v).sqrt();
                    
    
                    let dist_val = (1.0 - combined_score.powf(0.5)) * SCALE_FACTOR;
                    
                    let dist = dist_val.max(1.0) as u32;
                   
                    dist_matrix[u_idx][v_idx] = dist;
                    dist_matrix[v_idx][u_idx] = dist;
                }
            }
        }
    }
    // let mut max_weights = vec![1.0; num_contigs];
    // for (u, neighbors) in contacts.iter().enumerate() {
    //     if let Some(&u_idx) = id_to_idx.get(&u) {
    //         let max_w = neighbors.values().cloned().fold(0u32, |a, b| a.max(b as u32));
    //         max_weights[u_idx] = max_w as f64;
    //     }
    // }

    // for (u, neighbors) in contacts.iter().enumerate() {
    //     if let Some(&u_idx) = id_to_idx.get(&u) {
    //         let length_u = *contigsizes_idx.get(&u_idx).unwrap_or(&1) as f64;
            
    //         for (&v, &weight) in neighbors {
    //             if let Some(&v_idx) = id_to_idx.get(&v) {
    //                 if u_idx == v_idx { continue; }
    //                 let length_v = *contigsizes_idx.get(&v_idx).unwrap_or(&1) as f64;
                    
    //                 let w = weight as f64;
    //                 let score_u = w / max_weights[u_idx];
    //                 let score_v = w / max_weights[v_idx];
                    

    //                 let min_len = length_u.min(length_v);
    //                 let len_penalty = if min_len < 10_000.0 {
    //                     0.5 + 0.5 * (min_len / 10_000.0)
    //                 } else {
    //                     1.0
    //                 };

    //                 let combined_score = (score_u * score_v).sqrt() * len_penalty;
                    

    //                 let dist_val = (1.0 - combined_score.powf(0.5)) * SCALE_FACTOR;
                    
    //                 let dist = dist_val.max(1.0) as u32;
                   
    //                 dist_matrix[u_idx][v_idx] = dist;
    //                 dist_matrix[v_idx][u_idx] = dist;
    //             }
    //         }
    //     }
    // }

    for i in 0..num_contigs {
        dist_matrix[num_contigs][i] = 0;
        dist_matrix[i][num_contigs] = 0;
    }
    
    // To ensure symmetry and handle missing edges    
    for i in 0..num_contigs {
        for j in (i + 1)..num_contigs {
            if dist_matrix[i][j] == 0 {
                dist_matrix[i][j] = MAX_DIST;
                dist_matrix[j][i] = MAX_DIST;
            }
        }
    }

    let contig_dist_matrix = DistanceMatrix::new(dist_matrix);


    log::info!("Ordering with LKH...");
    let mut results = contig_dist_matrix.solve(iterations);
    
    let dummy_idx = num_contigs;

    let dummy_position = results.iter().position(|&x| x == dummy_idx).unwrap_or(0);
    log::info!("Reconnected cycle at dummy position.");

    results.rotate_left(dummy_position);

    let mut raw_indices: Vec<usize> = results.into_par_iter()
        .filter(|&x| x < num_contigs)
        .collect();

    let best_ordered_ids = raw_indices.iter().cloned().map(|idx| idx_to_id[idx]).collect::<Vec<usize>>();
    log::info!("Raw ordering LKH score: {}", calculate_fitness(&best_ordered_ids, matrix));

    let (best_ordered_ids, best_score) = (0..num_contigs).into_par_iter().flat_map(|i| {
        let mut candidates = Vec::with_capacity(4);

        let mut rot = raw_indices.clone();
        rot.rotate_left(i);
        candidates.push(rot);

        let mut global_rev = raw_indices.clone();
        global_rev.reverse();
        global_rev.rotate_left(i);
        candidates.push(global_rev);

        let mut rev_suffix = raw_indices.clone();
        rev_suffix[i..].reverse();
        rev_suffix.rotate_left(i);
        candidates.push(rev_suffix);

        let mut rev_prefix = raw_indices.clone();
        rev_prefix[..=i].reverse();
        if i + 1 < num_contigs {
            rev_prefix.rotate_left(i + 1);
            candidates.push(rev_prefix);
        } else {
            candidates.push(rev_prefix);
        }

        candidates.into_iter().map(|indices| {
            let ordered_ids: Vec<usize> = indices.iter().map(|&idx| idx_to_id[idx]).collect();
            let score = calculate_fitness(&ordered_ids, matrix);
            (ordered_ids, score)
        }).collect::<Vec<_>>() 
        
    }).min_by_key(|&(_, score)| score)
    .unwrap();

    log::info!("Best LKH score: {}", best_score);

    let mut optimized_order = best_ordered_ids;
    log::info!("LKH optimization finished.");
   
    Tour {
        contigs: optimized_order,
        signs: tour.signs.clone(),
    }
}

pub fn run_lkh_optimizer_dual_node(
    tour: &Tour,
    contigsizes: IndexMap<usize, usize>,
    matrix: &ContactMatrix, // 注意：这里需要更详细的 SplitContacts 信息，目前的 ContactMatrix 可能不够
    // 如果 matrix 只包含 contig 级别的互作，我们需要 split contacts 数据
    // 假设 matrix 已经被修改为包含 split 互作，或者我们需要传入 split_contacts
    // 为了演示，这里假设我们能从 matrix 获取端点互作，或者我们需要修改函数签名传入 split_contacts
    split_contacts: &crate::splitcontacts::SplitContacts, 
    contig2idx: &HashMap<String, usize>,
    iterations: usize,
    seed: u64,
) -> Tour {
    let num_contigs = contigsizes.len();
    log::info!("Starting LKH Dual-Node optimization with {} contigs ({} nodes)...", num_contigs, num_contigs * 2);

    // 1. 映射 ID
    let mut id_to_idx = HashMap::new();
    let mut idx_to_id = Vec::with_capacity(num_contigs);
    for (i, id) in tour.contigs.clone().into_iter().enumerate() {
        id_to_idx.insert(id, i);
        idx_to_id.push(id);
    }

    // 2. 构建距离矩阵 (2N + 1 维度，最后是 Dummy)
    let num_nodes = 2 * num_contigs;
    let dim = num_nodes + 1;
    let mut dist_matrix = vec![vec![0; dim]; dim];
    
    const SCALE_FACTOR: f64 = 100_000.0;
    const MAX_DIST: u32 = 10_000_000; 
    const INTERNAL_DIST: u32 = 0; // 强制连接 Head-Tail

    // 初始化为最大距离
    for i in 0..dim {
        for j in 0..dim {
            if i != j {
                dist_matrix[i][j] = MAX_DIST;
            }
        }
    }

    // 3. 设置 Contig 内部距离 (Head <-> Tail)
    for i in 0..num_contigs {
        let head = 2 * i;
        let tail = 2 * i + 1;
        dist_matrix[head][tail] = INTERNAL_DIST;
        dist_matrix[tail][head] = INTERNAL_DIST;
    }

    // 4. 填充 Contig 间距离
    // 我们需要遍历 SplitContacts 来获取端点间的互作
    // 假设 split_contacts.data 的 key 是 (ContigName, ContigName)，value 是 [0-0, 0-1, 1-0, 1-1]
    // 0 代表 Head (Start), 1 代表 Tail (End)
    
    // 预计算最大权重用于归一化 (可选，但推荐)
    let mut max_weights: Vec<f64> = vec![1.0; num_nodes];
    // 这里简化处理，直接用绝对值转换，或者您可以实现类似之前的相对强度逻辑
    for (pair, counts) in &split_contacts.data {
        let c1_name = &pair.Contig1;
        let c2_name = &pair.Contig2;

        if let (Some(&u_orig_idx), Some(&v_orig_idx)) = (contig2idx.get(c1_name), contig2idx.get(c2_name)) {
            if let (Some(&u_idx), Some(&v_idx)) = (id_to_idx.get(&u_orig_idx), id_to_idx.get(&v_orig_idx)) {
                if u_idx == v_idx { continue; }

                let u_head = 2 * u_idx;
                let u_tail = 2 * u_idx + 1;
                let v_head = 2 * v_idx;
                let v_tail = 2 * v_idx + 1;

                // counts: [0-0, 0-1, 1-0, 1-1] -> [H-H, H-T, T-H, T-T]
                let w_hh = counts[0];
                let w_ht = counts[1];
                let w_th = counts[2];
                let w_tt = counts[3];

                // 更新 u_head 的最大权重
                max_weights[u_head] = max_weights[u_head].max(w_hh).max(w_ht);
                // 更新 u_tail 的最大权重
                max_weights[u_tail] = max_weights[u_tail].max(w_th).max(w_tt);
                // 更新 v_head 的最大权重
                max_weights[v_head] = max_weights[v_head].max(w_hh).max(w_th);
                // 更新 v_tail 的最大权重
                max_weights[v_tail] = max_weights[v_tail].max(w_ht).max(w_tt);
            }
        }
    }

    for (pair, counts) in &split_contacts.data {
        let c1_name = &pair.Contig1;
        let c2_name = &pair.Contig2;

        if let (Some(&u_orig_idx), Some(&v_orig_idx)) = (contig2idx.get(c1_name), contig2idx.get(c2_name)) {
            // 转换到当前 tour 的索引
            if let (Some(&u_idx), Some(&v_idx)) = (id_to_idx.get(&u_orig_idx), id_to_idx.get(&v_orig_idx)) {
                if u_idx == v_idx { continue; }

                // u 的节点: 2*u_idx (Head), 2*u_idx+1 (Tail)
                let u_head = 2 * u_idx;
                let u_tail = 2 * u_idx + 1;
                let v_head = 2 * v_idx;
                let v_tail = 2 * v_idx + 1;

                // counts: [0-0, 0-1, 1-0, 1-1] -> [H-H, H-T, T-H, T-T]
                let w_hh = counts[0];
                let w_ht = counts[1];
                let w_th = counts[2];
                let w_tt = counts[3];
                let weight_to_dist = |w: f64| -> u32 {
                    if w <= 0.0 { return MAX_DIST; }
                    // 简单的对数转换，您可以替换为更复杂的归一化
                    let log_w = (w + 1.0).ln();
                    let max_log = 15.0; // 经验值
                    let norm = (log_w / max_log).min(1.0);
                    ((1.0 - norm) * SCALE_FACTOR) as u32
                };

                // 填充矩阵 (对称)
                // Head-Head
                let d_hh = weight_to_dist(w_hh);
                if d_hh < dist_matrix[u_head][v_head] {
                    dist_matrix[u_head][v_head] = d_hh;
                    dist_matrix[v_head][u_head] = d_hh;
                }

                // Head-Tail
                let d_ht = weight_to_dist(w_ht);
                if d_ht < dist_matrix[u_head][v_tail] {
                    dist_matrix[u_head][v_tail] = d_ht;
                    dist_matrix[v_tail][u_head] = d_ht;
                }

                // Tail-Head
                let d_th = weight_to_dist(w_th);
                if d_th < dist_matrix[u_tail][v_head] {
                    dist_matrix[u_tail][v_head] = d_th;
                    dist_matrix[v_head][u_tail] = d_th;
                }

                // Tail-Tail
                let d_tt = weight_to_dist(w_tt);
                if d_tt < dist_matrix[u_tail][v_tail] {
                    dist_matrix[u_tail][v_tail] = d_tt;
                    dist_matrix[v_tail][u_tail] = d_tt;
                }
                // 辅助函数：权重转距离
                // let calc_dist = |w: f64, node_a: usize, node_b: usize| -> u32 {
                //     if w <= 0.0 { return MAX_DIST; }
                    
                //     // 相对强度归一化
                //     let score_a = w / max_weights[node_a];
                //     let score_b = w / max_weights[node_b];
                    
                //     // 几何平均综合得分
                //     let combined_score = (score_a * score_b).sqrt();
                    
                //     // 转换为距离 (使用平方根拉伸差异)
                //     let dist_val = (1.0 - combined_score.powf(0.5)) * SCALE_FACTOR;
                //     dist_val.max(1.0) as u32
                // };


                // // 填充矩阵 (对称)
                // // Head-Head
                // let d_hh = calc_dist(w_hh, u_head, v_head);
                // if d_hh < dist_matrix[u_head][v_head] {
                //     dist_matrix[u_head][v_head] = d_hh;
                //     dist_matrix[v_head][u_head] = d_hh;
                // }

                // // Head-Tail
                // let d_ht = calc_dist(w_ht, u_head, v_tail);
                // if d_ht < dist_matrix[u_head][v_tail] {
                //     dist_matrix[u_head][v_tail] = d_ht;
                //     dist_matrix[v_tail][u_head] = d_ht;
                // }

                // // Tail-Head
                // let d_th = calc_dist(w_th, u_tail, v_head);
                // if d_th < dist_matrix[u_tail][v_head] {
                //     dist_matrix[u_tail][v_head] = d_th;
                //     dist_matrix[v_head][u_tail] = d_th;
                // }

                // // Tail-Tail
                // let d_tt = calc_dist(w_tt, u_tail, v_tail);
                // if d_tt < dist_matrix[u_tail][v_tail] {
                //     dist_matrix[u_tail][v_tail] = d_tt;
                //     dist_matrix[v_tail][u_tail] = d_tt;
                // }
            }
        }
    }

    // 5. Dummy 节点连接
    let dummy_idx = num_nodes;
    for i in 0..num_nodes {
        dist_matrix[dummy_idx][i] = 0;
        dist_matrix[i][dummy_idx] = 0;
    }

    // 6. 求解 TSP
    let contig_dist_matrix = DistanceMatrix::new(dist_matrix);
    log::info!("Solving Dual-Node TSP...");
    let mut results = contig_dist_matrix.solve(iterations);

    // 7. 解析结果
    // 找到 Dummy 的位置并旋转，使其成为起点/终点
    let dummy_pos = results.iter().position(|&x| x == dummy_idx).unwrap();
    results.rotate_left(dummy_pos);
    
    // 移除 Dummy
    let path: Vec<usize> = results.into_iter().filter(|&x| x != dummy_idx).collect();

    // 验证路径有效性并提取 Contig 顺序和方向
    let mut optimized_order = Vec::with_capacity(num_contigs);
    let mut optimized_signs = Vec::with_capacity(num_contigs);
    let mut visited = vec![false; num_contigs];

    // 遍历路径，每次处理两个节点 (一个 Contig)
    // 理想情况下，路径应该是 (u_node1, u_node2), (v_node1, v_node2)...
    // 其中 node1 和 node2 属于同一个 Contig
    
    let mut i = 0;
    while i < path.len() - 1 {
        let n1 = path[i];
        let n2 = path[i+1];
        
        let c1 = n1 / 2;
        let c2 = n2 / 2;

        if c1 == c2 {
            if !visited[c1] {
                optimized_order.push(idx_to_id[c1]);
                visited[c1] = true;

                if n1 % 2 == 0 {
                    optimized_signs.push(true);
                } else {
                    optimized_signs.push(false);
                }
            }
            i += 2;
        } else {

            if !visited[c1] {
                optimized_order.push(idx_to_id[c1]);
                visited[c1] = true;
      
                optimized_signs.push(true); 
            }
            i += 1;
        }
    }

    if i < path.len() {
        let n = path[i];
        let c = n / 2;
        if !visited[c] {
            optimized_order.push(idx_to_id[c]);
            optimized_signs.push(true);
        }
    }

    log::info!("Dual-Node optimization finished.");

    Tour {
        contigs: optimized_order,
        signs: optimized_signs,
    }
}



#[derive(Clone, Debug)]
pub struct MyMutation<G: EvolveGenotype> {
    _phantom: PhantomData<G>,
    pub initial_mutation_rate: f32,
    pub max_mutation_rate: f32,
    pub growth_step: f32, 
}

impl<G: EvolveGenotype> MyMutation<G> {
    pub fn new(initial_mutation_rate: f32, max_mutation_rate: f32, growth_step: f32) -> Self {
        Self {
            _phantom: PhantomData,
            initial_mutation_rate,
            max_mutation_rate,
            growth_step,
        }
    }
}


impl<G: EvolveGenotype> Mutate for MyMutation<G> {
    type Genotype = G;

    fn call<R: Rng, SR: StrategyReporter<Genotype = G>>(
        &mut self,
        _genotype: &G,
        state: &mut EvolveState<G>,
        _config: &EvolveConfig, 
        _reporter: &mut SR,    
        rng: &mut R,
    ) {
        let now = Instant::now();
        let generation = state.current_generation();
        let growth_step = self.growth_step;

        let current_prob = if growth_step <= 0.0 {
            self.initial_mutation_rate
        } else {
            let steps = generation / 10000; 
            let max_rate = self.max_mutation_rate;
            (self.initial_mutation_rate + (steps as f32 * growth_step))
                .min(max_rate)
        };

        let random_two_int = |rng: &mut R, len: usize| -> (usize, usize) {
            let mut i = rng.gen_range(0..len);
            let mut j = rng.gen_range(0..len);
            while i == j {
                j = rng.gen_range(0..len);
            }
            if i > j {
                std::mem::swap(&mut i, &mut j);
            }
            (i, j)
        };
   

        for chromosome in state.population.chromosomes.iter_mut().filter(|c| c.is_offspring()) {
            if !rng.gen_bool(current_prob as f64) {
                continue;
            }

            let genes = &mut chromosome.genes;
            let len = genes.len();
            if len < 2 { continue; }

            let strategy: f64 = rng.r#gen();
            

            let i = rng.gen_range(0..len);
            let mut j = rng.gen_range(0..len - 1);
            if j >= i { j += 1; } 
            

            let (min_idx, max_idx) = if i < j { (i, j) } else { (j, i) };

            if strategy < 0.2 {
                genes.swap(i, j);
            } else if strategy < 0.4 {
                let k = rng.gen_range(1..len);
                genes.rotate_left(k); 
            // } else if strategy < 0.55 {
            //     if len < 4 { 
            //         genes[min_idx..=max_idx].reverse();
            //         continue; 
            //     }
            //     let k = rng.gen_range(0..len);
            //     let mut cuts = [min_idx, max_idx, k];
            //     cuts.sort_unstable();

            //     if cuts[0] < cuts[1] && cuts[1] < cuts[2] {
            //         let (p1, p2, p3) = (cuts[0], cuts[1], cuts[2]);

            //         let mut new_genes = Vec::with_capacity(len);
            //         new_genes.extend_from_slice(&genes[0..p1]);   // A
            //         new_genes.extend_from_slice(&genes[p3..len]); // D
            //         new_genes.extend_from_slice(&genes[p2..p3]);  // C
            //         new_genes.extend_from_slice(&genes[p1..p2]);  // B
            //         *genes = new_genes;
            //     } else {
            //         genes[min_idx..=max_idx].reverse();
            //     }

            } else if strategy < 0.70 {
                if i < j {
                    genes[i..=j].rotate_left(1);
                } else {
                    genes[j..=i].rotate_right(1);
                }
            // } else if strategy < 0.9 {
            //     let block_len = max_idx - min_idx + 1;
            //     if block_len < len - 1 {
            //         let block: Vec<_> = genes.drain(min_idx..=max_idx).collect();
            //         let k = rng.gen_range(0..=(len - block_len));
            //         genes.splice(k..k, block);
            //     }
            } else {
                genes[min_idx..=max_idx].reverse();
            }
            
            chromosome.fitness_score = None;
        }
      
    }
}


#[derive(Clone, Debug)]
pub struct MyTournament<G: EvolveGenotype> {
    _phantom: PhantomData<G>,
    pub tournament_size: usize,
} 


impl<G: EvolveGenotype> Select for MyTournament<G> {
    type Genotype = G;

    fn call<R: Rng, SR: StrategyReporter<Genotype = G>>(
        &mut self,
        _genotype: &G,
        state: &mut EvolveState<G>,
        config: &EvolveConfig,
        _reporter: &mut SR,
        rng: &mut R,
    ) {
        let now = Instant::now();
        let target_population_size = config.target_population_size;
     
        let mut chromosomes = std::mem::take(&mut state.population.chromosomes);

        self.selection::<R>(
            &mut chromosomes,
            config.target_population_size,
            &mut state.population,
            config,
            rng,
        );

        state.population.chromosomes = chromosomes;

        // state.add_duration(StrategyAction::Select, now.elapsed());
    }
}



impl<G: EvolveGenotype> MyTournament<G> {
    pub fn new(tournament_size: usize) -> Self {
        Self {
            _phantom: PhantomData,
            tournament_size,
        }
    }


    pub fn selection<R: Rng>(
        &self,
        chromosomes: &mut Vec<Chromosome<G::Allele>>,
        selection_size: usize,
        population: &mut Population<G::Allele>,
        config: &EvolveConfig,
        rng: &mut R,
    ) {
        let mut working_population_size = chromosomes.len();
        let tournament_size = std::cmp::min(self.tournament_size, working_population_size);
        let selection_size = std::cmp::min(selection_size, working_population_size);

        let mut selected_chromosomes: Vec<Chromosome<G::Allele>> =
            Vec::with_capacity(selection_size);
        let mut sample_index: usize;
        let mut winning_index: usize;
        let mut sample_fitness_value: FitnessValue;
        let mut winning_fitness_value: FitnessValue;

        match config.fitness_ordering {
            FitnessOrdering::Maximize => {
                for _ in 0..selection_size {
                    winning_index = 0;
                    winning_fitness_value = FitnessValue::MIN;

                    for _ in 0..tournament_size {
                        sample_index = rng.gen_range(0..working_population_size);
                        sample_fitness_value = chromosomes[sample_index]
                            .fitness_score()
                            .unwrap_or(FitnessValue::MIN);

                        if sample_fitness_value >= winning_fitness_value {
                            winning_index = sample_index;
                            winning_fitness_value = sample_fitness_value;
                        }
                    }
                    let chromosome = chromosomes.swap_remove(winning_index);
                    selected_chromosomes.push(chromosome);
                    working_population_size -= 1;
                }
            }
            FitnessOrdering::Minimize => {
                for _ in 0..selection_size {
                    winning_index = 0;
                    winning_fitness_value = FitnessValue::MAX;

                    for _ in 0..tournament_size {
                        sample_index = rng.gen_range(0..working_population_size);
                        sample_fitness_value = chromosomes[sample_index]
                            .fitness_score()
                            .unwrap_or(FitnessValue::MAX);

                        if sample_fitness_value <= winning_fitness_value {
                            winning_index = sample_index;
                            winning_fitness_value = sample_fitness_value;
                        }
                    }
                    let chromosome = chromosomes.swap_remove(winning_index);
                    selected_chromosomes.push(chromosome);
                    working_population_size -= 1;
                }
            }
        };
        // Recycle all losing chromosomes to population's recycling bin
        population.truncate_external(chromosomes, 0);
        chromosomes.append(&mut selected_chromosomes);
    }
}


#[derive(Clone, Debug)]
pub enum Wrapper<G: EvolveGenotype> {
    Elite(SelectElite<G>),
    Tournament(SelectTournament<G>),
    MyTournament(MyTournament<G>),
}

impl<G: EvolveGenotype> Select for Wrapper<G> {
    type Genotype = G;

    fn before(&mut self, genotype: &G, state: &mut EvolveState<G>, config: &EvolveConfig) {
        match self {
            Wrapper::Elite(select) => select.before(genotype, state, config),
            Wrapper::Tournament(select) => select.before(genotype, state, config),
            Wrapper::MyTournament(select) => select.before(genotype, state, config),
        }
    }

    fn call<R: Rng, SR: StrategyReporter<Genotype = G>>(
        &mut self,
        genotype: &G,
        state: &mut EvolveState<G>,
        config: &EvolveConfig,
        reporter: &mut SR,
        rng: &mut R,
    ) {
        match self {
            Wrapper::Elite(select) => select.call(genotype, state, config, reporter, rng),
            Wrapper::Tournament(select) => select.call(genotype, state, config, reporter, rng),
            Wrapper::MyTournament(select) => select.call(genotype, state, config, reporter, rng),
        }
    }

    fn after(&mut self, genotype: &G, state: &mut EvolveState<G>, config: &EvolveConfig) {
        match self {
            Wrapper::Elite(select) => select.after(genotype, state, config),
            Wrapper::Tournament(select) => select.after(genotype, state, config),
            Wrapper::MyTournament(select) => select.after(genotype, state, config),
        }
    }
}

impl<G: EvolveGenotype> From<SelectElite<G>> for Wrapper<G> {
    fn from(select: SelectElite<G>) -> Self {
        Wrapper::Elite(select)
    }
}
impl<G: EvolveGenotype> From<SelectTournament<G>> for Wrapper<G> {
    fn from(select: SelectTournament<G>) -> Self {
        Wrapper::Tournament(select)
    }
}


impl<G: EvolveGenotype> From<MyTournament<G>> for Wrapper<G> {
    fn from(select: MyTournament<G>) -> Self {
        Wrapper::MyTournament(select)
    }
}




#[derive(Clone)]
pub struct Simple<G: EvolveGenotype> {
    pub buffer: Option<Vec<u8>>,
    pub period: usize,
    _phantom: PhantomData<G>,
}

impl<G: EvolveGenotype> Default for Simple<G> {
    fn default() -> Self {
        Self {
            buffer: None,
            period: 500,
            _phantom: PhantomData,
        }
    }
}

impl<G: EvolveGenotype> Simple<G> {
    pub fn new(period: usize) -> Self {
        Self {
            period,
            ..Default::default()
        }
    }

    pub fn new_with_buffer(period: usize) -> Self {
        Self {
            buffer: Some(Vec::new()),
            period,
            ..Default::default()
        }
    }

    fn writeln(&mut self, args: Arguments<'_>) {
        if let Some(buffer) = self.buffer.as_mut() {
            buffer.write_fmt(args).unwrap_or(());
            writeln!(buffer).unwrap_or(())
        } else {
            std::io::stdout().write_fmt(args).unwrap_or(());
            println!()
        }
    }
}
impl<G: EvolveGenotype> StrategyReporter for Simple<G> {
    type Genotype = G;

    fn flush(&mut self, output: &mut Vec<u8>) {
        if let Some(buffer) = self.buffer.as_mut() {
            output.append(buffer);
        }
    }

    fn on_enter<S: StrategyState<Self::Genotype>, C: StrategyConfig>(
        &mut self,
        genotype: &Self::Genotype,
        state: &S,
        config: &C,
    ) {
      
    }

    fn on_selection_complete<S: StrategyState<Self::Genotype>, C: StrategyConfig>(
        &mut self,
        genotype: &Self::Genotype,
        state: &S,
        config: &C,
    ) {
        if state.current_generation() % self.period == 0 {
            log::info!("Current generation {}; best fitness score: {:.6?}",
                state.current_generation(),
                state.best_fitness_score().unwrap()
            );

        
        }
    }


}


#[derive(Clone, Debug)]
pub struct MyOX1Crossover {
    pub mutation_rate: f64,
}

impl MyOX1Crossover {
    pub fn new(mutation_rate: f64) -> Self {
        Self { mutation_rate }
    }
}


impl Crossover for MyOX1Crossover {
    type Genotype = UniqueGenotype;

    fn call<R: Rng, SR: StrategyReporter<Genotype = Self::Genotype>>(
        &mut self,
        genotype: &Self::Genotype,
        state: &mut EvolveState<Self::Genotype>,
        _config: &EvolveConfig,
        _reporter: &mut SR,
        rng: &mut R,
    ) {

        let mut parent_indices: Vec<usize> = (0..state.population.size()).collect();

        let pop_size = state.population.size();
        let gene_size = genotype.genes_size();
        let mut new_chromosomes = Vec::with_capacity(pop_size);

        parent_indices.shuffle(rng);

        for chunk in parent_indices.chunks(2) {
            if chunk.len() < 2 {

                let p1 = &state.population.chromosomes[chunk[0]];
                new_chromosomes.push(p1.clone());
                continue;
            }

            if rng.r#gen::<f64>() > self.mutation_rate {
                let p1 = &state.population.chromosomes[chunk[0]];
                let p2 = &state.population.chromosomes[chunk[1]];
                new_chromosomes.push(p1.clone());
                new_chromosomes.push(p2.clone());
                continue;
            }

            let p1_idx = chunk[0];
            let p2_idx = chunk[1];
            let p1 = &state.population.chromosomes[p1_idx];
            let p2 = &state.population.chromosomes[p2_idx];


            let n = p1.genes.len();
 
            let mut i = rng.gen_range(0..n);
            let mut j = rng.gen_range(0..n);
            if i > j { std::mem::swap(&mut i, &mut j); }

            if i == j || (i == 0 && j == n - 1) {
                new_chromosomes.push(p1.clone());
                new_chromosomes.push(p2.clone());
                continue;
            }
            let create_child = |parent_a: &Vec<usize>, parent_b: &Vec<usize>| -> Vec<usize> {
                let mut child = vec![usize::MAX; n];
                let mut used = vec![false; n]; 
                for k in i..=j {
                    let gene = parent_a[k];
                    child[k] = gene;
                    if gene < n { used[gene] = true; }
                }

               
                let mut child_idx = (j + 1) % n;
                for k in 0..n {
                    let b_idx = (j + 1 + k) % n;
                    let gene = parent_b[b_idx];
                    if gene < n && !used[gene] {
                        child[child_idx] = gene;
                        used[gene] = true;
                        child_idx = (child_idx + 1) % n;
                    }
                }
                child
            };

            let c1_genes = create_child(&p1.genes, &p2.genes);
            let c2_genes = create_child(&p2.genes, &p1.genes);

            let mut c1 = p1.clone(); 
            c1.genes = c1_genes;
            c1.fitness_score = None; 

            let mut c2 = p2.clone();
            c2.genes = c2_genes;
            c2.fitness_score = None;

            new_chromosomes.push(c1);
            new_chromosomes.push(c2);

        }


        state.population.chromosomes = new_chromosomes;
    }
}


pub fn run_evolove_optimizer(
    tour: &Tour,
    contigsizes: IndexMap<usize, usize>,
    contacts: Vec<HashMap<usize, u32>>, 
    mutation_rate: f64,
    population_size: usize, 
    generations: usize,
    round: usize,
    resume: bool,
    seed: u64,
) -> Tour {

    
    let mut sizes: Vec<usize> = Vec::with_capacity(contigsizes.len());
    let mut initial_genes: Vec<usize> = Vec::with_capacity(contigsizes.len());
    for i in tour.clone().into_iter() {
        sizes.push(*contigsizes.get(&i).unwrap());
        initial_genes.push(i);
    }

    let matrix = ContactMatrix::new(&contigsizes, contacts);
    let num_contigs = matrix.num_contigs;
    
    let mut best_genes = initial_genes.clone();
    let mut best_fitness_score = 0.0;
    if !resume {
        initial_genes.shuffle(&mut SmallRng::seed_from_u64(seed));
    }
    best_fitness_score = calculate_fitness(&initial_genes, &matrix) as f64;
    for round_idx in 0..round {
        if round > 1{
            log::info!("Starting optimization round {} with score: {:.2}", round_idx + 1, best_fitness_score);
        } else {
            log::info!("Genetic Algorithm Optimization started.");
            log::info!("Generation-0, score: {:.2}", best_fitness_score);
        }
        let fitness = MyFitness { matrix: &matrix };
        let mut genotype = UniqueGenotype::builder()
            .with_allele_list((0..num_contigs).collect())
            .with_genes_hashing(false)
            .with_chromosome_recycling(false)
            .build()
            .unwrap();
  
        genotype.set_seed_genes_list(vec![initial_genes]);
        // let mut seed_population = vec![initial_genes.clone()];
        
        // let num_perturbed = (population_size as f64 * 0.5) as usize;
        // let mut rng = SmallRng::seed_from_u64(seed + round_idx as u64); // 确保每轮扰动不同

        // for _ in 0..num_perturbed {
        //     let mut perturbed_genes = initial_genes.clone();
        //     let len = perturbed_genes.len();
        //     if len > 2 {
        //         let num_swaps = (len as f64 * 0.1).max(5.0) as usize; // 5% 的差异
        //         for _ in 0..num_swaps {
        //             let i = rng.gen_range(0..len);
        //             let j = rng.gen_range(0..len);
        //             perturbed_genes.swap(i, j);
        //         }

        //         let block_start = rng.gen_range(0..len-1);
        //         let block_end = rng.gen_range(block_start+1..len).min(block_start + len/10);
        //         perturbed_genes[block_start..block_end].reverse();
        //     }
        //     seed_population.push(perturbed_genes);
        // }


        // genotype.set_seed_genes_list(seed_population);
    
    
        let mut evolve = Evolve::builder()
                .with_genotype(genotype)
                .with_target_population_size(population_size)
                .with_max_stale_generations(generations) 
                .with_max_generations(1000_000_000)      
                .with_par_fitness(true)    
                .with_target_fitness_score(isize::MIN) 
                .with_fitness_ordering(FitnessOrdering::Minimize)   
                .with_crossover(CrossoverClone::new(1.0))
                // .with_crossover(DirectlyClone::new())
                .with_mutate(MyMutation::new(mutation_rate as f32, 0.9, 0.0))
                .with_fitness(fitness)
                .with_select(MyTournament::new(3))
                .with_reporter(Simple::new(500))
                .with_rng_seed_from_u64(seed)   
                // .with_replace_on_equal_fitness(true)
                // .with_select(SelectTournament::new(0.4, 0.02, 3)) 
                // .with_extension(ExtensionMassExtinction::new(10, 0.05, 0.2))
                // .with_max_chromosome_age(10) 
                // .with_crossover(MyOX1Crossover::new(0.5))
                // .with_valid_fitness_score(10)             
                // .with_reporter(EvolveReporterSimple::new(500)) 
                // .with_target_fitness_score(isize::MIN) 
                .build()
                .unwrap();

        evolve.call();

        let current_generation = evolve.state.current_generation();
        
        let best_chromosome = evolve.best_chromosome().unwrap();
        best_genes = best_chromosome.genes;

        best_fitness_score = best_chromosome.fitness_score.unwrap() as f64;
        if round > 1{
            log::info!("Round {} finished. Best Score: {:.2}", round_idx + 1, best_fitness_score);
        } else {
            log::info!("Generation-{}, score: {:.2}", current_generation, best_fitness_score);

        }
        initial_genes = best_genes.clone();
    }

    log::info!("Optimization finished. Best Score: {:.4}", best_fitness_score);

    Tour {
        contigs: best_genes.clone(),
        signs: vec![true; num_contigs],
    }
}


pub fn run_hill_climbing(
    contigsizes: IndexMap<usize, usize>,
    contacts: Vec<HashMap<usize, u32>>, 
    mutation_rate: f64,
    population_size: usize, 
    generations: usize,
    round: usize,
    resume: bool,
    seed: u64,
) -> Tour {

   
    let mut initial_genes: Vec<usize> = Vec::with_capacity(contigsizes.len());
    for i in contigsizes.keys() {
        
        initial_genes.push(*i);
    }

    let matrix = ContactMatrix::new(&contigsizes, contacts);
    let num_contigs = matrix.num_contigs;
    
    
    let mut best_genes = initial_genes.clone();
    let mut best_fitness_score = 0.0;
    if !resume {
        initial_genes.shuffle(&mut SmallRng::seed_from_u64(seed));
    }
    best_fitness_score = calculate_fitness(&initial_genes, &matrix) as f64;
    
    for round_idx in 0..round {
        if round > 1 {
            log::info!("Starting optimization round {} with score: {:.2}", round_idx + 1, best_fitness_score);
        } else {
            log::info!("Hill Climbing Optimization started.");
            log::info!("Generation-0, score: {:.2}", best_fitness_score);
        }
        let fitness = MyFitness { matrix: &matrix };
        let mut genotype = UniqueGenotype::builder()
            .with_allele_list((0..num_contigs).collect())
            .with_genes_hashing(false)
            .with_chromosome_recycling(true)
            .build()
            .unwrap();

        genotype.set_seed_genes_list(vec![initial_genes]);

        let mut hill_climb_builder = HillClimb::builder()
            .with_rng_seed_from_u64(seed)
            .with_genotype(genotype)
            .with_fitness(fitness)
            .with_par_fitness(true)
            .with_variant(HillClimbVariant::Stochastic)
            .with_fitness_ordering(FitnessOrdering::Minimize)
            // .with_reporter(HillClimbReporterSimple::new(500))
            .with_max_stale_generations(generations)
            .with_max_generations(500_000_000) 
            .build()
            .unwrap();

        hill_climb_builder.call();

        let best_chromosome = hill_climb_builder.best_chromosome().unwrap();
        best_genes = best_chromosome.genes;
        best_fitness_score = best_chromosome.fitness_score.unwrap() as f64;
        if round > 1 {
            log::info!("Round {} finished. Best Score: {:.2}", round_idx + 1, best_fitness_score);
        } else {
            let current_generation = hill_climb_builder.state.current_generation();
            log::info!("Generation-{}, score: {:.2}", current_generation, best_fitness_score);
        }
        initial_genes = best_genes.clone();
    }

   
    log::info!("Optimization finished. Best Score: {:.4}", best_fitness_score);

    Tour {
        contigs: best_genes.clone(),
        signs: vec![true; num_contigs],
    }

}


pub fn run_hybrid(
    tour: &Tour,
    contigsizes: IndexMap<usize, usize>,
    contacts: Vec<HashMap<usize, u32>>, 
    mutation_rate: f64,
    population_size: usize, 
    generations: usize,
    round: usize,
    is_lkh: bool,
    is_ga: bool,
    resume: bool,
    seed: u64,
    threads: usize,
) -> Tour {
 
    let matrix = ContactMatrix::new(
        &contigsizes,
        contacts.clone(),
    );

    let best_tour = tour.clone();

    
        
    let best_tour = if is_lkh {
        log::info!("Ordering score before optimization: {:?}", calculate_fitness(&best_tour.contigs, &matrix));
        let tour = run_lkh_optimizer(
            &best_tour,
            contigsizes.clone(), 
            &matrix,
            10,
            seed,
        );
        log::info!("Ordering score after LKH optimization: {:?}", calculate_fitness(&tour.contigs, &matrix));
        tour
    } else {
        best_tour.clone()
        
    };

    // log::info!("Initial score before optimization: {:?}", calculate_fitness(&initial_contigs, &matrix));
    // two_opt_optimization(
    //     &mut initial_contigs,
    //     &matrix,
    //     num_contigs,
    // );
    
    // block_move_optimization(
    //     &mut initial_contigs,
    //     &matrix,
    //     100,
    // );
    // two_opt_optimization(
    //     &mut initial_contigs,
    //     &matrix,
    //     num_contigs,
    // );
    // let current_best_fitness = calculate_fitness(&initial_contigs, &matrix);
    // log::info!("Initial score after 2-Opt optimization: {:?}", current_best_fitness);
    // let mut new_contigsizes = IndexMap::new();
    // for i in initial_contigs.into_iter() {
    //     new_contigsizes.insert(i, *contigsizes.get(&i).unwrap());
    // }

    // let best_genes = best_tour.contigs;

    // let new_contigsizes = IndexMap::from_iter(
    //     best_genes.iter().map(|&i| (i, *contigsizes.get(&i).unwrap()))
    // );
    // let best_tour = run_hill_climbing(
    //     new_contigsizes,
    //     contacts.clone(),
    //     mutation_rate,
    //     population_size,
    //     1000,
    //     1,
    //     true,
    //     seed,
    // ); 


   
    let mut best_tour = if is_ga { 
        let tour = run_evolove_optimizer(
            &best_tour,
            contigsizes.clone(),
            contacts.clone(),
            mutation_rate,
            population_size,
            generations,
            1,
            true,
            seed,
        );


        let new_mutation_rate = if mutation_rate * 3.0 < 1.0 {
            mutation_rate * 3.0
        } else {
            1.0
        };
        let tour = run_evolove_optimizer(
                                &tour,
                                contigsizes.clone(),
                                contacts.clone(),
                                new_mutation_rate,
                                population_size,
                                generations,
                                1,
                                true,
                                seed,
                            );
        tour    
    } else {
        best_tour.clone()
    };

   

    log::info!("Final fitness: {}", calculate_fitness(&best_tour.contigs, &matrix));
    
    best_tour


}


