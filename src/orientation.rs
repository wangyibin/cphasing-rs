#![allow(unused)]
#![allow(unused_imports)]
use hashbrown::{HashMap};
use rand::Rng;
use rayon::prelude::*;
use super::order::Tour;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Orientation {
    Forward, // +
    Reverse, // -
}


#[derive(Debug, Clone)]
pub struct OrientedContig {
    pub id: usize,
    pub orientation: Orientation,
}

// pub fn orient_contigs(tour: &mut Tour<usize>, matrix: &DetailedContactMatrix) -> Vec<bool> {
//     let n = tour.contigs.len();
//     if n < 2 {
//         return vec![true; n];
//     }

//     let mut dp = vec![vec![0.0; 2]; n];
//     let mut path = vec![vec![0; 2]; n]; // 0: from +, 1: from -


//     dp[0][0] = 0.0;
//     dp[0][1] = 0.0;

//     for i in 1..n {
//         let u = tour.contigs[i-1];
//         let v = tour.contigs[i];

//         let score_pp = matrix.get(u, v, 1, 0);

//         let score_mp = matrix.get(u, v, 0, 0);


//         let score_pm = matrix.get(u, v, 1, 1);


//         let score_mm = matrix.get(u, v, 0, 1);

//         if dp[i-1][0] + score_pp >= dp[i-1][1] + score_mp {
//             dp[i][0] = dp[i-1][0] + score_pp;
//             path[i][0] = 0; 
//         } else {
//             dp[i][0] = dp[i-1][1] + score_mp;
//             path[i][0] = 1;
//         }

//         if dp[i-1][0] + score_pm >= dp[i-1][1] + score_mm {
//             dp[i][1] = dp[i-1][0] + score_pm;
//             path[i][1] = 0; 
//         } else {
//             dp[i][1] = dp[i-1][1] + score_mm;
//             path[i][1] = 1; 
//         }
//     }

//     let mut signs = vec![true; n];
//     let mut curr_dir = if dp[n-1][0] >= dp[n-1][1] { 0 } else { 1 };
    
//     for i in (0..n).rev() {
//         signs[i] = curr_dir == 0;
//         curr_dir = path[i][curr_dir];
//     }

//     tour.signs = signs.clone();
//     signs
// }

pub fn orient_contigs_sa(
    tour: &mut Tour<usize>, 
    matrix: &DetailedContactMatrix,
    max_iter: usize,
    start_temp: f64,
    cooling_rate: f64
) {
    let n = tour.contigs.len();
    let mut rng = rand::thread_rng();

    let mut current_score = calculate_all_pairs_score(tour, matrix);
    let mut best_score = current_score;
    let mut best_signs = tour.signs.clone();
    
    let mut temp = start_temp;

    let mut no_improve_count = 0;
    let max_no_improve = max_iter / 10; 

    println!("[SA] Start Score: {:.4}, Temp: {:.4}", current_score, temp);

    for i in 0..max_iter {

        let idx = rng.gen_range(0..n);

        let delta = calculate_flip_delta(tour, matrix, idx);

        let accept = if delta > 0.0 {
            true
        } else {
            let probability = (delta / temp).exp();
            rng.r#gen::<f64>() < probability
        };

        if accept {

            tour.signs[idx] = !tour.signs[idx];
            current_score += delta;

            if current_score > best_score {
                best_score = current_score;
                best_signs = tour.signs.clone();
                no_improve_count = 0;
            } else {
                no_improve_count += 1;
            }
        } else {
            no_improve_count += 1;
        }
        

        temp *= cooling_rate;

        if i % 5000 == 0 {
            // log::info!("Iter {}: Temp {:.4}, Score {:.2}, Best {:.2}", i, temp, current_score, best_score);
        }

        if temp < 1e-6 || no_improve_count > max_no_improve {
            break;
        }
    }
    
    println!("[SA] Final Best Score: {:.4}", best_score);

    tour.signs = best_signs;
}

fn calculate_flip_delta(tour: &Tour<usize>, matrix: &DetailedContactMatrix, idx: usize) -> f64 {
    let n = tour.contigs.len();
    let u = tour.contigs[idx];
    let u_sign_old = tour.signs[idx];
    let u_sign_new = !u_sign_old; 
    let mut delta = 0.0;


    let window_limit = 100; 
    let start = if window_limit > 0 && idx > window_limit { idx - window_limit } else { 0 };
    let end = if window_limit > 0 && idx + window_limit < n { idx + window_limit } else { n };

    for j in start..end {
        if idx == j { continue; }
        
        let v = tour.contigs[j];
        let v_sign = tour.signs[j];
        

        let v_end_for_u: usize;

        let score_old;
        let score_new;

        if idx < j {

            let u_end_old = if u_sign_old { 1 } else { 0 };
            let u_end_new = if u_sign_new { 1 } else { 0 };

            let v_end = if v_sign { 0 } else { 1 };
            
            score_old = matrix.get(u, v, u_end_old, v_end);
            score_new = matrix.get(u, v, u_end_new, v_end);
        } else {
            let v_end = if v_sign { 1 } else { 0 };
            

            let u_end_old = if u_sign_old { 0 } else { 1 };
            let u_end_new = if u_sign_new { 0 } else { 1 };
            
            score_old = matrix.get(v, u, v_end, u_end_old);
            score_new = matrix.get(v, u, v_end, u_end_new);
        }

        delta += score_new - score_old;
    }

    delta
}


pub fn robust_orient_contigs(
    tour: &mut Tour<usize>, 
    matrix: &DetailedContactMatrix, 
    iterations: usize
) {
    let n = tour.contigs.len();
    let mut forward_counts = vec![0; n];
    let mut rng = rand::thread_rng();


    // for _ in 0..iterations {
    //     let mut temp_tour = tour.clone();
    //     let signs = orient_contigs_with_noise(&mut temp_tour, matrix, &mut rng);
  
    //     for (i, &is_forward) in signs.iter().enumerate() {
    //         if is_forward {
    //             forward_counts[i] += 1;
    //         }
    //     }
    // }

    // let mut confidence = Vec::with_capacity(n);
    // let mut final_signs = Vec::with_capacity(n);

    // for &count in &forward_counts {
    //     let prob = count as f64 / iterations as f64;
    //     confidence.push(prob);
    //     final_signs.push(prob >= 0.5);
    // }

    let all_signs: Vec<Vec<bool>> = (0..iterations).into_par_iter().map(|_| {
        let mut rng = rand::thread_rng();
        let mut temp_tour = tour.clone();
        orient_contigs_with_noise(&mut temp_tour, matrix, &mut rng)
    }).collect();

    let mut forward_counts = vec![0; n];
    for signs in all_signs {
        for (i, &is_forward) in signs.iter().enumerate() {
            if is_forward {
                forward_counts[i] += 1;
            }
        }
    }

    let threshold = iterations / 2;
    for i in 0..n {
        tour.signs[i] = forward_counts[i] > threshold;
    }
    

    let final_signs = tour.signs.clone();
    let forward_score = calculate_orientation_score(tour, matrix);
    
    let rev_signs = tour.signs.iter().map(|&b| !b).collect::<Vec<bool>>();
    tour.signs = rev_signs;

    let reverse_score = calculate_orientation_score(tour, matrix);
    if reverse_score > forward_score {
       
        log::info!("Orientation flipped. Forward score: {}, Reverse score: {}", forward_score, reverse_score);
       

    } else {
        log::info!("Orientation kept. Forward score: {}, Reverse score: {}", forward_score, reverse_score);
        tour.signs = final_signs;
        
    }

}

fn orient_contigs_with_noise(
    tour: &mut Tour<usize>, 
    matrix: &DetailedContactMatrix,
    rng: &mut impl Rng
) -> Vec<bool> {
    let n = tour.contigs.len();
    let mut dp = vec![vec![0.0; 2]; n];
    let mut path = vec![vec![0; 2]; n];

    // Noise factor to control the amount of randomness
    let noise_factor = 0.05; 
    let consistency_bias = 1.0 ; 
    dp[0][0] = 0.0;
    dp[0][1] = 0.0;

    for i in 1..n {
        let u = tour.contigs[i-1];
        let v = tour.contigs[i];


        let raw_pp = matrix.get(u, v, 1, 0);
        let raw_mp = matrix.get(u, v, 0, 0);
        let raw_pm = matrix.get(u, v, 1, 1);
        let raw_mm = matrix.get(u, v, 0, 1);

        let mut add_noise = |score: f64| -> f64 {
            // if score <= 1e-6 { return 0.0; }
            let noise = rng.gen_range(-noise_factor..noise_factor);
            score * (1.0 + noise)
        };

        let score_pp = add_noise(raw_pp)  * consistency_bias;
        let score_mp = add_noise(raw_mp);
        let score_pm = add_noise(raw_pm);
        let score_mm = add_noise(raw_mm) * consistency_bias;
        // println!("Scores with noise for contigs {}-{}: pp: {}, mp: {}, pm: {}, mm: {}", u, v, score_pp, score_mp, score_pm, score_mm);
        if dp[i-1][0] + score_pp >= dp[i-1][1] + score_mp {
            dp[i][0] = dp[i-1][0] + score_pp;
            path[i][0] = 0;
        } else {
            dp[i][0] = dp[i-1][1] + score_mp;
            path[i][0] = 1;
        }

        if dp[i-1][0] + score_pm >= dp[i-1][1] + score_mm {
            dp[i][1] = dp[i-1][0] + score_pm;
            path[i][1] = 0;
        } else {
            dp[i][1] = dp[i-1][1] + score_mm;
            path[i][1] = 1;
        }

    }


    let mut signs = vec![true; n];
    let mut curr_dir = if dp[n-1][0] >= dp[n-1][1] { 0 } else { 1 };
  
    for i in (0..n).rev() {
        signs[i] = curr_dir == 0;
        curr_dir = path[i][curr_dir];
    }
    signs
}



pub fn optimize_orientations_flip_one(tour: &mut Tour<usize>, matrix: &DetailedContactMatrix) {
    let n = tour.contigs.len();
    if n < 2 { return; }

    let mut improved = true;
    let mut pass = 0;
    let max_passes = 20; 

    while improved && pass < max_passes {
        improved = false;
        pass += 1;

        for i in 0..n {
            let score_before = calculate_all_pairs_score(tour, matrix);

            tour.signs[i] = !tour.signs[i];

            let score_after = calculate_all_pairs_score(tour, matrix);

            if score_after > score_before {
                improved = true;
            } else {
                tour.signs[i] = !tour.signs[i];
            }
        }
    }
}


fn get_local_contribution(tour: &Tour<usize>, matrix: &DetailedContactMatrix, i: usize) -> f64 {
    let n = tour.contigs.len();
    let mut score = 0.0;

    if i > 0 {
        let u = tour.contigs[i - 1];
        let v = tour.contigs[i];
        
        let u_sign = tour.signs[i - 1];
        let v_sign = tour.signs[i];


        let u_end = if u_sign { 1 } else { 0 };
        // v(+) -> Head(0), v(-) -> Tail(1)
        let v_end = if v_sign { 0 } else { 1 };

        score += matrix.get(u, v, u_end, v_end);
    }


    if i < n - 1 {
        let u = tour.contigs[i];
        let v = tour.contigs[i + 1];
        
        let u_sign = tour.signs[i];
        let v_sign = tour.signs[i + 1];

        let u_end = if u_sign { 1 } else { 0 };
        let v_end = if v_sign { 0 } else { 1 };

        score += matrix.get(u, v, u_end, v_end);
    }

    score
}

#[derive(Debug, Clone)]
pub struct DetailedContactMatrix {
    // key: (u, v), value: (head_head, head_tail, tail_head, tail_tail)
    pub data: HashMap<(usize, usize), (f64, f64, f64, f64)>, 
}

impl DetailedContactMatrix {

    pub fn new(
        data: HashMap<(usize, usize), (f64, f64, f64, f64)>,
    ) -> Self {
        Self { data }
    }

    pub fn get(&self, u: usize, v: usize, u_end: usize, v_end: usize) -> f64 {
        let (k1, k2) = if u < v { (u, v) } else { (v, u) };
        
        let entry = self.data.get(&(k1, k2)).unwrap_or(&(0.0, 0.0, 0.0, 0.0));

        if u < v {
            match (u_end, v_end) {
                (0, 0) => entry.0, // u_Head - v_Head
                (0, 1) => entry.1, // u_Head - v_Tail
                (1, 0) => entry.2, // u_Tail - v_Head
                (1, 1) => entry.3, // u_Tail - v_Tail
                _ => 0.0,
            }
        } else {

            match (v_end, u_end) {
                (0, 0) => entry.0, // v_Head - u_Head
                (0, 1) => entry.1, // v_Head - u_Tail
                (1, 0) => entry.2, // v_Tail - u_Head
                (1, 1) => entry.3, // v_Tail - u_Tail
                _ => 0.0,

            }
        }
    }


}



pub fn calculate_orientation_score(tour: &Tour<usize>, matrix: &DetailedContactMatrix) -> f64 {
    tour.contigs
        .windows(2)
        .zip(tour.signs.windows(2))
        .map(|(contigs, signs)| {
            let u = contigs[0];
            let v = contigs[1];
            let u_sign = signs[0];
            let v_sign = signs[1];

            let u_end = if u_sign { 1 } else { 0 };
            let v_end = if v_sign { 0 } else { 1 };

           
            // println!("Calculating score for contigs {} (sign {}) and {} (sign {}): u_end {}, v_end {}", u, u_sign, v, v_sign, u_end, v_end);
            // println!("Score: {}", matrix.get(u, v, u_end, v_end));
            // println!("Matrix entry: {:?}", matrix.data.get(&(u.min(v), u.max(v))));
            matrix.get(u, v, u_end, v_end)

        })
        .sum()
}


pub fn calculate_all_pairs_score(tour: &Tour<usize>, matrix: &DetailedContactMatrix) -> f64 {
    let n = tour.contigs.len();
    let mut total_score = 0.0;

    for i in 0..n {
        let u = tour.contigs[i];
        let u_sign = tour.signs[i];

        let u_end = if u_sign { 1 } else { 0 };

        for j in (i + 1)..n {
            let v = tour.contigs[j];
            let v_sign = tour.signs[j];

            let v_end = if v_sign { 0 } else { 1 };

            total_score += matrix.get(u, v, u_end, v_end);
        }
    }

    total_score
}

pub fn calculate_adjacent_score(tour: &Tour<usize>, matrix: &DetailedContactMatrix) -> f64 {
    tour.contigs
        .windows(2)
        .zip(tour.signs.windows(2))
        .map(|(contigs, signs)| {
            let u = contigs[0];
            let v = contigs[1];
            let u_sign = signs[0];
            let v_sign = signs[1];

            let u_end = if u_sign { 1 } else { 0 };
            let v_end = if v_sign { 0 } else { 1 };

            matrix.get(u, v, u_end, v_end)
        })
        .sum()
}
