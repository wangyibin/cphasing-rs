// generate the minimizer sketch of a sequence
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use nthash::ntc64;

#[derive(Debug, Clone)]
pub struct MinimizerInfo {
    pub rid: u64,
    pub pos: u64,
    // pub rev: u32,
    // pub span: u32,
}

#[derive(Debug, Clone)]
pub struct MinimierData {
    pub minimizer: u64,
    pub info: MinimizerInfo,
}


// rolling hash algorithm
pub fn hash (seq: &str, k: usize) -> u64 {
    let mut h = 0;
    for i in 0..k {
        
        h = h << 2;
        h += match seq.chars().nth(i).unwrap() {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => 0,
        };
    }
    h
}

fn update_hash (h: u64, k: usize, c: char) -> u64 {
    let mut h = h;
    h = h << 2;
    h += match c {
        'A' => 0,
        'C' => 1,
        'G' => 2,
        'T' => 3,
        _ => 0,
    };
    h & ((1 << (2 * k)) - 1)
}

fn minimizer (seq: &str, rid: &u64, start: u64, k: usize) -> MinimierData {
    let mut h = hash(seq, k);
    let mut m = h;
    let mut pos = 0;
    for i in 1..(seq.len() - k + 1) {
        h = update_hash(h, k, seq.chars().nth(i + k - 1).unwrap());
        if h < m {
            m = h;
            pos = i as u64 ;
        }
    }

    let pos = start + pos;
    let m = MinimierData {
        minimizer: m,
        info: MinimizerInfo {
            rid: *rid,
            pos: pos,
            // rev: 0,
        }
    };
    println!("{}, {:?}", seq, m);
    m 
}

pub fn minimizer_nthash(seq: &[u8], rid: &u64, start: u64, k: usize) -> MinimierData {
    println!("seq: {:?}", seq);
    let mut h = ntc64(seq, 0, k);
    let mut m = h;
    let mut pos = 0;
    for i in 1..(seq.len() - k + 1) {
        h = ntc64(&seq[i..(i+k)], 0, k);
        if h < m {
            m = h;
            pos = i as u64 ;
        }
    }

    let pos = start + pos;
    let m = MinimierData {
        minimizer: m,
        info: MinimizerInfo {
            rid: *rid,
            pos: pos,
            // rev: 0,
        }
    };
    m 
}

pub fn sketch (seq: &String, rid: u64, k: usize, w: usize) -> Vec<MinimierData> {
    let mut sketch = Vec::new();
    let seq = seq.as_bytes();
    println!("seq len: {}", seq.len());
    for i in (0..seq.len() - k + w + 1 ).step_by(k + w - 1) {
        let m: MinimierData = if i + k + w > seq.len() {
            if seq.len() - i < k {
                break;
            }
            minimizer_nthash(&seq[i..], &rid, i.try_into().unwrap(), k)
        } else {
            minimizer_nthash(&seq[i..(i + k + w - 1)], &rid, i.try_into().unwrap(), k)
        };
        sketch.push(m);

    }
    sketch
}
