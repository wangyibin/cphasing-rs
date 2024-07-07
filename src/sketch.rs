// generate the minimizer sketch of a sequence
use anyhow::Result as AnyResult;
use std::cmp::Ordering;
use std::collections::{ HashMap, HashSet };
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};


// modified from nthash 
const MAXIMUM_K_SIZE: usize = u32::max_value() as usize;


const H_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'C' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'G' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'T' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'N' as usize] = 0;
    lookup
};

const RC_LOOKUP: [u64; 256] = {
    let mut lookup = [1; 256];
    lookup[b'A' as usize] = 0x2955_49f5_4be2_4456;
    lookup[b'C' as usize] = 0x2032_3ed0_8257_2324;
    lookup[b'G' as usize] = 0x3193_c185_62a0_2b4c;
    lookup[b'T' as usize] = 0x3c8b_fbb3_95c6_0474;
    lookup[b'N' as usize] = 0;
    lookup
};

#[inline(always)]
fn h(c: u8) -> u64 {
    let val = H_LOOKUP[c as usize];
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", c as char)
    }
    val
}

#[inline(always)]
fn rc(nt: u8) -> u64 {
    let val = RC_LOOKUP[nt as usize];
    if val == 1 {
        panic!("Non-ACGTN nucleotide encountered! {}", nt as char)
    }
    val
}

pub struct NtHashIterator<'a> {
    seq: &'a [u8],
    k: usize,
    fh: u64,
    rh: u64,
    current_idx: usize,
    max_idx: usize,
}


impl<'a> NtHashIterator<'a> {
    /// Creates a new NtHashIterator with internal state properly initialized.
    pub fn new(seq: &'a [u8], k: usize) -> AnyResult<NtHashIterator<'a>> {
       
        if k > seq.len() {
            return Err(anyhow::anyhow!("k must be less than or equal to the length of the sequence"));
        }
     
        assert!(k <= MAXIMUM_K_SIZE, "k must be less than or equal to {}", MAXIMUM_K_SIZE);
        let mut fh = 0;
        for (i, v) in seq[0..k].iter().enumerate() {
            fh ^= h(*v).rotate_left((k - i - 1) as u32);
        }

        let mut rh = 0;
        for (i, v) in seq[0..k].iter().rev().enumerate() {
            rh ^= rc(*v).rotate_left((k - i - 1) as u32);
        }

        Ok(NtHashIterator {
            seq,
            k,
            fh,
            rh,
            current_idx: 0,
            max_idx: seq.len() - k + 1,
        })
    }
}

impl<'a> Iterator for NtHashIterator<'a> {
    type Item = (u64, u8);

    fn next(&mut self) -> Option<(u64, u8)> {
        if self.current_idx == self.max_idx {
            return None;
        };

        if self.current_idx != 0 {
            let i = self.current_idx - 1;
            let seqi = self.seq[i];
            let seqk = self.seq[i + self.k];

            self.fh = self.fh.rotate_left(1) ^ h(seqi).rotate_left(self.k as u32) ^ h(seqk);

            self.rh = self.rh.rotate_right(1)
                ^ rc(seqi).rotate_right(1)
                ^ rc(seqk).rotate_left(self.k as u32 - 1);
        }

        self.current_idx += 1;
        if self.rh < self.fh {
            Some((self.rh, 1))
        } else {
            Some((self.fh, 0))
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.max_idx, Some(self.max_idx))
    }
}


#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MinimizerInfo {
    pub rid: u32,
    pub pos: u32,
    pub rev: u8,
    // pub span: u8,
}

impl Ord for MinimizerInfo {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rid.cmp(&other.rid).then(self.rev.cmp(&other.rev))
    }
}

impl PartialOrd for MinimizerInfo {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MinimizerData {
    pub minimizer: u64,
    pub info: MinimizerInfo,
}

impl Ord for MinimizerData {
    fn cmp(&self, other: &Self) -> Ordering {
        self.minimizer.cmp(&other.minimizer)
    }
}

impl PartialOrd for MinimizerData {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&x| match x {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    } ).collect()
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
            _ => 4,
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


fn minimizer (seq: &str, rid: &u32, start: u32, k: usize) -> MinimizerData {
    let mut h = hash(seq, k);
    let mut m = h;
    let mut pos = 0;
    for i in 1..(seq.len() - k + 1) {
        h = update_hash(h, k, seq.chars().nth(i + k - 1).unwrap());
        if h < m {
            m = h;
            pos = i as u32 ;
            
        }
    }

    let pos = start + pos;
    let m = MinimizerData {
        minimizer: m,
        info: MinimizerInfo {
            rid: *rid,
            pos: pos,
            rev: 0,
            // span: k.try_into().unwrap(),
        }
    };
    println!("{}, {:?}", seq, m);
    m 
}

pub fn minimizer_nthash(seq: &[u8], rid: &u32, start: u32, k: usize) -> MinimizerData {
    let seq_len =  seq.len();
    // let mut h = ntc64(seq, 0, k);
    let mut h = 0;
    let mut m = u64::MAX;
    let mut pos = 0;
    let mut rev: u8 = 0;

    let hash_iter = NtHashIterator::new(seq, k).unwrap();
    hash_iter.collect::<Vec<_>>().iter().enumerate().min_by_key(|&(_, x)| x).map(|(i, &x)| {
        m = x.0;
        rev = x.1;
        pos = i as u32;

    });
    

    let pos = start + pos;
    let m = MinimizerData {
        minimizer: m,
        info: MinimizerInfo {
            rid: *rid,
            pos: pos,
            rev: rev,
            // span: span,
        }
    };
    m 
}

pub fn sketch (seq: &Vec<u8>, rid: u32, k: usize, w: usize) -> Vec<MinimizerData> {
    let mut sketch = Vec::new();
    
    let seq = seq.iter()
                    .map(|&c| c.to_ascii_uppercase())
                    .map(|c| if c == b'A' || c == b'C' || c == b'G' || c == b'T' { c } else { b'N' }).collect::<Vec<u8>>();

    let seq_len = seq.len();

    let mut i: usize = 0;
    let seq_len_i32: i32 = seq_len.try_into().unwrap(); // 提前进行类型转换
    let end_cond = seq_len_i32 - k as i32;

    let mut exists_pos = HashSet::new();
    let count = 0;
    
    while i < end_cond as usize {

        let end_index = std::cmp::min(i + k + w - 1, seq_len);
        let m: MinimizerData = minimizer_nthash(&seq[i..end_index], &rid, i as u32, k);

        if exists_pos.insert(m.info.pos) {
            sketch.push(m);
        }

        i += 1;
    }

    sketch
}
