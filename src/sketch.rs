#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
use anyhow::Result as AnyResult;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
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

// #[inline(always)]
fn h(c: u8) -> u64 {
    unsafe { *H_LOOKUP.get_unchecked(c as usize) }
}

// #[inline(always)]
fn rc(nt: u8) -> u64 {
    unsafe { *RC_LOOKUP.get_unchecked(nt as usize) }
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
            return Err(anyhow::anyhow!(
                "k must be less than or equal to the length of the sequence"
            ));
        }

        assert!(
            k <= MAXIMUM_K_SIZE,
            "k must be less than or equal to {}",
            MAXIMUM_K_SIZE
        );
        let mut fh = 0;
        let mut rh = 0;

        for ((i, &v), &rv) in seq[0..k].iter().enumerate().zip(seq[0..k].iter().rev()) {
            let shift = (k - i - 1) as u32;
            fh ^= h(v).rotate_left(shift);
            rh ^= rc(rv).rotate_left(shift);
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

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct MinimizerInfo {
    pub pos: u64,
    pub rid: u32,
    pub rev: u8,
    pub span: u8,
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

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
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

/// Compact minimizer representation for contigs whose coordinates fit in 32
/// bits. The low bit of `rid_rev` stores the strand and the remaining bits
/// store the record id.
#[repr(C)]
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(crate) struct MinimizerData32 {
    minimizer: u64,
    pos: u32,
    rid_rev: u32,
}

impl MinimizerData32 {
    const MAX_RID: u32 = (1_u32 << 31) - 1;

    #[inline]
    fn new(minimizer: u64, rid: u32, pos: u64, rev: u8) -> Self {
        debug_assert!(rid <= Self::MAX_RID);
        debug_assert!(pos <= u64::from(u32::MAX));
        debug_assert!(rev <= 1);
        Self {
            minimizer,
            pos: pos as u32,
            rid_rev: (rid << 1) | u32::from(rev),
        }
    }
}

impl Ord for MinimizerData32 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.minimizer.cmp(&other.minimizer)
    }
}

impl PartialOrd for MinimizerData32 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub(crate) trait MinimizerRecord: Copy + Ord + Send + Sync {
    fn from_parts(minimizer: u64, rid: u32, pos: u64, rev: u8, span: u8) -> Self;
    fn dummy(rid: u32) -> Self;
    fn minimizer(&self) -> u64;
    fn rid(&self) -> u32;
    fn pos(&self) -> u64;
    fn rev(&self) -> u8;
}

impl MinimizerRecord for MinimizerData {
    #[inline]
    fn from_parts(minimizer: u64, rid: u32, pos: u64, rev: u8, span: u8) -> Self {
        Self {
            minimizer,
            info: MinimizerInfo {
                rid,
                pos,
                rev,
                span,
            },
        }
    }

    #[inline]
    fn dummy(rid: u32) -> Self {
        Self::from_parts(u64::MAX, rid, u64::MAX, 0, u8::MAX)
    }

    #[inline]
    fn minimizer(&self) -> u64 {
        self.minimizer
    }

    #[inline]
    fn rid(&self) -> u32 {
        self.info.rid
    }

    #[inline]
    fn pos(&self) -> u64 {
        self.info.pos
    }

    #[inline]
    fn rev(&self) -> u8 {
        self.info.rev
    }
}

impl MinimizerRecord for MinimizerData32 {
    #[inline]
    fn from_parts(minimizer: u64, rid: u32, pos: u64, rev: u8, _span: u8) -> Self {
        Self::new(minimizer, rid, pos, rev)
    }

    #[inline]
    fn dummy(rid: u32) -> Self {
        Self {
            minimizer: u64::MAX,
            pos: u32::MAX,
            rid_rev: rid << 1,
        }
    }

    #[inline]
    fn minimizer(&self) -> u64 {
        self.minimizer
    }

    #[inline]
    fn rid(&self) -> u32 {
        self.rid_rev >> 1
    }

    #[inline]
    fn pos(&self) -> u64 {
        u64::from(self.pos)
    }

    #[inline]
    fn rev(&self) -> u8 {
        (self.rid_rev & 1) as u8
    }
}

pub fn complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

// rolling hash algorithm
pub fn hash(seq: &str, k: usize) -> u64 {
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

fn update_hash(h: u64, k: usize, c: char) -> u64 {
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

fn minimizer(seq: &str, rid: &u32, start: u64, k: usize) -> MinimizerData {
    let mut h = hash(seq, k);
    let mut m = h;
    let mut pos = 0;
    for i in 1..(seq.len() - k + 1) {
        h = update_hash(h, k, seq.chars().nth(i + k - 1).unwrap());
        if h < m {
            m = h;
            pos = i as u32;
        }
    }

    let pos = start + pos as u64;
    let m = MinimizerData {
        minimizer: m,
        info: MinimizerInfo {
            rid: *rid,
            pos: pos,
            rev: 0,
            span: k as u8,
        },
    };
    println!("{}, {:?}", seq, m);
    m
}

pub fn minimizer_nthash(seq: &[u8], rid: &u32, start: u64, k: usize) -> AnyResult<MinimizerData> {
    let hash_iter = NtHashIterator::new(seq, k).unwrap();

    if let Some((i, (m, rev))) = hash_iter.enumerate().min_by_key(|&(_, x)| x.0) {
        return Ok(MinimizerData {
            minimizer: m,
            info: MinimizerInfo {
                rid: *rid,
                pos: start + i as u64,
                rev,
                span: k as u8,
            },
        });
    } else {
        return Err(anyhow::anyhow!("No minimizer found!"));
    };
}

const BASE_ENCODING: [u8; 256] = {
    let mut table = [4; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table[b'U' as usize] = 3;
    table[b'u' as usize] = 3;
    table
};

#[inline(always)]
fn yak_hash64(mut key: u64) -> u64 {
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8);
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4);
    key ^= key >> 28;
    key.wrapping_add(key << 31)
}

fn sketch_impl<M: MinimizerRecord>(seq: &[u8], rid: u32, k: usize, w: usize) -> Vec<M> {
    assert!((1..=63).contains(&k), "k must be between 1 and 63");
    assert!((1..256).contains(&w), "w must be between 1 and 255");
    if seq.is_empty() {
        return Vec::new();
    }

    let dummy = M::dummy(rid);
    let shift = k - 1;
    let mask = (1_u64 << k) - 1;
    let mut kmers = [0_u64; 4];
    let mut buffer = vec![dummy; w];
    let mut result = Vec::with_capacity(seq.len() / w + 1);
    let mut min = dummy;
    let mut min_pos = 0;
    let mut buffer_pos = 0;
    let mut valid_run = 0_usize;

    for (position, &base) in seq.iter().enumerate() {
        let code = BASE_ENCODING[base as usize];
        let mut info = dummy;
        if code < 4 {
            kmers[0] = ((kmers[0] << 1) | u64::from(code & 1)) & mask;
            kmers[1] = ((kmers[1] << 1) | u64::from(code >> 1)) & mask;
            kmers[2] = (kmers[2] >> 1) | (u64::from(1 - (code & 1)) << shift);
            kmers[3] = (kmers[3] >> 1) | (u64::from(1 - (code >> 1)) << shift);

            // The middle bit-plane identifies a reverse-complement palindrome,
            // whose strand is undefined in partig.
            if kmers[1] == kmers[3] {
                continue;
            }
            let rev = usize::from(kmers[1] >= kmers[3]);
            valid_run += 1;
            if valid_run >= k {
                info = M::from_parts(
                    yak_hash64(kmers[rev << 1]).wrapping_add(yak_hash64(kmers[(rev << 1) | 1])),
                    rid,
                    position as u64,
                    rev as u8,
                    k as u8,
                );
            }
        } else {
            valid_run = 0;
            kmers = [0; 4];
        }

        buffer[buffer_pos] = info;
        if valid_run == w + k - 1 && min.minimizer() != u64::MAX {
            for item in buffer[(buffer_pos + 1)..]
                .iter()
                .chain(buffer[..buffer_pos].iter())
            {
                if item.minimizer() == min.minimizer() && item.pos() != min.pos() {
                    result.push(*item);
                }
            }
        }

        if info.minimizer() <= min.minimizer() {
            if valid_run >= w + k && min.minimizer() != u64::MAX {
                result.push(min);
            }
            min = info;
            min_pos = buffer_pos;
        } else if buffer_pos == min_pos {
            if valid_run >= w + k - 1 && min.minimizer() != u64::MAX {
                result.push(min);
            }
            min = dummy;
            for (index, item) in buffer[(buffer_pos + 1)..]
                .iter()
                .enumerate()
                .map(|(i, item)| (buffer_pos + 1 + i, item))
                .chain(buffer[..=buffer_pos].iter().enumerate())
            {
                if item.minimizer() <= min.minimizer() {
                    min = *item;
                    min_pos = index;
                }
            }
            if valid_run >= w + k - 1 && min.minimizer() != u64::MAX {
                for item in buffer[(buffer_pos + 1)..]
                    .iter()
                    .chain(buffer[..=buffer_pos].iter())
                {
                    if item.minimizer() == min.minimizer() && item.pos() != min.pos() {
                        result.push(*item);
                    }
                }
            }
        }
        buffer_pos += 1;
        if buffer_pos == w {
            buffer_pos = 0;
        }
    }

    if min.minimizer() != u64::MAX {
        result.push(min);
    }
    result
}

/// Find symmetric `(w,k)` minimizers using the same hash and tie handling as
/// partig. The sequence is scanned once; no per-window k-mer rehashing or
/// normalized sequence copy is needed.
pub fn sketch(seq: &[u8], rid: u32, k: usize, w: usize) -> Vec<MinimizerData> {
    sketch_impl(seq, rid, k, w)
}

pub(crate) fn sketch32(seq: &[u8], rid: u32, k: usize, w: usize) -> Vec<MinimizerData32> {
    assert!(rid <= MinimizerData32::MAX_RID);
    assert!((seq.len() as u64) <= u64::from(u32::MAX) + 1);
    sketch_impl(seq, rid, k, w)
}

#[cfg(test)]
mod compact_tests {
    use super::*;

    #[test]
    fn compact_layout_and_packing_are_lossless() {
        assert_eq!(std::mem::size_of::<MinimizerData32>(), 16);
        let item = MinimizerData32::new(17, MinimizerData32::MAX_RID, u32::MAX.into(), 1);
        assert_eq!(item.minimizer(), 17);
        assert_eq!(item.rid(), MinimizerData32::MAX_RID);
        assert_eq!(item.pos(), u64::from(u32::MAX));
        assert_eq!(item.rev(), 1);
    }

    #[test]
    fn compact_and_wide_sketches_are_logically_identical() {
        let sequence = b"TTTCGACAGTTCTCCCTGGCACCTCTGAAAGCTTTCCTGGTTTNATTGTTGAAAGTCTTAGGGCTCAACTTGGTCAGCCCTTCTTCATGGAAATTGTTATGACCATGTGTTGGTCCATCTGGATGATGCGCAATGATGTCATTTTCAAAGGTTTAC";
        let wide = sketch(sequence, 7, 19, 19);
        let compact = sketch32(sequence, 7, 19, 19);
        assert_eq!(wide.len(), compact.len());
        for (wide, compact) in wide.iter().zip(&compact) {
            assert_eq!(wide.minimizer(), compact.minimizer());
            assert_eq!(wide.rid(), compact.rid());
            assert_eq!(wide.pos(), compact.pos());
            assert_eq!(wide.rev(), compact.rev());
        }
    }
}
