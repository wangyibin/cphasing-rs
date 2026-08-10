#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::{Context, Result as anyResult};
// use rdst::RadixSort;
use noodles::fasta;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use std::error::Error;
use std::io::{BufRead, BufReader, Write};
use std::mem::MaybeUninit;
use std::ops::Range;
use std::path::Path;
use std::process::exit;
// use petgraph::prelude::*;
// use petgraph::visit::NodeIndexable;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

use crate::core::{BaseTable, ContigPair, ContigPair2};
use crate::core::{common_reader, common_writer};
use crate::sketch::{MinimizerData, MinimizerData32, MinimizerRecord, sketch, sketch32};

// // maximal cliques
// pub fn bron_kerbosch(
//     graph: &Graph<(), u64, Undirected>,
//     mut r: HashSet<NodeIndex>,
//     mut p: HashSet<NodeIndex>,
//     mut x: HashSet<NodeIndex>,
//     cliques: &mut Vec<HashSet<NodeIndex>>,
// ) {
//     if p.is_empty() && x.is_empty() {
//         cliques.push(r.clone());
//         return;
//     }

//     let mut p_candidates = p.clone();
//     p_candidates.retain(|&v| {
//         let mut neighbors = graph.neighbors(v);
//         neighbors.all(|n| p.contains(&n))
//     });

//     let p_cloned = p.clone();
//     for v in p_cloned {
//         let mut neighbors = graph.neighbors(v);
//         let mut r_new = r.clone();
//         r_new.insert(v);

//         let p_new = p.intersection(&neighbors.clone().collect::<HashSet<_>>()).cloned().collect();
//         let x_new = x.intersection(&neighbors.collect::<HashSet<_>>()).cloned().collect();

//         bron_kerbosch(graph, r_new, p_new, x_new, cliques);

//         p.remove(&v);
//         x.insert(v);
//     }
// }

// pub fn find_cliques(graph: &Graph<(), u64, Undirected>) -> Vec<HashSet<NodeIndex>> {
//     let mut cliques = Vec::new();
//     let n = graph.node_count();
//     let mut all_nodes = (0..n).map(NodeIndex::new).collect::<HashSet<NodeIndex>>();
//     let r = HashSet::new();

//     let x = HashSet::new();

//     bron_kerbosch(graph, r, all_nodes, x, &mut cliques);
//     cliques
// }

#[derive(Debug, Clone)]
pub struct AlleleHeader {
    pub header: Vec<String>,
    pub contigsizes: HashMap<String, u64>,
    pub contigs: HashSet<String>,
    pub minimizer: HashMap<String, u64>,
    pub unique_minimizer: HashMap<String, u64>,
}

impl AlleleHeader {
    pub fn new() -> Self {
        Self {
            header: Vec::new(),
            contigsizes: HashMap::new(),
            contigs: HashSet::new(),
            minimizer: HashMap::new(),
            unique_minimizer: HashMap::new(),
        }
    }

    pub fn from_file(&mut self, file: &str) -> anyResult<()> {
        let input = common_reader(file);
        let reader = BufReader::new(input);
        let mut header: Vec<String> = Vec::new();
        let mut contigsizes: HashMap<String, u64> = HashMap::new();
        let mut contigs: HashSet<String> = HashSet::new();
        let mut minimizer: HashMap<String, u64> = HashMap::new();
        let mut unique_minimizer: HashMap<String, u64> = HashMap::new();

        for result in reader.lines() {
            let line = result?;
            if !line.starts_with("#") {
                break;
            }
            let line_vec = line
                .strip_prefix("#")
                .unwrap()
                .split(" ")
                .collect::<Vec<&str>>();
            let contig = line_vec[0].to_string();
            let size = line_vec[1].parse::<u64>().unwrap();
            let min = line_vec[2].parse::<u64>().unwrap();
            let unique_min = line_vec[3].parse::<u64>().unwrap();
            header.push(contig.clone());
            contigsizes.insert(contig.clone(), size);
            contigs.insert(contig.clone());
            minimizer.insert(contig.clone(), min);
            unique_minimizer.insert(contig.clone(), unique_min);
        }

        self.header = header;
        self.contigsizes = contigsizes;
        self.contigs = contigs;
        self.minimizer = minimizer;
        self.unique_minimizer = unique_minimizer;

        Ok(())
    }

    pub fn to_unique_minimizer_density(&mut self) -> HashMap<String, f64> {
        let mut density_hash: HashMap<String, f64> = HashMap::new();
        for (contig, min) in self.unique_minimizer.iter() {
            let minimizer = self.minimizer.get(contig).unwrap();
            let density = *min as f64 / *minimizer as f64;
            density_hash.insert(contig.clone(), density);
        }

        density_hash
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct AlleleRecord2 {
    pub idx1: u32,
    pub idx2: u32,
    pub contig1: String,
    pub contig2: String,
    pub mz1: u32,
    pub mz2: u32,
    pub mz_shared: u32,
    pub similarity: f64,
    pub strand: i8,
}

#[derive(Debug)]
pub struct AlleleTable2 {
    pub file: String,
    pub header: AlleleHeader,
    pub allele_records: Vec<AlleleRecord2>,
}

impl BaseTable for AlleleTable2 {
    fn new(name: &String) -> AlleleTable2 {
        let mut header = AlleleHeader::new();
        let _ = header.from_file(name);

        AlleleTable2 {
            file: name.clone(),
            header: header,
            allele_records: Vec::new(),
        }
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

impl AlleleTable2 {
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

    pub fn parse_header(&self) -> anyResult<AlleleHeader> {
        let mut header = AlleleHeader::new();
        let input = common_reader(&self.file);
        let reader = BufReader::new(input);

        for result in reader.lines() {
            let line = result?;
            if !line.starts_with("#") {
                break;
            }
            let line_vec = line
                .strip_prefix("#")
                .unwrap()
                .split(" ")
                .collect::<Vec<&str>>();
            let contig = line_vec[0].to_string();
            let size = line_vec[1].parse::<u64>().unwrap();
            let min = line_vec[2].parse::<u64>().unwrap();
            let unique_min = line_vec[3].parse::<u64>().unwrap();
            header.header.push(contig.clone());
            header.contigsizes.insert(contig.clone(), size);
            header.contigs.insert(contig.clone());
            header.minimizer.insert(contig.clone(), min);
            header.unique_minimizer.insert(contig.clone(), unique_min);
        }

        Ok(header)
    }

    pub fn allele_records(&self) -> anyResult<Vec<AlleleRecord2>> {
        let parse_result = self.parse();
        let mut records: Vec<AlleleRecord2> = Vec::new();

        match parse_result {
            Ok(mut rdr) => {
                for result in rdr.deserialize() {
                    let record: AlleleRecord2 = result?;
                    records.push(record);
                }
            }
            Err(e) => {
                eprintln!("Error parsing allele table: {}", e);
                exit(1);
            }
        }

        Ok(records)
    }

    pub fn get_allelic_contig_pairs(&self) -> HashSet<ContigPair2<'_>> {
        let mut contig_pairs: HashSet<ContigPair2> = HashSet::new();
        let mut records = &self.allele_records;
        // sort records by mz_shared in ascending order
        // records.sort_by(|a, b| a.mz_shared.partial_cmp(&b.mz_shared).unwrap());

        for record in records {
            let contig1 = &record.contig1;
            let contig2 = &record.contig2;
            let mut contig_pair = ContigPair2::new(&contig1, &contig2);
            contig_pair.order();
            contig_pairs.insert(contig_pair);
        }

        contig_pairs
    }

    pub fn get_allelic_record_by_contig_pairs(&self) -> HashMap<ContigPair2<'_>, &AlleleRecord2> {
        let mut data: HashMap<ContigPair2, &AlleleRecord2> = HashMap::new();
        let records = &self.allele_records;
        for record in records {
            let mut contig_pair = ContigPair2::new(&record.contig1, &record.contig2);
            contig_pair.order();
            data.insert(contig_pair, &record);
        }

        data
    }

    pub fn get_allelic_contigs_precise(
        &self,
        whitehash: &HashSet<&String>,
    ) -> HashMap<&String, Vec<Vec<&String>>> {
        let mut data: HashMap<&String, Vec<Vec<&String>>> = HashMap::new();
        let check_whitehash = !whitehash.is_empty();

        for record in &self.allele_records {
            if check_whitehash
                && (!whitehash.contains(&record.contig1) || !whitehash.contains(&record.contig2))
            {
                continue;
            }

            if !data.contains_key(&record.contig1) {
                data.insert(
                    &record.contig1,
                    vec![vec![&record.contig1, &record.contig2]],
                );
            } else {
                data.get_mut(&record.contig1).unwrap()[0].push(&record.contig2);
            }
        }

        data
    }

    pub fn get_allelic_contigs(
        &self,
        method: &str,
        whitehash: &HashSet<&String>,
    ) -> HashMap<&String, Vec<Vec<&String>>> {
        let mut data: HashMap<&String, Vec<Vec<&String>>> = HashMap::new();
        for record in &self.allele_records {
            if !whitehash.is_empty()
                && (!whitehash.contains(&record.contig1) || !whitehash.contains(&record.contig2))
            {
                continue;
            }

            if !data.contains_key(&record.contig1) {
                data.insert(
                    &record.contig1,
                    vec![vec![&record.contig1, &record.contig2]],
                );
            } else {
                if method == "fast" {
                    continue;
                }
                data.get_mut(&record.contig1)
                    .unwrap()
                    .push(vec![&record.contig1, &record.contig2]);
            }
        }

        data
    }

    pub fn write(&self, records: &Vec<AlleleRecord2>) -> anyResult<()> {
        let writer = common_writer(&self.file);
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);

        for record in records {
            wtr.serialize(record)?;
        }

        wtr.flush()?;
        Ok(())
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AlleleStrandRecord {
    pub idx1: u32,
    pub idx2: u32,
    pub contig1: String,
    pub contig2: String,
    pub length1: u32,
    pub length2: u32,
    pub matches: u32,
    pub identity: f32,
    pub strand: i8,
}

#[derive(Debug, Clone)]
pub struct AlleleStrandTable {
    pub file: String,
    pub allele_strand_records: Vec<AlleleStrandRecord>,
}

impl BaseTable for AlleleStrandTable {
    fn new(name: &String) -> AlleleStrandTable {
        AlleleStrandTable {
            file: name.clone(),
            allele_strand_records: Vec::new(),
        }
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

impl AlleleStrandTable {
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

    pub fn get_allele_strand_records(&mut self) {
        let parse_result = self.parse();

        match parse_result {
            Ok(mut rdr) => {
                for result in rdr.deserialize() {
                    let record: AlleleStrandRecord = result.unwrap();
                    self.allele_strand_records.push(record);
                }
            }
            Err(e) => {
                eprintln!("Error parsing allele table: {}", e);
                exit(1);
            }
        }
    }

    pub fn get_info(&self) -> HashMap<ContigPair2<'_>, &AlleleStrandRecord> {
        let mut data: HashMap<ContigPair2, &AlleleStrandRecord> = HashMap::new();

        for record in &self.allele_strand_records {
            if record.contig1 >= record.contig2 {
                continue;
            }

            let mut contig_pair = ContigPair2::new(&record.contig1, &record.contig2);

            data.insert(contig_pair, &record);
        }

        data
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AlleleRecord {
    pub data: Vec<String>,
}

#[derive(Debug)]
pub struct AlleleTable {
    pub file: String,
    pub allele_records: Vec<AlleleRecord>,
}

impl AlleleTable {
    pub fn new(name: &String) -> AlleleTable {
        AlleleTable {
            file: name.clone(),
            allele_records: Vec::new(),
        }
    }

    pub fn parse(&self) -> anyResult<Box<dyn BufRead + Send>> {
        let input = common_reader(&self.file);

        Ok(input)
    }

    pub fn allele_records(&self) -> anyResult<Vec<AlleleRecord>> {
        let parse_result = self.parse();
        let mut records: Vec<AlleleRecord> = Vec::new();

        let reader = BufReader::new(parse_result.unwrap());
        for result in reader.lines() {
            let line = result?;
            if line.starts_with("#") {
                continue;
            }
            //remove whitespace suffix
            let line = line.trim_end();
            let line_vec = line.split("\t").collect::<Vec<&str>>();

            // convert Vec<&str> to Vec<String>
            let record = line_vec[2..]
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>();
            if record.len() < 2 {
                continue;
            }
            records.push(AlleleRecord { data: record });
        }

        Ok(records)
    }

    pub fn to_contig_db(
        &self,
        method: &str,
        whitehash: &HashSet<&String>,
    ) -> HashMap<&String, Vec<&AlleleRecord>> {
        let mut data: HashMap<&String, Vec<&AlleleRecord>> = HashMap::new();
        for record in &self.allele_records {
            'inner: for contig in &record.data {
                if whitehash.len() > 0 {
                    if !whitehash.contains(contig) {
                        continue 'inner;
                    }
                }

                match method {
                    "fast" => {
                        if data.contains_key(contig) {
                            continue;
                        }
                        data.entry(contig).or_insert_with(Vec::new);
                        data.get_mut(contig).unwrap().push(record);
                    }
                    "precise" => {
                        data.entry(contig).or_insert_with(Vec::new);
                        data.get_mut(contig).unwrap().push(record);
                    }
                    _ => todo!(),
                }
            }
        }

        data
    }

    pub fn get_allelic_contig_pairs(
        &self,
        whitehash: &HashSet<&String>,
    ) -> HashSet<ContigPair2<'_>> {
        let mut contig_pairs: HashSet<ContigPair2> = HashSet::new();

        for record in &self.allele_records {
            for i in 0..record.data.len() - 1 {
                'inner: for j in (i + 1)..record.data.len() {
                    let contig1 = &record.data[i];
                    let contig2 = &record.data[j];
                    if whitehash.len() == 0 {
                        if !whitehash.contains(contig1) || !whitehash.contains(contig2) {
                            continue 'inner;
                        }
                    }
                    let contig_pair = ContigPair2::new(&contig1, &contig2);

                    contig_pairs.insert(contig_pair);
                }
            }
        }

        contig_pairs
    }
}

// https://github.com/lh3/partig

pub struct Uinfo {
    pub cnt1: u32,
    pub cnt2: u32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Anchor {
    key: u64,
    pos1: u64,
    pos2: u64,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct Anchor32 {
    key: u64,
    pos1: u32,
    pos2: u32,
}

#[inline]
fn anchor_key(rid1: u32, rid2: u32, rev: u8) -> u64 {
    debug_assert!(rid1 < (1 << 31) && rid2 < (1 << 31));
    (u64::from(rid1) << 33) | (u64::from(rid2) << 1) | u64::from(rev)
}

#[inline]
fn key_rid1(key: u64) -> u32 {
    (key >> 33) as u32
}

#[inline]
fn key_rid2(key: u64) -> u32 {
    (key as u32) >> 1
}

#[inline]
fn key_rev(key: u64) -> u8 {
    (key & 1) as u8
}

impl Anchor {
    #[inline]
    fn new(rid1: u32, pos1: u64, rid2: u32, pos2: u64, rev: u8) -> Self {
        Self {
            key: anchor_key(rid1, rid2, rev),
            pos1,
            pos2,
        }
    }

    #[inline]
    fn rid1(&self) -> u32 {
        key_rid1(self.key)
    }

    #[inline]
    fn rid2(&self) -> u32 {
        key_rid2(self.key)
    }

    #[inline]
    fn rev(&self) -> u8 {
        key_rev(self.key)
    }
}

impl Ord for Anchor {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

impl PartialOrd for Anchor {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Anchor32 {
    #[inline]
    fn new(rid1: u32, pos1: u64, rid2: u32, pos2: u64, rev: u8) -> Self {
        debug_assert!(pos1 <= u64::from(u32::MAX) && pos2 <= u64::from(u32::MAX));
        Self {
            key: anchor_key(rid1, rid2, rev),
            pos1: pos1 as u32,
            pos2: pos2 as u32,
        }
    }
}

impl Ord for Anchor32 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

impl PartialOrd for Anchor32 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

trait AnchorRecord: Copy + Ord + Send + Sync {
    fn key(&self) -> u64;
    fn pos1(&self) -> u64;
    fn pos2(&self) -> u64;
}

impl AnchorRecord for Anchor {
    #[inline]
    fn key(&self) -> u64 {
        self.key
    }

    #[inline]
    fn pos1(&self) -> u64 {
        self.pos1
    }

    #[inline]
    fn pos2(&self) -> u64 {
        self.pos2
    }
}

impl AnchorRecord for Anchor32 {
    #[inline]
    fn key(&self) -> u64 {
        self.key
    }

    #[inline]
    fn pos1(&self) -> u64 {
        u64::from(self.pos1)
    }

    #[inline]
    fn pos2(&self) -> u64 {
        u64::from(self.pos2)
    }
}

enum AnchorStorage {
    U32(Vec<Anchor32>),
    U64(Vec<Anchor>),
}

enum MinimizerStorage {
    U32(Vec<MinimizerData32>),
    U64(Vec<MinimizerData>),
}

struct TrimmedSequence {
    rid: u32,
    name: String,
    sequence: Vec<u8>,
    trimmed_range: Range<usize>,
}

impl TrimmedSequence {
    fn trimmed_sequence(&self) -> &[u8] {
        &self.sequence[self.trimmed_range.clone()]
    }

    fn trimmed_length(&self) -> usize {
        self.trimmed_range.len()
    }

    fn into_trimmed_sequence(mut self) -> Vec<u8> {
        if self.trimmed_range.start == 0 && self.trimmed_range.end == self.sequence.len() {
            return self.sequence;
        }
        let trimmed_length = self.trimmed_range.len();
        self.sequence.copy_within(self.trimmed_range, 0);
        self.sequence.truncate(trimmed_length);
        self.sequence
    }
}

struct SequenceMetadata {
    rid: u32,
    name: String,
    length: u64,
}

struct SketchedSequence {
    metadata: SequenceMetadata,
    minimizers: MinimizerStorage,
}

#[derive(Default)]
struct SketchAccumulator {
    metadata: Vec<SequenceMetadata>,
    minimizers: Option<MinimizerStorage>,
}

fn append_into_larger<T>(mut left: Vec<T>, mut right: Vec<T>) -> Vec<T> {
    if right.capacity() > left.capacity() {
        std::mem::swap(&mut left, &mut right);
    }
    left.append(&mut right);
    left
}

fn extend_wide_from_compact(
    wide: &mut Vec<MinimizerData>,
    compact: Vec<MinimizerData32>,
    k: usize,
) {
    wide.reserve(compact.len());
    wide.extend(compact.into_iter().map(|item| {
        MinimizerData::from_parts(
            item.minimizer(),
            item.rid(),
            item.pos(),
            item.rev(),
            k as u8,
        )
    }));
}

impl MinimizerStorage {
    fn len(&self) -> usize {
        match self {
            Self::U32(minimizers) => minimizers.len(),
            Self::U64(minimizers) => minimizers.len(),
        }
    }

    fn is_compact(&self) -> bool {
        matches!(self, Self::U32(_))
    }

    fn par_sort_unstable(&mut self) {
        match self {
            Self::U32(minimizers) => minimizers.par_sort_unstable(),
            Self::U64(minimizers) => minimizers.par_sort_unstable(),
        }
    }

    fn merge(self, other: Self, k: usize) -> Self {
        match (self, other) {
            (Self::U32(left), Self::U32(right)) => Self::U32(append_into_larger(left, right)),
            (Self::U64(left), Self::U64(right)) => Self::U64(append_into_larger(left, right)),
            (Self::U64(mut wide), Self::U32(compact)) => {
                extend_wide_from_compact(&mut wide, compact, k);
                Self::U64(wide)
            }
            (Self::U32(compact), Self::U64(mut wide)) => {
                extend_wide_from_compact(&mut wide, compact, k);
                Self::U64(wide)
            }
        }
    }
}

impl SketchAccumulator {
    fn push(&mut self, record: SketchedSequence, k: usize) {
        self.metadata.push(record.metadata);
        self.minimizers = Some(match self.minimizers.take() {
            Some(minimizers) => minimizers.merge(record.minimizers, k),
            None => record.minimizers,
        });
    }

    fn merge(mut self, mut other: Self, k: usize) -> Self {
        self.metadata = append_into_larger(self.metadata, other.metadata);
        self.minimizers = match (self.minimizers, other.minimizers) {
            (Some(left), Some(right)) => Some(left.merge(right, k)),
            (Some(minimizers), None) | (None, Some(minimizers)) => Some(minimizers),
            (None, None) => None,
        };
        self
    }
}

impl AnchorStorage {
    fn len(&self) -> usize {
        match self {
            Self::U32(anchor) => anchor.len(),
            Self::U64(anchor) => anchor.len(),
        }
    }

    fn radix_sort_by_key(&mut self) {
        match self {
            Self::U32(anchor) => parallel_radix_sort_by_key(anchor),
            Self::U64(anchor) => parallel_radix_sort_by_key(anchor),
        }
    }
}

const RADIX_BITS: u32 = 11;
const RADIX_BUCKETS: usize = 1 << RADIX_BITS;
const RADIX_MASK: u64 = (RADIX_BUCKETS as u64) - 1;
const RADIX_PASSES: u32 = 64_u32.div_ceil(RADIX_BITS);
const RADIX_MIN_LENGTH: usize = 1 << 16;

fn radix_pass<A: AnchorRecord>(source: &[A], destination: *mut A, shift: u32) {
    let workers = rayon::current_num_threads().max(1).min(source.len());
    let chunk_length = source.len().div_ceil(workers);
    let counts = source
        .par_chunks(chunk_length)
        .map(|chunk| {
            let mut local = vec![0_usize; RADIX_BUCKETS];
            for item in chunk {
                local[((item.key() >> shift) & RADIX_MASK) as usize] += 1;
            }
            local
        })
        .collect::<Vec<_>>();

    let mut offsets = vec![vec![0_usize; RADIX_BUCKETS]; counts.len()];
    let mut bucket_start = 0_usize;
    for bucket in 0..RADIX_BUCKETS {
        let mut cursor = bucket_start;
        for chunk in 0..counts.len() {
            offsets[chunk][bucket] = cursor;
            cursor += counts[chunk][bucket];
        }
        bucket_start = cursor;
    }
    debug_assert_eq!(bucket_start, source.len());

    // The prefix sums above assign every input chunk a disjoint output range
    // for every bucket. Each destination slot is therefore written exactly
    // once, and all writes complete before this function returns.
    let destination_address = destination as usize;
    source
        .par_chunks(chunk_length)
        .zip(offsets.into_par_iter())
        .for_each(|(chunk, mut local_offsets)| {
            let destination = destination_address as *mut A;
            for item in chunk {
                let bucket = ((item.key() >> shift) & RADIX_MASK) as usize;
                let index = local_offsets[bucket];
                local_offsets[bucket] += 1;
                // SAFETY: per-chunk prefix sums are disjoint, in bounds, and
                // cover the full destination exactly once.
                unsafe { destination.add(index).write(*item) };
            }
        });
}

fn parallel_radix_sort_by_key<A: AnchorRecord>(values: &mut Vec<A>) {
    if values.len() < RADIX_MIN_LENGTH {
        values.par_sort_unstable_by_key(AnchorRecord::key);
        return;
    }

    let mut scratch = Vec::<MaybeUninit<A>>::new();
    if let Err(error) = scratch.try_reserve_exact(values.len()) {
        log::warn!("Cannot allocate radix-sort scratch buffer ({error}); use comparison sort.");
        values.par_sort_unstable_by_key(AnchorRecord::key);
        return;
    }
    // Every pass writes every element before the scratch buffer is read.
    unsafe { scratch.set_len(values.len()) };

    for pass in 0..RADIX_PASSES {
        let shift = pass * RADIX_BITS;
        if pass % 2 == 0 {
            radix_pass(values, scratch.as_mut_ptr().cast::<A>(), shift);
        } else {
            // SAFETY: the preceding even pass initialized every scratch slot.
            let source =
                unsafe { std::slice::from_raw_parts(scratch.as_ptr().cast::<A>(), scratch.len()) };
            radix_pass(source, values.as_mut_ptr(), shift);
        }
    }

    // RADIX_PASSES is even, so the final pass always writes back to `values`.
    debug_assert_eq!(RADIX_PASSES % 2, 0);
}

#[inline]
fn use_32bit_anchor_coordinates(lengths: &[u64]) -> bool {
    lengths.iter().all(|&length| length <= u64::from(u32::MAX))
}

const MAX_PACKED_CONTIGS: usize = 1_usize << 31;

#[inline]
fn use_32bit_minimizer_coordinates(lengths: &[u64], contig_count: usize) -> bool {
    contig_count <= MAX_PACKED_CONTIGS
        && lengths
            .iter()
            .all(|&length| length <= u64::from(u32::MAX) + 1)
}

fn trimmed_sequence_range(sequence_length: usize, trim_length: usize) -> Range<usize> {
    if sequence_length > trim_length.saturating_mul(3) {
        trim_length..sequence_length - trim_length
    } else {
        0..sequence_length
    }
}

fn trimmed_fasta_records(
    file: &str,
    trim_length: usize,
) -> impl Iterator<Item = anyResult<TrimmedSequence>> + Send {
    let mut reader = fasta::Reader::new(common_reader(file));
    let mut raw_definition = String::new();
    let mut next_rid = 0_usize;
    let mut finished = false;

    std::iter::from_fn(move || {
        if finished {
            return None;
        }

        let result: anyResult<Option<TrimmedSequence>> = (|| {
            raw_definition.clear();
            let record_number = next_rid + 1;
            let definition_length =
                reader
                    .read_definition(&mut raw_definition)
                    .with_context(|| {
                        format!("failed to read FASTA definition for record {record_number}")
                    })?;
            if definition_length == 0 {
                return Ok(None);
            }
            let definition: fasta::record::Definition = raw_definition
                .parse()
                .with_context(|| format!("invalid FASTA definition for record {record_number}"))?;
            anyhow::ensure!(
                next_rid < MAX_PACKED_CONTIGS,
                "alleles supports at most {MAX_PACKED_CONTIGS} contigs"
            );

            let mut sequence = Vec::new();
            reader.read_sequence(&mut sequence).with_context(|| {
                format!("failed to read FASTA sequence for record {record_number}")
            })?;
            let trimmed_range = trimmed_sequence_range(sequence.len(), trim_length);

            let record = TrimmedSequence {
                rid: next_rid as u32,
                name: definition.name().to_owned(),
                sequence,
                trimmed_range,
            };
            next_rid += 1;
            Ok(Some(record))
        })();

        match result {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => {
                finished = true;
                None
            }
            Err(error) => {
                finished = true;
                Some(Err(error))
            }
        }
    })
}

#[derive(Debug)]
struct SplitRegion {
    name: String,
    start: usize,
    end: usize,
}

fn load_split_regions(file: &str) -> anyResult<HashMap<String, Vec<SplitRegion>>> {
    let reader = common_reader(file);
    let mut regions: HashMap<String, Vec<SplitRegion>> = HashMap::new();
    let mut split_names = HashSet::new();

    for (line_index, line) in reader.lines().enumerate() {
        let line_number = line_index + 1;
        let line = line.with_context(|| {
            format!("failed to read split-regions file `{file}` at line {line_number}")
        })?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields = line.split('\t').collect::<Vec<_>>();
        anyhow::ensure!(
            fields.len() == 4,
            "invalid split-regions line {line_number} in `{file}`: expected 4 tab-separated fields"
        );
        let split_name = fields[0];
        let source_contig = fields[1];
        anyhow::ensure!(
            !split_name.is_empty() && !source_contig.is_empty(),
            "invalid split-regions line {line_number} in `{file}`: contig names cannot be empty"
        );
        anyhow::ensure!(
            split_names.insert(split_name.to_owned()),
            "duplicate split sequence name `{split_name}` in `{file}`"
        );
        let start = fields[2].parse::<usize>().with_context(|| {
            format!("invalid start coordinate at line {line_number} in `{file}`")
        })?;
        let end = fields[3]
            .parse::<usize>()
            .with_context(|| format!("invalid end coordinate at line {line_number} in `{file}`"))?;
        anyhow::ensure!(
            start < end,
            "invalid interval {start}_{end} for `{split_name}`: start must be smaller than end"
        );
        regions
            .entry(source_contig.to_owned())
            .or_default()
            .push(SplitRegion {
                name: split_name.to_owned(),
                start,
                end,
            });
    }

    anyhow::ensure!(!regions.is_empty(), "split-regions file `{file}` is empty");
    Ok(regions)
}

fn split_fasta_records(
    fasta_file: &str,
    regions_file: &str,
    trim_length: usize,
) -> anyResult<impl Iterator<Item = anyResult<TrimmedSequence>> + Send> {
    let regions = load_split_regions(regions_file)?;
    let mut reader = fasta::Reader::new(common_reader(fasta_file));
    let mut raw_definition = String::new();
    let mut pending = VecDeque::new();
    let mut seen_sources = HashSet::new();
    let mut next_rid = 0_usize;
    let mut finished = false;

    Ok(std::iter::from_fn(move || {
        loop {
            if let Some(record) = pending.pop_front() {
                return Some(Ok(record));
            }
            if finished {
                return None;
            }

            raw_definition.clear();
            let result: anyResult<bool> = (|| {
                let definition_length = reader
                    .read_definition(&mut raw_definition)
                    .context("failed to read FASTA definition while applying split regions")?;
                if definition_length == 0 {
                    let mut missing = regions
                        .keys()
                        .filter(|name| !seen_sources.contains(*name))
                        .cloned()
                        .collect::<Vec<_>>();
                    missing.sort_unstable();
                    anyhow::ensure!(
                        missing.is_empty(),
                        "split-regions source contigs were not found in FASTA: {}",
                        missing.join(", ")
                    );
                    return Ok(false);
                }

                let definition: fasta::record::Definition = raw_definition
                    .parse()
                    .context("invalid FASTA definition while applying split regions")?;
                let mut sequence = Vec::new();
                reader.read_sequence(&mut sequence).with_context(|| {
                    format!("failed to read FASTA sequence `{}`", definition.name())
                })?;

                let Some(source_regions) = regions.get(definition.name()) else {
                    return Ok(true);
                };
                seen_sources.insert(definition.name().to_owned());
                for region in source_regions {
                    anyhow::ensure!(
                        region.end <= sequence.len(),
                        "split interval {}_{} for `{}` exceeds source contig `{}` length {}",
                        region.start,
                        region.end,
                        region.name,
                        definition.name(),
                        sequence.len()
                    );
                    anyhow::ensure!(
                        next_rid < MAX_PACKED_CONTIGS,
                        "alleles supports at most {MAX_PACKED_CONTIGS} contigs"
                    );
                    let split_sequence = sequence[region.start..region.end].to_vec();
                    let trimmed_range = trimmed_sequence_range(split_sequence.len(), trim_length);
                    pending.push_back(TrimmedSequence {
                        rid: next_rid as u32,
                        name: region.name.clone(),
                        sequence: split_sequence,
                        trimmed_range,
                    });
                    next_rid += 1;
                }
                Ok(true)
            })();

            match result {
                Ok(true) => continue,
                Ok(false) => {
                    finished = true;
                    return None;
                }
                Err(error) => {
                    finished = true;
                    return Some(Err(error));
                }
            }
        }
    }))
}

#[derive(Copy, Clone, Debug)]
pub struct MatchRecord {
    pub rid1: u32,
    pub rid2: u32,
    pub rev: u8,
    pub mz1: u32,
    pub mz2: u32,
    pub mz_shared: u32,
    pub similarity: f64,
}

impl PartialEq for MatchRecord {
    fn eq(&self, other: &Self) -> bool {
        self.rid1 == other.rid1 && self.rid2 == other.rid2 && self.rev == other.rev
    }
}

impl Eq for MatchRecord {}

impl Ord for MatchRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rid1
            .cmp(&other.rid1)
            .then(self.rid2.cmp(&other.rid2))
            .then(self.rev.cmp(&other.rev))
    }
}

impl PartialOrd for MatchRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub struct AllelesFasta {
    pub file: String,
    pub contigs: Vec<String>,
    pub contig_lengths: HashMap<String, u64>,
    pub split_regions: Option<String>,
}

#[derive(Debug, Clone, Copy)]
pub struct AllelesOptions {
    pub k: usize,
    pub w: usize,
    pub trim_length: usize,
    pub min_similarity: f64,
    pub max_occurrence: usize,
    pub min_chain: usize,
    pub diff_threshold: f64,
}

impl Default for AllelesOptions {
    fn default() -> Self {
        Self {
            k: 19,
            w: 19,
            trim_length: 0,
            min_similarity: 0.85,
            max_occurrence: 100,
            min_chain: 5,
            diff_threshold: 0.1,
        }
    }
}

impl BaseTable for AllelesFasta {
    fn new(name: &String) -> AllelesFasta {
        AllelesFasta {
            file: name.clone(),
            contigs: Vec::new(),
            contig_lengths: HashMap::new(),
            split_regions: None,
        }
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

impl AllelesFasta {
    pub fn set_split_regions(&mut self, split_regions: Option<String>) {
        self.split_regions = split_regions;
    }

    fn format_allele_table_similarity(record: &MatchRecord) -> String {
        // Match PartigRecords.to_alleletable(): partig first serializes the
        // core similarity with %.6f, then CPhasing raises it to either
        // directional shared-minimizer fraction (capped at one), and finally
        // writes four decimal places without insignificant trailing zeroes.
        let core_similarity = (record.similarity * 1_000_000.0).round() / 1_000_000.0;
        let similarity1 = f64::from(record.mz_shared) / f64::from(record.mz1);
        let similarity2 = f64::from(record.mz_shared) / f64::from(record.mz2);
        let similarity = core_similarity.max(similarity1).max(similarity2).min(1.0);
        let mut formatted = format!("{similarity:.4}");
        while formatted.ends_with('0') {
            formatted.pop();
        }
        if formatted.ends_with('.') {
            formatted.push('0');
        }
        formatted
    }

    fn commit_sequence_metadata(&mut self, mut metadata: Vec<SequenceMetadata>) -> anyResult<()> {
        metadata.sort_unstable_by_key(|record| record.rid);
        for (expected_rid, record) in metadata.iter().enumerate() {
            anyhow::ensure!(
                record.rid as usize == expected_rid,
                "internal FASTA record ordering error: expected rid {expected_rid}, found {}",
                record.rid
            );
        }

        let mut contigs = Vec::with_capacity(metadata.len());
        let mut contig_lengths = HashMap::with_capacity(metadata.len());
        for record in metadata {
            contig_lengths.insert(record.name.clone(), record.length);
            contigs.push(record.name);
        }
        self.contigs = contigs;
        self.contig_lengths = contig_lengths;
        Ok(())
    }

    pub fn seqs(&mut self, trim_length: usize) -> anyResult<Vec<Vec<u8>>> {
        let records = trimmed_fasta_records(&self.file, trim_length)
            .collect::<anyResult<Vec<TrimmedSequence>>>()?;
        let metadata = records
            .iter()
            .map(|record| SequenceMetadata {
                rid: record.rid,
                name: record.name.clone(),
                length: record.trimmed_length() as u64,
            })
            .collect();
        self.commit_sequence_metadata(metadata)?;
        Ok(records
            .into_iter()
            .map(TrimmedSequence::into_trimmed_sequence)
            .collect())
    }

    fn collect_minimizers(
        &mut self,
        k: usize,
        w: usize,
        trim_length: usize,
    ) -> anyResult<MinimizerStorage> {
        // `par_bridge` keeps FASTA parsing sequential while allowing the other
        // workers to sketch records already read. It uses the existing Rayon
        // pool, so this pipeline neither creates another worker pool nor queues
        // the complete FASTA in memory.
        let records: Box<dyn Iterator<Item = anyResult<TrimmedSequence>> + Send> =
            match self.split_regions.as_deref() {
                Some(regions_file) => {
                    Box::new(split_fasta_records(&self.file, regions_file, trim_length)?)
                }
                None => Box::new(trimmed_fasta_records(&self.file, trim_length)),
            };
        let accumulator = records
            .par_bridge()
            .map(|record| -> anyResult<SketchedSequence> {
                let record = record?;
                let length = record.trimmed_length() as u64;
                let minimizers = if use_32bit_minimizer_coordinates(&[length], 1) {
                    MinimizerStorage::U32(sketch32(record.trimmed_sequence(), record.rid, k, w))
                } else {
                    MinimizerStorage::U64(sketch(record.trimmed_sequence(), record.rid, k, w))
                };
                Ok(SketchedSequence {
                    metadata: SequenceMetadata {
                        rid: record.rid,
                        name: record.name,
                        length,
                    },
                    minimizers,
                })
            })
            .try_fold(SketchAccumulator::default, |mut accumulator, record| {
                accumulator.push(record?, k);
                Ok::<SketchAccumulator, anyhow::Error>(accumulator)
            })
            .try_reduce(SketchAccumulator::default, |left, right| {
                Ok(left.merge(right, k))
            })?;

        let sequence_count = accumulator.metadata.len();
        let mut minimizers = accumulator
            .minimizers
            .unwrap_or_else(|| MinimizerStorage::U32(Vec::new()));
        log::info!("Load `{sequence_count}` sequences.");
        if minimizers.is_compact() {
            log::info!("Use 16-byte minimizers with 32-bit coordinates.");
        } else {
            log::info!("Use 24-byte minimizers with 64-bit coordinates.");
        }
        log::info!("Collected {} minimizers.", minimizers.len());
        minimizers.par_sort_unstable();
        self.commit_sequence_metadata(accumulator.metadata)?;
        Ok(minimizers)
    }

    fn generate_canonical_anchors<M, A, F>(
        minimizers: &[M],
        ranges: &[Range<usize>],
        lengths: &[u64],
        k: usize,
        make_anchor: F,
    ) -> Vec<A>
    where
        M: MinimizerRecord,
        A: Send,
        F: Fn(u32, u64, u32, u64, u8) -> A + Send + Sync,
    {
        ranges
            .par_iter()
            .flat_map_iter(|range| {
                let group = &minimizers[range.clone()];
                let make_anchor = &make_anchor;
                (0..group.len()).flat_map(move |i| {
                    ((i + 1)..group.len()).filter_map(move |j| {
                        let left = &group[i];
                        let right = &group[j];
                        if left.rid() == right.rid() {
                            return None;
                        }

                        let rev = left.rev() ^ right.rev();
                        let (source, target) = if left.rid() < right.rid() {
                            (left, right)
                        } else {
                            (right, left)
                        };
                        let target_pos = if rev == 1 {
                            lengths[target.rid() as usize] - (target.pos() + 1 - k as u64) - 1
                        } else {
                            target.pos()
                        };
                        Some(make_anchor(
                            source.rid(),
                            source.pos(),
                            target.rid(),
                            target_pos,
                            rev,
                        ))
                    })
                })
            })
            .collect()
    }

    fn collect_anchors_for<M: MinimizerRecord>(
        &self,
        minimizers: Vec<M>,
        max_occurrence: usize,
        k: usize,
    ) -> anyResult<(AnchorStorage, Vec<usize>, Vec<usize>)> {
        let mut ranges = Vec::new();
        let mut unique_counts = vec![0; self.contigs.len()];
        let mut considered_counts = vec![0; self.contigs.len()];
        let mut start = 0;
        while start < minimizers.len() {
            let mut end = start + 1;
            while end < minimizers.len()
                && minimizers[end].minimizer() == minimizers[start].minimizer()
            {
                end += 1;
            }
            let count = end - start;
            if count == 1 {
                unique_counts[minimizers[start].rid() as usize] += 1;
            } else if count <= max_occurrence {
                for minimizer in &minimizers[start..end] {
                    considered_counts[minimizer.rid() as usize] += 1;
                }
                ranges.push(start..end);
            }
            start = end;
        }

        log::info!(
            "Collected {} non-repetitive minimizer groups.",
            ranges.len()
        );
        let lengths = self
            .contigs
            .iter()
            .map(|name| self.contig_lengths[name])
            .collect::<Vec<_>>();
        let mut anchor = if use_32bit_anchor_coordinates(&lengths) {
            log::info!("Use 32-bit anchor coordinates.");
            AnchorStorage::U32(Self::generate_canonical_anchors(
                &minimizers,
                &ranges,
                &lengths,
                k,
                Anchor32::new,
            ))
        } else {
            log::info!("Use 64-bit anchor coordinates.");
            AnchorStorage::U64(Self::generate_canonical_anchors(
                &minimizers,
                &ranges,
                &lengths,
                k,
                Anchor::new,
            ))
        };
        drop(minimizers);

        anchor.radix_sort_by_key();
        log::info!("Collected {} canonical anchors.", anchor.len());
        Ok((anchor, unique_counts, considered_counts))
    }

    fn collect_anchors(
        &self,
        minimizers: MinimizerStorage,
        max_occurrence: usize,
        k: usize,
    ) -> anyResult<(AnchorStorage, Vec<usize>, Vec<usize>)> {
        match minimizers {
            MinimizerStorage::U32(minimizers) => {
                self.collect_anchors_for(minimizers, max_occurrence, k)
            }
            MinimizerStorage::U64(minimizers) => {
                self.collect_anchors_for(minimizers, max_occurrence, k)
            }
        }
    }

    // longest increasing sequences
    pub fn lis<'a>(&'a self, anchor: &'a Vec<&Anchor>) -> usize {
        let n = anchor.len();
        let mut dp: Vec<usize> = vec![1; n];
        for i in 0..n {
            for j in 0..i {
                if anchor[i].pos2 > anchor[j].pos2 {
                    dp[i] = dp[i].max(dp[j] + 1);
                }
            }
        }

        dp.into_iter().max().unwrap()
    }

    // https://rosettacode.org/wiki/Longest_increasing_subsequence
    fn lis_optimized(positions: &[(u64, u64)], tails: &mut Vec<u64>) -> usize {
        tails.clear();
        if tails.capacity() < positions.len() {
            tails.reserve(positions.len() - tails.capacity());
        }
        for &(_, position) in positions {
            let mut lo = 0;
            let mut hi = tails.len();
            while lo < hi {
                let mid = (lo + hi) / 2;
                if tails[mid] < position {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            if lo == tails.len() {
                tails.push(position);
            } else {
                tails[lo] = position;
            }
        }
        tails.len()
    }

    fn calculate_simularity<A: AnchorRecord>(
        &self,
        anchor: &[A],
        contig_minimizer_counts: &[usize],
        k: usize,
        min_chain: usize,
        min_sim: f64,
    ) -> Vec<MatchRecord> {
        let mut ranges = Vec::new();
        let mut start = 0;
        while start < anchor.len() {
            let key = anchor[start].key();
            let mut end = start + 1;
            while end < anchor.len() && anchor[end].key() == key {
                end += 1;
            }
            if end - start >= min_chain {
                ranges.push(start..end);
            }
            start = end;
        }

        let lengths = self
            .contigs
            .iter()
            .map(|name| self.contig_lengths[name])
            .collect::<Vec<_>>();
        let res = ranges
            .par_iter()
            .map_init(
                || (Vec::<(u64, u64)>::new(), Vec::<u64>::new()),
                |(positions, tails), range| {
                    positions.clear();
                    positions.extend(
                        anchor[range.clone()]
                            .iter()
                            .map(|item| (item.pos1(), item.pos2())),
                    );
                    positions.sort_unstable();
                    let forward_m = Self::lis_optimized(positions, tails);
                    let first = &anchor[range.start];
                    let rid1 = key_rid1(first.key());
                    let rid2 = key_rid2(first.key());
                    let rev = key_rev(first.key());
                    let n1 = contig_minimizer_counts[rid1 as usize];
                    let n2 = contig_minimizer_counts[rid2 as usize];

                    let forward_similarity =
                        (2.0 * forward_m as f64 / (n1 + n2) as f64).powf(1.0 / k as f64);
                    let forward = if forward_m >= min_chain && forward_similarity >= min_sim {
                        MatchRecord {
                            rid1,
                            rid2,
                            rev,
                            mz1: n1 as u32,
                            mz2: n2 as u32,
                            mz_shared: forward_m as u32,
                            similarity: forward_similarity,
                        }
                    } else {
                        // The legacy symmetry pass discarded a directional
                        // match unless its reverse direction also passed.
                        return [None, None];
                    };

                    if rev == 0 {
                        for (source, target) in positions.iter_mut() {
                            std::mem::swap(source, target);
                        }
                    } else {
                        let length1 = lengths[rid1 as usize];
                        let length2 = lengths[rid2 as usize];
                        let span = k as u64;
                        for (source, reverse_target) in positions.iter_mut() {
                            let target = length2 - (*reverse_target + 1 - span) - 1;
                            let reverse_source = length1 - (*source + 1 - span) - 1;
                            *source = target;
                            *reverse_target = reverse_source;
                        }
                    }
                    positions.sort_unstable();
                    let reverse_m = Self::lis_optimized(positions, tails);
                    let reverse_similarity =
                        (2.0 * reverse_m as f64 / (n1 + n2) as f64).powf(1.0 / k as f64);
                    if reverse_m >= min_chain && reverse_similarity >= min_sim {
                        [
                            Some(forward),
                            Some(MatchRecord {
                                rid1: rid2,
                                rid2: rid1,
                                rev,
                                mz1: n2 as u32,
                                mz2: n1 as u32,
                                mz_shared: reverse_m as u32,
                                similarity: reverse_similarity,
                            }),
                        ]
                    } else {
                        [None, None]
                    }
                },
            )
            .flat_map_iter(|records| records.into_iter().flatten())
            .collect::<Vec<MatchRecord>>();

        log::info!("Collected {} matches.", res.len());

        res
    }

    fn filter_matches(matches: &mut Vec<MatchRecord>, min_chain: usize, diff_threshold: f64) {
        matches.par_sort_unstable();
        let mut keep = vec![false; matches.len()];
        let mut rid_ranges = matches
            .last()
            .map(|record| vec![0..0; record.rid1 as usize + 1])
            .unwrap_or_default();
        let mut start = 0;
        while start < matches.len() {
            let rid = matches[start].rid1;
            let mut end = start + 1;
            while end < matches.len() && matches[end].rid1 == rid {
                end += 1;
            }
            let max_shared = matches[start..end]
                .iter()
                .map(|record| record.mz_shared)
                .max()
                .unwrap_or(0);
            for index in start..end {
                let shared = matches[index].mz_shared;
                keep[index] = f64::from(shared) >= f64::from(max_shared) * diff_threshold
                    || shared.saturating_add(min_chain as u32) >= max_shared;
            }
            let mut index = start + 1;
            while index < end {
                if matches[index].rid2 == matches[index - 1].rid2 {
                    if matches[index].mz_shared > matches[index - 1].mz_shared {
                        keep[index - 1] = false;
                    } else {
                        keep[index] = false;
                    }
                }
                index += 1;
            }
            rid_ranges[rid as usize] = start..end;
            start = end;
        }

        for i in 0..matches.len() {
            if !keep[i] {
                continue;
            }
            let record = matches[i];
            let Some(range) = rid_ranges.get(record.rid2 as usize) else {
                keep[i] = false;
                continue;
            };
            let reverse = matches[range.clone()].binary_search_by(|candidate| {
                candidate
                    .rid2
                    .cmp(&record.rid1)
                    .then(candidate.rev.cmp(&record.rev))
            });
            if let Ok(offset) = reverse {
                keep[range.start + offset] = true;
            } else {
                // Preserve the old final symmetry pass defensively if a
                // caller ever provides an unpaired input.
                keep[i] = false;
            }
        }
        let mut cursor = 0;
        matches.retain(|_| {
            let retain = keep[cursor];
            cursor += 1;
            retain
        });
    }

    pub fn run(&mut self, k: usize, w: usize, m: f64, output: &String) {
        let options = AllelesOptions {
            k,
            w,
            min_similarity: m,
            ..AllelesOptions::default()
        };
        self.run_with_options(options, output).unwrap();
    }

    pub fn run_with_options(&mut self, options: AllelesOptions, output: &str) -> anyResult<()> {
        anyhow::ensure!(
            (1..=63).contains(&options.k),
            "k-mer size must be between 1 and 63"
        );
        anyhow::ensure!(
            (1..256).contains(&options.w),
            "window size must be between 1 and 255"
        );
        anyhow::ensure!(
            options.max_occurrence >= 2,
            "maximum occurrence must be at least 2"
        );
        anyhow::ensure!(
            options.min_chain > 0,
            "minimum chain length must be positive"
        );
        anyhow::ensure!(
            (0.0..=1.0).contains(&options.min_similarity),
            "minimum similarity must be between 0 and 1"
        );
        anyhow::ensure!(
            (0.0..=1.0).contains(&options.diff_threshold),
            "difference threshold must be between 0 and 1"
        );

        let start = std::time::Instant::now();
        let minimizers = self.collect_minimizers(options.k, options.w, options.trim_length)?;

        let (anchor, unique_counts, considered_counts) =
            self.collect_anchors(minimizers, options.max_occurrence, options.k)?;

        let mut matches = match &anchor {
            AnchorStorage::U32(anchor) => self.calculate_simularity(
                anchor,
                &considered_counts,
                options.k,
                options.min_chain,
                options.min_similarity,
            ),
            AnchorStorage::U64(anchor) => self.calculate_simularity(
                anchor,
                &considered_counts,
                options.k,
                options.min_chain,
                options.min_similarity,
            ),
        };
        drop(anchor);
        Self::filter_matches(&mut matches, options.min_chain, options.diff_threshold);
        matches.par_sort_unstable();

        let mut writer = common_writer(output);

        for (rid, contig) in self.contigs.iter().enumerate() {
            writeln!(
                writer,
                "#{} {} {} {}",
                contig, self.contig_lengths[contig], considered_counts[rid], unique_counts[rid]
            )?;
        }

        for (i, record) in matches.iter().enumerate() {
            let contig1 = self.contigs.get(record.rid1 as usize).unwrap();
            let contig2 = self.contigs.get(record.rid2 as usize).unwrap();
            let strand = if record.rev == 1 { "-1" } else { "1" };
            let similarity = Self::format_allele_table_similarity(record);
            writer.write_all(
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    i + 1,
                    i + 1,
                    contig1,
                    contig2,
                    record.mz1,
                    record.mz2,
                    record.mz_shared,
                    similarity,
                    strand
                )
                .as_bytes(),
            )?;
        }
        writer.flush()?;
        log::info!("Successful output allele table to `{}`.", output);
        log::info!(
            "Allele detection finished in {:.3}s.",
            start.elapsed().as_secs_f64()
        );
        Ok(())
    }
}

#[cfg(test)]
mod compact_alleles_tests {
    use super::*;

    fn alleles_with_lengths(lengths: &[u64]) -> AllelesFasta {
        let contigs = (0..lengths.len())
            .map(|index| format!("s{index}"))
            .collect::<Vec<_>>();
        let contig_lengths = contigs
            .iter()
            .cloned()
            .zip(lengths.iter().copied())
            .collect();
        AllelesFasta {
            file: String::new(),
            contigs,
            contig_lengths,
            split_regions: None,
        }
    }

    fn match_record(rid1: u32, rid2: u32, rev: u8, shared: u32) -> MatchRecord {
        MatchRecord {
            rid1,
            rid2,
            rev,
            mz1: 20,
            mz2: 20,
            mz_shared: shared,
            similarity: f64::from(shared) / 20.0,
        }
    }

    fn assert_symmetric(matches: &[MatchRecord]) {
        let keys = matches
            .iter()
            .map(|record| (record.rid1, record.rid2, record.rev))
            .collect::<HashSet<_>>();
        assert!(matches.iter().all(|record| keys.contains(&(
            record.rid2,
            record.rid1,
            record.rev
        ))));
    }

    #[test]
    fn compact_selector_and_wide_reverse_coordinate_preserve_boundaries() {
        let compact_limit = u64::from(u32::MAX) + 1;
        assert!(use_32bit_minimizer_coordinates(
            &[compact_limit],
            MAX_PACKED_CONTIGS
        ));
        assert!(!use_32bit_minimizer_coordinates(
            &[compact_limit],
            MAX_PACKED_CONTIGS + 1
        ));
        assert!(!use_32bit_minimizer_coordinates(&[compact_limit + 1], 1));

        let lengths = [100, compact_limit + 1];
        let alleles = alleles_with_lengths(&lengths);
        let minimizers = vec![
            MinimizerData::from_parts(7, 0, 18, 0, 19),
            MinimizerData::from_parts(7, 1, 18, 1, 19),
        ];
        let (anchors, _, considered) = alleles.collect_anchors_for(minimizers, 5, 19).unwrap();
        assert_eq!(considered, [1, 1]);
        let AnchorStorage::U64(anchors) = anchors else {
            panic!("a contig longer than 4 GiB must use wide anchors");
        };
        assert_eq!(anchors.len(), 1);
        assert_eq!(anchors[0].pos1, 18);
        assert_eq!(anchors[0].pos2, compact_limit);

        let boundary_alleles = alleles_with_lengths(&[100, compact_limit]);
        let compact = vec![
            MinimizerData32::from_parts(7, 0, 6, 0, 7),
            MinimizerData32::from_parts(7, 1, 6, 1, 7),
            MinimizerData32::from_parts(9, 0, 10, 0, 7),
            MinimizerData32::from_parts(9, 1, u32::MAX.into(), 1, 7),
        ];
        let (anchors, _, _) = boundary_alleles.collect_anchors_for(compact, 5, 7).unwrap();
        let AnchorStorage::U64(anchors) = anchors else {
            panic!("a 4 GiB contig requires wide anchors");
        };
        let positions = anchors
            .iter()
            .map(|anchor| (anchor.pos1, anchor.pos2))
            .collect::<HashSet<_>>();
        assert_eq!(
            positions,
            HashSet::from([(6, u64::from(u32::MAX)), (10, 6)])
        );
    }

    #[test]
    fn compact_batches_promote_losslessly_in_either_merge_order() {
        let compact = vec![
            MinimizerData32::from_parts(7, 0, 6, 0, 7),
            MinimizerData32::from_parts(11, 1, 9, 1, 7),
        ];
        let wide = vec![MinimizerData::from_parts(
            13,
            2,
            u64::from(u32::MAX) + 1,
            1,
            7,
        )];
        let snapshot = |storage: MinimizerStorage| {
            let MinimizerStorage::U64(values) = storage else {
                panic!("a mixed compact/wide merge must promote to wide storage");
            };
            let mut records = values
                .into_iter()
                .map(|item| {
                    assert_eq!(item.info.span, 7);
                    (item.minimizer(), item.rid(), item.pos(), item.rev())
                })
                .collect::<Vec<_>>();
            records.sort_unstable();
            records
        };
        let compact_first = snapshot(
            MinimizerStorage::U32(compact.clone()).merge(MinimizerStorage::U64(wide.clone()), 7),
        );
        let wide_first =
            snapshot(MinimizerStorage::U64(wide).merge(MinimizerStorage::U32(compact), 7));
        assert_eq!(compact_first, wide_first);
        assert_eq!(
            compact_first,
            [
                (7, 0, 6, 0),
                (11, 1, 9, 1),
                (13, 2, u64::from(u32::MAX) + 1, 1),
            ]
        );
    }

    #[test]
    fn compact_and_wide_minimizers_generate_identical_anchors() {
        let alleles = alleles_with_lengths(&[200, 200, 200]);
        let logical = [
            (7, 0, 20, 0),
            (7, 1, 22, 1),
            (9, 0, 40, 1),
            (9, 1, 42, 1),
            (9, 2, 44, 0),
        ];
        let wide = logical
            .iter()
            .map(|&(hash, rid, pos, rev)| MinimizerData::from_parts(hash, rid, pos, rev, 7))
            .collect();
        let compact = logical
            .iter()
            .map(|&(hash, rid, pos, rev)| MinimizerData32::from_parts(hash, rid, pos, rev, 7))
            .collect();

        let (wide_anchors, wide_unique, wide_considered) =
            alleles.collect_anchors_for(wide, 5, 7).unwrap();
        let (compact_anchors, compact_unique, compact_considered) =
            alleles.collect_anchors_for(compact, 5, 7).unwrap();
        assert_eq!(wide_unique, compact_unique);
        assert_eq!(wide_considered, compact_considered);

        let (AnchorStorage::U32(mut wide_anchors), AnchorStorage::U32(mut compact_anchors)) =
            (wide_anchors, compact_anchors)
        else {
            panic!("small coordinates must use compact anchors");
        };
        let order = |left: &Anchor32, right: &Anchor32| {
            left.key
                .cmp(&right.key)
                .then(left.pos1.cmp(&right.pos1))
                .then(left.pos2.cmp(&right.pos2))
        };
        wide_anchors.sort_unstable_by(order);
        compact_anchors.sort_unstable_by(order);
        assert_eq!(wide_anchors, compact_anchors);

        let snapshot = |matches: Vec<MatchRecord>| {
            matches
                .into_iter()
                .map(|record| {
                    (
                        record.rid1,
                        record.rid2,
                        record.rev,
                        record.mz1,
                        record.mz2,
                        record.mz_shared,
                        record.similarity.to_bits(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let wide_matches = alleles.calculate_simularity(&wide_anchors, &wide_considered, 7, 1, 0.0);
        let compact_matches =
            alleles.calculate_simularity(&compact_anchors, &compact_considered, 7, 1, 0.0);
        assert_eq!(snapshot(wide_matches), snapshot(compact_matches));
    }

    #[test]
    fn asymmetric_directional_candidate_is_not_emitted() {
        let alleles = alleles_with_lengths(&[100, 100]);
        let key = anchor_key(0, 1, 0);
        let one_way = [
            Anchor32 {
                key,
                pos1: 0,
                pos2: 0,
            },
            Anchor32 {
                key,
                pos1: 0,
                pos2: 1,
            },
        ];
        assert!(
            alleles
                .calculate_simularity(&one_way, &[2, 2], 19, 2, 0.0)
                .is_empty()
        );

        let both_ways = [
            Anchor32 {
                key,
                pos1: 0,
                pos2: 0,
            },
            Anchor32 {
                key,
                pos1: 1,
                pos2: 1,
            },
        ];
        let matches = alleles.calculate_simularity(&both_ways, &[2, 2], 19, 2, 0.0);
        assert_eq!(matches.len(), 2);
        assert_eq!((matches[0].rid1, matches[0].rid2), (0, 1));
        assert_eq!((matches[1].rid1, matches[1].rid2), (1, 0));
    }

    #[test]
    fn filtering_propagates_keep_decisions_to_reverse_matches() {
        let weak_pair = [match_record(0, 1, 0, 5), match_record(1, 0, 0, 5)];
        let strong_pair = [match_record(0, 2, 0, 10), match_record(2, 0, 0, 10)];
        let mut one_direction_passes = weak_pair.into_iter().chain(strong_pair).collect::<Vec<_>>();
        AllelesFasta::filter_matches(&mut one_direction_passes, 1, 0.8);
        assert_eq!(one_direction_passes.len(), 4);
        assert_symmetric(&one_direction_passes);

        let other_strong_pair = [match_record(1, 3, 0, 10), match_record(3, 1, 0, 10)];
        let mut neither_direction_passes = weak_pair
            .into_iter()
            .chain(strong_pair)
            .chain(other_strong_pair)
            .collect::<Vec<_>>();
        AllelesFasta::filter_matches(&mut neither_direction_passes, 1, 0.8);
        assert_eq!(neither_direction_passes.len(), 4);
        assert_symmetric(&neither_direction_passes);
        assert!(
            neither_direction_passes
                .iter()
                .all(|record| (record.rid1, record.rid2) != (0, 1)
                    && (record.rid1, record.rid2) != (1, 0))
        );
    }

    #[test]
    fn filtering_preserves_cross_direction_strand_winners_and_ties() {
        let mut orphan = vec![match_record(0, 1, 0, 10)];
        AllelesFasta::filter_matches(&mut orphan, 1, 0.0);
        assert!(orphan.is_empty());

        let mut wrong_strand_reverse = vec![match_record(0, 1, 0, 10), match_record(1, 0, 1, 10)];
        AllelesFasta::filter_matches(&mut wrong_strand_reverse, 1, 0.0);
        assert!(wrong_strand_reverse.is_empty());

        let mut cross_winners = vec![
            match_record(0, 1, 0, 100),
            match_record(0, 1, 1, 90),
            match_record(1, 0, 0, 80),
            match_record(1, 0, 1, 110),
        ];
        AllelesFasta::filter_matches(&mut cross_winners, 1, 0.0);
        assert_eq!(cross_winners.len(), 4);
        assert_symmetric(&cross_winners);

        let mut ties = vec![
            match_record(0, 1, 0, 100),
            match_record(0, 1, 1, 100),
            match_record(1, 0, 0, 100),
            match_record(1, 0, 1, 100),
        ];
        AllelesFasta::filter_matches(&mut ties, 1, 0.0);
        assert_eq!(ties.len(), 2);
        assert!(ties.iter().all(|record| record.rev == 0));
        assert_symmetric(&ties);
    }
}
