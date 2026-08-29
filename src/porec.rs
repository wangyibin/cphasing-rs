#![allow(unused)]
#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(unused_variables, unused_assignments)]
use anyhow::{Context, Result as anyResult, anyhow, bail};
use coitrees::{COITree, IntervalNode, IntervalTree};
use crossbeam_channel::{Receiver, Sender, bounded};
use itertools::{Combinations, Itertools};
use rand::prelude::*;
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering as AtomOrdering};
use std::sync::{Arc, Mutex};
use std::thread;

use crate::bed::{Bed3, Bed4};
use crate::core::{BaseTable, ChromSize, ChromSizeRecord, binify};
use crate::core::{common_reader, common_writer, coverage_integrals_at};
use crate::paf::PAFLine;
use crate::pairs::{PairHeader, PairRecord};
use crate::pqs::{
    copy_cn_info, merge_cn_info, write_broken_cn_info, write_materialized_cn_map,
    write_remapped_cn_info,
};

enum PosValue {
    U32(u32),
    U64(u64),
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PoreCRecordPlus {
    pub read_idx: u64,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub query_strand: char,
    pub target: String,
    #[serde(skip_serializing)]
    pub target_length: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub mapq: u8,
    pub identity: f32,
    pub filter_reason: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct PoreCRecord {
    pub read_idx: u64,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub query_strand: char,
    pub target: String,
    pub target_start: u64,
    pub target_end: u64,
    pub mapq: u8,
    pub identity: f32,
    pub filter_reason: String,
}

pub trait ConcatPqsRecord: Send + 'static {
    fn read_idx(&self) -> u64;
    fn query_length(&self) -> u32;
    fn query_start(&self) -> u32;
    fn query_end(&self) -> u32;
    fn query_strand(&self) -> char;
    fn target(&self) -> &str;
    fn target_start(&self) -> u64;
    fn target_end(&self) -> u64;
    fn mapq(&self) -> u8;
    fn identity(&self) -> f32;
    fn filter_reason(&self) -> &str;
}

impl ConcatPqsRecord for PoreCRecord {
    fn read_idx(&self) -> u64 {
        self.read_idx
    }

    fn query_length(&self) -> u32 {
        self.query_length
    }

    fn query_start(&self) -> u32 {
        self.query_start
    }

    fn query_end(&self) -> u32 {
        self.query_end
    }

    fn query_strand(&self) -> char {
        self.query_strand
    }

    fn target(&self) -> &str {
        &self.target
    }

    fn target_start(&self) -> u64 {
        self.target_start
    }

    fn target_end(&self) -> u64 {
        self.target_end
    }

    fn mapq(&self) -> u8 {
        self.mapq
    }

    fn identity(&self) -> f32 {
        self.identity
    }

    fn filter_reason(&self) -> &str {
        &self.filter_reason
    }
}

impl PartialOrd for PoreCRecordPlus {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.query_start.partial_cmp(&other.query_start)
    }
}

impl PartialEq for PoreCRecordPlus {
    fn eq(&self, other: &Self) -> bool {
        self.query_start == other.query_start
    }
}

impl PoreCRecordPlus {
    pub fn from_paf_record(
        record: PAFLine,
        read_idx: u64,
        identity: f32,
        filter_reason: String,
    ) -> Self {
        PoreCRecordPlus {
            read_idx,
            query_length: record.query_length,
            query_start: record.query_start,
            query_end: record.query_end,
            query_strand: record.query_strand,
            target: record.target,
            target_length: record.target_length,
            target_start: record.target_start,
            target_end: record.target_end,
            mapq: record.mapq,
            identity,
            filter_reason,
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.read_idx,
            self.query_length,
            self.query_start,
            self.query_end,
            self.query_strand,
            self.target,
            self.target_start,
            self.target_end,
            self.mapq,
            self.identity,
            self.filter_reason,
        )
    }

    pub fn write_to(&self, buf: &mut String) {
        use std::fmt::Write as _;

        let _ = write!(
            buf,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.read_idx,
            self.query_length,
            self.query_start,
            self.query_end,
            self.query_strand,
            self.target,
            self.target_start,
            self.target_end,
            self.mapq,
            self.identity,
            self.filter_reason
        );
    }

    pub fn is_in_regions(&self, interval_hash: &HashMap<String, Lapper<usize, u8>>) -> bool {
        let is_in_regions: bool = if let Some(interval) = interval_hash.get(&self.target) {
            let iv_start =
                interval.count((self.target_start - 1) as usize, self.target_start as usize);
            let iv_end = interval.count((self.target_end - 1) as usize, self.target_end as usize);
            if iv_start > 0 && iv_end > 0 {
                true
            } else {
                false
            }
        } else {
            false
        };
        is_in_regions
    }
}

impl PoreCRecord {
    pub fn from_paf_record(
        record: PAFLine,
        read_idx: u64,
        identity: f32,
        filter_reason: String,
    ) -> PoreCRecord {
        PoreCRecord {
            read_idx,
            query_length: record.query_length,
            query_start: record.query_start,
            query_end: record.query_end,
            query_strand: record.query_strand,
            target: record.target,
            target_start: record.target_start,
            target_end: record.target_end,
            mapq: record.mapq,
            identity,
            filter_reason,
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.read_idx,
            self.query_length,
            self.query_start,
            self.query_end,
            self.query_strand,
            self.target,
            self.target_start,
            self.target_end,
            self.mapq,
            self.identity,
            self.filter_reason,
        )
    }

    pub fn is_in_regions(&self, interval_hash: &HashMap<String, Lapper<usize, u8>>) -> bool {
        let is_in_regions: bool = if let Some(interval) = interval_hash.get(&self.target) {
            let iv_start =
                interval.count((self.target_start - 1) as usize, self.target_start as usize);
            let iv_end = interval.count((self.target_end - 1) as usize, self.target_end as usize);
            if iv_start > 0 && iv_end > 0 {
                true
            } else {
                false
            }
        } else {
            false
        };
        is_in_regions
    }
}

impl PartialOrd for PoreCRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.query_start.partial_cmp(&other.query_start)
    }
}

impl PartialEq for PoreCRecord {
    fn eq(&self, other: &Self) -> bool {
        self.query_start == other.query_start
    }
}

#[derive(Debug, Default)]
pub struct Concatemer {
    pub records: Vec<PoreCRecord>,
}

#[derive(Debug)]
struct PairAnchor {
    target: Arc<str>,
    target_start: u64,
    target_end: u64,
    query_strand: char,
    mapq: u8,
}

#[derive(Debug, Default)]
struct PairConcatemer {
    records: Vec<PairAnchor>,
}

impl PairConcatemer {
    fn clear(&mut self) {
        self.records.clear();
    }

    fn count(&self) -> usize {
        self.records.len()
    }

    fn sort(&mut self) {
        self.records
            .sort_unstable_by(|left, right| match left.target.cmp(&right.target) {
                Ordering::Equal => left.target_start.cmp(&right.target_start),
                other => other,
            });
    }
}

impl Concatemer {
    pub fn new() -> Concatemer {
        Concatemer {
            records: Vec::new(),
        }
    }

    pub fn push(&mut self, pcr: PoreCRecord) {
        self.records.push(pcr);
    }

    pub fn clear(&mut self) {
        self.records.clear();
    }

    pub fn count(&self) -> usize {
        self.records.len()
    }

    pub fn sort(&mut self) {
        // self.records.sort_by(| a, b | a.partial_cmp(&b).unwrap());
        // self.records.sort_unstable_by_key(|x| (x.target, x.target_start));
        self.records
            .sort_unstable_by(|a, b| match a.target.cmp(&b.target) {
                std::cmp::Ordering::Equal => a.target_start.cmp(&b.target_start),
                other => other,
            });
    }

    // pub fn decompose(&mut self) -> Combinations<std::vec::IntoIter<PoreCRecord>> {
    //     let r = self.records.clone();
    //     r.into_iter().combinations(2)

    // }
    pub fn decompose(&self) -> Combinations<std::slice::Iter<'_, PoreCRecord>> {
        self.records.iter().combinations(2)
    }
}

#[derive(Debug, Clone)]
pub struct ConcatemerSummary {
    pub summary: HashMap<u32, u64>,
}

impl ConcatemerSummary {
    pub fn new() -> ConcatemerSummary {
        ConcatemerSummary {
            summary: HashMap::<u32, u64>::new(),
        }
    }

    pub fn count(&mut self, concatemer: &Concatemer) {
        self.count_order(concatemer.count());
    }

    fn count_order(&mut self, order: usize) {
        let concatemer_count: u32 = order.try_into().unwrap();
        *self.summary.entry(concatemer_count).or_insert(0) += 1;
    }

    pub fn to_string(&self) -> String {
        let mut vec = self.summary.iter().collect::<Vec<_>>();
        vec.sort_by_key(|(key, _)| *key);

        vec.iter()
            .map(|(key, value)| format!("{}\t{}", key, value))
            .collect::<Vec<_>>()
            .join("\n")
    }

    pub fn save(&self, output: &String) {
        let mut wtr = common_writer(output);

        let result: String = self.to_string();

        wtr.write_all(result.as_bytes()).unwrap();
        log::info!("Successful output summary of concatemer `{}`", output);
    }
    // pub fn collapse(&self) -> HashMap<&str, u64> {

    // }
}

#[derive(Debug)]
pub struct PoreCTable {
    file: String,
}

const CONCAT_PQS_README: &str = r#"
# concat.pqs format
The _contigsizes file contains contig lengths.
The _metadata file declares the alignment-level concat schema.
Parquet shards contain complete concatemers and never split one read_idx
between files. q0 contains every alignment; q1 contains alignments with
mapping_quality >= 1.

concat.pqs/
|-- _contigsizes
|-- _metadata
|-- _metadata_counts
|-- _readme
|-- q0/*.parquet
|-- q1/*.parquet
"#;

fn is_concat_pqs(path: &Path) -> bool {
    if !path.is_dir() || !path.join("q0").is_dir() {
        return false;
    }
    std::fs::read_to_string(path.join("_metadata"))
        .map(|metadata| {
            metadata.contains("'format': 'concat'")
                || metadata.contains("'format': 'porec'")
                || metadata.contains("\"format\": \"concat\"")
                || metadata.contains("\"format\": \"porec\"")
        })
        .unwrap_or(false)
}

fn concat_pqs_files(path: &Path) -> anyResult<Vec<std::path::PathBuf>> {
    let q0 = path.join("q0");
    let mut files = std::fs::read_dir(&q0)
        .with_context(|| format!("failed to read concat PQS directory {}", q0.display()))?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "parquet")
        })
        .collect::<Vec<_>>();
    files.sort_by(|left, right| {
        let numeric_stem = |path: &Path| {
            path.file_stem()
                .and_then(|stem| stem.to_str())
                .and_then(|stem| stem.parse::<u64>().ok())
        };
        numeric_stem(left)
            .cmp(&numeric_stem(right))
            .then_with(|| left.cmp(right))
    });
    if files.is_empty() {
        bail!(
            "concat PQS input contains no q0 Parquet shards: {}",
            q0.display()
        );
    }
    Ok(files)
}

fn concat_frame_to_text(frame: &polars::prelude::DataFrame) -> anyResult<Vec<u8>> {
    use polars::prelude::*;

    let read_idx = frame.column("read_idx")?.as_materialized_series().u64()?;
    let read_length = frame
        .column("read_length")?
        .as_materialized_series()
        .u32()?;
    let read_start = frame.column("read_start")?.as_materialized_series().u32()?;
    let read_end = frame.column("read_end")?.as_materialized_series().u32()?;
    let strand_series = frame
        .column("strand")?
        .as_materialized_series()
        .cast(&DataType::String)?;
    let strand = strand_series.str()?;
    let chrom_series = frame
        .column("chrom")?
        .as_materialized_series()
        .cast(&DataType::String)?;
    let chrom = chrom_series.str()?;
    let start_series = frame
        .column("start")?
        .as_materialized_series()
        .cast(&DataType::UInt64)?;
    let start = start_series.u64()?;
    let end_series = frame
        .column("end")?
        .as_materialized_series()
        .cast(&DataType::UInt64)?;
    let end = end_series.u64()?;
    let mapping_quality = frame
        .column("mapping_quality")?
        .as_materialized_series()
        .u8()?;
    let identity = frame.column("identity")?.as_materialized_series().f32()?;
    let filter_series = frame
        .column("filter_reason")?
        .as_materialized_series()
        .cast(&DataType::String)?;
    let filter_reason = filter_series.str()?;

    let mut output = Vec::with_capacity(frame.height().saturating_mul(96));
    for row in 0..frame.height() {
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            read_idx
                .get(row)
                .ok_or_else(|| anyhow!("null read_idx at row {row}"))?,
            read_length
                .get(row)
                .ok_or_else(|| anyhow!("null read_length at row {row}"))?,
            read_start
                .get(row)
                .ok_or_else(|| anyhow!("null read_start at row {row}"))?,
            read_end
                .get(row)
                .ok_or_else(|| anyhow!("null read_end at row {row}"))?,
            strand
                .get(row)
                .ok_or_else(|| anyhow!("null strand at row {row}"))?,
            chrom
                .get(row)
                .ok_or_else(|| anyhow!("null chrom at row {row}"))?,
            start
                .get(row)
                .ok_or_else(|| anyhow!("null start at row {row}"))?,
            end.get(row)
                .ok_or_else(|| anyhow!("null end at row {row}"))?,
            mapping_quality
                .get(row)
                .ok_or_else(|| anyhow!("null mapping_quality at row {row}"))?,
            identity
                .get(row)
                .ok_or_else(|| anyhow!("null identity at row {row}"))?,
            filter_reason
                .get(row)
                .ok_or_else(|| anyhow!("null filter_reason at row {row}"))?,
        )?;
    }
    Ok(output)
}

fn concat_frame_to_records(frame: &polars::prelude::DataFrame) -> anyResult<Vec<PoreCRecord>> {
    use polars::prelude::*;

    let read_idx = frame.column("read_idx")?.as_materialized_series().u64()?;
    let read_length = frame
        .column("read_length")?
        .as_materialized_series()
        .u32()?;
    let read_start = frame.column("read_start")?.as_materialized_series().u32()?;
    let read_end = frame.column("read_end")?.as_materialized_series().u32()?;
    let strand_series = frame
        .column("strand")?
        .as_materialized_series()
        .cast(&DataType::String)?;
    let strand = strand_series.str()?;
    let chrom_series = frame
        .column("chrom")?
        .as_materialized_series()
        .cast(&DataType::String)?;
    let chrom = chrom_series.str()?;
    let start_series = frame
        .column("start")?
        .as_materialized_series()
        .cast(&DataType::UInt64)?;
    let start = start_series.u64()?;
    let end_series = frame
        .column("end")?
        .as_materialized_series()
        .cast(&DataType::UInt64)?;
    let end = end_series.u64()?;
    let mapping_quality = frame
        .column("mapping_quality")?
        .as_materialized_series()
        .u8()?;
    let identity = frame.column("identity")?.as_materialized_series().f32()?;
    let filter_series = frame
        .column("filter_reason")?
        .as_materialized_series()
        .cast(&DataType::String)?;
    let filter_reason = filter_series.str()?;

    let mut records = Vec::with_capacity(frame.height());
    for row in 0..frame.height() {
        let strand = strand
            .get(row)
            .ok_or_else(|| anyhow!("null strand at row {row}"))?;
        records.push(PoreCRecord {
            read_idx: read_idx
                .get(row)
                .ok_or_else(|| anyhow!("null read_idx at row {row}"))?,
            query_length: read_length
                .get(row)
                .ok_or_else(|| anyhow!("null read_length at row {row}"))?,
            query_start: read_start
                .get(row)
                .ok_or_else(|| anyhow!("null read_start at row {row}"))?,
            query_end: read_end
                .get(row)
                .ok_or_else(|| anyhow!("null read_end at row {row}"))?,
            query_strand: strand
                .chars()
                .next()
                .ok_or_else(|| anyhow!("empty strand at row {row}"))?,
            target: chrom
                .get(row)
                .ok_or_else(|| anyhow!("null chrom at row {row}"))?
                .to_string(),
            target_start: start
                .get(row)
                .ok_or_else(|| anyhow!("null start at row {row}"))?,
            target_end: end
                .get(row)
                .ok_or_else(|| anyhow!("null end at row {row}"))?,
            mapq: mapping_quality
                .get(row)
                .ok_or_else(|| anyhow!("null mapping_quality at row {row}"))?,
            identity: identity
                .get(row)
                .ok_or_else(|| anyhow!("null identity at row {row}"))?,
            filter_reason: filter_reason
                .get(row)
                .ok_or_else(|| anyhow!("null filter_reason at row {row}"))?
                .to_string(),
        });
    }
    Ok(records)
}

struct ConcatPqsTextReader {
    files: Vec<std::path::PathBuf>,
    next_file: usize,
    buffer: Cursor<Vec<u8>>,
    shard_scoped_read_idx: bool,
    next_logical_read_idx: u64,
}

impl ConcatPqsTextReader {
    fn new(path: &Path) -> anyResult<Self> {
        if !is_concat_pqs(path) {
            bail!(
                "PQS directory {} is not an alignment-level concat.pqs",
                path.display()
            );
        }
        let metadata = std::fs::read_to_string(path.join("_metadata"))?;
        let shard_scoped_read_idx = metadata.contains("'read_idx_scope': 'shard'")
            || metadata.contains("\"read_idx_scope\": \"shard\"");
        polars::enable_string_cache();
        Ok(Self {
            files: concat_pqs_files(path)?,
            next_file: 0,
            buffer: Cursor::new(Vec::new()),
            shard_scoped_read_idx,
            next_logical_read_idx: 1,
        })
    }

    fn load_next_file(&mut self) -> io::Result<bool> {
        use polars::prelude::*;

        while self.next_file < self.files.len() {
            let path = &self.files[self.next_file];
            self.next_file += 1;
            let mut frame = ParquetReader::new(File::open(path)?)
                .finish()
                .map_err(|error| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "failed to read concat PQS shard {}: {error}",
                            path.display()
                        ),
                    )
                })?;
            if self.shard_scoped_read_idx {
                let logical_ids = {
                    let raw_ids = frame
                        .column("read_idx")
                        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?
                        .as_materialized_series()
                        .u64()
                        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
                    let mut previous = None;
                    let mut logical_id = None;
                    let mut logical_ids = Vec::with_capacity(frame.height());
                    for raw_id in raw_ids.into_no_null_iter() {
                        if previous != Some(raw_id) {
                            previous = Some(raw_id);
                            logical_id = Some(self.next_logical_read_idx);
                            self.next_logical_read_idx =
                                self.next_logical_read_idx.checked_add(1).ok_or_else(|| {
                                    io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        "logical read_idx overflow while loading linked concat PQS",
                                    )
                                })?;
                        }
                        logical_ids.push(logical_id.expect("first row starts a logical read"));
                    }
                    logical_ids
                };
                frame
                    .replace("read_idx", Series::new("read_idx".into(), logical_ids))
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
            }
            let data = concat_frame_to_text(&frame).map_err(|error| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid concat PQS shard {}: {error}", path.display()),
                )
            })?;
            self.buffer = Cursor::new(data);
            if !self.buffer.get_ref().is_empty() {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

impl Read for ConcatPqsTextReader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let size = available.len().min(output.len());
        output[..size].copy_from_slice(&available[..size]);
        self.consume(size);
        Ok(size)
    }
}

impl BufRead for ConcatPqsTextReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.buffer.position() as usize >= self.buffer.get_ref().len()
            && !self.load_next_file()?
        {
            return Ok(&[]);
        }
        self.buffer.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        self.buffer.consume(amount);
    }
}

fn open_porec_reader(path: &str) -> anyResult<Box<dyn BufRead + Send + 'static>> {
    let path_obj = Path::new(path);
    if path_obj.is_dir() {
        log::info!("Load concat PQS `{}`", path);
        return Ok(Box::new(ConcatPqsTextReader::new(path_obj)?));
    }
    Ok(common_reader(path))
}

fn concat_pqs_metadata(chunksize: usize, position_type: &str) -> String {
    let creation_time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
    let position_dtype = position_type.to_lowercase();
    format!(
        "{{'format-version': '0.2.0',\n\
 'format-url': 'https://github.com/wangyibin/CPhasing',\n\
 'creation_time': '{creation_time}',\n\
 'is_pqs': True,\n\
 'format': 'concat',\n\
 'record_type': 'porec_alignment',\n\
 'chunksize': {chunksize},\n\
 'partition_key': 'read_idx',\n\
 'read_idx_scope': 'global',\n\
 'complete_groups': True,\n\
 'is_with_mapq': True,\n\
 'q1_min_mapq': 1,\n\
 'columns': ['read_idx', 'read_length', 'read_start', 'read_end',\n\
             'strand', 'chrom', 'start', 'end',\n\
             'mapping_quality', 'identity', 'filter_reason'],\n\
 'dtypes': {{'read_idx': 'uint64',\n\
             'read_length': 'uint32',\n\
             'read_start': 'uint32',\n\
             'read_end': 'uint32',\n\
             'strand': 'U',\n\
             'chrom': 'U',\n\
             'start': '{position_dtype}',\n\
             'end': '{position_dtype}',\n\
             'mapping_quality': 'uint8',\n\
             'identity': 'float32',\n\
             'filter_reason': 'U'}},\n\
 'schema': {{'read_idx': UInt64,\n\
            'read_length': UInt32,\n\
            'read_start': UInt32,\n\
            'read_end': UInt32,\n\
            'strand': Categorical,\n\
            'chrom': Categorical,\n\
            'start': {position_type},\n\
            'end': {position_type},\n\
            'mapping_quality': UInt8,\n\
            'identity': pl.Float32,\n\
            'filter_reason': Categorical}}}}\n"
    )
}

fn parse_concat_porec_line(line: &str, line_number: usize) -> anyResult<Option<PoreCRecord>> {
    let line = line.trim_end_matches(['\r', '\n']);
    if line.trim().is_empty() || line.starts_with('#') {
        return Ok(None);
    }
    let mut fields = [""; 11];
    let mut field_count = 0usize;
    for field in line.split('\t') {
        if field_count == fields.len() {
            field_count += 1;
            break;
        }
        fields[field_count] = field;
        field_count += 1;
    }
    if field_count != fields.len() {
        bail!(
            "malformed Pore-C input at line {line_number}: expected 11 tab-separated fields, found {}",
            field_count
        );
    }
    let parse = |index: usize, name: &str| -> anyResult<&str> {
        fields
            .get(index)
            .copied()
            .ok_or_else(|| anyhow!("missing {name} at Pore-C line {line_number}"))
    };
    let parse_number = |index: usize, name: &str| -> anyResult<&str> { parse(index, name) };
    let strand_text = parse(4, "strand")?;
    let mut strand_chars = strand_text.chars();
    let query_strand = strand_chars
        .next()
        .ok_or_else(|| anyhow!("empty strand at Pore-C line {line_number}"))?;
    if strand_chars.next().is_some() || !matches!(query_strand, '+' | '-') {
        bail!("invalid strand {strand_text:?} at Pore-C line {line_number}");
    }

    macro_rules! number {
        ($index:expr, $name:literal, $type:ty) => {{
            parse_number($index, $name)?
                .parse::<$type>()
                .with_context(|| {
                    format!(
                        "invalid {} {:?} at Pore-C line {}",
                        $name, fields[$index], line_number
                    )
                })?
        }};
    }

    Ok(Some(PoreCRecord {
        read_idx: number!(0, "read_idx", u64),
        query_length: number!(1, "read_length", u32),
        query_start: number!(2, "read_start", u32),
        query_end: number!(3, "read_end", u32),
        query_strand,
        target: parse(5, "chrom")?.to_string(),
        target_start: number!(6, "start", u64),
        target_end: number!(7, "end", u64),
        mapq: number!(8, "mapping_quality", u8),
        identity: number!(9, "identity", f32),
        filter_reason: parse(10, "filter_reason")?.to_string(),
    }))
}

fn parse_concat_porec_batch(batch: &str, first_line_number: usize) -> anyResult<Vec<PoreCRecord>> {
    let mut records = Vec::with_capacity(batch.len() / 96);
    for (offset, line) in batch.lines().enumerate() {
        if let Some(record) = parse_concat_porec_line(line, first_line_number + offset)? {
            records.push(record);
        }
    }
    Ok(records)
}

fn push_complete_concatemer<R>(
    sender: &crossbeam_channel::Sender<(usize, Vec<R>)>,
    shard: &mut Vec<R>,
    concatemer: &mut Vec<R>,
    chunksize: usize,
    shard_id: &mut usize,
) -> anyResult<()> {
    if concatemer.is_empty() {
        return Ok(());
    }
    if !shard.is_empty() && shard.len().saturating_add(concatemer.len()) > chunksize {
        sender
            .send((*shard_id, std::mem::take(shard)))
            .map_err(|_| anyhow!("concat PQS writer stopped before shard {shard_id}"))?;
        *shard_id += 1;
        *shard = Vec::with_capacity(chunksize.max(concatemer.len()));
    }
    shard.append(concatemer);
    Ok(())
}

fn count_complete_concatemers<R: ConcatPqsRecord>(records: &[R]) -> u64 {
    records
        .iter()
        .map(ConcatPqsRecord::read_idx)
        .fold((None, 0u64), |(previous, count), read_idx| {
            (
                Some(read_idx),
                count + u64::from(previous != Some(read_idx)),
            )
        })
        .1
}

fn count_frame_concatemers(frame: &polars::prelude::DataFrame) -> anyResult<u64> {
    let read_ids = frame.column("read_idx")?.as_materialized_series().u64()?;
    Ok(read_ids
        .into_no_null_iter()
        .fold((None, 0u64), |(previous, count), read_idx| {
            (
                Some(read_idx),
                count + u64::from(previous != Some(read_idx)),
            )
        })
        .1)
}

fn concat_records_to_frame<R: ConcatPqsRecord>(
    records: &[R],
    use_u32_positions: bool,
) -> anyResult<polars::prelude::DataFrame> {
    use polars::prelude::*;

    let read_idx = records
        .iter()
        .map(ConcatPqsRecord::read_idx)
        .collect::<Vec<_>>();
    let read_length = records
        .iter()
        .map(ConcatPqsRecord::query_length)
        .collect::<Vec<_>>();
    let read_start = records
        .iter()
        .map(ConcatPqsRecord::query_start)
        .collect::<Vec<_>>();
    let read_end = records
        .iter()
        .map(ConcatPqsRecord::query_end)
        .collect::<Vec<_>>();
    let strands = records
        .iter()
        .map(|record| {
            if record.query_strand() == '+' {
                "+"
            } else {
                "-"
            }
        })
        .collect::<Vec<_>>();
    let chroms = records
        .iter()
        .map(ConcatPqsRecord::target)
        .collect::<Vec<_>>();
    let mapping_quality = records
        .iter()
        .map(ConcatPqsRecord::mapq)
        .collect::<Vec<_>>();
    let identity = records
        .iter()
        .map(ConcatPqsRecord::identity)
        .collect::<Vec<_>>();
    let filter_reason = records
        .iter()
        .map(ConcatPqsRecord::filter_reason)
        .collect::<Vec<_>>();
    let start = if use_u32_positions {
        Series::new(
            "start".into(),
            records
                .iter()
                .map(|record| {
                    u32::try_from(record.target_start()).with_context(|| {
                        format!(
                            "Pore-C coordinate {} on {} exceeds UInt32 selected by chromsizes",
                            record.target_start(),
                            record.target()
                        )
                    })
                })
                .collect::<anyResult<Vec<_>>>()?,
        )
    } else {
        Series::new(
            "start".into(),
            records
                .iter()
                .map(ConcatPqsRecord::target_start)
                .collect::<Vec<_>>(),
        )
    };
    let end = if use_u32_positions {
        Series::new(
            "end".into(),
            records
                .iter()
                .map(|record| {
                    u32::try_from(record.target_end()).with_context(|| {
                        format!(
                            "Pore-C coordinate {} on {} exceeds UInt32 selected by chromsizes",
                            record.target_end(),
                            record.target()
                        )
                    })
                })
                .collect::<anyResult<Vec<_>>>()?,
        )
    } else {
        Series::new(
            "end".into(),
            records
                .iter()
                .map(ConcatPqsRecord::target_end)
                .collect::<Vec<_>>(),
        )
    };

    Ok(DataFrame::new(vec![
        Series::new("read_idx".into(), read_idx).into(),
        Series::new("read_length".into(), read_length).into(),
        Series::new("read_start".into(), read_start).into(),
        Series::new("read_end".into(), read_end).into(),
        Series::new("strand".into(), strands)
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?
            .into(),
        Series::new("chrom".into(), chroms)
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?
            .into(),
        start.into(),
        end.into(),
        Series::new("mapping_quality".into(), mapping_quality).into(),
        Series::new("identity".into(), identity).into(),
        Series::new("filter_reason".into(), filter_reason)
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?
            .into(),
    ])?)
}

/// Write complete Pore-C record batches directly to an alignment-level PQS.
///
/// Every input batch must end at a read boundary. Batch order itself may be
/// arbitrary (as with parallel PAF workers), because read IDs are globally
/// unique and each complete read remains inside one batch.
pub fn write_concat_pqs_record_batches<I, R>(
    batches: I,
    target_sizes: Arc<Mutex<HashMap<String, u64>>>,
    output: &str,
    chunksize: usize,
    threads: usize,
) -> anyResult<()>
where
    I: IntoIterator<Item = Vec<R>>,
    R: ConcatPqsRecord,
{
    use polars::prelude::*;

    if chunksize == 0 || threads == 0 {
        bail!("concat PQS chunksize and thread count must be at least 1");
    }
    let output_path = Path::new(output);
    if output_path.exists() {
        bail!(
            "concat PQS output already exists: {}",
            output_path.display()
        );
    }
    std::fs::create_dir_all(output_path.join("q0"))?;
    std::fs::create_dir_all(output_path.join("q1"))?;
    polars::enable_string_cache();

    // PAF target lengths are discovered concurrently with record production,
    // so coordinates use UInt64 without requiring a preliminary PAF scan.
    let use_u32_positions = false;
    let writer_threads = threads.saturating_div(3).max(1);
    let q0_records = Arc::new(AtomicU64::new(0));
    let q1_records = Arc::new(AtomicU64::new(0));
    let q0_concats = Arc::new(AtomicU64::new(0));
    let q1_concats = Arc::new(AtomicU64::new(0));
    let (sender, receiver) = bounded::<(usize, Vec<R>)>(writer_threads * 2);
    let mut handles = Vec::with_capacity(writer_threads);

    for _ in 0..writer_threads {
        let receiver = receiver.clone();
        let output = output.to_string();
        let q0_records = Arc::clone(&q0_records);
        let q1_records = Arc::clone(&q1_records);
        let q0_concats = Arc::clone(&q0_concats);
        let q1_concats = Arc::clone(&q1_concats);
        handles.push(thread::spawn(move || -> anyResult<()> {
            while let Ok((shard_id, records)) = receiver.recv() {
                let q0_concat_count = count_complete_concatemers(&records);
                let mut frame = concat_records_to_frame(&records, use_u32_positions)?;
                let mapq = frame
                    .column("mapping_quality")?
                    .as_materialized_series()
                    .u8()?;
                let mask = mapq.gt_eq(1);
                let mut q1_frame = frame.filter(&mask)?;

                q0_records.fetch_add(frame.height() as u64, AtomOrdering::Relaxed);
                q1_records.fetch_add(q1_frame.height() as u64, AtomOrdering::Relaxed);
                q0_concats.fetch_add(q0_concat_count, AtomOrdering::Relaxed);
                q1_concats.fetch_add(count_frame_concatemers(&q1_frame)?, AtomOrdering::Relaxed);

                let q0_path = format!("{output}/q0/{shard_id}.parquet");
                ParquetWriter::new(File::create(&q0_path)?)
                    .finish(&mut frame)
                    .with_context(|| format!("failed to write {q0_path}"))?;
                if q1_frame.height() > 0 {
                    let q1_path = format!("{output}/q1/{shard_id}.parquet");
                    ParquetWriter::new(File::create(&q1_path)?)
                        .finish(&mut q1_frame)
                        .with_context(|| format!("failed to write {q1_path}"))?;
                }
            }
            Ok(())
        }));
    }
    drop(receiver);

    let mut shard_id = 0usize;
    let mut shard = Vec::<R>::with_capacity(chunksize);
    let mut current_read_id = None::<u64>;
    let mut current_concatemer = Vec::<R>::new();
    let mut producer_result = Ok(());
    for batch in batches {
        if producer_result.is_err() {
            continue;
        }
        for record in batch {
            if let Some(read_id) = current_read_id {
                if record.read_idx() != read_id {
                    if let Err(error) = push_complete_concatemer(
                        &sender,
                        &mut shard,
                        &mut current_concatemer,
                        chunksize,
                        &mut shard_id,
                    ) {
                        producer_result = Err(error);
                        break;
                    }
                }
            }
            current_read_id = Some(record.read_idx());
            current_concatemer.push(record);
        }
    }
    if producer_result.is_ok() {
        producer_result = (|| -> anyResult<()> {
            push_complete_concatemer(
                &sender,
                &mut shard,
                &mut current_concatemer,
                chunksize,
                &mut shard_id,
            )?;
            if !shard.is_empty() {
                sender
                    .send((shard_id, shard))
                    .map_err(|_| anyhow!("concat PQS writer stopped before the final shard"))?;
            }
            Ok(())
        })();
    }
    drop(sender);

    let mut writer_result = Ok(());
    for handle in handles {
        let result = handle
            .join()
            .map_err(|_| anyhow!("concat PQS writer thread panicked"))
            .and_then(|result| result);
        if writer_result.is_ok() && result.is_err() {
            writer_result = result;
        }
    }
    producer_result?;
    writer_result?;

    let q0_n = q0_records.load(AtomOrdering::Relaxed);
    if q0_n == 0 {
        bail!("concat PQS conversion produced no Pore-C alignment records");
    }
    let q1_n = q1_records.load(AtomOrdering::Relaxed);
    let q0_concat_n = q0_concats.load(AtomOrdering::Relaxed);
    let q1_concat_n = q1_concats.load(AtomOrdering::Relaxed);

    let sizes = target_sizes.lock().unwrap();
    if sizes.is_empty() {
        bail!("concat PQS output requires target lengths");
    }
    let mut sorted_sizes = sizes.iter().collect::<Vec<_>>();
    sorted_sizes.sort_by(|left, right| left.0.cmp(right.0));
    let mut contigsizes =
        common_writer(output_path.join("_contigsizes").to_string_lossy().as_ref());
    for (target, length) in sorted_sizes {
        writeln!(contigsizes, "{target}\t{length}")?;
    }
    contigsizes.flush()?;

    let mut counts = common_writer(
        output_path
            .join("_metadata_counts")
            .to_string_lossy()
            .as_ref(),
    );
    writeln!(counts, "q0_records\t{q0_n}")?;
    writeln!(counts, "q1_records\t{q1_n}")?;
    writeln!(counts, "q0_concats\t{q0_concat_n}")?;
    writeln!(counts, "q1_concats\t{q1_concat_n}")?;
    counts.flush()?;
    std::fs::write(
        output_path.join("_metadata"),
        concat_pqs_metadata(chunksize, "UInt64"),
    )?;
    std::fs::write(output_path.join("_readme"), CONCAT_PQS_README)?;
    log::info!(
        "Successful concat PQS output `{}`: {} alignments in {} complete concatemers across {} shards",
        output,
        q0_n,
        q0_concat_n,
        shard_id + 1
    );
    Ok(())
}

impl BaseTable for PoreCTable {
    fn new(name: &String) -> PoreCTable {
        PoreCTable { file: name.clone() }
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

impl PoreCTable {
    pub fn parse(&self) -> anyResult<csv::Reader<Box<dyn BufRead + Send>>> {
        let input = open_porec_reader(&self.file)?;
        let rdr = csv::ReaderBuilder::new()
            .flexible(true)
            .has_headers(false)
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .from_reader(input);

        Ok(rdr)
    }

    pub fn parse2(&mut self) -> anyResult<Box<dyn BufRead + Send + 'static>> {
        let input = open_porec_reader(&self.file)?;

        Ok(input)
    }

    /// Convert an alignment-level Pore-C table to a sharded concatemer PQS.
    ///
    /// Shards target `chunksize` alignment rows, but boundaries are selected
    /// only between read IDs. A concatemer larger than `chunksize` therefore
    /// occupies one oversized shard rather than being split across files.
    pub fn to_concat_pqs(
        &mut self,
        chromsizes: &String,
        output: &String,
        chunksize: usize,
        threads: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;

        if chunksize == 0 {
            bail!("concat PQS chunksize must be at least 1");
        }
        if threads == 0 {
            bail!("concat PQS thread count must be at least 1");
        }
        let output_path = Path::new(output);
        if output_path.exists() {
            bail!(
                "concat PQS output already exists: {}",
                output_path.display()
            );
        }

        if is_concat_pqs(Path::new(&self.file)) {
            return clone_concat_pqs_native(
                Path::new(&self.file),
                output_path,
                Path::new(chromsizes),
                threads,
            );
        }

        let contigsizes = ChromSize::new(chromsizes);
        let contigsizes_data = contigsizes
            .to_vec()
            .map_err(|error| anyhow!("failed to read chromsizes file {chromsizes}: {error}"))?;
        let max_contig_size = contigsizes_data
            .iter()
            .map(|record| record.size)
            .max()
            .unwrap_or(0);
        let use_u32_positions = max_contig_size <= u32::MAX as u64;

        std::fs::create_dir_all(output_path.join("q0"))?;
        std::fs::create_dir_all(output_path.join("q1"))?;
        std::fs::copy(chromsizes, output_path.join("_contigsizes"))?;
        polars::enable_string_cache();

        let parser_threads = threads.saturating_mul(2).div_ceil(3).max(1);
        let writer_threads = threads.saturating_sub(parser_threads).max(1);
        log::info!(
            "porec2pqs pipeline: 1 reader, {} parsers, {} Parquet writers",
            parser_threads,
            writer_threads
        );

        let q0_records = Arc::new(AtomicU64::new(0));
        let q1_records = Arc::new(AtomicU64::new(0));
        let q0_concats = Arc::new(AtomicU64::new(0));
        let q1_concats = Arc::new(AtomicU64::new(0));
        let (sender, receiver) = bounded::<(usize, Vec<PoreCRecord>)>(writer_threads * 2);
        let mut handles = Vec::with_capacity(writer_threads);

        for _ in 0..writer_threads {
            let receiver = receiver.clone();
            let output = output.clone();
            let q0_records = Arc::clone(&q0_records);
            let q1_records = Arc::clone(&q1_records);
            let q0_concats = Arc::clone(&q0_concats);
            let q1_concats = Arc::clone(&q1_concats);
            handles.push(thread::spawn(move || -> anyResult<()> {
                while let Ok((shard_id, records)) = receiver.recv() {
                    let q0_concat_count = count_complete_concatemers(&records);
                    let mut frame = concat_records_to_frame(&records, use_u32_positions)?;
                    let mapq = frame
                        .column("mapping_quality")?
                        .as_materialized_series()
                        .u8()?;
                    let mask = mapq.gt_eq(1);
                    let mut q1_frame = frame.filter(&mask)?;

                    q0_records.fetch_add(frame.height() as u64, AtomOrdering::Relaxed);
                    q1_records.fetch_add(q1_frame.height() as u64, AtomOrdering::Relaxed);
                    q0_concats.fetch_add(q0_concat_count, AtomOrdering::Relaxed);
                    q1_concats
                        .fetch_add(count_frame_concatemers(&q1_frame)?, AtomOrdering::Relaxed);

                    let q0_path = format!("{output}/q0/{shard_id}.parquet");
                    ParquetWriter::new(File::create(&q0_path)?)
                        .finish(&mut frame)
                        .with_context(|| format!("failed to write {q0_path}"))?;
                    if q1_frame.height() > 0 {
                        let q1_path = format!("{output}/q1/{shard_id}.parquet");
                        ParquetWriter::new(File::create(&q1_path)?)
                            .finish(&mut q1_frame)
                            .with_context(|| format!("failed to write {q1_path}"))?;
                    }
                }
                Ok(())
            }));
        }
        drop(receiver);

        const PARSE_BATCH_BYTES: usize = 8 * 1024 * 1024;
        type ParseTask = (usize, usize, String);
        type ParseResult = (usize, anyResult<Vec<PoreCRecord>>);
        let (parse_sender, parse_receiver) = bounded::<ParseTask>(parser_threads * 2);
        let (parsed_sender, parsed_receiver) = bounded::<ParseResult>(parser_threads * 2);
        let reorder_window = parser_threads * 2;
        let (permit_sender, permit_receiver) = bounded::<()>(reorder_window);
        for _ in 0..reorder_window {
            permit_sender
                .send(())
                .expect("initialize porec2pqs reorder permits");
        }
        let mut parser_handles = Vec::with_capacity(parser_threads);
        for _ in 0..parser_threads {
            let receiver = parse_receiver.clone();
            let sender = parsed_sender.clone();
            parser_handles.push(thread::spawn(move || {
                while let Ok((batch_id, first_line_number, batch)) = receiver.recv() {
                    let records = parse_concat_porec_batch(&batch, first_line_number);
                    if sender.send((batch_id, records)).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(parse_receiver);
        drop(parsed_sender);

        let input = self.file.clone();
        let reader_handle = thread::spawn(move || -> anyResult<usize> {
            let mut reader = open_porec_reader(&input)?;
            let mut batch = String::with_capacity(PARSE_BATCH_BYTES + 1024);
            let mut batch_id = 0usize;
            let mut line_number = 0usize;
            let mut first_line_number = 1usize;

            loop {
                let bytes_read = reader.read_line(&mut batch)?;
                if bytes_read == 0 {
                    break;
                }
                line_number += 1;
                if batch.len() >= PARSE_BATCH_BYTES {
                    permit_receiver.recv().map_err(|_| {
                        anyhow!("porec2pqs ordered collector stopped while reading input")
                    })?;
                    parse_sender
                        .send((batch_id, first_line_number, std::mem::take(&mut batch)))
                        .map_err(|_| anyhow!("porec2pqs parsers stopped while reading input"))?;
                    batch_id += 1;
                    first_line_number = line_number + 1;
                    batch = String::with_capacity(PARSE_BATCH_BYTES + 1024);
                }
            }
            if !batch.is_empty() {
                permit_receiver.recv().map_err(|_| {
                    anyhow!("porec2pqs ordered collector stopped before the final batch")
                })?;
                parse_sender
                    .send((batch_id, first_line_number, batch))
                    .map_err(|_| anyhow!("porec2pqs parsers stopped before the final batch"))?;
                batch_id += 1;
            }
            Ok(batch_id)
        });

        let mut shard_id = 0usize;
        let mut pipeline_result = Ok(());
        let mut parsed_batches = 0usize;
        let mut next_batch_id = 0usize;
        let mut pending = BTreeMap::<usize, anyResult<Vec<PoreCRecord>>>::new();
        let mut current_read_id = None::<u64>;
        let mut current_concatemer = Vec::<PoreCRecord>::new();
        let mut shard = Vec::<PoreCRecord>::with_capacity(chunksize);

        for (batch_id, records) in parsed_receiver {
            parsed_batches += 1;
            pending.insert(batch_id, records);
            while let Some(records) = pending.remove(&next_batch_id) {
                next_batch_id += 1;
                // The reader may already have reached EOF and dropped its
                // receiver while the final parsed batches are being drained.
                let _ = permit_sender.send(());
                if pipeline_result.is_err() {
                    continue;
                }
                let records = match records {
                    Ok(records) => records,
                    Err(error) => {
                        pipeline_result = Err(error);
                        continue;
                    }
                };
                for record in records {
                    if let Some(read_id) = current_read_id {
                        if record.read_idx != read_id {
                            if let Err(error) = push_complete_concatemer(
                                &sender,
                                &mut shard,
                                &mut current_concatemer,
                                chunksize,
                                &mut shard_id,
                            ) {
                                pipeline_result = Err(error);
                                break;
                            }
                        }
                    }
                    current_read_id = Some(record.read_idx);
                    current_concatemer.push(record);
                }
            }
        }

        let reader_result = reader_handle
            .join()
            .map_err(|_| anyhow!("porec2pqs reader thread panicked"))
            .and_then(|result| result);
        let mut parser_result = Ok(());
        for handle in parser_handles {
            if handle.join().is_err() && parser_result.is_ok() {
                parser_result = Err(anyhow!("porec2pqs parser thread panicked"));
            }
        }
        if pipeline_result.is_ok() {
            if let Err(error) = &reader_result {
                pipeline_result = Err(anyhow!("{error}"));
            } else if let Err(error) = parser_result {
                pipeline_result = Err(error);
            } else if let Ok(expected_batches) = reader_result {
                if parsed_batches != expected_batches || !pending.is_empty() {
                    pipeline_result = Err(anyhow!(
                        "porec2pqs parsing pipeline lost batches: expected {expected_batches}, received {parsed_batches}"
                    ));
                }
            }
        }

        if pipeline_result.is_ok() {
            pipeline_result = (|| -> anyResult<()> {
                push_complete_concatemer(
                    &sender,
                    &mut shard,
                    &mut current_concatemer,
                    chunksize,
                    &mut shard_id,
                )?;
                if !shard.is_empty() {
                    sender
                        .send((shard_id, shard))
                        .map_err(|_| anyhow!("concat PQS writer stopped before the final shard"))?;
                }
                Ok(())
            })();
        }
        drop(sender);

        let mut writer_result = Ok(());
        for handle in handles {
            let result = handle
                .join()
                .map_err(|_| anyhow!("concat PQS writer thread panicked"))
                .and_then(|result| result);
            if writer_result.is_ok() && result.is_err() {
                writer_result = result;
            }
        }
        pipeline_result?;
        writer_result?;

        let q0_n = q0_records.load(AtomOrdering::Relaxed);
        let q1_n = q1_records.load(AtomOrdering::Relaxed);
        if q0_n == 0 {
            bail!("Pore-C input contains no alignment records");
        }
        let q0_concat_n = q0_concats.load(AtomOrdering::Relaxed);
        let q1_concat_n = q1_concats.load(AtomOrdering::Relaxed);
        let mut counts = common_writer(
            output_path
                .join("_metadata_counts")
                .to_string_lossy()
                .as_ref(),
        );
        writeln!(counts, "q0_records\t{q0_n}")?;
        writeln!(counts, "q1_records\t{q1_n}")?;
        writeln!(counts, "q0_concats\t{q0_concat_n}")?;
        writeln!(counts, "q1_concats\t{q1_concat_n}")?;
        counts.flush()?;

        let position_type = if use_u32_positions {
            "UInt32"
        } else {
            "UInt64"
        };
        let metadata = concat_pqs_metadata(chunksize, position_type);
        std::fs::write(output_path.join("_metadata"), metadata)?;
        std::fs::write(output_path.join("_readme"), CONCAT_PQS_README)?;
        log::info!(
            "Successful output concat PQS `{}`: {} alignments in {} complete concatemers across {} shards",
            output,
            q0_n,
            q0_concat_n,
            shard_id + 1
        );
        Ok(())
    }

    pub fn split_to_concat_pqs(
        &mut self,
        chromsizes: &String,
        output: &String,
        chunksize: usize,
        threads: usize,
    ) -> anyResult<()> {
        if is_concat_pqs(Path::new(&self.file)) {
            let mut input_sizes = BTreeMap::new();
            let mut requested_sizes = BTreeMap::new();
            read_concat_pqs_contigsizes(Path::new(&self.file), &mut input_sizes)?;
            read_contigsizes_file(Path::new(chromsizes), &mut requested_sizes)?;
            if input_sizes != requested_sizes {
                bail!("porec-split chromsizes differ from concat.pqs _contigsizes");
            }
            return rechunk_concat_pqs_native(
                Path::new(&self.file),
                Path::new(output.trim_end_matches('/')),
                chunksize,
                threads,
            );
        }
        self.to_concat_pqs(chromsizes, output, chunksize, threads)
    }

    pub fn to_pairs_pqs(
        &mut self,
        chromsizes: &String,
        output: &String,
        chunksize: usize,
        min_quality: u8,
        min_order: usize,
        max_order: usize,
        threads: usize,
    ) -> anyResult<()> {
        use crate::pqs::_METADATA;
        use crate::pqs::_README as _readme;
        use polars::prelude::*;

        if threads == 0 {
            bail!("porec2pairs thread count must be at least 1");
        }
        polars::enable_string_cache();

        if is_concat_pqs(Path::new(&self.file)) {
            return self.to_pairs_pqs_native(
                chromsizes,
                output,
                chunksize,
                min_quality,
                min_order,
                max_order,
                threads,
            );
        }

        let mut rdr = self.parse2().expect("Failed to open input file");
        log::info!(
            "Only retain concatemer that order in the range of [{}, {})",
            min_order,
            max_order
        );

        let _ = std::fs::create_dir_all(output);
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));
        std::fs::copy(chromsizes, format!("{}/_contigsizes", output))?;

        let contigsizes = ChromSize::new(chromsizes);
        let contigsizes_data = contigsizes.to_vec().unwrap();
        let max_contig_size = contigsizes_data.iter().map(|x| x.size).max().unwrap_or(0);
        let pos_type = if max_contig_size < 4294967295 {
            DataType::UInt32
        } else {
            DataType::UInt64
        };

        let pos_type_string = match pos_type {
            DataType::UInt32 => "UInt32",
            DataType::UInt64 => "UInt64",
            _ => "UInt32",
        };
        let mut wtr = common_writer(format!("{}/_readme", output).as_str());
        wtr.write_all(_readme.as_bytes()).unwrap();
        wtr.flush().unwrap();

        let create_date_time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string();

        let mut wtr = common_writer(format!("{}/_metadata", output).as_str());
        let mut _metadata = _METADATA.to_string();
        _metadata = _metadata.replace("REPLACE", &create_date_time);
        _metadata = _metadata.replace("CHUNKSIZE", &chunksize.to_string());
        _metadata = _metadata.replace("pos_type_lower", pos_type_string.to_lowercase().as_str());
        _metadata = _metadata.replace("pos_type", pos_type_string);
        wtr.write_all(_metadata.as_bytes()).unwrap();
        wtr.flush().unwrap();

        let q0_total = Arc::new(AtomicU64::new(0));
        let q1_total = Arc::new(AtomicU64::new(0));
        let concatemer_summary = Arc::new(Mutex::new(ConcatemerSummary::new()));

        let num_workers = threads;
        let (sender, receiver) = bounded::<(usize, Vec<u8>, u64)>(2);
        let mut handles = vec![];

        for _ in 0..num_workers {
            let rx = receiver.clone();
            let output = output.clone();
            let q0_total = Arc::clone(&q0_total);
            let q1_total = Arc::clone(&q1_total);
            let summary_lock = Arc::clone(&concatemer_summary);

            handles.push(thread::spawn(move || {
                while let Ok((chunk_id, chunk, mut current_id)) = rx.recv() {
                    let mut records = Vec::<PairConcatemer>::with_capacity(chunk.len() / 400);
                    let mut concatemer = PairConcatemer::default();
                    let mut target_interner = HashSet::<Arc<str>>::new();
                    let mut previous_read_idx = u64::MAX;

                    for raw_line in chunk.split(|byte| *byte == b'\n') {
                        let trimmed = std::str::from_utf8(raw_line).unwrap().trim_end();
                        if trimmed.is_empty() || trimmed.starts_with('#') {
                            continue;
                        }

                        let mut parts = trimmed.split('\t');
                        let read_idx = parts
                            .next()
                            .and_then(|value| value.parse::<u64>().ok())
                            .unwrap_or(0);
                        if previous_read_idx != u64::MAX && read_idx != previous_read_idx {
                            let order = concatemer.count();
                            if order >= min_order && order < max_order {
                                records.push(std::mem::take(&mut concatemer));
                            } else {
                                concatemer.clear();
                            }
                        }

                        let _ = parts.next();
                        let _ = parts.next();
                        let _ = parts.next();
                        let query_strand = parts
                            .next()
                            .and_then(|value| value.chars().next())
                            .unwrap_or('+');
                        let target_str = parts.next().unwrap_or("");
                        let target_start = parts
                            .next()
                            .and_then(|value| value.parse::<u64>().ok())
                            .unwrap_or(0);
                        let target_end = parts
                            .next()
                            .and_then(|value| value.parse::<u64>().ok())
                            .unwrap_or(0);
                        let mapq = parts
                            .next()
                            .and_then(|value| value.parse::<u8>().ok())
                            .unwrap_or(0);
                        let target = if let Some(target) = target_interner.get(target_str) {
                            Arc::clone(target)
                        } else {
                            let target = Arc::<str>::from(target_str);
                            target_interner.insert(Arc::clone(&target));
                            target
                        };
                        concatemer.records.push(PairAnchor {
                            target,
                            target_start,
                            target_end,
                            query_strand,
                            mapq,
                        });
                        previous_read_idx = read_idx;
                    }

                    let order = concatemer.count();
                    if order >= min_order && order < max_order {
                        records.push(concatemer);
                    }
                    drop(target_interner);
                    drop(chunk);

                    let mut local_sum = HashMap::new();
                    for c in &records {
                        *local_sum.entry(c.count() as u32).or_insert(0) += 1;
                    }
                    {
                        let mut global = summary_lock.lock().unwrap();
                        for (k, v) in local_sum {
                            *global.summary.entry(k).or_insert(0) += v;
                        }
                    }

                    let est_capacity = records
                        .iter()
                        .map(|c| c.count() * (c.count() - 1) / 2)
                        .sum();
                    if est_capacity == 0 {
                        continue;
                    }

                    let mut read_idx_vec = Vec::with_capacity(est_capacity);
                    let mut chrom1_vec = Vec::with_capacity(est_capacity);
                    let mut chrom2_vec = Vec::with_capacity(est_capacity);

                    let use_u32 = max_contig_size < 4294967295;
                    let mut pos1_u32 = Vec::with_capacity(if use_u32 { est_capacity } else { 0 });
                    let mut pos2_u32 = Vec::with_capacity(if use_u32 { est_capacity } else { 0 });
                    let mut pos1_u64 = Vec::with_capacity(if use_u32 { 0 } else { est_capacity });
                    let mut pos2_u64 = Vec::with_capacity(if use_u32 { 0 } else { est_capacity });

                    let mut strand1_vec = Vec::with_capacity(est_capacity);
                    let mut strand2_vec = Vec::with_capacity(est_capacity);
                    let mut mapq_vec = Vec::with_capacity(est_capacity);

                    for concatemer in &mut records {
                        concatemer.sort();
                        let recs = &concatemer.records;
                        let n = recs.len();
                        for i in 0..n {
                            let r1 = &recs[i];
                            let r1_pos = (r1.target_start + r1.target_end) >> 1;
                            let r1_strand = if r1.query_strand == '+' { "+" } else { "-" };

                            for j in i + 1..n {
                                let r2 = &recs[j];
                                let r2_pos = (r2.target_start + r2.target_end) >> 1;

                                current_id += 1;
                                read_idx_vec.push(current_id);

                                chrom1_vec.push(r1.target.as_ref());
                                chrom2_vec.push(r2.target.as_ref());

                                if use_u32 {
                                    pos1_u32.push(r1_pos as u32);
                                    pos2_u32.push(r2_pos as u32);
                                } else {
                                    pos1_u64.push(r1_pos as u64);
                                    pos2_u64.push(r2_pos as u64);
                                }

                                strand1_vec.push(r1_strand);
                                strand2_vec.push(if r2.query_strand == '+' { "+" } else { "-" });

                                mapq_vec.push(std::cmp::min(r1.mapq, r2.mapq));
                            }
                        }
                    }

                    let pos1_series = if use_u32 {
                        Series::new("pos1".into(), pos1_u32)
                    } else {
                        Series::new("pos1".into(), pos1_u64)
                    };
                    let pos2_series = if use_u32 {
                        Series::new("pos2".into(), pos2_u32)
                    } else {
                        Series::new("pos2".into(), pos2_u64)
                    };

                    let s_read_idx = Series::new("read_idx".into(), read_idx_vec)
                        .cast(&DataType::String)
                        .unwrap();
                    let s_chrom1 = Series::new("chrom1".into(), chrom1_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_chrom2 = Series::new("chrom2".into(), chrom2_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_strand1 = Series::new("strand1".into(), strand1_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_strand2 = Series::new("strand2".into(), strand2_vec)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
                        .unwrap();
                    let s_mapq = Series::new("mapq".into(), mapq_vec);

                    let mut df = DataFrame::new(vec![
                        s_read_idx.into(),
                        s_chrom1.into(),
                        pos1_series.into(),
                        s_chrom2.into(),
                        pos2_series.into(),
                        s_strand1.into(),
                        s_strand2.into(),
                        s_mapq.into(),
                    ])
                    .unwrap();

                    q0_total.fetch_add(df.height() as u64, AtomOrdering::Relaxed);
                    let path0 = format!("{}/q0/{}.parquet", output, chunk_id);
                    ParquetWriter::new(File::create(path0).unwrap())
                        .finish(&mut df)
                        .unwrap();

                    let mapq_col = df
                        .column("mapq".into())
                        .unwrap()
                        .as_materialized_series()
                        .u8()
                        .unwrap();
                    let mask = mapq_col.gt_eq(1);
                    let mut df_q1 = df.filter(&mask).unwrap();
                    drop(df);

                    q1_total.fetch_add(df_q1.height() as u64, AtomOrdering::Relaxed);

                    if df_q1.height() > 0 {
                        let path1 = format!("{}/q1/{}.parquet", output, chunk_id);
                        ParquetWriter::new(File::create(path1).unwrap())
                            .finish(&mut df_q1)
                            .unwrap();
                    }
                }
            }));
        }

        let mut line_buf = String::new();
        let mut chunk = Vec::<u8>::with_capacity(4 * 1024 * 1024);
        let mut previous_read_idx: u64 = 0;
        let mut current_read_order = 0usize;
        let mut current_chunk_lines_len = 0usize;
        let mut current_chunk_pair_count = 0u64;
        let mut global_pair_offset = 0u64;
        let mut chunk_id = 0;
        let mut first_iteration = true;

        while rdr.read_line(&mut line_buf)? > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let mut parts = trimmed.split('\t');
            let read_idx = parts
                .next()
                .and_then(|value| value.parse::<u64>().ok())
                .unwrap_or(0);

            if !first_iteration && read_idx != previous_read_idx {
                if current_read_order >= min_order && current_read_order < max_order {
                    let order = current_read_order as u64;
                    current_chunk_pair_count += order * (order - 1) / 2;
                }
                current_read_order = 0;

                if current_chunk_lines_len >= chunksize {
                    sender
                        .send((chunk_id, std::mem::take(&mut chunk), global_pair_offset))
                        .unwrap();
                    global_pair_offset += current_chunk_pair_count;
                    chunk_id += 1;
                    current_chunk_lines_len = 0;
                    current_chunk_pair_count = 0;
                    chunk = Vec::with_capacity(4 * 1024 * 1024);
                }
            }

            first_iteration = false;
            previous_read_idx = read_idx;
            current_read_order += 1;
            current_chunk_lines_len += 1;
            chunk.extend_from_slice(line_buf.as_bytes());
            line_buf.clear();
        }

        if current_chunk_lines_len > 0 {
            sender.send((chunk_id, chunk, global_pair_offset)).unwrap();
        }

        drop(sender);
        for handle in handles {
            handle.join().unwrap();
        }

        {
            let q0_n = q0_total.load(AtomOrdering::Relaxed);
            let q1_n = q1_total.load(AtomOrdering::Relaxed);

            let mut wtr = common_writer(format!("{}/_metadata_counts", output).as_str());
            writeln!(wtr, "q0_records\t{}", q0_n).unwrap();
            writeln!(wtr, "q1_records\t{}", q1_n).unwrap();
            wtr.flush().unwrap();
        }

        let output_prefix = if output == "-" {
            Path::new(&self.file)
                .with_extension("")
                .to_str()
                .unwrap()
                .to_string()
        } else {
            Path::new(&output)
                .with_extension("")
                .to_str()
                .unwrap()
                .to_string()
        };

        let final_summary = concatemer_summary.lock().unwrap();
        final_summary.save(&format!("{}.concatemer.summary", output_prefix));

        Ok(())
    }

    fn to_pairs_pqs_native(
        &self,
        chromsizes: &String,
        output: &String,
        chunksize: usize,
        min_quality: u8,
        min_order: usize,
        max_order: usize,
        threads: usize,
    ) -> anyResult<()> {
        use crate::pqs::_METADATA;
        use crate::pqs::_README as PAIRS_PQS_README;
        use polars::prelude::*;

        let input_path = Path::new(&self.file);
        let output_path = Path::new(output.trim_end_matches('/'));
        if output_path.exists() {
            bail!("pairs PQS output already exists: {}", output_path.display());
        }
        let contigsizes_data = ChromSize::new(chromsizes)
            .to_vec()
            .map_err(|error| anyhow!("failed to read chromsizes {chromsizes}: {error}"))?;
        let max_contig_size = contigsizes_data
            .iter()
            .map(|record| record.size)
            .max()
            .unwrap_or(0);
        let use_u32_positions = max_contig_size < u32::MAX as u64;
        let position_type = if use_u32_positions {
            "UInt32"
        } else {
            "UInt64"
        };

        std::fs::create_dir_all(output_path.join("q0"))?;
        std::fs::create_dir_all(output_path.join("q1"))?;
        std::fs::copy(chromsizes, output_path.join("_contigsizes"))?;
        let creation_time = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
        let metadata = _METADATA
            .replace("REPLACE", &creation_time.to_string())
            .replace("CHUNKSIZE", &chunksize.to_string())
            .replace("pos_type_lower", &position_type.to_lowercase())
            .replace("pos_type", position_type);
        std::fs::write(output_path.join("_metadata"), metadata)?;
        std::fs::write(output_path.join("_readme"), PAIRS_PQS_README)?;
        copy_cn_info(input_path, output_path)?;
        polars::enable_string_cache();

        let shards = concat_pqs_files(input_path)?;
        let q0_total = Arc::new(AtomicU64::new(0));
        let q1_total = Arc::new(AtomicU64::new(0));
        let concatemer_summary = Arc::new(Mutex::new(ConcatemerSummary::new()));
        let writer_workers = (threads / 2).clamp(1, 4);
        let transform_workers = threads.saturating_sub(writer_workers).max(1);
        let (raw_sender, raw_receiver) = bounded::<(usize, DataFrame)>(2);
        let (pairs_sender, pairs_receiver) =
            bounded::<(usize, DataFrame, HashMap<u32, u64>)>(writer_workers);
        log::info!(
            "Native concat PQS to pairs PQS: {} shards with one projected reader, {} pair workers, and {} bounded writer workers",
            shards.len(),
            transform_workers,
            writer_workers,
        );

        let mut writer_handles = Vec::with_capacity(writer_workers);
        for _ in 0..writer_workers {
            let receiver = pairs_receiver.clone();
            let output_path = output_path.to_path_buf();
            let q0_total = Arc::clone(&q0_total);
            let q1_total = Arc::clone(&q1_total);
            let concatemer_summary = Arc::clone(&concatemer_summary);
            writer_handles.push(thread::spawn(move || -> anyResult<()> {
                while let Ok((shard_id, mut pairs, local_summary)) = receiver.recv() {
                    let q0_count = pairs.height();
                    let q0_path = output_path.join(format!("q0/{shard_id}.parquet"));
                    ParquetWriter::new(File::create(&q0_path)?)
                        .finish(&mut pairs)
                        .with_context(|| format!("failed to write {}", q0_path.display()))?;
                    q0_total.fetch_add(q0_count as u64, AtomOrdering::Relaxed);

                    // Materialize q1 only after q0 compression has finished, then release q0
                    // before compressing q1. This shortens the lifetime of the duplicated frame.
                    let pair_mapq = pairs.column("mapq")?.as_materialized_series().u8()?;
                    let mut q1_pairs = pairs.filter(&pair_mapq.gt_eq(1))?;
                    drop(pairs);
                    let q1_count = q1_pairs.height();
                    if q1_count > 0 {
                        let q1_path = output_path.join(format!("q1/{shard_id}.parquet"));
                        ParquetWriter::new(File::create(&q1_path)?)
                            .finish(&mut q1_pairs)
                            .with_context(|| format!("failed to write {}", q1_path.display()))?;
                    }
                    q1_total.fetch_add(q1_count as u64, AtomOrdering::Relaxed);

                    let mut summary = concatemer_summary.lock().unwrap();
                    for (order, count) in local_summary {
                        *summary.summary.entry(order).or_insert(0) += count;
                    }
                }
                Ok(())
            }));
        }
        drop(pairs_receiver);

        let mut transform_handles = Vec::with_capacity(transform_workers);
        for _ in 0..transform_workers {
            let receiver = raw_receiver.clone();
            let sender = pairs_sender.clone();
            transform_handles.push(thread::spawn(move || -> anyResult<()> {
                while let Ok((shard_id, frame)) = receiver.recv() {
                    let read_idx = frame.column("read_idx")?.as_materialized_series().u64()?;
                    let read_start = frame.column("read_start")?.as_materialized_series().u32()?;
                    let strand_series = frame
                        .column("strand")?
                        .as_materialized_series()
                        .cast(&DataType::String)?;
                    let strand = strand_series.str()?;
                    let chrom_series = frame
                        .column("chrom")?
                        .as_materialized_series()
                        .cast(&DataType::String)?;
                    let chrom = chrom_series.str()?;
                    let start_source = frame.column("start")?.as_materialized_series();
                    let end_source = frame.column("end")?.as_materialized_series();
                    let native_u32_positions = use_u32_positions
                        && start_source.dtype() == &DataType::UInt32
                        && end_source.dtype() == &DataType::UInt32;
                    let start_u64_series = if native_u32_positions {
                        None
                    } else {
                        Some(start_source.cast(&DataType::UInt64)?)
                    };
                    let end_u64_series = if native_u32_positions {
                        None
                    } else {
                        Some(end_source.cast(&DataType::UInt64)?)
                    };
                    let start_u32 = native_u32_positions
                        .then(|| start_source.u32())
                        .transpose()?;
                    let end_u32 = native_u32_positions.then(|| end_source.u32()).transpose()?;
                    let start_u64 = start_u64_series.as_ref().map(|s| s.u64()).transpose()?;
                    let end_u64 = end_u64_series.as_ref().map(|s| s.u64()).transpose()?;
                    let mapq = frame
                        .column("mapping_quality")?
                        .as_materialized_series()
                        .u8()?;

                    let estimated_pairs = frame.height().saturating_mul(2);
                    let mut pair_ids = Vec::<u64>::with_capacity(estimated_pairs);
                    let mut chrom1 = Vec::<&str>::with_capacity(estimated_pairs);
                    let mut chrom2 = Vec::<&str>::with_capacity(estimated_pairs);
                    let mut pos1_u32 = Vec::<u32>::with_capacity(
                        native_u32_positions.then_some(estimated_pairs).unwrap_or(0),
                    );
                    let mut pos2_u32 = Vec::<u32>::with_capacity(
                        native_u32_positions.then_some(estimated_pairs).unwrap_or(0),
                    );
                    let mut pos1_u64 = Vec::<u64>::with_capacity(
                        (!native_u32_positions)
                            .then_some(estimated_pairs)
                            .unwrap_or(0),
                    );
                    let mut pos2_u64 = Vec::<u64>::with_capacity(
                        (!native_u32_positions)
                            .then_some(estimated_pairs)
                            .unwrap_or(0),
                    );
                    let mut strand1 = Vec::<&str>::with_capacity(estimated_pairs);
                    let mut strand2 = Vec::<&str>::with_capacity(estimated_pairs);
                    let mut pair_mapq = Vec::<u8>::with_capacity(estimated_pairs);
                    let mut local_pair_id = 0u64;
                    let shard_id_prefix = (shard_id as u64)
                        .checked_shl(48)
                        .ok_or_else(|| anyhow!("pair shard id overflow: {shard_id}"))?;
                    let mut local_summary = HashMap::<u32, u64>::new();

                    let mut group_start = 0usize;
                    while group_start < frame.height() {
                        let group_read_idx = read_idx
                            .get(group_start)
                            .ok_or_else(|| anyhow!("null read_idx at row {group_start}"))?;
                        let mut group_end = group_start + 1;
                        while group_end < frame.height()
                            && read_idx.get(group_end) == Some(group_read_idx)
                        {
                            group_end += 1;
                        }
                        let mut rows = (group_start..group_end).collect::<Vec<_>>();
                        let order = rows.len();
                        if order >= min_order && order < max_order {
                            *local_summary.entry(order as u32).or_insert(0) += 1;
                            rows.sort_unstable_by_key(|&row| read_start.get(row).unwrap_or(0));
                            for left in 0..order {
                                for right in left + 1..order {
                                    let row1 = rows[left];
                                    let row2 = rows[right];
                                    local_pair_id =
                                        local_pair_id.checked_add(1).ok_or_else(|| {
                                            anyhow!("pair id overflow in shard {shard_id}")
                                        })?;
                                    if local_pair_id >= (1u64 << 48) {
                                        bail!("too many pairs in shard {shard_id}");
                                    }
                                    pair_ids.push(shard_id_prefix | local_pair_id);
                                    chrom1.push(
                                        chrom
                                            .get(row1)
                                            .ok_or_else(|| anyhow!("null chrom at row {row1}"))?,
                                    );
                                    chrom2.push(
                                        chrom
                                            .get(row2)
                                            .ok_or_else(|| anyhow!("null chrom at row {row2}"))?,
                                    );
                                    if native_u32_positions {
                                        let start = start_u32.as_ref().unwrap();
                                        let end = end_u32.as_ref().unwrap();
                                        pos1_u32.push(
                                            ((start.get(row1).unwrap_or(0) as u64
                                                + end.get(row1).unwrap_or(0) as u64)
                                                >> 1)
                                                as u32,
                                        );
                                        pos2_u32.push(
                                            ((start.get(row2).unwrap_or(0) as u64
                                                + end.get(row2).unwrap_or(0) as u64)
                                                >> 1)
                                                as u32,
                                        );
                                    } else {
                                        let start = start_u64.as_ref().unwrap();
                                        let end = end_u64.as_ref().unwrap();
                                        pos1_u64.push(
                                            (start.get(row1).unwrap_or(0)
                                                + end.get(row1).unwrap_or(0))
                                                >> 1,
                                        );
                                        pos2_u64.push(
                                            (start.get(row2).unwrap_or(0)
                                                + end.get(row2).unwrap_or(0))
                                                >> 1,
                                        );
                                    }
                                    strand1.push(strand.get(row1).unwrap_or("+"));
                                    strand2.push(strand.get(row2).unwrap_or("+"));
                                    pair_mapq.push(
                                        mapq.get(row1)
                                            .unwrap_or(0)
                                            .min(mapq.get(row2).unwrap_or(0)),
                                    );
                                }
                            }
                        }
                        group_start = group_end;
                    }

                    if pair_ids.is_empty() {
                        continue;
                    }
                    let pair_ids =
                        Series::new("read_idx".into(), pair_ids).cast(&DataType::String)?;
                    let chrom1 = Series::new("chrom1".into(), chrom1)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;
                    let chrom2 = Series::new("chrom2".into(), chrom2)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;
                    let pos1 = if native_u32_positions {
                        Series::new("pos1".into(), pos1_u32)
                    } else {
                        Series::new("pos1".into(), pos1_u64).cast(if use_u32_positions {
                            &DataType::UInt32
                        } else {
                            &DataType::UInt64
                        })?
                    };
                    let pos2 = if native_u32_positions {
                        Series::new("pos2".into(), pos2_u32)
                    } else {
                        Series::new("pos2".into(), pos2_u64).cast(if use_u32_positions {
                            &DataType::UInt32
                        } else {
                            &DataType::UInt64
                        })?
                    };
                    let strand1 = Series::new("strand1".into(), strand1)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;
                    let strand2 = Series::new("strand2".into(), strand2)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;
                    let pairs = DataFrame::new(vec![
                        pair_ids.into(),
                        chrom1.into(),
                        pos1.into(),
                        chrom2.into(),
                        pos2.into(),
                        strand1.into(),
                        strand2.into(),
                        Series::new("mapq".into(), pair_mapq).into(),
                    ])?;
                    sender
                        .send((shard_id, pairs, local_summary))
                        .map_err(|_| anyhow!("native porec2pairs writer pipeline disconnected"))?;
                }
                Ok(())
            }));
        }
        drop(raw_receiver);
        drop(pairs_sender);

        // Keep input I/O sequential and projected; pair expansion and output compression run
        // concurrently behind bounded queues so at most a few whole shards are resident.
        let projected_columns = [
            "read_idx",
            "read_start",
            "strand",
            "chrom",
            "start",
            "end",
            "mapping_quality",
        ]
        .into_iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
        let read_result =
            shards
                .iter()
                .enumerate()
                .try_for_each(|(shard_id, input_shard)| -> anyResult<()> {
                    let frame = ParquetReader::new(File::open(input_shard)?)
                        .with_columns(Some(projected_columns.clone()))
                        .finish()
                        .with_context(|| {
                            format!("failed to read concat PQS shard {}", input_shard.display())
                        })?;
                    raw_sender
                        .send((shard_id, frame))
                        .map_err(|_| anyhow!("native porec2pairs transform pipeline disconnected"))
                });
        drop(raw_sender);

        let mut transform_result = Ok(());
        for handle in transform_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(error)) if transform_result.is_ok() => transform_result = Err(error),
                Err(_) if transform_result.is_ok() => {
                    transform_result = Err(anyhow!("native porec2pairs transform worker panicked"))
                }
                _ => {}
            }
        }
        let mut writer_result = Ok(());
        for handle in writer_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(error)) if writer_result.is_ok() => writer_result = Err(error),
                Err(_) if writer_result.is_ok() => {
                    writer_result = Err(anyhow!("native porec2pairs writer worker panicked"))
                }
                _ => {}
            }
        }
        read_result?;
        transform_result?;
        writer_result?;

        let mut counts = common_writer(
            output_path
                .join("_metadata_counts")
                .to_string_lossy()
                .as_ref(),
        );
        writeln!(
            counts,
            "q0_records\t{}",
            q0_total.load(AtomOrdering::Relaxed)
        )?;
        writeln!(
            counts,
            "q1_records\t{}",
            q1_total.load(AtomOrdering::Relaxed)
        )?;
        counts.flush()?;
        let output_prefix = Path::new(output).with_extension("");
        concatemer_summary.lock().unwrap().save(&format!(
            "{}.concatemer.summary",
            output_prefix.to_string_lossy()
        ));
        log::info!(
            "Successful native pairs PQS output `{}`: {} q0 pairs across {} input shards",
            output,
            q0_total.load(AtomOrdering::Relaxed),
            shards.len()
        );
        Ok(())
    }

    pub fn to_pairs(
        &mut self,
        chromsizes: &String,
        output: &String,
        min_quality: u8,
        min_order: usize,
        max_order: usize,
        threads: usize,
    ) -> Result<(), Box<dyn Error>> {
        if threads == 0 {
            return Err("porec2pairs thread count must be at least 1".into());
        }
        let mut rdr = self.parse2().expect("Failed to open input file");
        log::info!(
            "Only retain concatemer that order in the range of [{}, {})",
            min_order,
            max_order
        );

        let chromsizes_obj = ChromSize::new(chromsizes);
        let chromsizes_data = chromsizes_obj.to_vec().unwrap();
        let mut ph = PairHeader::new();
        ph.from_chromsizes(chromsizes_data);

        let writer = common_writer(output);
        let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, writer);
        writer.write_all(ph.to_string().as_bytes()).unwrap();

        let (sender, receiver) = bounded::<(usize, Vec<PairConcatemer>, u64)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let num_workers = threads;
        let buffer_pool_size = num_workers.saturating_mul(2).max(1);
        let (buffer_sender, buffer_receiver) = bounded::<Vec<u8>>(buffer_pool_size);
        for _ in 0..buffer_pool_size {
            buffer_sender.send(Vec::with_capacity(1024 * 1024)).unwrap();
        }

        for _ in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let buffers = buffer_receiver.clone();
            thread::spawn(move || {
                while let Ok((chunk_id, batch, mut current_id)) = rx.recv() {
                    let mut local_buf = buffers.recv().unwrap();
                    local_buf.clear();
                    for mut concatemer in batch {
                        concatemer.sort();
                        let records = &concatemer.records;
                        let n = records.len();
                        for i in 0..n {
                            for j in i + 1..n {
                                current_id += 1;
                                let r1 = &records[i];
                                let r2 = &records[j];
                                let mapq = std::cmp::min(r1.mapq, r2.mapq);
                                let pos1 = (r1.target_start + r1.target_end) / 2;
                                let pos2 = (r2.target_start + r2.target_end) / 2;

                                writeln!(
                                    local_buf,
                                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                    current_id,
                                    r1.target.as_ref(),
                                    pos1,
                                    r2.target.as_ref(),
                                    pos2,
                                    r1.query_strand,
                                    r2.query_strand,
                                    mapq
                                )
                                .unwrap();
                            }
                        }
                    }
                    tx.send((chunk_id, local_buf)).unwrap();
                }
            });
        }
        drop(buffer_receiver);
        drop(out_sender);

        let mut global_pair_counter: u64 = 0;
        let write_handle = thread::spawn(move || {
            let mut pending = BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(mut data) = pending.remove(&next_chunk) {
                    writer.write_all(&data).unwrap();
                    data.clear();
                    let _ = buffer_sender.send(data);
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut line_buf = String::new();
        let mut concatemer = PairConcatemer::default();
        let mut batch = Vec::with_capacity(5000);
        let mut target_interner = HashSet::<Arc<str>>::new();
        let mut old_read_idx: u64 = 0;
        let mut first_iteration = true;
        let mut chunk_id = 0;
        let mut concatemer_summary = ConcatemerSummary::new();

        let mut global_pair_counter: u64 = 0;
        while rdr.read_line(&mut line_buf)? > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let mut parts = trimmed.split('\t');

            let read_idx = match parts.next() {
                Some(s) => s.parse::<u64>().unwrap_or(0),
                None => {
                    line_buf.clear();
                    continue;
                }
            };

            if !first_iteration && read_idx != old_read_idx {
                let order = concatemer.count();
                if order >= min_order && order < max_order {
                    concatemer_summary.count_order(order);
                    batch.push(std::mem::take(&mut concatemer));
                    if batch.len() >= 5000 {
                        let pairs_in_batch: u64 = batch
                            .iter()
                            .map(|c| {
                                let n = c.count() as u64;
                                n * (n - 1) / 2
                            })
                            .sum();
                        sender
                            .send((chunk_id, std::mem::take(&mut batch), global_pair_counter))
                            .unwrap();

                        global_pair_counter += pairs_in_batch;
                        chunk_id += 1;
                    }
                } else {
                    concatemer.clear();
                }
            }

            let _ = parts.next();
            let _ = parts.next();
            let _ = parts.next();
            let q_strand = parts.next().and_then(|s| s.chars().next()).unwrap_or('+');
            let target_str = parts.next().unwrap_or("");
            let t_start = parts
                .next()
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(0);
            let t_end = parts
                .next()
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(0);
            let mapq = parts.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);

            first_iteration = false;
            old_read_idx = read_idx;

            if mapq >= min_quality {
                let target = if let Some(target) = target_interner.get(target_str) {
                    Arc::clone(target)
                } else {
                    let target = Arc::<str>::from(target_str);
                    target_interner.insert(Arc::clone(&target));
                    target
                };
                concatemer.records.push(PairAnchor {
                    target,
                    target_start: t_start,
                    target_end: t_end,
                    query_strand: q_strand,
                    mapq,
                });
            }

            line_buf.clear();
        }

        if !batch.is_empty() {
            sender.send((chunk_id, batch, global_pair_counter)).unwrap();
        }

        drop(sender);
        write_handle.join().unwrap();

        let output_prefix = Path::new(output).with_extension("");
        concatemer_summary.save(&format!(
            "{}.concatemer.summary",
            output_prefix.to_str().unwrap()
        ));

        log::info!("Successful output pairs `{}`", output);
        Ok(())
    }

    pub fn intersect(&mut self, hcr_bed: &String, invert: bool, output: &String) {
        type IvU8 = Interval<usize, u8>;
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        let writer = common_writer(output);
        let mut wtr = csv::WriterBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_writer(writer);

        for (i, line) in self.parse().unwrap().records().enumerate() {
            let record = match line {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}", i + 1);
                    continue;
                }
            };

            let target_start = record[6].parse::<usize>().unwrap();
            let target_end = record[7].parse::<usize>().unwrap();

            let is_in_regions = interval_hash.get(&record[5]).map_or(false, |interval| {
                interval.count(target_start, target_end) > 0
            });

            if is_in_regions ^ invert {
                let _ = wtr.write_record(&record);
            }
        }

        log::info!(
            "Successful output intersection porec table into `{}`",
            output
        );
    }

    pub fn intersect_concat_pqs_native(
        &self,
        hcr_bed: &String,
        invert: bool,
        output: &String,
        threads: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;

        let input = Path::new(&self.file);
        if !is_concat_pqs(input) {
            bail!("native concat PQS intersect requires concat.pqs input");
        }
        let intervals = Arc::new(Bed3::new(hcr_bed).to_interval_hash());
        let mut sizes = BTreeMap::new();
        read_concat_pqs_contigsizes(input, &mut sizes)?;
        transform_concat_pqs_native(
            input,
            Path::new(output.trim_end_matches('/')),
            sizes.into_iter().collect(),
            threads,
            true,
            move |_shard_id, frame| {
                let chrom_series = frame
                    .column("chrom")?
                    .as_materialized_series()
                    .cast(&DataType::String)?;
                let chrom = chrom_series.str()?;
                let start_series = frame
                    .column("start")?
                    .as_materialized_series()
                    .cast(&DataType::UInt64)?;
                let start = start_series.u64()?;
                let end_series = frame
                    .column("end")?
                    .as_materialized_series()
                    .cast(&DataType::UInt64)?;
                let end = end_series.u64()?;
                let mut keep = Vec::with_capacity(frame.height());
                for row in 0..frame.height() {
                    let target = chrom
                        .get(row)
                        .ok_or_else(|| anyhow!("null chrom at row {row}"))?;
                    let target_start = usize::try_from(
                        start
                            .get(row)
                            .ok_or_else(|| anyhow!("null start at row {row}"))?,
                    )?;
                    let target_end = usize::try_from(
                        end.get(row)
                            .ok_or_else(|| anyhow!("null end at row {row}"))?,
                    )?;
                    let overlaps = intervals
                        .get(target)
                        .is_some_and(|regions| regions.count(target_start, target_end) > 0);
                    keep.push(overlaps ^ invert);
                }
                Ok(frame.filter(&BooleanChunked::from_slice("keep".into(), &keep))?)
            },
        )?;
        log::info!(
            "Successful native concat PQS intersection into `{}`",
            output
        );
        Ok(())
    }

    pub fn intersect_multi_threads(&mut self, hcr_bed: &String, invert: bool, output: &String) {
        type IvU8 = Interval<usize, u8>;
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        let wtr = common_writer(output);

        let (sender, receiver) = bounded::<Vec<String>>(1000);

        let mut handles = vec![];
        let wtr = Arc::new(Mutex::new(wtr));

        for _ in 0..10 {
            let interval_hash = interval_hash.clone();
            let wtr = Arc::clone(&wtr);
            let receiver = receiver.clone();
            handles.push(thread::spawn(move || {
                while let Ok(records) = receiver.recv() {
                    let data = records
                        .par_iter()
                        .filter_map(|record| {
                            let record = record.split("\t").collect::<Vec<_>>();
                            let target_start = record[6].parse::<usize>().unwrap();
                            let target_end = record[7].parse::<usize>().unwrap();

                            let is_in_regions =
                                interval_hash.get(record[5]).map_or(false, |interval| {
                                    interval.count(target_start, target_end) > 0
                                });

                            if is_in_regions ^ invert {
                                Some(record.iter().join("\t"))
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>();

                    if !data.is_empty() {
                        let mut wtr = wtr.lock().unwrap();
                        let data = data.join("\n") + "\n";
                        wtr.write_all(data.as_bytes()).unwrap();
                    }
                }
            }));
        }

        let batch_size = 10_000;
        let mut batch = Vec::with_capacity(batch_size);
        for (idx, record) in self.parse2().unwrap().lines().enumerate() {
            let record = match record {
                Ok(v) => v,
                Err(error) => {
                    log::warn!("Could not parse line {}", idx + 1);
                    continue;
                }
            };
            batch.push(record);
            if batch.len() == batch_size {
                sender.send(std::mem::take(&mut batch)).unwrap();
            }
        }

        drop(sender);

        for handle in handles {
            handle.join().unwrap();
        }

        log::info!(
            "Successful output intersection porec table into `{}`",
            output
        );
    }

    pub fn intersect_multi_threads_coitree(
        &mut self,
        hcr_bed: &String,
        invert: bool,
        output: &String,
    ) {
        log::info!("Building COI-Trees from BED...");
        let bed = Bed3::new(hcr_bed);
        let lapper_map = bed.to_interval_hash();

        let mut tree_map: HashMap<String, COITree<u8, u32>> = HashMap::new();
        for (chrom, lapper) in lapper_map {
            let nodes: Vec<IntervalNode<u8, u32>> = lapper
                .intervals
                .into_iter()
                .map(|iv| IntervalNode::new(iv.start as i32, iv.stop as i32, iv.val))
                .collect();

            tree_map.insert(chrom, COITree::new(&nodes));
        }
        let interval_hash = Arc::new(tree_map);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let num_workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8);
        log::info!(
            "Intersecting Bed in parallel using {} workers...",
            num_workers
        );

        let mut handles = vec![];
        for _ in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let ih = Arc::clone(&interval_hash);

            handles.push(thread::spawn(move || {
                let mut local_buf = Vec::with_capacity(2 * 1024 * 1024);
                let mut tab_indices = [0usize; 8];

                while let Ok((chunk_id, records)) = rx.recv() {
                    local_buf.clear();
                    for record in records {
                        let bytes = record.as_bytes();
                        let len = bytes.len();

                        let mut tab_count = 0;
                        for idx in 0..len {
                            if bytes[idx] == b'\t' {
                                if tab_count < 8 {
                                    tab_indices[tab_count] = idx;
                                    tab_count += 1;
                                } else {
                                    break;
                                }
                            }
                        }

                        if tab_count < 8 {
                            continue;
                        }

                        let target = &record[tab_indices[4] + 1..tab_indices[5]];
                        let start_str = &record[tab_indices[5] + 1..tab_indices[6]];
                        let end_str = &record[tab_indices[6] + 1..tab_indices[7]];

                        let target_start = match start_str.parse::<i32>() {
                            Ok(val) => val,
                            Err(_) => continue,
                        };
                        let target_end = match end_str.parse::<i32>() {
                            Ok(val) => val,
                            Err(_) => continue,
                        };

                        let is_in_regions = ih.get(target).map_or(false, |tree| {
                            let mut has_overlap = false;
                            tree.query(target_start, target_end, |_| {
                                has_overlap = true;
                            });
                            has_overlap
                        });

                        if is_in_regions ^ invert {
                            local_buf.extend_from_slice(record.as_bytes());
                            local_buf.push(b'\n');
                        }
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender);

        let out_path = output.clone();
        let write_handle = thread::spawn(move || {
            let wtr_file = common_writer(&out_path);
            let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, wtr_file);
            let mut pending = BTreeMap::new();
            let mut next_chunk = 0;

            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    if !data.is_empty() {
                        writer.write_all(&data).unwrap();
                    }
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut rdr = self
            .parse2()
            .expect("Failed to open porec table for reading");
        let batch_size = 10_000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;
        let mut line_buf = String::new();

        while rdr.read_line(&mut line_buf).unwrap_or(0) > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }
            batch.push(std::mem::take(&mut line_buf));
            if batch.len() >= batch_size {
                sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                batch = Vec::with_capacity(batch_size);
                chunk_id += 1;
            }
        }
        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles {
            h.join().unwrap();
        }
        write_handle.join().unwrap();

        log::info!(
            "Successful output intersection porec table into `{}`",
            output
        );
    }

    fn target_sizes_for_pqs(&self, chromsizes: Option<&String>) -> anyResult<HashMap<String, u64>> {
        let mut sizes = BTreeMap::new();
        if let Some(chromsizes) = chromsizes {
            read_contigsizes_file(Path::new(chromsizes), &mut sizes)?;
        } else if Path::new(&self.file).is_dir() {
            read_concat_pqs_contigsizes(Path::new(&self.file), &mut sizes)?;
        } else {
            bail!(
                "PQS output from a text Pore-C table requires --chromsizes; input was {}",
                self.file
            );
        }
        Ok(sizes.into_iter().collect())
    }

    fn break_contigs_pqs(
        &mut self,
        break_bed: &String,
        output: &String,
        threads: usize,
        chunksize: usize,
        chromsizes: Option<&String>,
    ) -> anyResult<()> {
        type IvString = Interval<usize, String>;
        let mut target_sizes = self.target_sizes_for_pqs(chromsizes)?;
        for record in Bed4::new(break_bed) {
            let length = record.end.saturating_sub(record.start).saturating_add(1) as u64;
            if let Some(previous) = target_sizes.insert(record.gene.clone(), length) {
                if previous != length {
                    bail!(
                        "conflicting lengths for corrected contig {}: {} and {}",
                        record.gene,
                        previous,
                        length
                    );
                }
            }
        }
        let interval_hash = Bed4::new(break_bed).to_interval_hash();
        if is_concat_pqs(Path::new(&self.file)) {
            return self.break_contigs_pqs_native(interval_hash, target_sizes, output, threads);
        }
        let mut lines = self.parse2()?.lines().enumerate();
        let batches = std::iter::from_fn(move || {
            let mut records = Vec::with_capacity(10_000);
            while records.len() < 10_000 {
                let (line_number, line_result) = match lines.next() {
                    Some(line) => line,
                    None => return (!records.is_empty()).then_some(records),
                };
                let line = match line_result {
                    Ok(line) => line,
                    Err(error) => {
                        log::warn!("Failed to read Pore-C line {}: {}", line_number + 1, error);
                        continue;
                    }
                };
                let mut record = match parse_concat_porec_line(&line, line_number + 1) {
                    Ok(Some(record)) => record,
                    Ok(None) => continue,
                    Err(error) => {
                        log::warn!("{}", error);
                        continue;
                    }
                };
                if let Some(intervals) = interval_hash.get(&record.target) {
                    if let Some(iv) = intervals
                        .find(record.target_start as usize, record.target_end as usize)
                        .next()
                    {
                        let break_length = iv.stop.saturating_sub(iv.start).saturating_add(1);
                        if record.target_start as usize >= iv.start {
                            let new_start = record.target_start - iv.start as u64 + 1;
                            let new_end = record.target_end - iv.start as u64 + 1;
                            if new_end <= break_length as u64 {
                                record.target = iv.val.clone();
                                record.target_start = new_start;
                                record.target_end = new_end;
                            }
                        }
                    }
                }
                records.push(record);
            }
            Some(records)
        });
        write_concat_pqs_record_batches(
            batches,
            Arc::new(Mutex::new(target_sizes)),
            output.trim_end_matches('/'),
            chunksize,
            threads,
        )
    }

    fn break_contigs_pqs_native(
        &self,
        interval_hash: HashMap<String, Lapper<usize, String>>,
        target_sizes: HashMap<String, u64>,
        output: &String,
        threads: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;

        let input_path = Path::new(&self.file);
        let output_path = Path::new(output.trim_end_matches('/'));
        if output_path.exists() {
            bail!(
                "concat PQS output already exists: {}",
                output_path.display()
            );
        }
        let shards = concat_pqs_files(input_path)?;
        std::fs::create_dir_all(output_path.join("q0"))?;
        std::fs::create_dir_all(output_path.join("q1"))?;
        polars::enable_string_cache();
        let broken_cn = interval_hash
            .iter()
            .map(|(source, intervals)| {
                (
                    source.clone(),
                    intervals
                        .iter()
                        .map(|interval| interval.val.clone())
                        .collect(),
                )
            })
            .collect::<BTreeMap<String, std::collections::HashSet<String>>>();
        let intervals = Arc::new(interval_hash);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("failed to create native porec-break thread pool")?;

        log::info!(
            "Native concat PQS break: {} q0 shards with {} workers (no text conversion)",
            shards.len(),
            threads
        );
        pool.install(|| {
            shards
                .par_iter()
                .enumerate()
                .try_for_each(|(shard_id, input_shard)| -> anyResult<()> {
                    let mut frame = ParquetReader::new(File::open(input_shard)?)
                        .finish()
                        .with_context(|| {
                            format!("failed to read concat PQS shard {}", input_shard.display())
                        })?;
                    let chrom = frame.column("chrom")?.categorical()?;
                    if chrom.null_count() != 0 {
                        bail!("porec-break does not support null contig identifiers");
                    }
                    let chrom_codes = chrom.physical();
                    let chrom_rev_map = chrom.get_rev_map();
                    let max_input_code = chrom_codes.max().unwrap_or(0) as usize;

                    let mut output_names = Vec::new();
                    for code in 0..=max_input_code {
                        if let Some(name) = chrom_rev_map.get_optional(code as u32) {
                            output_names.push(name.to_string());
                        }
                    }
                    for source_intervals in intervals.values() {
                        output_names
                            .extend(source_intervals.iter().map(|interval| interval.val.clone()));
                    }
                    output_names.sort_unstable();
                    output_names.dedup();

                    let mut builder = CategoricalChunkedBuilder::new(
                        "chrom".into(),
                        output_names.len(),
                        CategoricalOrdering::Physical,
                    );
                    for name in &output_names {
                        builder.append_value(name);
                    }
                    let categorical_template = builder.finish();
                    let registered_codes: Vec<u32> = categorical_template
                        .physical()
                        .into_no_null_iter()
                        .collect();
                    let output_codes: HashMap<String, u32> = output_names
                        .iter()
                        .cloned()
                        .zip(registered_codes.into_iter())
                        .collect();

                    let mut unchanged_codes = vec![None; max_input_code + 1];
                    let mut encoded_intervals = vec![None; max_input_code + 1];
                    for code in 0..=max_input_code {
                        let Some(name) = chrom_rev_map.get_optional(code as u32) else {
                            continue;
                        };
                        unchanged_codes[code] = output_codes.get(name).copied();
                        if let Some(source_intervals) = intervals.get(name) {
                            encoded_intervals[code] = Some(
                                source_intervals
                                    .iter()
                                    .map(|interval| {
                                        (
                                            interval.start,
                                            interval.stop,
                                            *output_codes
                                                .get(&interval.val)
                                                .expect("break fragment was not registered"),
                                        )
                                    })
                                    .collect::<Vec<_>>(),
                            );
                        }
                    }
                    let original_position_dtype = frame.column("start")?.dtype().clone();
                    if frame.column("end")?.dtype() != &original_position_dtype {
                        bail!("concat PQS start/end columns use different dtypes");
                    }

                    macro_rules! remap_typed_positions {
                        ($position_type:ty, $start:expr, $end:expr) => {{
                            let start = $start;
                            let end = $end;
                            let mut output_chroms = Vec::<u32>::with_capacity(frame.height());
                            let mut output_starts =
                                Vec::<$position_type>::with_capacity(frame.height());
                            let mut output_ends =
                                Vec::<$position_type>::with_capacity(frame.height());
                            for ((source_code, source_start), source_end) in chrom_codes
                                .into_no_null_iter()
                                .zip(start.into_no_null_iter())
                                .zip(end.into_no_null_iter())
                            {
                                let source_start_usize = source_start as usize;
                                let source_end_usize = source_end as usize;
                                let source_code_index = source_code as usize;
                                let mut target_code = unchanged_codes[source_code_index]
                                    .expect("input contig was not registered");
                                let mut target_start = source_start;
                                let mut target_end = source_end;
                                if let Some(source_intervals) =
                                    encoded_intervals[source_code_index].as_ref()
                                {
                                    if let Some((interval_start, interval_stop, output_code)) =
                                        source_intervals.iter().find(|(start, stop, _)| {
                                            source_start_usize >= *start
                                                && source_start_usize < *stop
                                                && source_end_usize > *start
                                                && source_end_usize <= *stop
                                        })
                                    {
                                        target_code = *output_code;
                                        target_start = (source_start_usize - *interval_start + 1)
                                            as $position_type;
                                        target_end = (source_end_usize - *interval_start + 1)
                                            as $position_type;
                                    }
                                }
                                output_chroms.push(target_code);
                                output_starts.push(target_start);
                                output_ends.push(target_end);
                            }

                            let output_physical =
                                UInt32Chunked::from_vec("chrom".into(), output_chroms);
                            let output_chroms = unsafe {
                                CategoricalChunked::from_cats_and_rev_map_unchecked(
                                    output_physical,
                                    categorical_template.get_rev_map().clone(),
                                    false,
                                    CategoricalOrdering::Physical,
                                )
                            };
                            (
                                output_chroms.into_series(),
                                Series::new("start".into(), output_starts),
                                Series::new("end".into(), output_ends),
                            )
                        }};
                    }

                    let (output_chroms, output_starts, output_ends) = match original_position_dtype
                    {
                        DataType::UInt32 => remap_typed_positions!(
                            u32,
                            frame.column("start")?.as_materialized_series().u32()?,
                            frame.column("end")?.as_materialized_series().u32()?
                        ),
                        DataType::UInt64 => remap_typed_positions!(
                            u64,
                            frame.column("start")?.as_materialized_series().u64()?,
                            frame.column("end")?.as_materialized_series().u64()?
                        ),
                        ref dtype => {
                            bail!("unsupported concat PQS position dtype for porec-break: {dtype}")
                        }
                    };

                    frame.replace("chrom", output_chroms)?;
                    frame.replace("start", output_starts)?;
                    frame.replace("end", output_ends)?;

                    let mapq = frame
                        .column("mapping_quality")?
                        .as_materialized_series()
                        .u8()?;
                    let mut q1_frame = frame.filter(&mapq.gt_eq(1))?;
                    let q0_path = output_path.join(format!("q0/{shard_id}.parquet"));
                    ParquetWriter::new(File::create(&q0_path)?)
                        .finish(&mut frame)
                        .with_context(|| format!("failed to write {}", q0_path.display()))?;
                    if q1_frame.height() > 0 {
                        let q1_path = output_path.join(format!("q1/{shard_id}.parquet"));
                        ParquetWriter::new(File::create(&q1_path)?)
                            .finish(&mut q1_frame)
                            .with_context(|| format!("failed to write {}", q1_path.display()))?;
                    }
                    Ok(())
                })
        })?;

        write_concat_pqs_sidecars(
            input_path,
            output_path,
            target_sizes,
            Some(input_path.join("_metadata_counts")),
            true,
        )?;
        write_broken_cn_info(input_path, &broken_cn, output_path)?;
        log::info!(
            "Successful native concat PQS break into `{}` across {} shards",
            output,
            shards.len()
        );
        Ok(())
    }

    pub fn break_contigs(
        &mut self,
        break_bed: &String,
        output: &String,
        threads: usize,
        chunksize: usize,
        chromsizes: Option<&String>,
    ) -> anyResult<()> {
        if is_concat_pqs_path(output) {
            return self.break_contigs_pqs(break_bed, output, threads, chunksize, chromsizes);
        }
        type IvString = Interval<usize, String>;
        let bed = Bed4::new(break_bed);

        let interval_hash = Arc::new(bed.to_interval_hash());

        let wtr_file = common_writer(output);
        let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, wtr_file);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let num_workers = threads;
        let mut handles = vec![];

        for _ in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let ih = Arc::clone(&interval_hash);

            handles.push(thread::spawn(move || {
                let mut local_buf = Vec::with_capacity(1024 * 1024);
                while let Ok((chunk_id, batch)) = rx.recv() {
                    local_buf.clear();
                    for line in batch {
                        let trimmed = line.trim_end();
                        let fields: Vec<&str> = trimmed.split('\t').collect();

                        if fields.len() < 8 {
                            local_buf.extend_from_slice(line.as_bytes());
                            local_buf.push(b'\n');
                            continue;
                        }

                        let target = fields[5];
                        let mut processed = false;

                        if let Some(intervals) = ih.get(target) {
                            if let (Ok(t_start), Ok(t_end)) =
                                (fields[6].parse::<usize>(), fields[7].parse::<usize>())
                            {
                                if let Some(iv) = intervals.find(t_start, t_end).next() {
                                    let break_contig_length = iv.stop - iv.start + 1;

                                    if t_start >= iv.start {
                                        let new_target_start = t_start - iv.start + 1;
                                        let new_target_end = t_end - iv.start + 1;

                                        if new_target_end <= break_contig_length {
                                            for i in 0..5 {
                                                local_buf.extend_from_slice(fields[i].as_bytes());
                                                local_buf.push(b'\t');
                                            }
                                            local_buf.extend_from_slice(iv.val.as_bytes());
                                            local_buf.push(b'\t');

                                            let _ = write!(
                                                local_buf,
                                                "{}\t{}\t",
                                                new_target_start, new_target_end
                                            );
                                            for i in 8..=10 {
                                                if i < fields.len() {
                                                    local_buf
                                                        .extend_from_slice(fields[i].as_bytes());
                                                }
                                                if i < 10 {
                                                    local_buf.push(b'\t');
                                                }
                                            }
                                            local_buf.push(b'\n');
                                            processed = true;
                                        }
                                    }
                                }
                            }
                        }

                        if !processed {
                            local_buf.extend_from_slice(line.as_bytes());
                            local_buf.push(b'\n');
                        }
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender);

        let write_handle = thread::spawn(move || {
            let mut pending = std::collections::BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    writer.write_all(&data).unwrap();
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut rdr = self.parse2().expect("Failed to open table for reading");
        let batch_size = 5000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;

        for line_res in rdr.lines() {
            if let Ok(line) = line_res {
                batch.push(line);
                if batch.len() >= batch_size {
                    sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                    batch = Vec::with_capacity(batch_size);
                    chunk_id += 1;
                }
            }
        }
        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles {
            h.join().unwrap();
        }
        write_handle.join().unwrap();

        log::info!(
            "Successful output contigs corrected porec table into `{}`",
            output
        );
        Ok(())
    }

    fn chr_porec_to_contig_pqs_native(
        &self,
        contig_bed: &String,
        output: &String,
        threads: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;

        #[derive(Clone)]
        struct ContigInterval {
            start: u64,
            end: u64,
            contig: String,
        }

        let input = Path::new(&self.file);
        if !is_concat_pqs(input) {
            bail!("native porec-chr2ctg PQS output requires concat.pqs input");
        }
        let reader = BufReader::new(common_reader(contig_bed));
        let mut mapping = HashMap::<String, Vec<ContigInterval>>::new();
        let mut cn_mapping = BTreeMap::<String, std::collections::HashSet<String>>::new();
        let mut target_sizes = HashMap::<String, u64>::new();
        for (line_number, line) in reader.lines().enumerate() {
            let line = line?;
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.is_empty() || fields[0].starts_with('#') {
                continue;
            }
            if fields.len() < 4 {
                bail!(
                    "chr2ctg BED line {} has fewer than four columns",
                    line_number + 1
                );
            }
            let start = fields[1].parse::<u64>()?;
            let end = fields[2].parse::<u64>()?;
            if end <= start {
                bail!("chr2ctg BED line {} has an empty interval", line_number + 1);
            }
            let contig = fields[3].to_string();
            let length = end - start;
            if let Some(previous) = target_sizes.insert(contig.clone(), length) {
                if previous != length {
                    bail!("conflicting chr2ctg lengths for {contig}: {previous} and {length}");
                }
            }
            mapping
                .entry(fields[0].to_string())
                .or_default()
                .push(ContigInterval { start, end, contig });
            cn_mapping
                .entry(fields[0].to_string())
                .or_default()
                .insert(fields[3].to_string());
        }
        for intervals in mapping.values_mut() {
            intervals.sort_by_key(|interval| interval.start);
        }
        let mapping = Arc::new(mapping);
        transform_concat_pqs_native(
            input,
            Path::new(output.trim_end_matches('/')),
            target_sizes,
            threads,
            false,
            move |_shard_id, mut frame| {
                let chrom_series = frame
                    .column("chrom")?
                    .as_materialized_series()
                    .cast(&DataType::String)?;
                let chrom = chrom_series.str()?;
                let start_dtype = frame.column("start")?.dtype().clone();
                let start_series = frame
                    .column("start")?
                    .as_materialized_series()
                    .cast(&DataType::UInt64)?;
                let start = start_series.u64()?;
                let end_series = frame
                    .column("end")?
                    .as_materialized_series()
                    .cast(&DataType::UInt64)?;
                let end = end_series.u64()?;
                let mut keep = Vec::with_capacity(frame.height());
                let mut output_chrom = Vec::with_capacity(frame.height());
                let mut output_start = Vec::with_capacity(frame.height());
                let mut output_end = Vec::with_capacity(frame.height());
                for row in 0..frame.height() {
                    let source = chrom
                        .get(row)
                        .ok_or_else(|| anyhow!("null chrom at row {row}"))?;
                    let source_start = start
                        .get(row)
                        .ok_or_else(|| anyhow!("null start at row {row}"))?;
                    let source_end = end
                        .get(row)
                        .ok_or_else(|| anyhow!("null end at row {row}"))?;
                    let interval = (source_start < source_end)
                        .then(|| mapping.get(source))
                        .flatten()
                        .and_then(|intervals| {
                            intervals.iter().find(|interval| {
                                source_start >= interval.start && source_end <= interval.end
                            })
                        });
                    if let Some(interval) = interval {
                        keep.push(true);
                        output_chrom.push(interval.contig.clone());
                        output_start.push(source_start - interval.start);
                        output_end.push(source_end - interval.start);
                    } else {
                        keep.push(false);
                        output_chrom.push(String::new());
                        output_start.push(0);
                        output_end.push(0);
                    }
                }
                frame.replace(
                    "chrom",
                    Series::new("chrom".into(), output_chrom)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?,
                )?;
                frame.replace(
                    "start",
                    Series::new("start".into(), output_start).cast(&start_dtype)?,
                )?;
                frame.replace(
                    "end",
                    Series::new("end".into(), output_end).cast(&start_dtype)?,
                )?;
                Ok(frame.filter(&BooleanChunked::from_slice("keep".into(), &keep))?)
            },
        )?;
        write_remapped_cn_info(input, &cn_mapping, Path::new(output.trim_end_matches('/')))?;
        log::info!(
            "Successful native chromosome-to-contig concat PQS into `{}`",
            output
        );
        Ok(())
    }

    pub fn chr_porec_to_contig_porec(
        &mut self,
        contig_bed: &String,
        output: &String,
        threads: usize,
    ) -> anyResult<()> {
        if is_concat_pqs_path(output) {
            return self.chr_porec_to_contig_pqs_native(contig_bed, output, threads);
        }
        #[derive(Debug, Clone)]
        struct ContigInterval {
            start: u64,
            end: u64,
            contig: String,
        }

        log::info!("Loading contig mapping BED from `{}`...", contig_bed);
        let f = common_reader(contig_bed);
        let rdr = std::io::BufReader::new(f);
        let mut chrom_to_intervals: HashMap<String, Vec<ContigInterval>> = HashMap::new();

        for line in rdr.lines().flatten() {
            let s = line.trim();
            if s.is_empty() || s.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = s.split_whitespace().collect();
            if fields.len() < 4 {
                continue;
            }
            let chrom = fields[0].to_string();
            let start: u64 = fields[1].parse().unwrap_or(0);
            let end: u64 = fields[2].parse().unwrap_or(0);
            let contig = fields[3].to_string();

            chrom_to_intervals
                .entry(chrom)
                .or_default()
                .push(ContigInterval { start, end, contig });
        }

        for ivs in chrom_to_intervals.values_mut() {
            ivs.sort_by_key(|iv| iv.start);
        }

        let chrom_intervals = Arc::new(chrom_to_intervals);
        let wtr_file = common_writer(output);
        let mut writer = std::io::BufWriter::with_capacity(1024 * 1024, wtr_file);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        log::info!("Converting chromosome-level Pore-C alignments to contig-level in parallel...");
        let mut handles = vec![];

        for _ in 0..threads {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let ci = Arc::clone(&chrom_intervals);

            handles.push(thread::spawn(move || {
                let mut local_buf = Vec::with_capacity(1024 * 1024);
                while let Ok((chunk_id, batch)) = rx.recv() {
                    local_buf.clear();
                    for line in batch {
                        let trimmed = line.trim_end();
                        let fields: Vec<&str> = trimmed.split('\t').collect();
                        if fields.len() < 8 {
                            continue;
                        }

                        let target = fields[5];
                        let t_start: u64 = fields[6].parse().unwrap_or(0);
                        let t_end: u64 = fields[7].parse().unwrap_or(0);

                        if t_start >= t_end {
                            continue;
                        }

                        let mid_0based = t_start + (t_end - t_start) / 2;

                        let mut mapped = false;
                        if let Some(ivs) = ci.get(target) {
                            let idx_res = ivs.binary_search_by(|iv| {
                                if mid_0based < iv.start {
                                    std::cmp::Ordering::Greater
                                } else if mid_0based >= iv.end {
                                    std::cmp::Ordering::Less
                                } else {
                                    std::cmp::Ordering::Equal
                                }
                            });

                            if let Ok(idx) = idx_res {
                                let iv = &ivs[idx];
                                if t_start >= iv.start && t_end <= iv.end {
                                    let new_start = t_start - iv.start;
                                    let new_end = t_end - iv.start;

                                    for i in 0..5 {
                                        local_buf.extend_from_slice(fields[i].as_bytes());
                                        local_buf.push(b'\t');
                                    }
                                    local_buf.extend_from_slice(iv.contig.as_bytes());
                                    local_buf.push(b'\t');
                                    let _ = write!(local_buf, "{}\t{}\t", new_start, new_end);

                                    for i in 8..fields.len() {
                                        local_buf.extend_from_slice(fields[i].as_bytes());
                                        if i < fields.len() - 1 {
                                            local_buf.push(b'\t');
                                        }
                                    }
                                    local_buf.push(b'\n');
                                    mapped = true;
                                }
                            }
                        }
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender);

        let write_handle = thread::spawn(move || {
            let mut pending = std::collections::BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    writer.write_all(&data).unwrap();
                    next_chunk += 1;
                }
            }
            writer.flush().unwrap();
        });

        let mut rdr = self.parse2().expect("Failed to open table for reading");
        let batch_size = 5000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;

        for line_res in rdr.lines() {
            if let Ok(line) = line_res {
                batch.push(line);
                if batch.len() >= batch_size {
                    sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                    batch = Vec::with_capacity(batch_size);
                    chunk_id += 1;
                }
            }
        }
        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles {
            h.join().unwrap();
        }
        write_handle.join().unwrap();

        log::info!(
            "Successfully generated contig-level porec table into `{}`",
            output
        );
        Ok(())
    }

    pub fn dup_v0_3_0(&mut self, collapsed_list: &String, seed: usize, output: &String) {
        let reader = common_reader(collapsed_list);
        let mut collapsed_contigs: HashMap<String, Vec<String>> = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let s: Vec<&str> = record.split("\t").collect();
            if s.len() != 2 {
                log::warn!("Invalid record: {}", record);
                continue;
            }
            let contig1 = s[0].to_string();
            let contig2 = s[1].to_string();
            collapsed_contigs
                .entry(contig1.clone())
                .or_insert(vec![contig1])
                .push(contig2);
        }

        let mut seed_array = [0u8; 32];
        let seed_bytes = seed.to_ne_bytes();
        seed_array[..8].copy_from_slice(&seed_bytes);
        let mut rng = StdRng::from_seed(seed_array);

        let reader = open_porec_reader(&self.file).expect("Failed to open Pore-C input");
        let mut wtr = common_writer(output);

        let mut current_read_idx = String::new();
        let mut read_decisions: HashMap<String, String> = HashMap::new();

        for record in reader.lines() {
            let record = record.unwrap();
            let trimmed = record.trim();
            if trimmed.is_empty() {
                continue;
            }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 6 {
                writeln!(wtr, "{}", trimmed).unwrap();
                continue;
            }

            let read_idx = fields[0];
            let contig = fields[5];

            if read_idx != current_read_idx {
                current_read_idx = read_idx.to_string();
                read_decisions.clear();
            }

            let mut final_contig = contig;

            if let Some(candidates) = collapsed_contigs.get(contig) {
                final_contig = read_decisions.entry(contig.to_string()).or_insert_with(|| {
                    let idx = rng.gen_range(0..candidates.len());
                    candidates[idx].clone()
                });
            }

            for (i, field) in fields.iter().enumerate() {
                if i > 0 {
                    write!(wtr, "\t").unwrap();
                }
                if i == 5 {
                    write!(wtr, "{}", final_contig).unwrap();
                } else {
                    write!(wtr, "{}", field).unwrap();
                }
            }
            writeln!(wtr).unwrap();
        }
    }

    fn dup_pqs(
        &mut self,
        collapsed_list: &String,
        seed: usize,
        output: &String,
        threads: usize,
        chunksize: usize,
        chromsizes: Option<&String>,
    ) -> anyResult<()> {
        let reader = common_reader(collapsed_list);
        let mut collapsed_contigs: HashMap<String, Vec<String>> = HashMap::new();
        let mut target_sizes = self.target_sizes_for_pqs(chromsizes)?;
        for (line_number, line) in reader.lines().enumerate() {
            let line = line?;
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() != 2 {
                log::warn!(
                    "Invalid collapsed-contig record at line {}: {}",
                    line_number + 1,
                    line
                );
                continue;
            }
            let source = fields[0].to_string();
            let duplicate = fields[1].to_string();
            let source_length = target_sizes.get(&source).copied().ok_or_else(|| {
                anyhow!("collapsed source contig {source:?} is missing from chromsizes")
            })?;
            if let Some(previous) = target_sizes.insert(duplicate.clone(), source_length) {
                if previous != source_length {
                    bail!(
                        "conflicting length for duplicated contig {}: {} and {}",
                        duplicate,
                        previous,
                        source_length
                    );
                }
            }
            collapsed_contigs
                .entry(source.clone())
                .or_insert_with(|| vec![source])
                .push(duplicate);
        }

        if is_concat_pqs(Path::new(&self.file)) {
            return self.dup_concat_pqs_native(
                collapsed_contigs,
                target_sizes,
                seed,
                output,
                threads,
            );
        }

        let mut seed_array = [0u8; 32];
        let seed_bytes = seed.to_ne_bytes();
        seed_array[..seed_bytes.len()].copy_from_slice(&seed_bytes);
        let mut rng = StdRng::from_seed(seed_array);
        let mut lines = self.parse2()?.lines().enumerate();
        let mut current_read_idx = None;
        let mut read_decisions = HashMap::<String, String>::new();
        let batches = std::iter::from_fn(move || {
            let mut records = Vec::with_capacity(10_000);
            while records.len() < 10_000 {
                let (line_number, line_result) = match lines.next() {
                    Some(line) => line,
                    None => return (!records.is_empty()).then_some(records),
                };
                let line = match line_result {
                    Ok(line) => line,
                    Err(error) => {
                        log::warn!("Failed to read Pore-C line {}: {}", line_number + 1, error);
                        continue;
                    }
                };
                let mut record = match parse_concat_porec_line(&line, line_number + 1) {
                    Ok(Some(record)) => record,
                    Ok(None) => continue,
                    Err(error) => {
                        log::warn!("{}", error);
                        continue;
                    }
                };
                if current_read_idx != Some(record.read_idx) {
                    current_read_idx = Some(record.read_idx);
                    read_decisions.clear();
                }
                if let Some(candidates) = collapsed_contigs.get(&record.target) {
                    let selected = read_decisions
                        .entry(record.target.clone())
                        .or_insert_with(|| candidates[rng.gen_range(0..candidates.len())].clone());
                    record.target = selected.clone();
                }
                records.push(record);
            }
            Some(records)
        });
        write_concat_pqs_record_batches(
            batches,
            Arc::new(Mutex::new(target_sizes)),
            output.trim_end_matches('/'),
            chunksize,
            threads,
        )
    }

    fn dup_concat_pqs_native(
        &self,
        collapsed_contigs: HashMap<String, Vec<String>>,
        target_sizes: HashMap<String, u64>,
        seed: usize,
        output: &String,
        threads: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;

        if threads == 0 {
            bail!("porec-dup thread count must be at least 1");
        }
        let input_path = Path::new(&self.file);
        let output_path = Path::new(output.trim_end_matches('/'));
        if output_path.exists() {
            bail!(
                "concat PQS output already exists: {}",
                output_path.display()
            );
        }

        let shards = concat_pqs_files(input_path)?;
        std::fs::create_dir_all(output_path.join("q0"))?;
        std::fs::create_dir_all(output_path.join("q1"))?;
        polars::enable_string_cache();

        let collapsed_contigs = Arc::new(collapsed_contigs);
        let q0_records = AtomicU64::new(0);
        let q1_records = AtomicU64::new(0);
        let q0_concats = AtomicU64::new(0);
        let q1_concats = AtomicU64::new(0);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("failed to create native porec-dup thread pool")?;

        log::info!(
            "Native concat PQS duplication: {} q0 shards with {} workers (no text conversion)",
            shards.len(),
            threads
        );
        pool.install(|| {
            shards
                .par_iter()
                .enumerate()
                .try_for_each(|(shard_id, input_shard)| -> anyResult<()> {
                    let mut frame = ParquetReader::new(File::open(input_shard)?)
                        .finish()
                        .with_context(|| {
                            format!("failed to read concat PQS shard {}", input_shard.display())
                        })?;

                    let read_idx = frame.column("read_idx")?.as_materialized_series().u64()?;
                    let chrom_series = frame
                        .column("chrom")?
                        .as_materialized_series()
                        .cast(&DataType::String)?;
                    let chrom = chrom_series.str()?;

                    let mut seed_array = [0u8; 32];
                    let seed_bytes = seed.to_ne_bytes();
                    seed_array[..seed_bytes.len()].copy_from_slice(&seed_bytes);
                    let shard_bytes = (shard_id as u64).to_ne_bytes();
                    seed_array[8..16].copy_from_slice(&shard_bytes);
                    let mut rng = StdRng::from_seed(seed_array);
                    let mut current_read_idx = None;
                    let mut read_decisions = HashMap::<String, String>::new();
                    let mut output_chroms = Vec::<String>::with_capacity(frame.height());

                    for row in 0..frame.height() {
                        let row_read_idx = read_idx
                            .get(row)
                            .ok_or_else(|| anyhow!("null read_idx at row {row}"))?;
                        if current_read_idx != Some(row_read_idx) {
                            current_read_idx = Some(row_read_idx);
                            read_decisions.clear();
                        }
                        let source = chrom
                            .get(row)
                            .ok_or_else(|| anyhow!("null chrom at row {row}"))?;
                        let selected = if let Some(candidates) = collapsed_contigs.get(source) {
                            read_decisions
                                .entry(source.to_string())
                                .or_insert_with(|| {
                                    candidates[rng.gen_range(0..candidates.len())].clone()
                                })
                                .clone()
                        } else {
                            source.to_string()
                        };
                        output_chroms.push(selected);
                    }

                    let chrom = Series::new("chrom".into(), output_chroms)
                        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))?;
                    frame.replace("chrom", chrom)?;

                    let mapq = frame
                        .column("mapping_quality")?
                        .as_materialized_series()
                        .u8()?;
                    let mask = mapq.gt_eq(1);
                    let mut q1_frame = frame.filter(&mask)?;

                    q0_records.fetch_add(frame.height() as u64, AtomOrdering::Relaxed);
                    q1_records.fetch_add(q1_frame.height() as u64, AtomOrdering::Relaxed);
                    q0_concats.fetch_add(count_frame_concatemers(&frame)?, AtomOrdering::Relaxed);
                    q1_concats
                        .fetch_add(count_frame_concatemers(&q1_frame)?, AtomOrdering::Relaxed);

                    let q0_path = output_path.join(format!("q0/{shard_id}.parquet"));
                    ParquetWriter::new(File::create(&q0_path)?)
                        .finish(&mut frame)
                        .with_context(|| format!("failed to write {}", q0_path.display()))?;
                    if q1_frame.height() > 0 {
                        let q1_path = output_path.join(format!("q1/{shard_id}.parquet"));
                        ParquetWriter::new(File::create(&q1_path)?)
                            .finish(&mut q1_frame)
                            .with_context(|| format!("failed to write {}", q1_path.display()))?;
                    }
                    Ok(())
                })
        })?;

        let mut sizes = target_sizes.into_iter().collect::<Vec<_>>();
        sizes.sort_by(|left, right| left.0.cmp(&right.0));
        let mut contigsizes =
            common_writer(output_path.join("_contigsizes").to_string_lossy().as_ref());
        for (target, length) in sizes {
            writeln!(contigsizes, "{target}\t{length}")?;
        }
        contigsizes.flush()?;

        let mut counts = common_writer(
            output_path
                .join("_metadata_counts")
                .to_string_lossy()
                .as_ref(),
        );
        writeln!(
            counts,
            "q0_records\t{}",
            q0_records.load(AtomOrdering::Relaxed)
        )?;
        writeln!(
            counts,
            "q1_records\t{}",
            q1_records.load(AtomOrdering::Relaxed)
        )?;
        writeln!(
            counts,
            "q0_concats\t{}",
            q0_concats.load(AtomOrdering::Relaxed)
        )?;
        writeln!(
            counts,
            "q1_concats\t{}",
            q1_concats.load(AtomOrdering::Relaxed)
        )?;
        counts.flush()?;

        std::fs::copy(input_path.join("_metadata"), output_path.join("_metadata"))?;
        std::fs::write(output_path.join("_readme"), CONCAT_PQS_README)?;
        let materialized = collapsed_contigs
            .iter()
            .map(|(source, copies)| {
                (
                    source.clone(),
                    copies
                        .iter()
                        .filter(|copy| *copy != source)
                        .cloned()
                        .collect::<std::collections::HashSet<_>>(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        write_materialized_cn_map(input_path, &materialized, output_path)?;
        log::info!(
            "Successful native concat PQS duplication into `{}`: {} q0 records across {} shards",
            output,
            q0_records.load(AtomOrdering::Relaxed),
            shards.len()
        );
        Ok(())
    }

    pub fn dup(
        &mut self,
        collapsed_list: &String,
        seed: usize,
        output: &String,
        threads: usize,
        chunksize: usize,
        chromsizes: Option<&String>,
    ) -> anyResult<()> {
        if is_concat_pqs_path(output) {
            return self.dup_pqs(collapsed_list, seed, output, threads, chunksize, chromsizes);
        }
        let reader = common_reader(collapsed_list);
        let mut collapsed_contigs: HashMap<String, Vec<String>> = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let s: Vec<&str> = record.split("\t").collect();
            if s.len() != 2 {
                log::warn!("Invalid record: {}", record);
                continue;
            }
            let contig1 = s[0].to_string();
            let contig2 = s[1].to_string();
            collapsed_contigs
                .entry(contig1.clone())
                .or_insert(vec![contig1])
                .push(contig2);
        }

        let mut seed_array = [0u8; 32];
        let seed_bytes = seed.to_ne_bytes();
        seed_array[..8].copy_from_slice(&seed_bytes);

        let (sender, receiver) = bounded::<(usize, Vec<String>)>(200);
        let (out_sender, out_receiver) = bounded::<(usize, Vec<u8>)>(200);

        let collapsed_contigs = Arc::new(collapsed_contigs);
        let num_workers = threads;
        let mut handles = vec![];

        log::info!(
            "Duplicating collapsed contigs in parallel using {} workers...",
            num_workers
        );

        for worker_id in 0..num_workers {
            let rx = receiver.clone();
            let tx = out_sender.clone();
            let collapsed = Arc::clone(&collapsed_contigs);

            let mut thread_seed = seed_array;
            let thread_offset = (worker_id as u64).to_ne_bytes();
            thread_seed[8..16].copy_from_slice(&thread_offset);

            handles.push(thread::spawn(move || {
                let mut rng = StdRng::from_seed(thread_seed);
                let mut current_read_idx = String::new();
                let mut read_decisions: HashMap<String, String> = HashMap::new();
                let mut local_buf = Vec::with_capacity(1024 * 1024);

                while let Ok((chunk_id, batch)) = rx.recv() {
                    local_buf.clear();
                    for line in batch {
                        let trimmed = line.trim_end();
                        if trimmed.is_empty() {
                            continue;
                        }

                        let mut tab_indices = [0usize; 6];
                        let mut tab_count = 0;
                        for (idx, &b) in trimmed.as_bytes().iter().enumerate() {
                            if b == b'\t' {
                                tab_indices[tab_count] = idx;
                                tab_count += 1;
                                if tab_count == 6 {
                                    break;
                                }
                            }
                        }

                        if tab_count < 6 {
                            local_buf.extend_from_slice(trimmed.as_bytes());
                            local_buf.push(b'\n');
                            continue;
                        }

                        let read_idx = &trimmed[0..tab_indices[0]];
                        let contig = &trimmed[tab_indices[4] + 1..tab_indices[5]];

                        if read_idx != current_read_idx {
                            current_read_idx.clear();
                            current_read_idx.push_str(read_idx);
                            read_decisions.clear();
                        }

                        let mut final_contig = contig;
                        if let Some(candidates) = collapsed.get(contig) {
                            final_contig =
                                read_decisions.entry(contig.to_string()).or_insert_with(|| {
                                    let idx = rng.gen_range(0..candidates.len());
                                    candidates[idx].clone()
                                });
                        }

                        local_buf.extend_from_slice(trimmed[0..tab_indices[4] + 1].as_bytes());
                        local_buf.extend_from_slice(final_contig.as_bytes());
                        local_buf.extend_from_slice(trimmed[tab_indices[5]..].as_bytes());
                        local_buf.push(b'\n');
                    }
                    tx.send((chunk_id, local_buf.clone())).unwrap();
                }
            }));
        }
        drop(out_sender);

        let out_path = output.clone();
        let write_handle = thread::spawn(move || {
            let mut wtr = common_writer(&out_path);
            let mut pending = BTreeMap::new();
            let mut next_chunk = 0;
            while let Ok((chunk_id, data)) = out_receiver.recv() {
                pending.insert(chunk_id, data);
                while let Some(data) = pending.remove(&next_chunk) {
                    wtr.write_all(&data).unwrap();
                    next_chunk += 1;
                }
            }
            wtr.flush().unwrap();
        });

        let mut rdr = self
            .parse2()
            .expect("Failed to open input file for reading");
        let batch_size = 10000;
        let mut batch = Vec::with_capacity(batch_size);
        let mut chunk_id = 0;
        let mut line_buf = String::new();
        let mut previous_read_idx = String::new();
        let mut first_iteration = true;

        while rdr.read_line(&mut line_buf).unwrap_or(0) > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let read_idx = trimmed.split('\t').next().unwrap_or("");
            if !first_iteration && read_idx != previous_read_idx {
                if batch.len() >= batch_size {
                    sender.send((chunk_id, std::mem::take(&mut batch))).unwrap();
                    batch = Vec::with_capacity(batch_size);
                    chunk_id += 1;
                }
            }

            first_iteration = false;
            previous_read_idx.clear();
            previous_read_idx.push_str(read_idx);
            batch.push(std::mem::take(&mut line_buf));
        }

        if !batch.is_empty() {
            sender.send((chunk_id, batch)).unwrap();
        }
        drop(sender);

        for h in handles {
            h.join().unwrap();
        }
        write_handle.join().unwrap();

        log::info!("Successfully executed parallelized dup to `{}`", output);
        Ok(())
    }

    pub fn split(&mut self, output: &String) {
        let reader = open_porec_reader(&self.file).expect("Failed to open Pore-C input");
        let mut wtr = common_writer(output);

        let (sender, receiver) = bounded::<Vec<(u32, String)>>(100);

        log::info!("Successful output split porec table into `{}`", output);
    }

    pub fn to_pe_pair_reads_single(
        &mut self,
        fasta_path: &str,
        output_prefix: &str,
        read_len: Option<usize>,
    ) -> anyResult<()> {
        log::info!("Loading reference genome from `{}`...", fasta_path);
        let mut genome: HashMap<String, Vec<u8>> = HashMap::new();
        let mut reader = common_reader(&fasta_path.to_string());
        let mut line = String::new();
        let mut current_ctg = String::new();
        let mut current_seq = Vec::new();

        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim();
            if trimmed.starts_with('>') {
                if !current_ctg.is_empty() {
                    genome.insert(current_ctg.clone(), current_seq.clone());
                    current_seq.clear();
                }
                current_ctg = trimmed[1..].split_whitespace().next().unwrap().to_string();
            } else {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
            line.clear();
        }
        if !current_ctg.is_empty() {
            genome.insert(current_ctg, current_seq);
        }

        let out_r1 = format!("{}_R1.fa.gz", output_prefix);
        let out_r2 = format!("{}_R2.fa.gz", output_prefix);
        let mut wtr_r1 = common_writer(&out_r1);
        let mut wtr_r2 = common_writer(&out_r2);

        let mut rdr = self
            .parse2()
            .expect("Failed to open porec table for reading");
        let mut concatemer = Concatemer::new();
        let mut old_read_idx: u64 = u64::MAX;
        let mut line_buf = String::new();
        let mut pair_id = 0;

        let revcomp = |seq: &[u8]| -> Vec<u8> {
            seq.iter()
                .rev()
                .map(|&c| match c {
                    b'A' | b'a' => b'T',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    b'T' | b't' => b'A',
                    b'N' | b'n' => b'N',
                    _ => b'N',
                })
                .collect()
        };

        log::info!("Generating paired-end reads...");

        while rdr.read_line(&mut line_buf)? > 0 {
            let trimmed = line_buf.trim_end();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                line_buf.clear();
                continue;
            }

            let mut parts = trimmed.split('\t');
            let read_idx = parts.next().unwrap_or("0").parse::<u64>().unwrap_or(0);

            if read_idx != old_read_idx && old_read_idx != u64::MAX {
                concatemer.sort();
                for pair in concatemer.decompose() {
                    let r1 = &pair[0];
                    let r2 = &pair[1];
                    pair_id += 1;

                    let process_read = |r: &PoreCRecord| -> Option<Vec<u8>> {
                        let seq = genome.get(&r.target)?;
                        // 1-based to 0-based
                        let start = r.target_start.saturating_sub(1) as usize;
                        let end = (r.target_end as usize).min(seq.len());
                        if start >= end {
                            return None;
                        }

                        let mut fragment = seq[start..end].to_vec();
                        if r.query_strand == '-' {
                            fragment = revcomp(&fragment);
                        }

                        if let Some(len) = read_len {
                            fragment.truncate(len);
                        }
                        Some(fragment)
                    };

                    if let (Some(seq1), Some(seq2)) = (process_read(&r1), process_read(&r2)) {
                        writeln!(
                            wtr_r1,
                            ">read_{} 1\n{}",
                            pair_id,
                            String::from_utf8_lossy(&seq1)
                        )
                        .unwrap();
                        writeln!(
                            wtr_r2,
                            ">read_{} 2\n{}",
                            pair_id,
                            String::from_utf8_lossy(&seq2)
                        )
                        .unwrap();
                    }
                }
                concatemer.clear();
            }

            let q_len = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_start = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_end = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let q_strand = parts.next().and_then(|s| s.chars().next()).unwrap_or('+');
            let target = parts.next().unwrap_or("").to_string();
            let t_start = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let t_end = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            let mapq = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0);

            concatemer.push(PoreCRecord {
                read_idx,
                query_length: q_len,
                query_start: q_start,
                query_end: q_end,
                query_strand: q_strand,
                target,
                target_start: t_start,
                target_end: t_end,
                mapq,
                identity: 0.0,
                filter_reason: "".to_string(),
            });

            old_read_idx = read_idx;
            line_buf.clear();
        }

        concatemer.sort();
        for pair in concatemer.decompose() {
            let r1 = &pair[0];
            let r2 = &pair[1];
            pair_id += 1;
            let process_read = |r: &PoreCRecord| -> Option<Vec<u8>> {
                let seq = genome.get(&r.target)?;
                let start = r.target_start.saturating_sub(1) as usize;
                let end = (r.target_end as usize).min(seq.len());
                if start >= end {
                    return None;
                }

                let mut fragment = seq[start..end].to_vec();
                if r.query_strand == '-' {
                    fragment = revcomp(&fragment);
                }

                if let Some(len) = read_len {
                    fragment.truncate(len);
                }
                Some(fragment)
            };

            if let (Some(seq1), Some(seq2)) = (process_read(&r1), process_read(&r2)) {
                writeln!(
                    wtr_r1,
                    ">read_{} 1\n{}",
                    pair_id,
                    String::from_utf8_lossy(&seq1)
                )
                .unwrap();
                writeln!(
                    wtr_r2,
                    ">read_{} 2\n{}",
                    pair_id,
                    String::from_utf8_lossy(&seq2)
                )
                .unwrap();
            }
        }

        log::info!(
            "Successfully generated PE reads to {}_R[12].fa",
            output_prefix
        );
        Ok(())
    }

    pub fn to_pe_pair_reads(
        &mut self,
        fasta_path: &str,
        output_prefix: &str,
        read_len: Option<usize>,
    ) -> anyResult<()> {
        log::info!("Loading reference genome from `{}`...", fasta_path);
        let mut genome_map: HashMap<String, Vec<u8>> = HashMap::new();
        let mut reader = common_reader(&fasta_path.to_string());
        let mut line = String::new();
        let mut current_ctg = String::new();
        let mut current_seq = Vec::new();

        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim();
            if trimmed.starts_with('>') {
                if !current_ctg.is_empty() {
                    genome_map.insert(current_ctg.clone(), current_seq.clone());
                    current_seq.clear();
                }
                current_ctg = trimmed[1..].split_whitespace().next().unwrap().to_string();
            } else {
                current_seq.extend_from_slice(trimmed.as_bytes());
            }
            line.clear();
        }
        if !current_ctg.is_empty() {
            genome_map.insert(current_ctg, current_seq);
        }

        let genome = Arc::new(genome_map);

        let out_r1 = format!("{}_R1.fa.gz", output_prefix);
        let out_r2 = format!("{}_R2.fa.gz", output_prefix);

        log::info!("Generating paired-end reads in parallel...");

        let (tx_work, rx_work) = bounded::<Vec<Concatemer>>(200);
        let (tx_write, rx_write) = bounded::<(Vec<u8>, Vec<u8>)>(200);

        let pair_id_counter = Arc::new(AtomicU64::new(1));
        let num_workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8);

        let mut handles = Vec::with_capacity(num_workers);

        for _ in 0..num_workers {
            let rx = rx_work.clone();
            let tx = tx_write.clone();
            let r#gen = Arc::clone(&genome);
            let counter = Arc::clone(&pair_id_counter);

            handles.push(thread::spawn(move || {
                let mut buf_r1 = Vec::with_capacity(2 * 1024 * 1024);
                let mut buf_r2 = Vec::with_capacity(2 * 1024 * 1024);

                let revcomp = |seq: &[u8]| -> Vec<u8> {
                    seq.iter()
                        .rev()
                        .map(|&c| match c {
                            b'A' | b'a' => b'T',
                            b'C' | b'c' => b'G',
                            b'G' | b'g' => b'C',
                            b'T' | b't' => b'A',
                            b'N' | b'n' => b'N',
                            _ => b'N',
                        })
                        .collect()
                };

                while let Ok(batch) = rx.recv() {
                    for mut concatemer in batch {
                        concatemer.sort();
                        for pair in concatemer.decompose() {
                            let r1 = &pair[0];
                            let r2 = &pair[1];

                            let process_read = |r: &PoreCRecord| -> Option<Vec<u8>> {
                                let seq = r#gen.get(&r.target)?;
                                let start = r.target_start.saturating_sub(1) as usize;
                                let end = (r.target_end as usize).min(seq.len());
                                if start >= end {
                                    return None;
                                }

                                let mut fragment = seq[start..end].to_vec();
                                if r.query_strand == '-' {
                                    fragment = revcomp(&fragment);
                                }

                                if let Some(len) = read_len {
                                    fragment.truncate(len);
                                }
                                Some(fragment)
                            };

                            if let (Some(seq1), Some(seq2)) = (process_read(&r1), process_read(&r2))
                            {
                                let pid = counter.fetch_add(1, AtomOrdering::Relaxed);
                                use std::io::Write;
                                let _ = writeln!(&mut buf_r1, ">read_{}/1\n{}", pid, unsafe {
                                    std::str::from_utf8_unchecked(&seq1)
                                });
                                let _ = writeln!(&mut buf_r2, ">read_{}/2\n{}", pid, unsafe {
                                    std::str::from_utf8_unchecked(&seq2)
                                });
                            }
                        }
                    }
                    if buf_r1.len() >= 1024 * 1024 {
                        tx.send((buf_r1.clone(), buf_r2.clone())).unwrap();
                        buf_r1.clear();
                        buf_r2.clear();
                    }
                }

                if !buf_r1.is_empty() {
                    tx.send((buf_r1, buf_r2)).unwrap();
                }
            }));
        }
        drop(tx_write); // Drop the original sender so writer knows when to stop

        // Writer thread
        let write_handle = thread::spawn(move || {
            let mut wtr_r1 = common_writer(&out_r1);
            let mut wtr_r2 = common_writer(&out_r2);
            while let Ok((chunk_r1, chunk_r2)) = rx_write.recv() {
                wtr_r1.write_all(&chunk_r1).unwrap();
                wtr_r2.write_all(&chunk_r2).unwrap();
            }
            wtr_r1.flush().unwrap();
            wtr_r2.flush().unwrap();
        });

        let batch_size = 2000;
        let mut batch = Vec::with_capacity(batch_size);
        if is_concat_pqs(Path::new(&self.file)) {
            use polars::prelude::*;
            let shards = concat_pqs_files(Path::new(&self.file))?;
            log::info!(
                "Reading {} concat PQS shards directly for porec2reads",
                shards.len()
            );
            for shard in shards {
                let frame = ParquetReader::new(File::open(&shard)?)
                    .finish()
                    .with_context(|| {
                        format!("failed to read concat PQS shard {}", shard.display())
                    })?;
                let records = concat_frame_to_records(&frame)?;
                let mut concatemer = Concatemer::new();
                let mut previous = None;
                for record in records {
                    if previous.is_some() && previous != Some(record.read_idx) {
                        batch.push(std::mem::take(&mut concatemer));
                        if batch.len() >= batch_size {
                            tx_work.send(std::mem::take(&mut batch)).unwrap();
                        }
                    }
                    previous = Some(record.read_idx);
                    concatemer.push(record);
                }
                if concatemer.count() > 0 {
                    batch.push(concatemer);
                }
            }
        } else {
            let mut rdr = self
                .parse2()
                .expect("Failed to open porec table for reading");
            let mut concatemer = Concatemer::new();
            let mut old_read_idx: u64 = u64::MAX;
            let mut line_buf = String::new();
            while rdr.read_line(&mut line_buf)? > 0 {
                let trimmed = line_buf.trim_end();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    line_buf.clear();
                    continue;
                }
                let mut parts = trimmed.split('\t');
                let read_idx = parts.next().unwrap_or("0").parse::<u64>().unwrap_or(0);
                if read_idx != old_read_idx && old_read_idx != u64::MAX {
                    batch.push(std::mem::take(&mut concatemer));
                    if batch.len() >= batch_size {
                        tx_work.send(std::mem::take(&mut batch)).unwrap();
                    }
                }
                concatemer.push(PoreCRecord {
                    read_idx,
                    query_length: parts.next().and_then(|s| s.parse().ok()).unwrap_or(0),
                    query_start: parts.next().and_then(|s| s.parse().ok()).unwrap_or(0),
                    query_end: parts.next().and_then(|s| s.parse().ok()).unwrap_or(0),
                    query_strand: parts.next().and_then(|s| s.chars().next()).unwrap_or('+'),
                    target: parts.next().unwrap_or("").to_string(),
                    target_start: parts.next().and_then(|s| s.parse().ok()).unwrap_or(0),
                    target_end: parts.next().and_then(|s| s.parse().ok()).unwrap_or(0),
                    mapq: parts.next().and_then(|s| s.parse().ok()).unwrap_or(0),
                    identity: parts.next().and_then(|s| s.parse().ok()).unwrap_or(0.0),
                    filter_reason: parts.next().unwrap_or("").to_string(),
                });
                old_read_idx = read_idx;
                line_buf.clear();
            }
            if concatemer.count() > 0 {
                batch.push(concatemer);
            }
        }
        if !batch.is_empty() {
            tx_work.send(batch).unwrap();
        }
        drop(tx_work);

        // Wait for workers and writer
        for h in handles {
            h.join().unwrap();
        }
        write_handle.join().unwrap();

        log::info!(
            "Successfully generated PE reads to {}_R[12].fa.gz",
            output_prefix
        );
        Ok(())
    }

    fn downsample_concat_pqs_native(
        &self,
        output: &String,
        seed: u64,
        min_quality: u8,
        min_order: usize,
        max_order: usize,
        target_reads: Option<u64>,
        target_pairs: Option<u64>,
        frac: Option<f64>,
        frac_by_pairs: bool,
        threads: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;

        let input = Path::new(&self.file);
        let shards = concat_pqs_files(input)?;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("failed to create native concat PQS downsample thread pool")?;
        let per_shard = pool.install(|| {
            shards
                .par_iter()
                .enumerate()
                .map(
                    |(shard_id, shard)| -> anyResult<Vec<((usize, u64), u32, u64)>> {
                        let frame = ParquetReader::new(File::open(shard)?)
                            .finish()
                            .with_context(|| {
                                format!("failed to read concat PQS shard {}", shard.display())
                            })?;
                        let read_idx = frame.column("read_idx")?.as_materialized_series().u64()?;
                        let mapq = frame
                            .column("mapping_quality")?
                            .as_materialized_series()
                            .u8()?;
                        let mut result = Vec::new();
                        let mut previous = None;
                        let mut order = 0usize;
                        for row in 0..frame.height() {
                            let read = read_idx
                                .get(row)
                                .ok_or_else(|| anyhow!("null read_idx at row {row}"))?;
                            if previous.is_some() && previous != Some(read) {
                                if order >= min_order && order < max_order {
                                    let pairs = (order as u64)
                                        .saturating_mul(order.saturating_sub(1) as u64)
                                        / 2;
                                    result.push((
                                        (shard_id, previous.unwrap()),
                                        order as u32,
                                        pairs,
                                    ));
                                }
                                order = 0;
                            }
                            if mapq.get(row).unwrap_or(0) >= min_quality {
                                order += 1;
                            }
                            previous = Some(read);
                        }
                        if let Some(read) = previous {
                            if order >= min_order && order < max_order {
                                let pairs = (order as u64)
                                    .saturating_mul(order.saturating_sub(1) as u64)
                                    / 2;
                                result.push(((shard_id, read), order as u32, pairs));
                            }
                        }
                        Ok(result)
                    },
                )
                .collect::<anyResult<Vec<_>>>()
        })?;
        let mut per_read = per_shard.into_iter().flatten().collect::<Vec<_>>();
        let total_reads = per_read.len() as u64;
        let total_pairs = per_read.iter().map(|entry| entry.2).sum::<u64>();
        let (want_reads, want_pairs) = if let Some(reads) = target_reads {
            (Some(reads.min(total_reads)), None)
        } else if let Some(pairs) = target_pairs {
            (None, Some(pairs.min(total_pairs)))
        } else {
            let fraction = frac.expect("validated fraction");
            if frac_by_pairs {
                (None, Some(((total_pairs as f64) * fraction).round() as u64))
            } else {
                (Some(((total_reads as f64) * fraction).round() as u64), None)
            }
        };
        per_read.shuffle(&mut StdRng::seed_from_u64(seed));
        let mut chosen = std::collections::HashSet::<(usize, u64)>::new();
        if let Some(reads) = want_reads {
            chosen.extend(per_read.iter().take(reads as usize).map(|entry| entry.0));
        } else if let Some(target) = want_pairs {
            let mut accumulated = 0u64;
            for entry in &per_read {
                if accumulated >= target {
                    break;
                }
                chosen.insert(entry.0);
                accumulated = accumulated.saturating_add(entry.2);
            }
        }
        let selected_count = chosen.len();
        let chosen = Arc::new(chosen);
        let mut sizes = BTreeMap::new();
        read_concat_pqs_contigsizes(input, &mut sizes)?;
        transform_concat_pqs_native(
            input,
            Path::new(output.trim_end_matches('/')),
            sizes.into_iter().collect(),
            threads,
            true,
            move |shard_id, frame| {
                let read_idx = frame.column("read_idx")?.as_materialized_series().u64()?;
                let keep = read_idx
                    .into_iter()
                    .map(|read| read.is_some_and(|read| chosen.contains(&(shard_id, read))))
                    .collect::<Vec<_>>();
                Ok(frame.filter(&BooleanChunked::from_slice("keep".into(), &keep))?)
            },
        )?;
        log::info!(
            "Native concat PQS downsample selected {} of {} eligible reads into `{}`",
            selected_count,
            total_reads,
            output
        );
        Ok(())
    }

    pub fn downsample(
        &mut self,
        output: &String,
        seed: u64,
        min_quality: u8,
        min_order: usize,
        max_order: usize,
        target_reads: Option<u64>,
        target_pairs: Option<u64>,
        frac: Option<f64>,
        frac_by_pairs: bool,
        threads: usize,
    ) -> anyResult<()> {
        if [
            target_reads.is_some(),
            target_pairs.is_some(),
            frac.is_some(),
        ]
        .into_iter()
        .filter(|x| *x)
        .count()
            != 1
        {
            anyhow::bail!("Specify exactly one of --reads, --pairs, or --frac");
        }
        if let Some(f) = frac {
            if !(0.0 < f && f <= 1.0) {
                anyhow::bail!("--frac must be in (0, 1]");
            }
        }
        if is_concat_pqs_path(output) {
            if !is_concat_pqs(Path::new(&self.file)) {
                bail!("concat PQS downsample output requires concat.pqs input");
            }
            return self.downsample_concat_pqs_native(
                output,
                seed,
                min_quality,
                min_order,
                max_order,
                target_reads,
                target_pairs,
                frac,
                frac_by_pairs,
                threads,
            );
        }

        log::info!(
            "Downsample pass1 scanning... (min_mapq={}, order in [{}, {}))",
            min_quality,
            min_order,
            max_order
        );

        let mut rdr = self.parse2().expect("Failed to open input file");
        let mut line = String::new();

        let mut prev_read: Option<u64> = None;
        let mut cur_order: usize = 0;

        let mut per_read: Vec<(u64, u32, u64)> = Vec::new();

        #[inline]
        fn pairs_of(order: usize) -> u64 {
            let n = order as u64;
            n.saturating_mul(n.saturating_sub(1)) / 2
        }

        while rdr.read_line(&mut line)? > 0 {
            let t = line.trim_end();
            if t.is_empty() || t.starts_with('#') {
                line.clear();
                continue;
            }

            let mut it = t.split('\t');
            let read_idx = it.next().and_then(|s| s.parse::<u64>().ok()).unwrap_or(0);

            let mut mapq: u8 = 0;
            for _ in 0..7 {
                if it.next().is_none() {
                    break;
                }
            }
            if let Some(mq_s) = it.next() {
                mapq = mq_s.parse::<u8>().unwrap_or(0);
            }

            let is_eligible_piece = mapq >= min_quality;

            match prev_read {
                None => {
                    prev_read = Some(read_idx);
                    cur_order = if is_eligible_piece { 1 } else { 0 };
                }
                Some(pr) if pr == read_idx => {
                    if is_eligible_piece {
                        cur_order += 1;
                    }
                }
                Some(pr) => {
                    // finalize previous read
                    if cur_order >= min_order && cur_order < max_order {
                        let p = pairs_of(cur_order);
                        per_read.push((pr, cur_order as u32, p));
                    }
                    // start new
                    prev_read = Some(read_idx);
                    cur_order = if is_eligible_piece { 1 } else { 0 };
                }
            }

            line.clear();
        }
        // finalize last
        if let Some(pr) = prev_read {
            if cur_order >= min_order && cur_order < max_order {
                let p = pairs_of(cur_order);
                per_read.push((pr, cur_order as u32, p));
            }
        }

        let total_reads = per_read.len() as u64;
        let total_pairs: u64 = per_read.iter().map(|x| x.2).sum();

        log::info!(
            "Eligible reads: {} ; eligible virtual pairs: {}",
            total_reads,
            total_pairs
        );

        let (want_reads, want_pairs): (Option<u64>, Option<u64>) = if let Some(n) = target_reads {
            (Some(n.min(total_reads)), None)
        } else if let Some(p) = target_pairs {
            (None, Some(p.min(total_pairs)))
        } else {
            // frac
            let f = frac.unwrap();
            if frac_by_pairs {
                (None, Some(((total_pairs as f64) * f).round() as u64))
            } else {
                (Some(((total_reads as f64) * f).round() as u64), None)
            }
        };

        let mut rng = StdRng::seed_from_u64(seed);
        per_read.shuffle(&mut rng);

        let mut chosen: std::collections::HashSet<u64> = std::collections::HashSet::new();

        if let Some(k) = want_reads {
            for (rid, _ord, _p) in per_read.iter().take(k as usize) {
                chosen.insert(*rid);
            }
            log::info!("Selected reads: {}", chosen.len());
        } else if let Some(tp) = want_pairs {
            let mut acc: u64 = 0;
            for (rid, _ord, p) in per_read.iter() {
                if acc >= tp {
                    break;
                }
                chosen.insert(*rid);
                acc = acc.saturating_add(*p);
            }
            log::info!(
                "Selected reads: {} (approx pairs: {} / target {})",
                chosen.len(),
                acc,
                tp
            );
        }

        log::info!("Downsample pass2 writing output to `{}`...", output);
        let mut rdr2 = self.parse2().expect("Failed to open input file (pass2)");
        let mut wtr = common_writer(output);

        let mut line2 = String::new();
        while rdr2.read_line(&mut line2)? > 0 {
            let t = line2.trim_end();
            if t.is_empty() {
                line2.clear();
                continue;
            }
            if t.starts_with('#') {
                wtr.write_all(line2.as_bytes())?;
                line2.clear();
                continue;
            }

            let read_idx = t
                .split('\t')
                .next()
                .and_then(|s| s.parse::<u64>().ok())
                .unwrap_or(u64::MAX);

            if chosen.contains(&read_idx) {
                wtr.write_all(line2.as_bytes())?;
            }

            line2.clear();
        }

        log::info!(
            "Downsample finished. Output reads: {} (seed={})",
            chosen.len(),
            seed
        );
        Ok(())
    }

    pub fn to_depth(
        &mut self,
        chromsize: &String,
        window_size: usize,
        step_size: usize,
        min_mapq: u8,
        output: &String,
    ) -> Result<(), Box<dyn Error>> {
        let mut chrom_names = Vec::new();
        let mut chrom_sizes = Vec::new();
        let mut chrom_map = HashMap::new();

        log::info!("Loading chromsizes from `{}`...", chromsize);
        {
            let input = common_reader(chromsize);
            let buf = std::io::BufReader::new(input);
            for line in buf.lines() {
                let line = line?;
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }

                let mut spl = line.split_whitespace();
                let chr = spl.next().unwrap_or("").to_string();
                let size_str = spl.next().unwrap_or("0");
                let size = size_str.parse::<usize>().unwrap_or(0);
                chrom_map.insert(chr.clone(), chrom_names.len());
                chrom_names.push(chr);
                chrom_sizes.push(size);
            }
        }

        let num_chroms = chrom_names.len();
        let chrom_map = Arc::new(chrom_map);
        let mut global_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];
        if is_concat_pqs(Path::new(&self.file)) {
            use polars::prelude::*;

            polars::enable_string_cache();
            let shards = concat_pqs_files(Path::new(&self.file))?;
            log::info!(
                "Calculating coverage directly from {} concat PQS shards...",
                shards.len()
            );
            let shared_events = Arc::new(
                (0..num_chroms)
                    .map(|_| Mutex::new(Vec::<(u32, i32)>::new()))
                    .collect::<Vec<_>>(),
            );
            let projected_columns = ["chrom", "start", "end", "mapping_quality"]
                .into_iter()
                .map(str::to_string)
                .collect::<Vec<_>>();
            shards.par_iter().try_for_each(|shard| -> anyResult<()> {
                let frame = ParquetReader::new(File::open(shard)?)
                    .with_columns(Some(projected_columns.clone()))
                    .finish()
                    .with_context(|| {
                        format!("failed to read concat PQS shard {}", shard.display())
                    })?;
                let chrom = frame.column("chrom")?.categorical()?;
                let chrom_codes = chrom.physical();
                let chrom_rev_map = chrom.get_rev_map();
                let max_chrom_code = chrom_codes.max().unwrap_or(0) as usize;
                let mut chrom_indices = vec![None; max_chrom_code + 1];
                for code in 0..=max_chrom_code {
                    let Some(name) = chrom_rev_map.get_optional(code as u32) else {
                        continue;
                    };
                    chrom_indices[code] = chrom_map.get(name).copied();
                }
                let start_series = frame
                    .column("start")?
                    .as_materialized_series()
                    .cast(&DataType::UInt64)?;
                let start = start_series.u64()?;
                let end_series = frame
                    .column("end")?
                    .as_materialized_series()
                    .cast(&DataType::UInt64)?;
                let end = end_series.u64()?;
                let mapq = frame
                    .column("mapping_quality")?
                    .as_materialized_series()
                    .u8()?;
                let mut local_events = vec![Vec::<(u32, i32)>::new(); num_chroms];
                for row in 0..frame.height() {
                    if mapq.get(row).unwrap_or(0) < min_mapq {
                        continue;
                    }
                    let Some(chrom_code) = chrom_codes.get(row) else {
                        continue;
                    };
                    let Some(idx) = chrom_indices.get(chrom_code as usize).copied().flatten()
                    else {
                        continue;
                    };
                    let source_start = start.get(row).unwrap_or(0);
                    let source_end = end.get(row).unwrap_or(0);
                    let start_u32 = u32::try_from(source_start)
                        .with_context(|| format!("coverage start {source_start} exceeds UInt32"))?;
                    let end_u32 = u32::try_from(source_end)
                        .with_context(|| format!("coverage end {source_end} exceeds UInt32"))?;
                    local_events[idx].push((start_u32, 1));
                    local_events[idx].push((end_u32, -1));
                }
                for (idx, mut events) in local_events.into_iter().enumerate() {
                    if !events.is_empty() {
                        shared_events[idx].lock().unwrap().append(&mut events);
                    }
                }
                Ok(())
            })?;
            let shared_events = Arc::try_unwrap(shared_events)
                .map_err(|_| anyhow!("porec2depth event collectors are still shared"))?;
            global_events = shared_events
                .into_iter()
                .map(|events| {
                    events
                        .into_inner()
                        .map_err(|_| anyhow!("porec2depth event collector was poisoned"))
                })
                .collect::<anyResult<Vec<_>>>()?;
        } else {
            let (sender, receiver) = bounded::<Vec<String>>(200);
            let num_workers = std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(8);
            let mut handles = Vec::new();
            for _ in 0..num_workers {
                let receiver = receiver.clone();
                let chrom_map = chrom_map.clone();
                handles.push(thread::spawn(move || {
                    let mut local_events: Vec<Vec<(u32, i32)>> = vec![Vec::new(); num_chroms];
                    while let Ok(lines) = receiver.recv() {
                        for line in lines {
                            if line.is_empty() || line.starts_with('#') {
                                continue;
                            }
                            let mut parts = line.split('\t');
                            let target = parts.nth(5).unwrap_or("");
                            let t_start_str = parts.next().unwrap_or("0");
                            let t_end_str = parts.next().unwrap_or("0");
                            let mapq_str = parts.next().unwrap_or("0");
                            let mapq = mapq_str.parse::<u8>().unwrap_or(0);
                            if mapq < min_mapq {
                                continue;
                            }
                            if let Some(&idx) = chrom_map.get(target) {
                                let start = t_start_str.parse::<u32>().unwrap_or(0);
                                let end = t_end_str.parse::<u32>().unwrap_or(0);
                                local_events[idx].push((start, 1));
                                local_events[idx].push((end, -1));
                            }
                        }
                    }
                    local_events
                }));
            }

            log::info!("Parsing Pore-C table to calculate coverage depth...");
            let mut rdr = self
                .parse2()
                .expect("Failed to open porec table for reading");
            let mut line_buf = String::new();
            let batch_size = 10_000;
            let mut batch = Vec::with_capacity(batch_size);
            while rdr.read_line(&mut line_buf)? > 0 {
                if !line_buf.trim().is_empty() {
                    batch.push(line_buf.clone());
                }
                if batch.len() >= batch_size {
                    sender.send(std::mem::take(&mut batch)).unwrap();
                }
                line_buf.clear();
            }
            if !batch.is_empty() {
                sender.send(batch).unwrap();
            }
            drop(sender);
            for handle in handles {
                let local_events = handle.join().unwrap();
                for (i, events) in local_events.into_iter().enumerate() {
                    global_events[i].extend(events);
                }
            }
        }

        log::info!("Aggregating coverage across parallel windows...");

        let output_clone = output.clone();
        // Using Rayon Scope to handle parallel generation and sequential output file writing
        rayon::scope(|s| {
            let (tx, rx) = bounded::<(usize, String)>(1024);

            s.spawn(move |_| {
                use std::io::Write;
                let writer = common_writer(&output_clone);
                let mut buf_writer = std::io::BufWriter::with_capacity(1024 * 1024, writer);
                let mut pending = std::collections::BTreeMap::new();
                let mut next_idx = 0;

                while let Ok((idx, data)) = rx.recv() {
                    pending.insert(idx, data);
                    while let Some(d) = pending.remove(&next_idx) {
                        buf_writer.write_all(d.as_bytes()).unwrap();
                        next_idx += 1;
                    }
                }
                buf_writer.flush().unwrap();
            });

            global_events
                .into_par_iter()
                .enumerate()
                .for_each_with(tx, |tx, (i, mut events)| {
                    if events.is_empty() {
                        tx.send((i, String::new())).unwrap();
                        return;
                    }

                    let chrom_len = chrom_sizes[i];
                    let chrom_name = &chrom_names[i];

                    events.par_sort_unstable_by_key(|e| e.0);

                    if window_size == 0 || step_size == 0 {
                        tx.send((i, String::new())).unwrap();
                        return;
                    }
                    let window_starts = (0..chrom_len).step_by(step_size).collect::<Vec<_>>();
                    let window_ends = window_starts
                        .iter()
                        .map(|start| start.saturating_add(window_size).min(chrom_len))
                        .collect::<Vec<_>>();
                    let mut query_positions = Vec::with_capacity(window_starts.len() * 2);
                    query_positions.extend_from_slice(&window_starts);
                    query_positions.extend_from_slice(&window_ends);
                    query_positions.sort_unstable();
                    query_positions.dedup();
                    let integrals = coverage_integrals_at(&events, &query_positions, chrom_len);

                    let mut chunk_out = String::with_capacity(4096);
                    use std::fmt::Write;

                    for (start, end) in window_starts.into_iter().zip(window_ends) {
                        let len = end - start;
                        if len == 0 {
                            continue;
                        }
                        let start_index = query_positions.binary_search(&start).unwrap();
                        let end_index = query_positions.binary_search(&end).unwrap();
                        let sum = integrals[end_index] - integrals[start_index];
                        let mean = sum as f64 / len as f64;

                        writeln!(
                            &mut chunk_out,
                            "{}\t{}\t{}\t{:.3}",
                            chrom_name, start, end, mean
                        )
                        .unwrap();
                    }

                    tx.send((i, chunk_out)).unwrap();
                });
        });

        log::info!(
            "Successful output coverage (window {}bp, step {}bp) to `{}`",
            window_size,
            step_size,
            output
        );
        Ok(())
    }
}

fn is_concat_pqs_path(path: &str) -> bool {
    let path = path.trim_end_matches('/');
    path.ends_with(".concat.pqs") || path.ends_with(".porec.pqs")
}

fn read_contigsizes_file(path: &Path, merged: &mut BTreeMap<String, u64>) -> anyResult<()> {
    let reader = common_reader(path.to_string_lossy().as_ref());
    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        let chrom = fields.next().ok_or_else(|| {
            anyhow!(
                "missing contig name at {}:{}",
                path.display(),
                line_number + 1
            )
        })?;
        let length = fields
            .next()
            .ok_or_else(|| {
                anyhow!(
                    "missing contig length at {}:{}",
                    path.display(),
                    line_number + 1
                )
            })?
            .parse::<u64>()
            .with_context(|| {
                format!(
                    "invalid contig length at {}:{}",
                    path.display(),
                    line_number + 1
                )
            })?;
        if let Some(previous) = merged.insert(chrom.to_string(), length) {
            if previous != length {
                bail!(
                    "conflicting lengths for contig {chrom}: {previous} and {length} in {}",
                    path.display()
                );
            }
        }
    }
    Ok(())
}

fn read_concat_pqs_contigsizes(input: &Path, merged: &mut BTreeMap<String, u64>) -> anyResult<()> {
    read_contigsizes_file(&input.join("_contigsizes"), merged)
}

#[derive(Clone, Copy, Default)]
struct ConcatPqsCounts {
    q0_records: u64,
    q1_records: u64,
    q0_concats: u64,
    q1_concats: u64,
}

impl ConcatPqsCounts {
    fn checked_add(self, other: Self, input: &Path) -> anyResult<Self> {
        let add = |left: u64, right: u64, field: &str| {
            left.checked_add(right).ok_or_else(|| {
                anyhow!(
                    "{field} overflow while merging metadata counts from {}",
                    input.display()
                )
            })
        };
        Ok(Self {
            q0_records: add(self.q0_records, other.q0_records, "q0_records")?,
            q1_records: add(self.q1_records, other.q1_records, "q1_records")?,
            q0_concats: add(self.q0_concats, other.q0_concats, "q0_concats")?,
            q1_concats: add(self.q1_concats, other.q1_concats, "q1_concats")?,
        })
    }
}

fn read_concat_pqs_counts(input: &Path) -> anyResult<ConcatPqsCounts> {
    let path = input.join("_metadata_counts");
    let contents = std::fs::read_to_string(&path)
        .with_context(|| format!("zero-read PQS merge requires {}", path.display()))?;
    let mut values = BTreeMap::<&str, u64>::new();
    for (line_number, line) in contents.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        let key = fields.next().ok_or_else(|| {
            anyhow!(
                "missing count name at {}:{}",
                path.display(),
                line_number + 1
            )
        })?;
        let value = fields
            .next()
            .ok_or_else(|| {
                anyhow!(
                    "missing count value at {}:{}",
                    path.display(),
                    line_number + 1
                )
            })?
            .parse::<u64>()
            .with_context(|| format!("invalid count at {}:{}", path.display(), line_number + 1))?;
        values.insert(key, value);
    }
    let get = |key| {
        values
            .get(key)
            .copied()
            .ok_or_else(|| anyhow!("zero-read PQS merge requires `{key}` in {}", path.display()))
    };
    Ok(ConcatPqsCounts {
        q0_records: get("q0_records")?,
        q1_records: get("q1_records")?,
        q0_concats: get("q0_concats")?,
        q1_concats: get("q1_concats")?,
    })
}

fn concat_pqs_metadata_layout(input: &Path) -> anyResult<(&'static str, usize)> {
    let path = input.join("_metadata");
    let metadata = std::fs::read_to_string(&path)?;
    let position_type = if metadata.contains("'start': UInt32") {
        "UInt32"
    } else if metadata.contains("'start': UInt64") {
        "UInt64"
    } else {
        bail!(
            "cannot determine start/end position type from {}; zero-read merge cannot normalize incompatible Parquet schemas",
            path.display()
        );
    };
    let chunksize = metadata
        .lines()
        .find_map(|line| {
            line.trim()
                .strip_prefix("'chunksize':")
                .and_then(|value| value.trim().trim_end_matches(',').parse::<usize>().ok())
        })
        .unwrap_or(1_000_000);
    Ok((position_type, chunksize))
}

fn link_concat_pqs_shard(source: &Path, destination: &Path) -> anyResult<()> {
    match std::fs::hard_link(source, destination) {
        Ok(()) => Ok(()),
        Err(hard_link_error) => {
            #[cfg(unix)]
            {
                let absolute_source = std::fs::canonicalize(source)?;
                std::os::unix::fs::symlink(&absolute_source, destination).with_context(|| {
                    format!(
                        "failed to hard-link {} to {} ({hard_link_error}); symlink fallback also failed",
                        source.display(),
                        destination.display()
                    )
                })?;
                Ok(())
            }
            #[cfg(not(unix))]
            {
                bail!(
                    "failed to hard-link {} to {}: {hard_link_error}; zero-read merge has no fallback on this platform",
                    source.display(),
                    destination.display()
                )
            }
        }
    }
}

fn write_concat_pqs_sidecars(
    metadata_source: &Path,
    output: &Path,
    target_sizes: HashMap<String, u64>,
    counts_source: Option<std::path::PathBuf>,
    preserve_cn_info: bool,
) -> anyResult<()> {
    let mut sizes = target_sizes.into_iter().collect::<Vec<_>>();
    sizes.sort_by(|left, right| left.0.cmp(&right.0));
    let mut contigsizes = common_writer(output.join("_contigsizes").to_string_lossy().as_ref());
    for (target, length) in sizes {
        writeln!(contigsizes, "{target}\t{length}")?;
    }
    contigsizes.flush()?;
    std::fs::copy(metadata_source.join("_metadata"), output.join("_metadata"))?;
    if let Some(counts_source) = counts_source {
        std::fs::copy(counts_source, output.join("_metadata_counts"))?;
    }
    std::fs::write(output.join("_readme"), CONCAT_PQS_README)?;
    if preserve_cn_info {
        copy_cn_info(metadata_source, output)?;
    }
    Ok(())
}

fn write_concat_pqs_counts_file(output: &Path, counts: ConcatPqsCounts) -> anyResult<()> {
    let mut writer = common_writer(output.join("_metadata_counts").to_string_lossy().as_ref());
    writeln!(writer, "q0_records\t{}", counts.q0_records)?;
    writeln!(writer, "q1_records\t{}", counts.q1_records)?;
    writeln!(writer, "q0_concats\t{}", counts.q0_concats)?;
    writeln!(writer, "q1_concats\t{}", counts.q1_concats)?;
    writer.flush()?;
    Ok(())
}

fn transform_concat_pqs_native<F>(
    input: &Path,
    output: &Path,
    target_sizes: HashMap<String, u64>,
    threads: usize,
    preserve_cn_info: bool,
    transform: F,
) -> anyResult<()>
where
    F: Fn(usize, polars::prelude::DataFrame) -> anyResult<polars::prelude::DataFrame> + Sync + Send,
{
    use polars::prelude::*;

    if threads == 0 {
        bail!("concat PQS transform requires at least one thread");
    }
    if output.exists() {
        bail!("concat PQS output already exists: {}", output.display());
    }
    let shards = concat_pqs_files(input)?;
    std::fs::create_dir_all(output.join("q0"))?;
    std::fs::create_dir_all(output.join("q1"))?;
    polars::enable_string_cache();
    let q0_records = AtomicU64::new(0);
    let q1_records = AtomicU64::new(0);
    let q0_concats = AtomicU64::new(0);
    let q1_concats = AtomicU64::new(0);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("failed to create native concat PQS transform thread pool")?;

    pool.install(|| {
        shards
            .par_iter()
            .enumerate()
            .try_for_each(|(shard_id, shard)| -> anyResult<()> {
                let frame = ParquetReader::new(File::open(shard)?)
                    .finish()
                    .with_context(|| {
                        format!("failed to read concat PQS shard {}", shard.display())
                    })?;
                let mut q0 = transform(shard_id, frame)?;
                let mapq = q0
                    .column("mapping_quality")?
                    .as_materialized_series()
                    .u8()?;
                let mut q1 = q0.filter(&mapq.gt_eq(1))?;

                q0_records.fetch_add(q0.height() as u64, AtomOrdering::Relaxed);
                q1_records.fetch_add(q1.height() as u64, AtomOrdering::Relaxed);
                q0_concats.fetch_add(count_frame_concatemers(&q0)?, AtomOrdering::Relaxed);
                q1_concats.fetch_add(count_frame_concatemers(&q1)?, AtomOrdering::Relaxed);

                let q0_path = output.join(format!("q0/{shard_id}.parquet"));
                ParquetWriter::new(File::create(&q0_path)?)
                    .finish(&mut q0)
                    .with_context(|| format!("failed to write {}", q0_path.display()))?;
                if q1.height() > 0 {
                    let q1_path = output.join(format!("q1/{shard_id}.parquet"));
                    ParquetWriter::new(File::create(&q1_path)?)
                        .finish(&mut q1)
                        .with_context(|| format!("failed to write {}", q1_path.display()))?;
                }
                Ok(())
            })
    })?;

    let counts = ConcatPqsCounts {
        q0_records: q0_records.load(AtomOrdering::Relaxed),
        q1_records: q1_records.load(AtomOrdering::Relaxed),
        q0_concats: q0_concats.load(AtomOrdering::Relaxed),
        q1_concats: q1_concats.load(AtomOrdering::Relaxed),
    };
    write_concat_pqs_sidecars(input, output, target_sizes, None, preserve_cn_info)?;
    write_concat_pqs_counts_file(output, counts)?;
    Ok(())
}

fn rechunk_concat_pqs_native(
    input: &Path,
    output: &Path,
    chunksize: usize,
    threads: usize,
) -> anyResult<()> {
    use polars::prelude::*;

    if chunksize == 0 || threads == 0 {
        bail!("porec-split chunksize and thread count must be at least 1");
    }
    if output.exists() {
        bail!("concat PQS output already exists: {}", output.display());
    }
    let input_shards = concat_pqs_files(input)?;
    std::fs::create_dir_all(output.join("q0"))?;
    std::fs::create_dir_all(output.join("q1"))?;
    polars::enable_string_cache();

    let writer_threads = threads.saturating_div(3).clamp(1, 4);
    let (sender, receiver) = bounded::<(usize, DataFrame)>(writer_threads * 2);
    let q0_records = Arc::new(AtomicU64::new(0));
    let q1_records = Arc::new(AtomicU64::new(0));
    let q0_concats = Arc::new(AtomicU64::new(0));
    let q1_concats = Arc::new(AtomicU64::new(0));
    let mut handles = Vec::with_capacity(writer_threads);
    for _ in 0..writer_threads {
        let receiver = receiver.clone();
        let output = output.to_path_buf();
        let q0_records = Arc::clone(&q0_records);
        let q1_records = Arc::clone(&q1_records);
        let q0_concats = Arc::clone(&q0_concats);
        let q1_concats = Arc::clone(&q1_concats);
        handles.push(thread::spawn(move || -> anyResult<()> {
            while let Ok((shard_id, mut q0)) = receiver.recv() {
                let mapq = q0
                    .column("mapping_quality")?
                    .as_materialized_series()
                    .u8()?;
                let mut q1 = q0.filter(&mapq.gt_eq(1))?;
                q0_records.fetch_add(q0.height() as u64, AtomOrdering::Relaxed);
                q1_records.fetch_add(q1.height() as u64, AtomOrdering::Relaxed);
                q0_concats.fetch_add(count_frame_concatemers(&q0)?, AtomOrdering::Relaxed);
                q1_concats.fetch_add(count_frame_concatemers(&q1)?, AtomOrdering::Relaxed);
                let q0_path = output.join(format!("q0/{shard_id}.parquet"));
                ParquetWriter::new(File::create(&q0_path)?).finish(&mut q0)?;
                if q1.height() > 0 {
                    let q1_path = output.join(format!("q1/{shard_id}.parquet"));
                    ParquetWriter::new(File::create(&q1_path)?).finish(&mut q1)?;
                }
            }
            Ok(())
        }));
    }
    drop(receiver);

    let mut output_shard = 0usize;
    let mut accumulated = None::<DataFrame>;
    let mut accumulated_rows = 0usize;
    let mut accumulated_read_ids = std::collections::HashSet::<u64>::new();
    for input_shard in input_shards {
        let frame = ParquetReader::new(File::open(&input_shard)?)
            .finish()
            .with_context(|| {
                format!("failed to read concat PQS shard {}", input_shard.display())
            })?;
        let ids = frame.column("read_idx")?.as_materialized_series().u64()?;
        let mut group_start = 0usize;
        while group_start < frame.height() {
            let read_id = ids
                .get(group_start)
                .ok_or_else(|| anyhow!("null read_idx in {}", input_shard.display()))?;
            let mut group_end = group_start + 1;
            while group_end < frame.height() && ids.get(group_end) == Some(read_id) {
                group_end += 1;
            }
            let group_len = group_end - group_start;
            if accumulated_rows > 0
                && (accumulated_rows.saturating_add(group_len) > chunksize
                    || accumulated_read_ids.contains(&read_id))
            {
                sender
                    .send((output_shard, accumulated.take().expect("non-empty shard")))
                    .map_err(|_| anyhow!("porec-split writer stopped"))?;
                output_shard += 1;
                accumulated_rows = 0;
                accumulated_read_ids.clear();
            }
            let group = frame.slice(group_start as i64, group_len);
            if let Some(shard) = accumulated.as_mut() {
                shard.vstack_mut(&group)?;
            } else {
                accumulated = Some(group);
            }
            accumulated_rows += group_len;
            accumulated_read_ids.insert(read_id);
            group_start = group_end;
        }
    }
    if let Some(shard) = accumulated {
        sender
            .send((output_shard, shard))
            .map_err(|_| anyhow!("porec-split writer stopped before final shard"))?;
        output_shard += 1;
    }
    drop(sender);
    for handle in handles {
        handle
            .join()
            .map_err(|_| anyhow!("porec-split writer thread panicked"))??;
    }

    let counts = ConcatPqsCounts {
        q0_records: q0_records.load(AtomOrdering::Relaxed),
        q1_records: q1_records.load(AtomOrdering::Relaxed),
        q0_concats: q0_concats.load(AtomOrdering::Relaxed),
        q1_concats: q1_concats.load(AtomOrdering::Relaxed),
    };
    let mut sizes = BTreeMap::new();
    read_concat_pqs_contigsizes(input, &mut sizes)?;
    write_concat_pqs_sidecars(input, output, sizes.into_iter().collect(), None, true)?;
    let metadata_path = output.join("_metadata");
    let metadata = std::fs::read_to_string(&metadata_path)?
        .lines()
        .map(|line| {
            if line.trim_start().starts_with("'chunksize':") {
                format!(" 'chunksize': {chunksize},")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    std::fs::write(metadata_path, metadata)?;
    write_concat_pqs_counts_file(output, counts)?;
    log::info!(
        "Successful native concat PQS split into `{}`: {} complete-read shards",
        output.display(),
        output_shard
    );
    Ok(())
}

fn clone_concat_pqs_native(
    input: &Path,
    output: &Path,
    chromsizes: &Path,
    threads: usize,
) -> anyResult<()> {
    let mut input_sizes = BTreeMap::<String, u64>::new();
    let mut requested_sizes = BTreeMap::<String, u64>::new();
    read_concat_pqs_contigsizes(input, &mut input_sizes)?;
    read_contigsizes_file(chromsizes, &mut requested_sizes)?;
    if input_sizes != requested_sizes {
        bail!(
            "porec2pqs cannot zero-read clone {} because {} differs from its _contigsizes",
            input.display(),
            chromsizes.display()
        );
    }

    let q0_files = concat_pqs_files(input)?;
    std::fs::create_dir_all(output.join("q0"))?;
    std::fs::create_dir_all(output.join("q1"))?;
    let mut links = Vec::<(std::path::PathBuf, std::path::PathBuf)>::new();
    for (shard_id, q0_path) in q0_files.iter().enumerate() {
        links.push((
            q0_path.clone(),
            output.join(format!("q0/{shard_id}.parquet")),
        ));
        let q1_path = input
            .join("q1")
            .join(q0_path.file_name().expect("q0 shard filename"));
        if q1_path.is_file() {
            links.push((q1_path, output.join(format!("q1/{shard_id}.parquet"))));
        }
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("failed to create native porec2pqs clone thread pool")?;
    pool.install(|| {
        links
            .par_iter()
            .try_for_each(|(source, destination)| link_concat_pqs_shard(source, destination))
    })?;
    write_concat_pqs_sidecars(
        input,
        output,
        input_sizes.into_iter().collect(),
        Some(input.join("_metadata_counts")),
        true,
    )?;
    log::info!(
        "Successful zero-read linked porec2pqs clone of {} into `{}`: {} shards",
        input.display(),
        output.display(),
        q0_files.len()
    );
    Ok(())
}

fn merge_concat_pqs_native(input: &[&String], output: &str, threads: usize) -> anyResult<()> {
    if input.is_empty() {
        bail!("porec-merge requires at least one input");
    }
    if threads == 0 {
        bail!("porec-merge thread count must be at least 1");
    }
    let output_path = Path::new(output);
    if output_path.exists() {
        bail!(
            "native concat PQS output already exists: {}",
            output_path.display()
        );
    }

    let mut contigsizes = BTreeMap::<String, u64>::new();
    let mut counts = ConcatPqsCounts::default();
    let mut expected_position_type = None;
    let mut chunksize = 1usize;
    let mut output_shard_id = 0usize;
    let mut links = Vec::<(std::path::PathBuf, std::path::PathBuf)>::new();
    for input_path in input {
        let input_path = Path::new(input_path.as_str());
        if !is_concat_pqs(input_path) {
            bail!(
                "native PQS merge requires concat.pqs/porec.pqs inputs; `{}` is not one",
                input_path.display()
            );
        }
        read_concat_pqs_contigsizes(input_path, &mut contigsizes)?;
        counts = counts.checked_add(read_concat_pqs_counts(input_path)?, input_path)?;
        let (position_type, input_chunksize) = concat_pqs_metadata_layout(input_path)?;
        if let Some(expected) = expected_position_type {
            if expected != position_type {
                bail!(
                    "cannot zero-read merge incompatible position types {expected} and {position_type} in {}",
                    input_path.display()
                );
            }
        } else {
            expected_position_type = Some(position_type);
        }
        chunksize = chunksize.max(input_chunksize);
        let q0_files = concat_pqs_files(input_path)?;
        for (local_index, q0_path) in q0_files.iter().enumerate() {
            let output_id = output_shard_id + local_index;
            links.push((
                q0_path.clone(),
                output_path.join(format!("q0/{output_id}.parquet")),
            ));
            let q1_path = input_path
                .join("q1")
                .join(q0_path.file_name().expect("q0 shard filename"));
            if q1_path.is_file() {
                links.push((q1_path, output_path.join(format!("q1/{output_id}.parquet"))));
            }
        }
        output_shard_id += q0_files.len();
    }
    if counts.q0_records == 0 {
        bail!("native concat PQS merge produced no q0 records");
    }

    std::fs::create_dir_all(output_path.join("q0"))?;
    std::fs::create_dir_all(output_path.join("q1"))?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("failed to create native porec-merge thread pool")?;
    pool.install(|| {
        links
            .par_iter()
            .try_for_each(|(source, destination)| link_concat_pqs_shard(source, destination))
    })?;

    let mut contigsizes_writer =
        common_writer(output_path.join("_contigsizes").to_string_lossy().as_ref());
    for (chrom, length) in contigsizes {
        writeln!(contigsizes_writer, "{chrom}\t{length}")?;
    }
    contigsizes_writer.flush()?;
    let mut counts_writer = common_writer(
        output_path
            .join("_metadata_counts")
            .to_string_lossy()
            .as_ref(),
    );
    writeln!(counts_writer, "q0_records\t{}", counts.q0_records)?;
    writeln!(counts_writer, "q1_records\t{}", counts.q1_records)?;
    writeln!(counts_writer, "q0_concats\t{}", counts.q0_concats)?;
    writeln!(counts_writer, "q1_concats\t{}", counts.q1_concats)?;
    counts_writer.flush()?;
    let metadata = concat_pqs_metadata(
        chunksize,
        expected_position_type.expect("validated at least one PQS input"),
    )
    .replace("'read_idx_scope': 'global'", "'read_idx_scope': 'shard'");
    std::fs::write(output_path.join("_metadata"), metadata)?;
    std::fs::write(output_path.join("_readme"), CONCAT_PQS_README)?;
    merge_cn_info(input, output_path)?;
    log::info!(
        "Successful zero-read linked merge of {} concat PQS inputs into `{}`: {} q0 records, {} complete concatemers, {} shards",
        input.len(),
        output,
        counts.q0_records,
        counts.q0_concats,
        output_shard_id
    );
    Ok(())
}

pub fn merge_porec_tables(input: Vec<&String>, output: &String, threads: usize) -> anyResult<()> {
    let output_path = output.trim_end_matches('/');
    if output_path.ends_with(".concat.pqs") || output_path.ends_with(".porec.pqs") {
        let result = merge_concat_pqs_native(&input, output_path, threads);
        if result.is_err() && Path::new(output_path).exists() {
            let _ = std::fs::remove_dir_all(output_path);
        }
        return result;
    }

    let mut wtr = common_writer(output);
    let mut idx = 0;

    for file in input {
        let reader = open_porec_reader(file).expect("Failed to open Pore-C input");
        let mut max_idx: u64 = 0;
        for (i, line) in reader.lines().enumerate() {
            let line = line.unwrap();
            let mut line = line.split("\t");
            let mut read_idx: u64 = line.next().unwrap().parse().unwrap();

            if read_idx > max_idx {
                max_idx = read_idx;
            }

            read_idx += idx as u64;
            let mut record = line.collect::<Vec<&str>>().join("\t");
            record = format!("{}\t{}", read_idx, record);
            wtr.write_all(record.to_string().as_bytes()).unwrap();
            wtr.write_all(b"\n").unwrap();
        }
        idx += max_idx + 1;
    }
    log::info!("Successful output merge porec tables into `{}`", output);
    Ok(())
}
