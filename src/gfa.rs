use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use hashbrown::HashMap;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::BTreeSet;
use std::env;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Instant;
use tempfile::{Builder, NamedTempFile, tempdir_in, tempfile_in};

use crate::clm::{CLMB_DEFAULT_BLOCK_SIZE, ClmbReader, ClmbWriter, encode_endpoint, is_clmb_file};
use crate::core::common_reader;

const IO_BUFFER_SIZE: usize = 256 * 1024;
const CONTACT_CHUNK_SIZE: usize = 8 * 1024 * 1024;
const CLUSTER_CLMB_BLOCK_SIZE: usize = 256 * 1024;

type EndContactKey = (u32, u8, u32, u8);
type SegmentEndKey = (u32, u8);
type RankedContact = (EndContactKey, f64);

fn canonical_contact_key(left: SegmentEndKey, right: SegmentEndKey) -> EndContactKey {
    let (left, right) = if left <= right {
        (left, right)
    } else {
        (right, left)
    };
    (left.0, left.1, right.0, right.1)
}

fn contact_ends(key: EndContactKey) -> (SegmentEndKey, SegmentEndKey) {
    ((key.0, key.1), (key.2, key.3))
}

struct SegmentTable {
    ids: HashMap<String, u32>,
    names: Vec<String>,
    lengths: Vec<u64>,
}

struct ContactChunkReader {
    reader: Box<dyn BufRead + Send + 'static>,
    next_line: usize,
    finished: bool,
}

impl ContactChunkReader {
    fn new(path: &str) -> Self {
        Self {
            reader: common_reader(path),
            next_line: 1,
            finished: false,
        }
    }
}

impl Iterator for ContactChunkReader {
    type Item = Result<(usize, Vec<u8>)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let start_line = self.next_line;
        let mut chunk = Vec::with_capacity(CONTACT_CHUNK_SIZE + IO_BUFFER_SIZE);
        let mut lines = 0usize;
        while chunk.len() < CONTACT_CHUNK_SIZE {
            match self.reader.read_until(b'\n', &mut chunk) {
                Ok(0) => {
                    self.finished = true;
                    break;
                }
                Ok(_) => lines += 1,
                Err(error) => {
                    self.finished = true;
                    return Some(Err(error.into()));
                }
            }
        }
        if chunk.is_empty() {
            None
        } else {
            self.next_line += lines;
            Some(Ok((start_line, chunk)))
        }
    }
}

fn parse_contact_end(value: &str, line_number: usize) -> Result<(&str, u8)> {
    let (segment, suffix) = value.rsplit_once('_').with_context(|| {
        format!(
            "invalid split-contact end label `{}` on line {}",
            value, line_number
        )
    })?;
    if segment.is_empty() {
        bail!(
            "empty segment in split-contact end label on line {}",
            line_number
        );
    }
    let side = match suffix {
        "0" => 0,
        "1" => 1,
        _ => bail!(
            "invalid split-contact end label `{}` on line {}; expected '<segment>_0' or '<segment>_1'",
            value,
            line_number
        ),
    };
    Ok((segment, side))
}

fn aggregate_contact_chunk(
    segments: &SegmentTable,
    start_line: usize,
    chunk: Vec<u8>,
) -> Result<HashMap<EndContactKey, f64>> {
    let mut counts = HashMap::<EndContactKey, f64>::new();
    for (offset, raw_line) in chunk.split(|value| *value == b'\n').enumerate() {
        let line_number = start_line + offset;
        let raw_line = raw_line.strip_suffix(b"\r").unwrap_or(raw_line);
        if raw_line.is_empty() || raw_line.starts_with(b"#") {
            continue;
        }
        let record = std::str::from_utf8(raw_line)
            .with_context(|| format!("non-UTF-8 split-contact line {}", line_number))?;
        let mut fields = record.split('\t');
        let raw1 = fields
            .next()
            .with_context(|| format!("missing first endpoint on line {}", line_number))?;
        let raw2 = fields
            .next()
            .with_context(|| format!("missing second endpoint on line {}", line_number))?;
        let raw_count = fields
            .next()
            .with_context(|| format!("missing contact count on line {}", line_number))?;
        let (segment1, side1) = parse_contact_end(raw1, line_number)?;
        let (segment2, side2) = parse_contact_end(raw2, line_number)?;
        let (Some(&segment1), Some(&segment2)) =
            (segments.ids.get(segment1), segments.ids.get(segment2))
        else {
            continue;
        };
        if segment1 == segment2 {
            continue;
        }
        let count = raw_count
            .parse::<f64>()
            .with_context(|| format!("invalid contact count on line {}", line_number))?;
        if !count.is_finite() || count < 0.0 {
            bail!(
                "contact count on line {} must be finite and non-negative",
                line_number
            );
        }
        let key = canonical_contact_key((segment1, side1), (segment2, side2));
        *counts.entry(key).or_insert(0.0) += count;
    }
    Ok(counts)
}

fn parse_gfa_links(path: &str, segments: &SegmentTable) -> Result<Vec<EndContactKey>> {
    let mut links = Vec::new();
    let mut reader = common_reader(path);
    let mut line = String::new();
    let mut line_number = 0usize;
    while reader.read_line(&mut line)? != 0 {
        line_number += 1;
        let record = line.trim_end_matches(&['\r', '\n'][..]);
        if record.is_empty() || record.starts_with('#') {
            line.clear();
            continue;
        }
        let fields: Vec<_> = record.split('\t').collect();
        if fields.len() != 4 {
            bail!(
                "malformed GFA link line {}: expected 4 tab-separated fields",
                line_number
            );
        }
        let side = |value: &str| -> Result<u8> {
            match value {
                "L" => Ok(0),
                "R" => Ok(1),
                _ => bail!(
                    "invalid GFA segment end `{}` on line {}",
                    value,
                    line_number
                ),
            }
        };
        let segment1 = *segments.ids.get(fields[0]).with_context(|| {
            format!(
                "unknown GFA-link segment `{}` on line {}",
                fields[0], line_number
            )
        })?;
        let segment2 = *segments.ids.get(fields[2]).with_context(|| {
            format!(
                "unknown GFA-link segment `{}` on line {}",
                fields[2], line_number
            )
        })?;
        links.push(canonical_contact_key(
            (segment1, side(fields[1])?),
            (segment2, side(fields[3])?),
        ));
        line.clear();
    }
    Ok(links)
}

fn contact_score(key: EndContactKey, count: f64, segments: &SegmentTable) -> f64 {
    let half_length1 = (segments.lengths[key.0 as usize] as f64 / 2.0).max(1.0);
    let half_length2 = (segments.lengths[key.2 as usize] as f64 / 2.0).max(1.0);
    count / (half_length1 * half_length2)
}

fn update_top_contacts(
    top_contacts: &mut HashMap<SegmentEndKey, [Option<RankedContact>; 2]>,
    endpoint: SegmentEndKey,
    key: EndContactKey,
    score: f64,
) {
    let top = top_contacts.entry(endpoint).or_insert([None, None]);
    if top[0].is_none_or(|(_, current)| score > current) {
        top[1] = top[0];
        top[0] = Some((key, score));
    } else if top[0].is_some_and(|(current_key, _)| current_key != key)
        && top[1].is_none_or(|(_, current)| score > current)
    {
        top[1] = Some((key, score));
    }
}

fn competitor_score(
    top_contacts: &HashMap<SegmentEndKey, [Option<RankedContact>; 2]>,
    endpoint: SegmentEndKey,
    excluded: EndContactKey,
) -> f64 {
    let Some(top) = top_contacts.get(&endpoint) else {
        return 0.0;
    };
    for candidate in top.iter().flatten() {
        if candidate.0 != excluded {
            return candidate.1;
        }
    }
    0.0
}

fn merge_contact_counts(
    mut left: HashMap<EndContactKey, f64>,
    right: HashMap<EndContactKey, f64>,
) -> Result<HashMap<EndContactKey, f64>> {
    if left.len() < right.len() {
        return merge_contact_counts(right, left);
    }
    for (key, count) in right {
        *left.entry(key).or_insert(0.0) += count;
    }
    Ok(left)
}

/// Aggregate and normalize split-contact evidence between allowed GFA ends.
///
/// The output is deliberately simple TSV so Python can retain the existing
/// topology matching and audit-table implementation without scanning the
/// complete contacts file itself.
pub fn aggregate_gfa_end_contacts(
    segments_path: &str,
    contacts_path: &str,
    gfa_links_path: &str,
    output_path: &str,
    threads: usize,
) -> Result<()> {
    let mut segment_lengths = HashMap::<String, u64>::new();
    let mut reader = common_reader(segments_path);
    let mut line = String::new();
    let mut line_number = 0usize;
    while reader.read_line(&mut line)? != 0 {
        line_number += 1;
        let record = line.trim_end_matches(&['\r', '\n'][..]);
        if record.is_empty() || record.starts_with('#') {
            line.clear();
            continue;
        }
        let mut fields = record.split('\t');
        let segment = fields.next().context("missing segment name")?;
        let length = fields
            .next()
            .with_context(|| format!("missing segment length on line {}", line_number))?
            .parse::<u64>()
            .with_context(|| format!("invalid segment length on line {}", line_number))?;
        segment_lengths.insert(segment.to_string(), length);
        line.clear();
    }

    let mut names: Vec<_> = segment_lengths.keys().cloned().collect();
    names.sort_unstable();
    let lengths = names.iter().map(|name| segment_lengths[name]).collect();
    let ids = names
        .iter()
        .enumerate()
        .map(|(index, name)| (name.clone(), index as u32))
        .collect();
    let segments = SegmentTable {
        ids,
        names,
        lengths,
    };
    let links = parse_gfa_links(gfa_links_path, &segments)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("cannot create GFA contact aggregation thread pool")?;
    let counts = pool.install(|| {
        ContactChunkReader::new(contacts_path)
            .par_bridge()
            .map(|chunk| {
                let (start_line, chunk) = chunk?;
                aggregate_contact_chunk(&segments, start_line, chunk)
            })
            .try_reduce(HashMap::new, merge_contact_counts)
    })?;

    let mut top_contacts = HashMap::<SegmentEndKey, [Option<RankedContact>; 2]>::new();
    for (&key, &count) in &counts {
        let score = contact_score(key, count, &segments);
        let (end1, end2) = contact_ends(key);
        update_top_contacts(&mut top_contacts, end1, key, score);
        update_top_contacts(&mut top_contacts, end2, key, score);
    }

    let mut writer = fast_writer(output_path)?;
    writeln!(
        writer,
        "#Segment1\tEnd1\tSegment2\tEnd2\tCount\tScore\tCompetitor1\tCompetitor2"
    )?;
    for key in links {
        let count = counts.get(&key).copied().unwrap_or(0.0);
        let score = contact_score(key, count, &segments);
        let (end1, end2) = contact_ends(key);
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            segments.names[key.0 as usize],
            if key.1 == 0 { 'L' } else { 'R' },
            segments.names[key.2 as usize],
            if key.3 == 0 { 'L' } else { 'R' },
            count,
            score,
            competitor_score(&top_contacts, end1, key),
            competitor_score(&top_contacts, end2, key),
        )?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Clone, Debug)]
struct OrientedMapping {
    orientation: char,
    left_offset: u64,
    right_offset: u64,
}

#[derive(Clone, Debug)]
struct SegmentMapping {
    unit: String,
    halves: [String; 2],
    plus: OrientedMapping,
    minus: OrientedMapping,
}

#[derive(Debug)]
struct MappingTable {
    segments: HashMap<String, SegmentMapping>,
    first_unit: Option<String>,
}

impl MappingTable {
    fn from_path(path: &str) -> Result<Self> {
        let mut segments = HashMap::new();
        let mut first_unit = None;
        let mut reader = common_reader(path);
        let mut line = String::new();
        let mut line_number = 0usize;
        while reader.read_line(&mut line)? != 0 {
            line_number += 1;
            let record = line.trim_end_matches(&['\r', '\n'][..]);
            if record.is_empty() || record.starts_with('#') {
                line.clear();
                continue;
            }
            let fields: Vec<&str> = record.split('\t').collect();
            if fields.len() != 10 {
                bail!(
                    "Malformed GFA remap table at line {}: expected 10 tab-separated fields, found {}",
                    line_number,
                    fields.len()
                );
            }
            let parse_orientation = |value: &str| -> Result<char> {
                match value {
                    "+" => Ok('+'),
                    "-" => Ok('-'),
                    _ => bail!("invalid mapped orientation `{}`", value),
                }
            };
            let parse_offset = |value: &str| -> Result<u64> {
                value
                    .parse::<u64>()
                    .with_context(|| format!("invalid block offset `{}`", value))
            };
            let mapping = SegmentMapping {
                unit: fields[1].to_string(),
                halves: [fields[2].to_string(), fields[3].to_string()],
                plus: OrientedMapping {
                    orientation: parse_orientation(fields[4])?,
                    left_offset: parse_offset(fields[5])?,
                    right_offset: parse_offset(fields[6])?,
                },
                minus: OrientedMapping {
                    orientation: parse_orientation(fields[7])?,
                    left_offset: parse_offset(fields[8])?,
                    right_offset: parse_offset(fields[9])?,
                },
            };
            if first_unit.is_none() {
                first_unit = Some(mapping.unit.clone());
            }
            segments.insert(fields[0].to_string(), mapping);
            line.clear();
        }
        Ok(Self {
            segments,
            first_unit,
        })
    }

    fn map_clm_token(&self, token: &str, left_role: bool) -> Option<(&str, char, u64)> {
        let (segment, orientation) = if let Some(segment) = token.strip_suffix('+') {
            (segment, '+')
        } else if let Some(segment) = token.strip_suffix('-') {
            (segment, '-')
        } else {
            return None;
        };
        self.map_clm_parts(segment, orientation, left_role)
    }

    fn map_clm_parts(
        &self,
        segment: &str,
        orientation: char,
        left_role: bool,
    ) -> Option<(&str, char, u64)> {
        let mapping = self.segments.get(segment)?;
        let oriented = if orientation == '+' {
            &mapping.plus
        } else if orientation == '-' {
            &mapping.minus
        } else {
            return None;
        };
        let offset = if left_role {
            oriented.left_offset
        } else {
            oriented.right_offset
        };
        Some((&mapping.unit, oriented.orientation, offset))
    }

    fn units(&self) -> Vec<String> {
        self.segments
            .values()
            .map(|mapping| mapping.unit.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }
}

#[derive(Clone, Copy)]
struct ContractedHalfMapping {
    unit: u32,
    halves: [u32; 2],
}

struct ContactContractionTable {
    segments: HashMap<String, ContractedHalfMapping>,
    labels: Vec<String>,
    first_unit: Option<String>,
}

impl ContactContractionTable {
    fn from_mapping(mapping: &MappingTable) -> Self {
        let mut units = mapping
            .segments
            .values()
            .map(|item| item.unit.clone())
            .collect::<Vec<_>>();
        units.sort_unstable();
        units.dedup();
        let unit_ids = units
            .into_iter()
            .enumerate()
            .map(|(index, unit)| (unit, index as u32))
            .collect::<HashMap<_, _>>();

        let mut labels = mapping
            .segments
            .values()
            .flat_map(|item| item.halves.iter().cloned())
            .collect::<Vec<_>>();
        labels.sort_unstable();
        labels.dedup();
        let label_ids = labels
            .iter()
            .enumerate()
            .map(|(index, label)| (label.clone(), index as u32))
            .collect::<HashMap<_, _>>();

        let segments = mapping
            .segments
            .iter()
            .map(|(segment, item)| {
                (
                    segment.clone(),
                    ContractedHalfMapping {
                        unit: unit_ids[&item.unit],
                        halves: [label_ids[&item.halves[0]], label_ids[&item.halves[1]]],
                    },
                )
            })
            .collect();
        Self {
            segments,
            labels,
            first_unit: mapping.first_unit.clone(),
        }
    }
}

type ContractedContactKey = (u32, u32);

fn parse_half_label(label: &str) -> Option<(&str, usize)> {
    let (segment, suffix) = label.rsplit_once('_')?;
    let index = match suffix {
        "0" => 0,
        "1" => 1,
        _ => return None,
    };
    Some((segment, index))
}

fn contract_contact_chunk(
    mapping: &ContactContractionTable,
    start_line: usize,
    chunk: Vec<u8>,
) -> Result<HashMap<ContractedContactKey, f64>> {
    let mut counts = HashMap::<ContractedContactKey, f64>::new();
    // split.contacts is ordered by segment pair.  Cache both segment lookups
    // within each chunk: the left segment commonly spans thousands of rows,
    // and the right segment commonly spans all four half combinations.
    let mut cached1: Option<(&str, Option<ContractedHalfMapping>)> = None;
    let mut cached2: Option<(&str, Option<ContractedHalfMapping>)> = None;
    for (offset, raw_line) in chunk.split(|value| *value == b'\n').enumerate() {
        let line_number = start_line + offset;
        let raw_line = raw_line.strip_suffix(b"\r").unwrap_or(raw_line);
        if raw_line.is_empty() || raw_line.starts_with(b"#") {
            continue;
        }
        let record = std::str::from_utf8(raw_line)
            .with_context(|| format!("non-UTF-8 split-contact line {}", line_number))?;
        let mut fields = record.splitn(3, '\t');
        let raw1 = fields
            .next()
            .with_context(|| format!("missing first contacts endpoint on line {}", line_number))?;
        let raw2 = fields
            .next()
            .with_context(|| format!("missing second contacts endpoint on line {}", line_number))?;
        let raw_count = fields
            .next()
            .with_context(|| format!("missing contacts count on line {}", line_number))?;
        let (Some((segment1, half1)), Some((segment2, half2))) =
            (parse_half_label(raw1), parse_half_label(raw2))
        else {
            continue;
        };
        let mapping1 = if cached1.is_some_and(|(segment, _)| segment == segment1) {
            cached1.unwrap().1
        } else {
            let value = mapping.segments.get(segment1).copied();
            cached1 = Some((segment1, value));
            value
        };
        let mapping2 = if cached2.is_some_and(|(segment, _)| segment == segment2) {
            cached2.unwrap().1
        } else {
            let value = mapping.segments.get(segment2).copied();
            cached2 = Some((segment2, value));
            value
        };
        let (Some(mapping1), Some(mapping2)) = (mapping1, mapping2) else {
            continue;
        };
        let mapped1 = (mapping1.unit, mapping1.halves[half1]);
        let mapped2 = (mapping2.unit, mapping2.halves[half2]);
        if mapped1.0 == mapped2.0 {
            continue;
        }
        let count = raw_count.parse::<f64>().with_context(|| {
            format!(
                "invalid contacts count `{}` on line {}",
                raw_count, line_number
            )
        })?;
        let key = if mapped1.1 <= mapped2.1 {
            (mapped1.1, mapped2.1)
        } else {
            (mapped2.1, mapped1.1)
        };
        *counts.entry(key).or_insert(0.0) += count;
    }
    Ok(counts)
}

fn merge_contracted_contacts(
    mut left: HashMap<ContractedContactKey, f64>,
    right: HashMap<ContractedContactKey, f64>,
) -> Result<HashMap<ContractedContactKey, f64>> {
    if left.len() < right.len() {
        return merge_contracted_contacts(right, left);
    }
    for (key, count) in right {
        *left.entry(key).or_insert(0.0) += count;
    }
    Ok(left)
}

fn fast_writer(path: &str) -> Result<Box<dyn Write>> {
    let file = File::create(path).with_context(|| format!("cannot create `{}`", path))?;
    let buffered = BufWriter::with_capacity(IO_BUFFER_SIZE, file);
    if Path::new(path).extension().and_then(|value| value.to_str()) == Some("gz") {
        Ok(Box::new(GzEncoder::new(buffered, Compression::fast())))
    } else {
        Ok(Box::new(buffered))
    }
}

fn run_sort(
    source: &Path,
    output: &Path,
    temporary_directory: &Path,
    keys: &[&str],
    buffer_size: &str,
    threads: usize,
) -> Result<()> {
    let mut command = Command::new("sort");
    command
        .arg("--stable")
        .arg("--field-separator=\t")
        .arg(format!("--buffer-size={}", buffer_size))
        .arg(format!("--parallel={}", threads.max(1)))
        .arg(format!(
            "--temporary-directory={}",
            temporary_directory.display()
        ));
    for key in keys {
        command.arg(format!("--key={}", key));
    }
    let result = command
        .arg(format!("--output={}", output.display()))
        .arg(source)
        .env("LC_ALL", "C")
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .output()
        .context("failed to start `sort`; install GNU coreutils or ensure `sort` is on PATH")?;
    if !result.status.success() {
        bail!(
            "external sort failed: {}",
            String::from_utf8_lossy(&result.stderr).trim()
        );
    }
    Ok(())
}

fn temporary_root(requested: Option<&str>) -> Result<PathBuf> {
    if let Some(path) = requested {
        let path = PathBuf::from(path);
        if !path.is_dir() {
            bail!("temporary directory `{}` does not exist", path.display());
        }
        Ok(path)
    } else {
        Ok(env::temp_dir())
    }
}

fn remap_contacts(mapping: &MappingTable, input: &str, output: &str, threads: usize) -> Result<()> {
    let contraction = ContactContractionTable::from_mapping(mapping);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("cannot create contact contraction thread pool")?;
    let counts = pool.install(|| {
        ContactChunkReader::new(input)
            .par_bridge()
            .map(|chunk| {
                let (start_line, chunk) = chunk?;
                contract_contact_chunk(&contraction, start_line, chunk)
            })
            .try_reduce(HashMap::new, merge_contracted_contacts)
    })?;

    let mut contacts = counts.into_iter().collect::<Vec<_>>();
    let sort_started = Instant::now();
    contacts.sort_unstable_by_key(|item| item.0);
    log::info!(
        "GFA preprocessing timing: output sorting: {:.3} s",
        sort_started.elapsed().as_secs_f64()
    );
    let mut writer = fast_writer(output)?;
    for ((end1, end2), count) in &contacts {
        write_contact_record(
            &mut writer,
            &contraction.labels[*end1 as usize],
            &contraction.labels[*end2 as usize],
            *count,
        )?;
    }
    if contacts.is_empty() {
        if let Some(unit) = &contraction.first_unit {
            writeln!(writer, "{}_0\t{}_1\t0", unit, unit)?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_contact_record(writer: &mut dyn Write, end1: &str, end2: &str, count: f64) -> Result<()> {
    if count.fract() == 0.0 {
        writeln!(writer, "{}\t{}\t{:.0}", end1, end2, count)?;
    } else {
        writeln!(writer, "{}\t{}\t{}", end1, end2, count)?;
    }
    Ok(())
}

fn flush_clm_group(spool: &mut File, writer: &mut dyn Write, key: &str, count: u64) -> Result<()> {
    spool.flush()?;
    spool.seek(SeekFrom::Start(0))?;
    write!(writer, "{}\t{}\t", key, count)?;
    io::copy(spool, writer)?;
    writer.write_all(b"\n")?;
    spool.set_len(0)?;
    spool.seek(SeekFrom::Start(0))?;
    Ok(())
}

fn clm_key(record: &str) -> Result<&str> {
    record
        .split_once('\t')
        .map(|fields| fields.0)
        .context("missing sorted CLM key")
}

fn clm_chunk_boundaries(path: &Path, threads: usize) -> Result<Vec<u64>> {
    let length = path.metadata()?.len();
    if length == 0 {
        return Ok(vec![0]);
    }
    let mut boundaries = vec![0];
    for index in 1..threads.max(1) {
        let target = length.saturating_mul(index as u64) / threads.max(1) as u64;
        let mut reader = BufReader::with_capacity(IO_BUFFER_SIZE, File::open(path)?);
        reader.seek(SeekFrom::Start(target))?;
        if target > 0 {
            let mut partial = String::new();
            reader.read_line(&mut partial)?;
        }
        let mut first_key = String::new();
        let mut line = String::new();
        loop {
            let position = reader.stream_position()?;
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let key = clm_key(line.trim_end_matches(&['\r', '\n'][..]))?;
            if first_key.is_empty() {
                first_key.push_str(key);
            } else if first_key != key {
                boundaries.push(position);
                break;
            }
        }
    }
    boundaries.push(length);
    boundaries.sort_unstable();
    boundaries.dedup();
    Ok(boundaries)
}

fn process_clm_chunk(
    sorted_path: &Path,
    output_path: &Path,
    temporary_directory: &Path,
    start: u64,
    end: u64,
) -> Result<()> {
    let mut file = File::open(sorted_path)?;
    file.seek(SeekFrom::Start(start))?;
    let mut reader = BufReader::with_capacity(IO_BUFFER_SIZE, file.take(end - start));
    let mut writer = fast_writer(
        output_path
            .to_str()
            .context("temporary CLM output path is not valid UTF-8")?,
    )?;
    let mut spool = tempfile_in(temporary_directory)?;
    let mut line = String::new();
    let mut current_key = String::new();
    let mut distance_count = 0u64;
    while reader.read_line(&mut line)? != 0 {
        let record = line.trim_end_matches(&['\r', '\n'][..]);
        let mut fields = record.splitn(4, '\t');
        let key = fields.next().context("missing sorted CLM key")?;
        fields.next().context("missing sorted CLM sequence")?;
        let shift = fields
            .next()
            .context("missing sorted CLM shift")?
            .parse::<u128>()?;
        let distances = fields.next().context("missing sorted CLM distances")?;
        if !current_key.is_empty() && current_key != key {
            flush_clm_group(&mut spool, &mut writer, &current_key, distance_count)?;
            distance_count = 0;
        }
        if current_key.is_empty() || current_key != key {
            current_key.clear();
            current_key.push_str(key);
        }
        let mut adjusted = String::with_capacity(distances.len() + 32);
        for value in distances.split_whitespace() {
            let original = value
                .parse::<i128>()
                .with_context(|| format!("invalid CLM distance `{}`", value))?;
            let shifted = (original + shift as i128).max(2);
            if distance_count > 0 || !adjusted.is_empty() {
                adjusted.push(' ');
            }
            write!(adjusted, "{}", shifted)?;
            distance_count += 1;
        }
        spool.write_all(adjusted.as_bytes())?;
        line.clear();
    }
    if !current_key.is_empty() {
        flush_clm_group(&mut spool, &mut writer, &current_key, distance_count)?;
    }
    writer.flush()?;
    Ok(())
}

fn process_clm_chunks(
    sorted_path: &Path,
    output: &str,
    temporary_directory: &Path,
    threads: usize,
) -> Result<()> {
    let boundaries = clm_chunk_boundaries(sorted_path, threads)?;
    if boundaries.len() < 2 {
        let mut writer = fast_writer(output)?;
        writer.flush()?;
        return Ok(());
    }
    let compressed = Path::new(output)
        .extension()
        .and_then(|value| value.to_str())
        == Some("gz");
    let suffix = if compressed { ".clm.gz" } else { ".clm" };
    let shards = boundaries
        .windows(2)
        .map(|_| {
            Builder::new()
                .suffix(suffix)
                .tempfile_in(temporary_directory)
        })
        .collect::<io::Result<Vec<_>>>()?;

    std::thread::scope(|scope| -> Result<()> {
        let handles = boundaries
            .windows(2)
            .zip(shards.iter())
            .map(|(range, shard)| {
                let sorted_path = sorted_path.to_path_buf();
                let shard_path = shard.path().to_path_buf();
                let temporary_directory = temporary_directory.to_path_buf();
                let start = range[0];
                let end = range[1];
                scope.spawn(move || {
                    process_clm_chunk(&sorted_path, &shard_path, &temporary_directory, start, end)
                })
            })
            .collect::<Vec<_>>();
        for handle in handles {
            handle
                .join()
                .map_err(|_| anyhow::anyhow!("CLM block worker panicked"))??;
        }
        Ok(())
    })?;

    let mut writer = BufWriter::with_capacity(IO_BUFFER_SIZE, File::create(output)?);
    for shard in shards {
        let mut reader = BufReader::with_capacity(IO_BUFFER_SIZE, File::open(shard.path())?);
        io::copy(&mut reader, &mut writer)?;
    }
    writer.flush()?;
    Ok(())
}

#[allow(dead_code)]
fn parse_sort_buffer(value: &str) -> Option<u64> {
    let (digits, multiplier) = match value.as_bytes().last().copied() {
        Some(b'K') | Some(b'k') => (&value[..value.len() - 1], 1024u64),
        Some(b'M') | Some(b'm') => (&value[..value.len() - 1], 1024u64.pow(2)),
        Some(b'G') | Some(b'g') => (&value[..value.len() - 1], 1024u64.pow(3)),
        _ => (value, 1u64),
    };
    digits.parse::<u64>().ok()?.checked_mul(multiplier)
}

#[allow(dead_code)]
fn divided_sort_buffer(value: &str, divisor: usize) -> String {
    parse_sort_buffer(value)
        .map(|bytes| (bytes / divisor.max(1) as u64).max(1024 * 1024).to_string())
        .unwrap_or_else(|| value.to_string())
}

fn parse_clm_endpoint_token(token: &str) -> Result<(&str, char)> {
    if let Some(segment) = token.strip_suffix('+') {
        Ok((segment, '+'))
    } else if let Some(segment) = token.strip_suffix('-') {
        Ok((segment, '-'))
    } else {
        bail!("CLM endpoint `{token}` has no orientation")
    }
}

fn adjusted_clm_distances(distances: &[u64], shift: u64) -> Result<Cow<'_, [u64]>> {
    if shift == 0 && distances.iter().all(|distance| *distance >= 2) {
        return Ok(Cow::Borrowed(distances));
    }
    distances
        .iter()
        .map(|distance| {
            distance
                .checked_add(shift)
                .map(|value| value.max(2))
                .context("contracted CLM distance exceeds u64")
        })
        .collect::<Result<Vec<_>>>()
        .map(Cow::Owned)
}

fn stream_clm_records(
    input: &str,
    mut consume: impl FnMut(&str, char, &str, char, &[u64]) -> Result<()>,
) -> Result<()> {
    if is_clmb_file(input)? {
        let mut reader = ClmbReader::open(input)?;
        let contigs = reader.header.contigs.clone();
        while let Some(block) = reader.next_block()? {
            for record in block {
                let orientation1 = if record.orientation1() == 0 { '+' } else { '-' };
                let orientation2 = if record.orientation2() == 0 { '+' } else { '-' };
                consume(
                    &contigs[record.contig1() as usize],
                    orientation1,
                    &contigs[record.contig2() as usize],
                    orientation2,
                    &record.distances,
                )?;
            }
        }
        return Ok(());
    }

    let mut reader = common_reader(input);
    let mut line = String::new();
    let mut line_number = 0usize;
    while reader.read_line(&mut line)? != 0 {
        line_number += 1;
        let record = line.trim_end_matches(&['\r', '\n'][..]);
        if record.is_empty() || record.starts_with('#') {
            line.clear();
            continue;
        }
        let mut fields = record.splitn(3, '\t');
        let pair = fields.next().context("missing CLM contig pair")?;
        let declared_count = fields
            .next()
            .context("missing CLM count")?
            .parse::<usize>()
            .with_context(|| format!("invalid CLM count on line {line_number}"))?;
        let distance_field = fields.next().context("missing CLM distances")?;
        let mut tokens = pair.split_whitespace();
        let token1 = tokens.next().context("missing first CLM token")?;
        let token2 = tokens.next().context("missing second CLM token")?;
        if tokens.next().is_some() {
            bail!("too many CLM endpoints on line {line_number}");
        }
        let (segment1, orientation1) = parse_clm_endpoint_token(token1)?;
        let (segment2, orientation2) = parse_clm_endpoint_token(token2)?;
        let distances = distance_field
            .split_whitespace()
            .map(|value| {
                value.parse::<u64>().with_context(|| {
                    format!("invalid CLM distance `{value}` on line {line_number}")
                })
            })
            .collect::<Result<Vec<_>>>()?;
        if distances.len() != declared_count {
            bail!(
                "CLM count mismatch on line {line_number}: declared {declared_count}, observed {}",
                distances.len()
            );
        }
        consume(segment1, orientation1, segment2, orientation2, &distances)?;
        line.clear();
    }
    Ok(())
}

fn remap_clm_to_clmb(mapping: &MappingTable, input: &str, output: &str) -> Result<()> {
    let units = mapping.units();
    let unit_ids = units
        .iter()
        .enumerate()
        .map(|(id, unit)| Ok((unit.clone(), u32::try_from(id)?)))
        .collect::<Result<HashMap<_, _>>>()?;
    let mut writer = ClmbWriter::create(output, &units, CLMB_DEFAULT_BLOCK_SIZE, None, None)?;
    stream_clm_records(
        input,
        |segment1, orientation1, segment2, orientation2, distances| {
            let (Some(mapped1), Some(mapped2)) = (
                mapping.map_clm_parts(segment1, orientation1, true),
                mapping.map_clm_parts(segment2, orientation2, false),
            ) else {
                return Ok(());
            };
            if mapped1.0 == mapped2.0 {
                return Ok(());
            }
            let shift = mapped1
                .2
                .checked_add(mapped2.2)
                .context("contracted CLM shift exceeds u64")?;
            let adjusted = adjusted_clm_distances(distances, shift)?;
            let endpoint1 =
                encode_endpoint(unit_ids[mapped1.0], if mapped1.1 == '+' { 0 } else { 1 })?;
            let endpoint2 =
                encode_endpoint(unit_ids[mapped2.0], if mapped2.1 == '+' { 0 } else { 1 })?;
            writer.write_record(endpoint1, endpoint2, adjusted.as_ref())
        },
    )?;
    writer.finish()
}

fn remap_clm_by_cluster_to_clmb(
    mapping: &MappingTable,
    cluster_path: &str,
    input: &str,
    output_directory: &str,
) -> Result<()> {
    let output_directory = PathBuf::from(output_directory);
    std::fs::create_dir_all(&output_directory)?;
    let units = mapping.units();
    let unit_ids = units
        .iter()
        .enumerate()
        .map(|(id, unit)| Ok((unit.clone(), u32::try_from(id)?)))
        .collect::<Result<HashMap<_, _>>>()?;

    let mut groups = Vec::new();
    let mut unit_groups: HashMap<String, Vec<usize>> = HashMap::new();
    let mut cluster_reader = common_reader(cluster_path);
    let mut line = String::new();
    while cluster_reader.read_line(&mut line)? != 0 {
        let mut fields = line.split_whitespace();
        if let Some(group) = fields.next() {
            fields.next();
            let group_index = groups.len();
            groups.push(group.to_string());
            for unit in fields {
                let memberships = unit_groups.entry(unit.to_string()).or_default();
                if memberships.last() != Some(&group_index) {
                    memberships.push(group_index);
                }
            }
        }
        line.clear();
    }
    let mut writers = groups
        .iter()
        .map(|group| {
            ClmbWriter::create_synchronous(
                output_directory.join(format!("{group}.clmb")),
                &units,
                CLUSTER_CLMB_BLOCK_SIZE,
                None,
                None,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    stream_clm_records(
        input,
        |segment1, orientation1, segment2, orientation2, distances| {
            let (Some(mapped1), Some(mapped2)) = (
                mapping.map_clm_parts(segment1, orientation1, true),
                mapping.map_clm_parts(segment2, orientation2, false),
            ) else {
                return Ok(());
            };
            if mapped1.0 == mapped2.0 {
                return Ok(());
            }
            let (Some(groups1), Some(groups2)) =
                (unit_groups.get(mapped1.0), unit_groups.get(mapped2.0))
            else {
                return Ok(());
            };
            let shift = mapped1
                .2
                .checked_add(mapped2.2)
                .context("contracted CLM shift exceeds u64")?;
            let adjusted = adjusted_clm_distances(distances, shift)?;
            let endpoint1 =
                encode_endpoint(unit_ids[mapped1.0], if mapped1.1 == '+' { 0 } else { 1 })?;
            let endpoint2 =
                encode_endpoint(unit_ids[mapped2.0], if mapped2.1 == '+' { 0 } else { 1 })?;
            for group in groups1 {
                if groups2.binary_search(group).is_ok() {
                    writers[*group].write_record(endpoint1, endpoint2, adjusted.as_ref())?;
                }
            }
            Ok(())
        },
    )?;
    for writer in writers {
        writer.finish()?;
    }
    Ok(())
}

#[allow(dead_code)]
fn remap_clm_by_cluster(
    mapping: &MappingTable,
    cluster_path: &str,
    input: &str,
    output_directory: &str,
    temporary_directory: &Path,
    sort_buffer: &str,
    threads: usize,
) -> Result<()> {
    let workspace = tempdir_in(temporary_directory)
        .context("cannot create temporary directory for clustered remapped CLM")?;
    let output_directory = PathBuf::from(output_directory);
    std::fs::create_dir_all(&output_directory)?;

    let mut groups = Vec::new();
    let mut unit_groups: HashMap<String, Vec<usize>> = HashMap::new();
    let mut cluster_reader = common_reader(cluster_path);
    let mut line = String::new();
    while cluster_reader.read_line(&mut line)? != 0 {
        let mut fields = line.split_whitespace();
        if let Some(group) = fields.next() {
            fields.next();
            let group_index = groups.len();
            groups.push(group.to_string());
            for unit in fields {
                let memberships = unit_groups.entry(unit.to_string()).or_default();
                if memberships.last() != Some(&group_index) {
                    memberships.push(group_index);
                }
            }
        }
        line.clear();
    }

    let mapped_paths = (0..groups.len())
        .map(|index| workspace.path().join(format!("group-{}.mapped.tsv", index)))
        .collect::<Vec<_>>();
    let sorted_paths = (0..groups.len())
        .map(|index| workspace.path().join(format!("group-{}.sorted.tsv", index)))
        .collect::<Vec<_>>();
    let mut writers = mapped_paths
        .iter()
        .map(|path| File::create(path).map(|file| BufWriter::with_capacity(IO_BUFFER_SIZE, file)))
        .collect::<io::Result<Vec<_>>>()?;

    let mut reader = common_reader(input);
    let mut sequence = 0u64;
    while reader.read_line(&mut line)? != 0 {
        let record = line.trim_end_matches(&['\r', '\n'][..]);
        if record.is_empty() || record.starts_with('#') {
            line.clear();
            continue;
        }
        let mut fields = record.splitn(3, '\t');
        let pair = fields.next().context("missing CLM contig pair")?;
        fields.next().context("missing CLM count")?;
        let distances = fields.next().context("missing CLM distances")?;
        let mut tokens = pair.split_whitespace();
        let token1 = tokens.next().context("missing first CLM token")?;
        let token2 = tokens.next().context("missing second CLM token")?;
        if let (Some(mapped1), Some(mapped2)) = (
            mapping.map_clm_token(token1, true),
            mapping.map_clm_token(token2, false),
        ) {
            if mapped1.0 != mapped2.0 {
                let groups1 = unit_groups.get(mapped1.0);
                let groups2 = unit_groups.get(mapped2.0);
                if let (Some(groups1), Some(groups2)) = (groups1, groups2) {
                    let shift = mapped1.2 as u128 + mapped2.2 as u128;
                    for group in groups1 {
                        if groups2.binary_search(group).is_ok() {
                            writeln!(
                                writers[*group],
                                "{}{} {}{}\t{}\t{}\t{}",
                                mapped1.0,
                                mapped1.1,
                                mapped2.0,
                                mapped2.1,
                                sequence,
                                shift,
                                distances
                            )?;
                            sequence += 1;
                        }
                    }
                }
            }
        }
        line.clear();
    }
    for writer in &mut writers {
        writer.flush()?;
    }
    drop(writers);

    let nonempty = mapped_paths
        .iter()
        .enumerate()
        .filter_map(|(index, path)| {
            path.metadata()
                .ok()
                .filter(|metadata| metadata.len() > 0)
                .map(|_| index)
        })
        .collect::<Vec<_>>();
    let worker_count = threads.max(1).min(nonempty.len().max(1));
    let threads_per_group = (threads.max(1) / worker_count).max(1);
    let group_sort_buffer = divided_sort_buffer(sort_buffer, worker_count);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(worker_count)
        .build()?;
    pool.install(|| {
        nonempty.par_iter().try_for_each(|index| -> Result<()> {
            run_sort(
                &mapped_paths[*index],
                &sorted_paths[*index],
                workspace.path(),
                &["1,1", "2,2n"],
                &group_sort_buffer,
                threads_per_group,
            )?;
            let output = output_directory.join(format!("{}.clm", groups[*index]));
            process_clm_chunks(
                &sorted_paths[*index],
                output
                    .to_str()
                    .context("clustered CLM output path is not valid UTF-8")?,
                workspace.path(),
                threads_per_group,
            )
        })
    })?;

    for (index, group) in groups.iter().enumerate() {
        let output = output_directory.join(format!("{}.clm", group));
        if !output.exists() || mapped_paths[index].metadata()?.len() == 0 {
            File::create(output)?;
        }
    }
    Ok(())
}

fn remap_clm(
    mapping: &MappingTable,
    input: &str,
    output: &str,
    temporary_directory: &Path,
    sort_buffer: &str,
    threads: usize,
) -> Result<()> {
    let workspace = tempdir_in(temporary_directory)
        .context("cannot create temporary directory for remapped CLM")?;
    let mut mapped = NamedTempFile::new_in(workspace.path())?;
    let sorted = NamedTempFile::new_in(workspace.path())?;
    {
        let mut writer = BufWriter::with_capacity(IO_BUFFER_SIZE, mapped.as_file_mut());
        let mut reader = common_reader(input);
        let mut line = String::new();
        let mut sequence = 0u64;
        while reader.read_line(&mut line)? != 0 {
            let record = line.trim_end_matches(&['\r', '\n'][..]);
            if record.is_empty() || record.starts_with('#') {
                line.clear();
                continue;
            }
            let mut fields = record.splitn(3, '\t');
            let pair = fields.next().context("missing CLM contig pair")?;
            fields.next().context("missing CLM count")?;
            let distances = fields.next().context("missing CLM distances")?;
            let mut tokens = pair.split_whitespace();
            let token1 = tokens.next().context("missing first CLM token")?;
            let token2 = tokens.next().context("missing second CLM token")?;
            let mapped1 = mapping.map_clm_token(token1, true);
            let mapped2 = mapping.map_clm_token(token2, false);
            if let (Some(mapped1), Some(mapped2)) = (mapped1, mapped2) {
                if mapped1.0 != mapped2.0 {
                    let shift = mapped1.2 as u128 + mapped2.2 as u128;
                    writeln!(
                        writer,
                        "{}{} {}{}\t{}\t{}\t{}",
                        mapped1.0, mapped1.1, mapped2.0, mapped2.1, sequence, shift, distances
                    )?;
                    sequence += 1;
                }
            }
            line.clear();
        }
        writer.flush()?;
    }
    run_sort(
        mapped.path(),
        sorted.path(),
        workspace.path(),
        &["1,1", "2,2n"],
        sort_buffer,
        threads,
    )?;

    process_clm_chunks(sorted.path(), output, workspace.path(), threads)
}

pub fn contract_gfa_scaffolding_inputs(
    mapping_path: &str,
    contacts_input: &str,
    contacts_output: &str,
    clm_input: Option<&str>,
    clm_output: Option<&str>,
    cluster_path: Option<&str>,
    clm_output_directory: Option<&str>,
    temporary_directory: Option<&str>,
    sort_buffer: &str,
    threads: usize,
) -> Result<()> {
    let full_clm_mode = clm_input.is_some() && clm_output.is_some();
    let split_clm_mode =
        clm_input.is_some() && cluster_path.is_some() && clm_output_directory.is_some();
    if clm_input.is_some() && full_clm_mode == split_clm_mode {
        bail!("supply either --clm-output or both --clusters and --clm-output-dir");
    }
    let mapping = MappingTable::from_path(mapping_path)?;
    let temporary_root = temporary_root(temporary_directory)?;
    let threads = threads.max(1);
    let concurrent_contact_threads = (threads / 2).max(1);
    let concurrent_clm_threads = threads.saturating_sub(concurrent_contact_threads).max(1);
    if split_clm_mode {
        let input = clm_input.unwrap();
        let cluster_path = cluster_path.unwrap();
        let output_directory = clm_output_directory.unwrap();
        std::thread::scope(|scope| -> Result<()> {
            let contacts = scope.spawn(|| {
                remap_contacts(
                    &mapping,
                    contacts_input,
                    contacts_output,
                    concurrent_contact_threads,
                )
            });
            let clm = scope.spawn(|| {
                remap_clm_by_cluster_to_clmb(&mapping, cluster_path, input, output_directory)
            });
            contacts
                .join()
                .map_err(|_| anyhow::anyhow!("contacts remapping thread panicked"))??;
            clm.join()
                .map_err(|_| anyhow::anyhow!("clustered CLM remapping thread panicked"))??;
            Ok(())
        })?;
    } else if let (Some(input), Some(output)) = (clm_input, clm_output) {
        std::thread::scope(|scope| -> Result<()> {
            let contacts = scope.spawn(|| {
                remap_contacts(
                    &mapping,
                    contacts_input,
                    contacts_output,
                    concurrent_contact_threads,
                )
            });
            let clm = scope.spawn(|| {
                if output.ends_with(".clmb") {
                    remap_clm_to_clmb(&mapping, input, output)
                } else {
                    remap_clm(
                        &mapping,
                        input,
                        output,
                        &temporary_root,
                        sort_buffer,
                        concurrent_clm_threads,
                    )
                }
            });
            contacts
                .join()
                .map_err(|_| anyhow::anyhow!("contacts remapping thread panicked"))??;
            clm.join()
                .map_err(|_| anyhow::anyhow!("CLM remapping thread panicked"))??;
            Ok(())
        })?;
    } else {
        remap_contacts(&mapping, contacts_input, contacts_output, threads)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn borrows_unchanged_clm_distances() {
        let unchanged = [2, 10, 100];
        assert!(matches!(
            adjusted_clm_distances(&unchanged, 0).unwrap(),
            Cow::Borrowed(_)
        ));
        assert_eq!(adjusted_clm_distances(&[1, 2], 0).unwrap().as_ref(), [2, 2]);
        assert_eq!(
            adjusted_clm_distances(&[2, 10], 5).unwrap().as_ref(),
            [7, 15]
        );
    }

    #[test]
    fn aggregates_normalized_gfa_end_contacts() {
        let directory = tempfile::tempdir().unwrap();
        let segments = directory.path().join("segments.tsv");
        let contacts = directory.path().join("contacts.tsv");
        let links = directory.path().join("links.tsv");
        let output = directory.path().join("output.tsv");
        fs::write(&segments, "#Segment\tLength\nA\t100\nB\t200\n").unwrap();
        fs::write(
            &contacts,
            "A_1\tB_0\t2.5\nB_0\tA_1\t3.5\nA_0\tA_1\t10\nA_1\tmissing_0\t20\n",
        )
        .unwrap();
        fs::write(&links, "A\tR\tB\tL\n").unwrap();

        aggregate_gfa_end_contacts(
            segments.to_str().unwrap(),
            contacts.to_str().unwrap(),
            links.to_str().unwrap(),
            output.to_str().unwrap(),
            4,
        )
        .unwrap();

        assert_eq!(
            fs::read_to_string(output).unwrap(),
            "#Segment1\tEnd1\tSegment2\tEnd2\tCount\tScore\tCompetitor1\tCompetitor2\nA\tR\tB\tL\t6\t0.0012\t0\t0\n"
        );
    }

    #[test]
    fn remaps_contacts_and_updates_all_clm_distances() {
        let directory = tempfile::tempdir().unwrap();
        let mapping = directory.path().join("mapping.tsv");
        let contacts = directory.path().join("input.contacts");
        let clm = directory.path().join("input.clm");
        let output_contacts = directory.path().join("output.contacts");
        let output_clm = directory.path().join("output.clm");
        fs::write(
            &mapping,
            "#Segment\tUnit\tHalf0\tHalf1\tPlusOrientation\tPlusLeftOffset\tPlusRightOffset\tMinusOrientation\tMinusLeftOffset\tMinusRightOffset\n\
A\tblock\tblock_0\tblock_0\t+\t135\t100\t-\t100\t135\n\
B\tblock\tblock_1\tblock_1\t+\t0\t95\t-\t95\t0\n\
C\tC\tC_0\tC_1\t+\t0\t7\t-\t0\t11\n",
        )
        .unwrap();
        fs::write(&contacts, "A_0\tC_1\t2.5\nA_0\tC_1\t3.5\nA_1\tB_0\t20\n").unwrap();
        fs::write(
            &clm,
            "A+ C+\t2\t10 20\nA+ C-\t1\t10\nB- C+\t1\t30\nB- C-\t1\t30\nA+ B+\t1\t5\n",
        )
        .unwrap();

        contract_gfa_scaffolding_inputs(
            mapping.to_str().unwrap(),
            contacts.to_str().unwrap(),
            output_contacts.to_str().unwrap(),
            Some(clm.to_str().unwrap()),
            Some(output_clm.to_str().unwrap()),
            None,
            None,
            Some(directory.path().to_str().unwrap()),
            "1M",
            4,
        )
        .unwrap();

        assert_eq!(
            fs::read_to_string(output_contacts).unwrap(),
            "C_1\tblock_0\t6\n"
        );
        assert_eq!(
            fs::read_to_string(output_clm).unwrap(),
            "block+ C+\t2\t152 162\nblock+ C-\t1\t156\nblock- C+\t1\t132\nblock- C-\t1\t136\n"
        );
    }
}
