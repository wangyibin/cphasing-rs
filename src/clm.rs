use anyhow::{Context, Result as anyResult, bail};
// use indicatif::ProgressBar;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{Receiver as BlockReceiver, SyncSender, TryRecvError};

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

use crossbeam_channel::{Receiver, Sender, bounded, unbounded};
use indexmap::IndexMap;
use memmap2::Mmap;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};
use std::thread;

use crate::core::BaseTable;
use crate::core::{ContigPair2, ContigPair3};
use crate::core::{common_reader, common_writer};

pub const CLMB_MAGIC: [u8; 8] = *b"CPCLMB1\0";
pub const CLMB_FLAG_GZIP_BLOCKS: u32 = 1;
pub const CLMB_BLOCK_FLAG_U64_DISTANCES: u32 = 1;
pub const CLMB_DEFAULT_BLOCK_SIZE: usize = 8 * 1024 * 1024;
pub const CLMB_UNKNOWN_COUNT: u64 = u64::MAX;
const SPLIT_CLMB_BLOCK_SIZE: usize = 256 * 1024;
const SPLIT_CLM_WORKERS: usize = 8;

pub fn clm_output_prefix(output: &str) -> &str {
    if let Some(prefix) = output.strip_suffix(".clmb") {
        return prefix;
    }
    let prefix = output.strip_suffix(".gz").unwrap_or(output);
    prefix.strip_suffix(".clm").unwrap_or(prefix)
}

#[derive(Debug, Clone)]
pub struct Clm {
    file: String,
}

impl BaseTable for Clm {
    fn new(name: &String) -> Self {
        Clm { file: name.clone() }
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

impl Clm {
    // pub fn split_clm(&self, cluster_file: &String, output_dir: &String) -> anyResult<()> {
    //     let mut cluster_map: HashMap<String, Vec<String>> = HashMap::new();

    //     std::fs::create_dir_all(output_dir)?;

    //     let mut cluster_file = BufReader::new(common_reader(cluster_file));
    //     let mut line = String::new();
    //     while cluster_file.read_line(&mut line)? > 0 {
    //         let mut iter = line.split_whitespace();
    //         let cluster = iter.next().unwrap();
    //         let cluster_set = cluster_map.entry(cluster.to_string()).or_insert(Vec::new());
    //         for item in iter {
    //             cluster_set.push(item.to_string());
    //         }
    //         line.clear();
    //     }

    //     let mut cluster_paired_map: HashMap<ContigPair3, &String> = HashMap::new();
    //     for (cluster, contigs) in &cluster_map {
    //         for i in 0..contigs.len() {
    //             for j in i+1..contigs.len() {
    //                 let (a, b) = if contigs[i] <= contigs[j] {
    //                     (&contigs[i], &contigs[j])
    //                 } else {
    //                     (&contigs[j], &contigs[i])
    //                 };
    //                 let contig_pair = ContigPair3::new(a, b);

    //                 cluster_paired_map.insert(contig_pair, cluster);
    //             }
    //         }

    //     }

    //     let reader = common_reader(&self.file);

    //     let mut writer_map: HashMap<&String, Box<dyn Write + Send>> = HashMap::with_capacity(cluster_map.len());
    //     for (cluster, _) in &cluster_map {
    //         let file_name = format!("{}/{}.clm", output_dir, cluster);
    //         let writer = common_writer(&file_name);
    //         writer_map.insert(cluster, writer);
    //     }
    //     fn strip_orient(s: &str) -> &str {
    //         if let Some(last) = s.as_bytes().last() {
    //             if *last == b'+' || *last == b'-' {
    //                 return &s[..s.len() - 1];
    //             }
    //         }
    //         s
    //     }
    //     for line in reader.lines() {
    //         let line = line?;
    //         let mut iter = line.split_whitespace();
    //         let c1_raw = match iter.next() {
    //             Some(x) => x,
    //             None => continue,
    //         };
    //         let c2_raw = match iter.next() {
    //             Some(x) => x,
    //             None => continue,
    //         };

    //         let contig1 = strip_orient(c1_raw);
    //         let contig2 = strip_orient(c2_raw);
    //         let (a, b) = if contig1 <= contig2 { (contig1, contig2) } else { (contig2, contig1) };
    //         let contig_pair = ContigPair3::new(a, b);

    //         let cluster = match cluster_paired_map.get(&contig_pair) {
    //             Some(cluster) => cluster,
    //             None => continue
    //         };

    //         let writer = writer_map.get_mut(cluster).unwrap();
    //         writeln!(writer, "{}", line)?;

    //     }

    //     Ok(())
    // }
    pub fn split_clm(&self, cluster_file: &String, output_dir: &String) -> anyResult<()> {
        let mut cluster_map: HashMap<String, Vec<String>> = HashMap::new();
        std::fs::create_dir_all(output_dir)?;

        let mut cluster_reader = BufReader::new(common_reader(cluster_file));
        let mut line = String::new();
        while cluster_reader.read_line(&mut line)? > 0 {
            let mut iter = line.split_whitespace();
            if let Some(cluster) = iter.next() {
                let cluster_set = cluster_map.entry(cluster.to_string()).or_default();
                for item in iter {
                    cluster_set.push(item.to_string());
                }
            }
            line.clear();
        }

        if is_clmb_file(&self.file)? {
            return self.split_clmb(cluster_map, output_dir);
        }

        let mut contig_to_cids: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut writers = Vec::new();

        for (cluster, contigs) in cluster_map {
            let cid = writers.len();
            for c in contigs {
                contig_to_cids.entry(c).or_default().push(cid);
            }
            let file_path = format!("{}/{}.clm", output_dir, cluster);
            let writer = Arc::new(Mutex::new(BufWriter::new(common_writer(&file_path))));
            writers.push(writer);
        }

        let contig_to_cids = Arc::new(contig_to_cids);
        let writers = Arc::new(writers);
        let (tx, rx) = bounded::<Vec<String>>(200);

        let mut worker_handles = Vec::new();
        for _ in 0..SPLIT_CLM_WORKERS {
            let rx = rx.clone();
            let contig_to_cids = Arc::clone(&contig_to_cids);
            let writers = Arc::clone(&writers);

            worker_handles.push(thread::spawn(move || {
                #[inline]
                fn strip_orient(s: &str) -> &str {
                    let b = s.as_bytes();
                    if !b.is_empty() && (b[b.len() - 1] == b'+' || b[b.len() - 1] == b'-') {
                        &s[..s.len() - 1]
                    } else {
                        s
                    }
                }

                #[inline]
                fn write_to_common_clusters(
                    cids1: &[usize],
                    cids2: &[usize],
                    writers: &[Arc<Mutex<BufWriter<Box<dyn Write + Send>>>>],
                    line: &str,
                ) {
                    for &cid in cids1 {
                        if cids2.iter().any(|&x| x == cid) {
                            let mut w = writers[cid].lock().unwrap();
                            writeln!(w, "{}", line).unwrap();
                        }
                    }
                }

                for batch in rx {
                    for line in batch {
                        let mut parts = line.splitn(3, |c: char| c == '\t' || c == ' ');
                        let c1_raw = match parts.next() {
                            Some(x) => x,
                            None => continue,
                        };
                        let c2_raw = match parts.next() {
                            Some(x) => x,
                            None => continue,
                        };

                        let contig1 = strip_orient(c1_raw);
                        let contig2 = strip_orient(c2_raw);

                        let cids1 = contig_to_cids.get(contig1);
                        let cids2 = contig_to_cids.get(contig2);
                        if let (Some(c1s), Some(c2s)) = (cids1, cids2) {
                            write_to_common_clusters(c1s, c2s, &writers, &line);
                        }
                    }
                }
            }));
        }

        let reader = common_reader(&self.file);
        let mut batch = Vec::with_capacity(2000);
        for line in reader.lines() {
            batch.push(line?);
            if batch.len() >= 2000 {
                tx.send(batch).unwrap();
                batch = Vec::with_capacity(2000);
            }
        }
        if !batch.is_empty() {
            tx.send(batch).unwrap();
        }
        drop(tx);
        for h in worker_handles {
            h.join().unwrap();
        }

        for w in writers.iter() {
            w.lock().unwrap().flush()?;
        }

        Ok(())
    }

    fn split_clmb(
        &self,
        cluster_map: HashMap<String, Vec<String>>,
        output_dir: &str,
    ) -> anyResult<()> {
        let mut reader = ClmbReader::open(&self.file)?;
        let contigs = reader.header.contigs.clone();
        let mut contig_to_cids: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut writers = Vec::with_capacity(cluster_map.len());

        for (cluster, cluster_contigs) in cluster_map {
            let cid = writers.len();
            for contig in cluster_contigs {
                contig_to_cids.entry(contig).or_default().push(cid);
            }
            writers.push(Mutex::new(ClmbWriter::create_synchronous(
                Path::new(output_dir).join(format!("{cluster}.clmb")),
                &contigs,
                SPLIT_CLMB_BLOCK_SIZE,
                None,
                None,
            )?));
        }

        let contigs = Arc::new(contigs);
        let contig_to_cids = Arc::new(contig_to_cids);
        let writers = Arc::new(writers);
        let (tx, rx) = bounded::<Vec<ClmbRecord>>(SPLIT_CLM_WORKERS * 2);
        let mut worker_handles = Vec::with_capacity(SPLIT_CLM_WORKERS);
        for _ in 0..SPLIT_CLM_WORKERS {
            let rx = rx.clone();
            let contigs = Arc::clone(&contigs);
            let contig_to_cids = Arc::clone(&contig_to_cids);
            let writers = Arc::clone(&writers);
            worker_handles.push(thread::spawn(move || -> anyResult<()> {
                for block in rx {
                    for record in block {
                        let contig1 = &contigs[record.contig1() as usize];
                        let contig2 = &contigs[record.contig2() as usize];
                        let (Some(cids1), Some(cids2)) =
                            (contig_to_cids.get(contig1), contig_to_cids.get(contig2))
                        else {
                            continue;
                        };
                        for &cid in cids1 {
                            if cids2.contains(&cid) {
                                writers[cid]
                                    .lock()
                                    .map_err(|_| {
                                        anyhow::anyhow!("CLMB split writer lock poisoned")
                                    })?
                                    .write_record(
                                        record.endpoint1,
                                        record.endpoint2,
                                        &record.distances,
                                    )?;
                            }
                        }
                    }
                }
                Ok(())
            }));
        }

        let read_result = (|| -> anyResult<()> {
            while let Some(block) = reader.next_block()? {
                tx.send(block)
                    .map_err(|_| anyhow::anyhow!("CLMB split workers stopped unexpectedly"))?;
            }
            Ok(())
        })();
        drop(tx);

        let mut worker_error = None;
        for handle in worker_handles {
            match handle.join() {
                Ok(Ok(())) => {}
                Ok(Err(error)) if worker_error.is_none() => worker_error = Some(error),
                Ok(Err(_)) => {}
                Err(_) if worker_error.is_none() => {
                    worker_error = Some(anyhow::anyhow!("CLMB split worker panicked"));
                }
                Err(_) => {}
            }
        }
        if let Some(error) = worker_error {
            return Err(error);
        }
        read_result?;

        let writers = Arc::try_unwrap(writers)
            .map_err(|_| anyhow::anyhow!("CLMB split writers are still in use"))?;
        for writer in writers {
            writer
                .into_inner()
                .map_err(|_| anyhow::anyhow!("CLMB split writer lock poisoned"))?
                .finish()?;
        }
        Ok(())
    }

    pub fn parse_clm(&self) -> anyResult<IndexMap<String, Vec<u32>>> {
        let reader = common_reader(&self.file);
        let mut data: IndexMap<String, Vec<u32>> = IndexMap::new();

        for line in reader.lines() {
            let line = line?;
            let mut iter = line.split('\t');
            let contigs = iter.next().unwrap();
            let _ = iter.next().unwrap();

            let values = iter.next().unwrap();
            let values = values
                .split(' ')
                .map(|x| x.parse::<u32>().unwrap())
                .collect::<Vec<u32>>();
            data.entry(contigs.to_string())
                .or_insert_with(Vec::new)
                .extend(values);
        }

        Ok(data)
    }
}

pub fn merge_clm(clm_files: Vec<String>, output: &String) -> anyResult<()> {
    let mut wtr = BufWriter::new(common_writer(output));

    let (sender, receiver) = bounded(100);
    let mut data: IndexMap<String, Vec<u32>> = IndexMap::new();

    let producer_handle = thread::spawn(move || {
        for clm_file in clm_files {
            let reader = common_reader(&clm_file);
            for line in reader.lines() {
                let line = line.unwrap();
                sender.send(line).unwrap();
            }
        }
    });

    let consumer_handles: Vec<_> = (0..8)
        .map(|_| {
            let receiver = receiver.clone();
            thread::spawn(move || {
                let mut local_data: IndexMap<String, Vec<u32>> = IndexMap::new();
                for line in receiver.iter() {
                    let mut iter = line.split('\t');
                    let contigs = iter.next().unwrap();
                    let _ = iter.next().unwrap();

                    let values = iter.next().unwrap();
                    let values = values
                        .split(' ')
                        .map(|x| x.parse::<u32>().unwrap())
                        .collect::<Vec<u32>>();
                    local_data
                        .entry(contigs.to_string())
                        .or_insert_with(Vec::new)
                        .extend(values);
                }
                local_data
            })
        })
        .collect();

    producer_handle.join().unwrap();

    for handle in consumer_handles {
        let local_data = handle.join().unwrap();
        for (contigs, values) in local_data {
            data.entry(contigs).or_insert_with(Vec::new).extend(values);
        }
    }

    let (sender, receiver) = bounded(100);

    let producer_handle = thread::spawn(move || {
        for (contigs, values) in data {
            sender.send((contigs, values)).unwrap();
        }
    });

    let consumer_handles: Vec<_> = (0..8)
        .map(|_| {
            let receiver = receiver.clone();
            thread::spawn(move || {
                let mut local_buffer = Vec::new();
                for (contigs, values) in receiver.iter() {
                    let count = values.len();
                    let values_str = values
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<String>>()
                        .join(" ");
                    local_buffer.push(format!("{}\t{}\t{}", contigs, count, values_str));
                }
                local_buffer
            })
        })
        .collect();

    producer_handle.join().unwrap();

    for handle in consumer_handles {
        let local_buffer = handle.join().unwrap();
        for line in local_buffer {
            writeln!(wtr, "{}", line).unwrap();
        }
    }

    wtr.flush().unwrap();

    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClmbRecord {
    pub endpoint1: u32,
    pub endpoint2: u32,
    pub distances: Vec<u64>,
}

impl ClmbRecord {
    #[inline]
    pub fn contig1(&self) -> u32 {
        self.endpoint1 >> 1
    }

    #[inline]
    pub fn contig2(&self) -> u32 {
        self.endpoint2 >> 1
    }

    #[inline]
    pub fn orientation1(&self) -> u8 {
        (self.endpoint1 & 1) as u8
    }

    #[inline]
    pub fn orientation2(&self) -> u8 {
        (self.endpoint2 & 1) as u8
    }
}

#[derive(Debug, Clone)]
pub struct ClmbHeader {
    pub flags: u32,
    pub target_block_size: u32,
    pub record_count: u64,
    pub distance_count: u64,
    pub contigs: Vec<String>,
}

pub struct ClmbWriter {
    path: PathBuf,
    writer: BufWriter<File>,
    contig_count: u32,
    target_block_size: usize,
    raw_block: Vec<u8>,
    block_records: u32,
    block_distances: u32,
    block_flags: u32,
    observed_records: u64,
    observed_distances: u64,
    declared_records: u64,
    declared_distances: u64,
    compression_jobs: Option<SyncSender<PendingClmbBlock>>,
    compressed_blocks: BlockReceiver<Result<CompressedClmbBlock, String>>,
    compression_workers: Vec<thread::JoinHandle<()>>,
    pending_compressed: BTreeMap<u64, CompressedClmbBlock>,
    submitted_blocks: u64,
    written_blocks: u64,
    max_in_flight: u64,
    compression_level: u32,
}

pub struct ClmbBatchWriter<'a> {
    writer: &'a Mutex<ClmbWriter>,
    records: Vec<ClmbRecord>,
    encoded_bytes: usize,
    target_bytes: usize,
}

impl<'a> ClmbBatchWriter<'a> {
    pub fn new(writer: &'a Mutex<ClmbWriter>, target_bytes: usize) -> Self {
        Self {
            writer,
            records: Vec::new(),
            encoded_bytes: 0,
            target_bytes: target_bytes.max(1),
        }
    }

    pub fn write_record(
        &mut self,
        endpoint1: u32,
        endpoint2: u32,
        distances: Vec<u64>,
    ) -> anyResult<()> {
        let distance_width = if distances.iter().any(|&value| value > u32::MAX as u64) {
            8
        } else {
            4
        };
        self.encoded_bytes = self
            .encoded_bytes
            .saturating_add(12usize.saturating_add(distances.len().saturating_mul(distance_width)));
        self.records.push(ClmbRecord {
            endpoint1,
            endpoint2,
            distances,
        });
        if self.encoded_bytes >= self.target_bytes {
            self.flush()?;
        }
        Ok(())
    }

    pub fn flush(&mut self) -> anyResult<()> {
        if self.records.is_empty() {
            return Ok(());
        }
        let mut writer = self.writer.lock().unwrap();
        for record in self.records.drain(..) {
            writer.write_record(record.endpoint1, record.endpoint2, &record.distances)?;
        }
        self.encoded_bytes = 0;
        Ok(())
    }
}

impl Drop for ClmbBatchWriter<'_> {
    fn drop(&mut self) {
        self.flush().expect("failed to flush CLMB record batch");
    }
}

struct PendingClmbBlock {
    id: u64,
    raw: Vec<u8>,
    record_count: u32,
    distance_count: u32,
    flags: u32,
}

struct CompressedClmbBlock {
    id: u64,
    raw_length: u32,
    compressed: Vec<u8>,
    record_count: u32,
    distance_count: u32,
    flags: u32,
}

impl ClmbWriter {
    pub fn create(
        path: impl AsRef<Path>,
        contigs: &[String],
        target_block_size: usize,
        record_count: Option<u64>,
        distance_count: Option<u64>,
    ) -> anyResult<Self> {
        let compression_threads = std::env::var("CPHASING_IO_THREADS")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or_else(|| rayon::current_num_threads().min(16))
            .max(1);
        Self::create_configured(
            path,
            contigs,
            target_block_size,
            record_count,
            distance_count,
            compression_threads,
        )
    }

    pub fn create_synchronous(
        path: impl AsRef<Path>,
        contigs: &[String],
        target_block_size: usize,
        record_count: Option<u64>,
        distance_count: Option<u64>,
    ) -> anyResult<Self> {
        Self::create_configured(
            path,
            contigs,
            target_block_size,
            record_count,
            distance_count,
            0,
        )
    }

    fn create_configured(
        path: impl AsRef<Path>,
        contigs: &[String],
        target_block_size: usize,
        record_count: Option<u64>,
        distance_count: Option<u64>,
        compression_threads: usize,
    ) -> anyResult<Self> {
        if target_block_size == 0 || target_block_size > u32::MAX as usize {
            bail!(
                "CLMB target block size must be between 1 and {} bytes",
                u32::MAX
            );
        }
        let contig_count = u32::try_from(contigs.len()).context("too many CLMB contigs")?;
        let path = path.as_ref().to_path_buf();
        let file = File::create(&path)
            .with_context(|| format!("cannot create CLMB file `{}`", path.display()))?;
        let mut writer =
            BufWriter::with_capacity(target_block_size.min(CLMB_DEFAULT_BLOCK_SIZE), file);
        writer.write_all(&CLMB_MAGIC)?;
        write_u32(&mut writer, CLMB_FLAG_GZIP_BLOCKS)?;
        write_u32(&mut writer, target_block_size as u32)?;
        write_u64(&mut writer, record_count.unwrap_or(CLMB_UNKNOWN_COUNT))?;
        write_u64(&mut writer, distance_count.unwrap_or(CLMB_UNKNOWN_COUNT))?;
        write_u32(&mut writer, contig_count)?;
        let mut unique = HashSet::with_capacity(contigs.len());
        for contig in contigs {
            if !unique.insert(contig.as_str()) {
                bail!("duplicate contig `{}` in CLMB dictionary", contig);
            }
            let bytes = contig.as_bytes();
            write_u32(
                &mut writer,
                u32::try_from(bytes.len()).context("CLMB contig name is too long")?,
            )?;
            writer.write_all(bytes)?;
        }
        let compression_level = std::env::var("CPHASING_CLMB_COMPRESSION_LEVEL")
            .ok()
            .and_then(|value| value.parse::<u32>().ok())
            .filter(|level| *level <= 9)
            .unwrap_or(0);
        if compression_threads == 0 {
            log::debug!(
                "CLMB compression is synchronous at level {}.",
                compression_level
            );
        } else {
            log::info!(
                "CLMB compression uses {} workers at level {}.",
                compression_threads,
                compression_level
            );
        }
        let max_in_flight = compression_threads.saturating_mul(2).max(2);
        let (job_sender, job_receiver) = std::sync::mpsc::sync_channel(max_in_flight);
        let (result_sender, result_receiver) = std::sync::mpsc::sync_channel(max_in_flight);
        let mut compression_workers = Vec::with_capacity(compression_threads);
        let shared_jobs = Arc::new(Mutex::new(job_receiver));
        for worker_id in 0..compression_threads {
            let jobs = Arc::clone(&shared_jobs);
            let results = result_sender.clone();
            compression_workers.push(
                thread::Builder::new()
                    .name(format!("clmb-gzip-{worker_id}"))
                    .spawn(move || {
                        loop {
                            let job = {
                                let receiver = jobs.lock().unwrap();
                                receiver.recv()
                            };
                            let Ok(job) = job else {
                                break;
                            };
                            if results
                                .send(compress_clmb_block(job, compression_level))
                                .is_err()
                            {
                                break;
                            }
                        }
                    })?,
            );
        }
        drop(result_sender);
        Ok(Self {
            path,
            writer,
            contig_count,
            target_block_size,
            raw_block: Vec::with_capacity(target_block_size),
            block_records: 0,
            block_distances: 0,
            block_flags: 0,
            observed_records: 0,
            observed_distances: 0,
            declared_records: record_count.unwrap_or(CLMB_UNKNOWN_COUNT),
            declared_distances: distance_count.unwrap_or(CLMB_UNKNOWN_COUNT),
            compression_jobs: Some(job_sender),
            compressed_blocks: result_receiver,
            compression_workers,
            pending_compressed: BTreeMap::new(),
            submitted_blocks: 0,
            written_blocks: 0,
            max_in_flight: max_in_flight as u64,
            compression_level,
        })
    }

    pub fn create_with_compression_threads(
        path: impl AsRef<Path>,
        contigs: &[String],
        target_block_size: usize,
        record_count: Option<u64>,
        distance_count: Option<u64>,
        compression_threads: usize,
    ) -> anyResult<Self> {
        Self::create_configured(
            path,
            contigs,
            target_block_size,
            record_count,
            distance_count,
            compression_threads,
        )
    }

    pub fn write_record(
        &mut self,
        endpoint1: u32,
        endpoint2: u32,
        distances: &[u64],
    ) -> anyResult<()> {
        if endpoint1 >> 1 >= self.contig_count || endpoint2 >> 1 >= self.contig_count {
            bail!("CLMB endpoint references a contig outside the dictionary");
        }
        let distance_count = u32::try_from(distances.len())
            .context("one CLMB record contains more than u32::MAX distances")?;
        let block_flags = if distances.iter().any(|&value| value > u32::MAX as u64) {
            CLMB_BLOCK_FLAG_U64_DISTANCES
        } else {
            0
        };
        if !self.raw_block.is_empty() && self.block_flags != block_flags {
            self.flush_block()?;
        }
        self.block_flags = block_flags;
        let distance_width = if block_flags == CLMB_BLOCK_FLAG_U64_DISTANCES {
            8
        } else {
            4
        };
        let record_size = 12usize
            .checked_add(
                distances
                    .len()
                    .checked_mul(distance_width)
                    .context("CLMB record is too large")?,
            )
            .context("CLMB record is too large")?;
        if !self.raw_block.is_empty()
            && self.raw_block.len().saturating_add(record_size) > self.target_block_size
        {
            self.flush_block()?;
        }
        append_u32(&mut self.raw_block, endpoint1);
        append_u32(&mut self.raw_block, endpoint2);
        append_u32(&mut self.raw_block, distance_count);
        self.raw_block.reserve(distances.len() * distance_width);
        for &distance in distances {
            if block_flags == CLMB_BLOCK_FLAG_U64_DISTANCES {
                append_u64(&mut self.raw_block, distance);
            } else {
                append_u32(&mut self.raw_block, distance as u32);
            }
        }
        self.block_records = self
            .block_records
            .checked_add(1)
            .context("too many records in one CLMB block")?;
        self.block_distances = self
            .block_distances
            .checked_add(distance_count)
            .context("too many distances in one CLMB block")?;
        self.observed_records += 1;
        self.observed_distances = self
            .observed_distances
            .checked_add(distance_count as u64)
            .context("too many distances in CLMB output")?;
        Ok(())
    }

    fn flush_block(&mut self) -> anyResult<()> {
        if self.block_records == 0 {
            return Ok(());
        }
        while self.submitted_blocks - self.written_blocks >= self.max_in_flight {
            self.receive_compressed_block(true)?;
        }
        let raw = std::mem::replace(
            &mut self.raw_block,
            Vec::with_capacity(self.target_block_size),
        );
        let job = PendingClmbBlock {
            id: self.submitted_blocks,
            raw,
            record_count: self.block_records,
            distance_count: self.block_distances,
            flags: self.block_flags,
        };
        self.submitted_blocks += 1;
        if self.compression_workers.is_empty() {
            let block =
                compress_clmb_block(job, self.compression_level).map_err(anyhow::Error::msg)?;
            self.accept_compressed_block(block)?;
            self.block_records = 0;
            self.block_distances = 0;
            self.block_flags = 0;
            return Ok(());
        }
        self.compression_jobs
            .as_ref()
            .context("CLMB compression workers already stopped")?
            .send(job)
            .map_err(|_| anyhow::anyhow!("CLMB compression workers stopped unexpectedly"))?;
        self.block_records = 0;
        self.block_distances = 0;
        self.block_flags = 0;
        while self.receive_compressed_block(false)? {}
        Ok(())
    }

    fn receive_compressed_block(&mut self, blocking: bool) -> anyResult<bool> {
        let result =
            if blocking {
                Some(self.compressed_blocks.recv().map_err(|_| {
                    anyhow::anyhow!("CLMB compression workers stopped unexpectedly")
                })?)
            } else {
                match self.compressed_blocks.try_recv() {
                    Ok(result) => Some(result),
                    Err(TryRecvError::Empty) => None,
                    Err(TryRecvError::Disconnected) => None,
                }
            };
        let Some(result) = result else {
            return Ok(false);
        };
        let block = result.map_err(anyhow::Error::msg)?;
        self.accept_compressed_block(block)?;
        Ok(true)
    }

    fn accept_compressed_block(&mut self, block: CompressedClmbBlock) -> anyResult<()> {
        self.pending_compressed.insert(block.id, block);
        while let Some(block) = self.pending_compressed.remove(&self.written_blocks) {
            write_u32(&mut self.writer, block.raw_length)?;
            write_u32(
                &mut self.writer,
                u32::try_from(block.compressed.len())
                    .context("compressed CLMB block is too large")?,
            )?;
            write_u32(&mut self.writer, block.record_count)?;
            write_u32(&mut self.writer, block.distance_count)?;
            write_u32(&mut self.writer, block.flags)?;
            self.writer.write_all(&block.compressed)?;
            self.written_blocks += 1;
        }
        Ok(())
    }

    pub fn finish(mut self) -> anyResult<()> {
        self.flush_block()?;
        drop(self.compression_jobs.take());
        while self.written_blocks < self.submitted_blocks {
            self.receive_compressed_block(true)?;
        }
        for worker in self.compression_workers.drain(..) {
            if worker.join().is_err() {
                bail!("CLMB compression worker panicked");
            }
        }
        if self.declared_records != CLMB_UNKNOWN_COUNT
            && self.declared_records != self.observed_records
        {
            bail!(
                "CLMB record count mismatch while writing `{}`: declared {}, wrote {}",
                self.path.display(),
                self.declared_records,
                self.observed_records
            );
        }
        if self.declared_distances != CLMB_UNKNOWN_COUNT
            && self.declared_distances != self.observed_distances
        {
            bail!(
                "CLMB distance count mismatch while writing `{}`: declared {}, wrote {}",
                self.path.display(),
                self.declared_distances,
                self.observed_distances
            );
        }
        self.writer.flush()?;
        Ok(())
    }
}

fn compress_clmb_block(
    job: PendingClmbBlock,
    compression_level: u32,
) -> Result<CompressedClmbBlock, String> {
    let raw_length = u32::try_from(job.raw.len())
        .map_err(|_| "CLMB block exceeds u32::MAX bytes".to_string())?;
    let mut encoder = GzEncoder::new(Vec::new(), Compression::new(compression_level));
    encoder
        .write_all(&job.raw)
        .map_err(|error| error.to_string())?;
    let compressed = encoder.finish().map_err(|error| error.to_string())?;
    Ok(CompressedClmbBlock {
        id: job.id,
        raw_length,
        compressed,
        record_count: job.record_count,
        distance_count: job.distance_count,
        flags: job.flags,
    })
}

pub struct ClmbReader {
    path: PathBuf,
    reader: BufReader<File>,
    pub header: ClmbHeader,
    observed_records: u64,
    observed_distances: u64,
    finished: bool,
}

impl ClmbReader {
    pub fn open(path: impl AsRef<Path>) -> anyResult<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path)
            .with_context(|| format!("cannot open CLMB file `{}`", path.display()))?;
        let mut reader = BufReader::with_capacity(CLMB_DEFAULT_BLOCK_SIZE, file);
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if magic != CLMB_MAGIC {
            bail!("`{}` is not a CLMB v1 file", path.display());
        }
        let flags = read_u32(&mut reader)?;
        if flags != CLMB_FLAG_GZIP_BLOCKS {
            bail!("unsupported CLMB flags {flags:#x} in `{}`", path.display());
        }
        let target_block_size = read_u32(&mut reader)?;
        if target_block_size == 0 {
            bail!(
                "invalid zero CLMB target block size in `{}`",
                path.display()
            );
        }
        let record_count = read_u64(&mut reader)?;
        let distance_count = read_u64(&mut reader)?;
        let contig_count = read_u32(&mut reader)?;
        let mut contigs = Vec::with_capacity(contig_count as usize);
        let mut unique = HashSet::with_capacity(contig_count as usize);
        for _ in 0..contig_count {
            let length = read_u32(&mut reader)? as usize;
            if length > 16 * 1024 * 1024 {
                bail!("unreasonable CLMB contig-name length {length}");
            }
            let mut bytes = vec![0u8; length];
            reader.read_exact(&mut bytes)?;
            let contig = String::from_utf8(bytes).context("CLMB contig name is not UTF-8")?;
            if !unique.insert(contig.clone()) {
                bail!("duplicate contig `{contig}` in CLMB dictionary");
            }
            contigs.push(contig);
        }
        Ok(Self {
            path,
            reader,
            header: ClmbHeader {
                flags,
                target_block_size,
                record_count,
                distance_count,
                contigs,
            },
            observed_records: 0,
            observed_distances: 0,
            finished: false,
        })
    }

    pub fn next_block(&mut self) -> anyResult<Option<Vec<ClmbRecord>>> {
        if self.finished {
            return Ok(None);
        }
        let Some(raw_length) = read_u32_or_eof(&mut self.reader)? else {
            self.finished = true;
            self.validate_totals()?;
            return Ok(None);
        };
        let compressed_length = read_u32(&mut self.reader)?;
        let record_count = read_u32(&mut self.reader)?;
        let distance_count = read_u32(&mut self.reader)?;
        let block_flags = read_u32(&mut self.reader)?;
        if block_flags & !CLMB_BLOCK_FLAG_U64_DISTANCES != 0 {
            bail!(
                "unsupported CLMB block flags {block_flags:#x} in `{}`",
                self.path.display()
            );
        }
        let distance_width = if block_flags & CLMB_BLOCK_FLAG_U64_DISTANCES != 0 {
            8usize
        } else {
            4usize
        };
        let mut compressed = vec![0u8; compressed_length as usize];
        self.reader.read_exact(&mut compressed)?;
        let mut decoder = GzDecoder::new(compressed.as_slice());
        let mut raw = Vec::with_capacity(raw_length as usize);
        decoder.read_to_end(&mut raw)?;
        if raw.len() != raw_length as usize {
            bail!(
                "CLMB block length mismatch in `{}`: expected {}, decoded {}",
                self.path.display(),
                raw_length,
                raw.len()
            );
        }
        let mut records = Vec::with_capacity(record_count as usize);
        let mut offset = 0usize;
        let mut observed_distances = 0u32;
        for _ in 0..record_count {
            if raw.len().saturating_sub(offset) < 12 {
                bail!("truncated CLMB record header in `{}`", self.path.display());
            }
            let endpoint1 = take_u32(&raw, &mut offset)?;
            let endpoint2 = take_u32(&raw, &mut offset)?;
            if endpoint1 >> 1 >= self.header.contigs.len() as u32
                || endpoint2 >> 1 >= self.header.contigs.len() as u32
            {
                bail!(
                    "CLMB endpoint outside dictionary in `{}`",
                    self.path.display()
                );
            }
            let count = take_u32(&raw, &mut offset)?;
            let required = (count as usize)
                .checked_mul(distance_width)
                .context("CLMB distance array is too large")?;
            if raw.len().saturating_sub(offset) < required {
                bail!("truncated CLMB distance array in `{}`", self.path.display());
            }
            let mut distances = Vec::with_capacity(count as usize);
            for _ in 0..count {
                if distance_width == 8 {
                    distances.push(take_u64(&raw, &mut offset)?);
                } else {
                    distances.push(take_u32(&raw, &mut offset)? as u64);
                }
            }
            observed_distances = observed_distances
                .checked_add(count)
                .context("too many distances in CLMB block")?;
            records.push(ClmbRecord {
                endpoint1,
                endpoint2,
                distances,
            });
        }
        if offset != raw.len() {
            bail!("trailing bytes in CLMB block in `{}`", self.path.display());
        }
        if observed_distances != distance_count {
            bail!(
                "CLMB block distance count mismatch in `{}`: declared {}, decoded {}",
                self.path.display(),
                distance_count,
                observed_distances
            );
        }
        self.observed_records += record_count as u64;
        self.observed_distances += observed_distances as u64;
        Ok(Some(records))
    }

    fn validate_totals(&self) -> anyResult<()> {
        if self.header.record_count != CLMB_UNKNOWN_COUNT
            && self.header.record_count != self.observed_records
        {
            bail!(
                "CLMB record total mismatch in `{}`: declared {}, decoded {}",
                self.path.display(),
                self.header.record_count,
                self.observed_records
            );
        }
        if self.header.distance_count != CLMB_UNKNOWN_COUNT
            && self.header.distance_count != self.observed_distances
        {
            bail!(
                "CLMB distance total mismatch in `{}`: declared {}, decoded {}",
                self.path.display(),
                self.header.distance_count,
                self.observed_distances
            );
        }
        Ok(())
    }
}

pub fn is_clmb_file(path: impl AsRef<Path>) -> anyResult<bool> {
    let mut file = File::open(path.as_ref())?;
    let mut magic = [0u8; 8];
    match file.read_exact(&mut magic) {
        Ok(()) => Ok(magic == CLMB_MAGIC),
        Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(error) => Err(error.into()),
    }
}

#[derive(Debug)]
struct TextClmRecord {
    contig1: String,
    orientation1: u8,
    contig2: String,
    orientation2: u8,
    distances: Vec<u64>,
}

fn parse_text_clm_record(line: &str) -> anyResult<Option<TextClmRecord>> {
    let line = line.trim_end_matches(['\r', '\n']);
    if line.is_empty() || line.starts_with('#') {
        return Ok(None);
    }
    let mut fields = line.splitn(3, '\t');
    let pair = fields.next().context("missing CLM contig pair")?;
    let declared_count = fields
        .next()
        .context("missing CLM distance count")?
        .parse::<usize>()
        .context("invalid CLM distance count")?;
    let distance_field = fields.next().context("missing CLM distances")?;
    let mut endpoints = pair.split_whitespace();
    let endpoint1 = endpoints.next().context("missing first CLM endpoint")?;
    let endpoint2 = endpoints.next().context("missing second CLM endpoint")?;
    if endpoints.next().is_some() {
        bail!("too many CLM endpoints in `{pair}`");
    }
    let (contig1, orientation1) = parse_text_endpoint(endpoint1)?;
    let (contig2, orientation2) = parse_text_endpoint(endpoint2)?;
    let distances = distance_field
        .split_whitespace()
        .map(|value| {
            value
                .parse::<u64>()
                .with_context(|| format!("invalid CLM distance `{value}`"))
        })
        .collect::<anyResult<Vec<_>>>()?;
    if distances.len() != declared_count {
        bail!(
            "CLM count mismatch for `{pair}`: declared {declared_count}, observed {}",
            distances.len()
        );
    }
    Ok(Some(TextClmRecord {
        contig1: contig1.to_string(),
        orientation1,
        contig2: contig2.to_string(),
        orientation2,
        distances,
    }))
}

fn parse_text_endpoint(endpoint: &str) -> anyResult<(&str, u8)> {
    if let Some(contig) = endpoint.strip_suffix('+') {
        Ok((contig, 0))
    } else if let Some(contig) = endpoint.strip_suffix('-') {
        Ok((contig, 1))
    } else {
        bail!("CLM endpoint `{endpoint}` has no + or - orientation")
    }
}

fn write_text_clm_record(
    writer: &mut dyn Write,
    contig1: &str,
    orientation1: u8,
    contig2: &str,
    orientation2: u8,
    distances: &[u64],
) -> anyResult<()> {
    let orientation1 = if orientation1 == 0 { '+' } else { '-' };
    let orientation2 = if orientation2 == 0 { '+' } else { '-' };
    write!(
        writer,
        "{contig1}{orientation1} {contig2}{orientation2}\t{}\t",
        distances.len()
    )?;
    for (index, distance) in distances.iter().enumerate() {
        if index > 0 {
            writer.write_all(b" ")?;
        }
        write!(writer, "{distance}")?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn convert_clm(input: &str, output: &str, block_size: usize) -> anyResult<()> {
    let input_is_clmb = is_clmb_file(input)?;
    let output_is_clmb = output.ends_with(".clmb");
    if output_is_clmb {
        if input_is_clmb {
            let mut reader = ClmbReader::open(input)?;
            let mut writer = ClmbWriter::create(
                output,
                &reader.header.contigs,
                block_size,
                (reader.header.record_count != CLMB_UNKNOWN_COUNT)
                    .then_some(reader.header.record_count),
                (reader.header.distance_count != CLMB_UNKNOWN_COUNT)
                    .then_some(reader.header.distance_count),
            )?;
            while let Some(block) = reader.next_block()? {
                for record in block {
                    writer.write_record(record.endpoint1, record.endpoint2, &record.distances)?;
                }
            }
            return writer.finish();
        }

        let mut contigs = Vec::new();
        let mut contig_ids: FxHashMap<String, u32> = FxHashMap::default();
        let mut record_count = 0u64;
        let mut distance_count = 0u64;
        let mut reader = common_reader(input);
        let mut line = String::new();
        while reader.read_line(&mut line)? != 0 {
            if let Some(record) = parse_text_clm_record(&line)? {
                for contig in [&record.contig1, &record.contig2] {
                    if !contig_ids.contains_key(contig) {
                        let id = u32::try_from(contigs.len()).context("too many CLM contigs")?;
                        contig_ids.insert(contig.clone(), id);
                        contigs.push(contig.clone());
                    }
                }
                record_count += 1;
                distance_count = distance_count
                    .checked_add(record.distances.len() as u64)
                    .context("too many CLM distances")?;
            }
            line.clear();
        }

        let mut writer = ClmbWriter::create(
            output,
            &contigs,
            block_size,
            Some(record_count),
            Some(distance_count),
        )?;
        let mut reader = common_reader(input);
        while reader.read_line(&mut line)? != 0 {
            if let Some(record) = parse_text_clm_record(&line)? {
                let endpoint1 = encode_endpoint(contig_ids[&record.contig1], record.orientation1)?;
                let endpoint2 = encode_endpoint(contig_ids[&record.contig2], record.orientation2)?;
                writer.write_record(endpoint1, endpoint2, &record.distances)?;
            }
            line.clear();
        }
        writer.finish()
    } else {
        let mut writer = common_writer(output);
        if input_is_clmb {
            let mut reader = ClmbReader::open(input)?;
            let contigs = reader.header.contigs.clone();
            while let Some(block) = reader.next_block()? {
                for record in block {
                    write_text_clm_record(
                        writer.as_mut(),
                        &contigs[record.contig1() as usize],
                        record.orientation1(),
                        &contigs[record.contig2() as usize],
                        record.orientation2(),
                        &record.distances,
                    )?;
                }
            }
        } else {
            let mut reader = common_reader(input);
            std::io::copy(&mut reader, &mut writer)?;
        }
        writer.flush()?;
        Ok(())
    }
}

#[inline]
pub fn encode_endpoint(contig_id: u32, orientation: u8) -> anyResult<u32> {
    if orientation > 1 {
        bail!("CLMB orientation must be 0 (+) or 1 (-)");
    }
    contig_id
        .checked_mul(2)
        .and_then(|value| value.checked_add(orientation as u32))
        .context("CLMB contig ID is too large")
}

fn append_u32(buffer: &mut Vec<u8>, value: u32) {
    buffer.extend_from_slice(&value.to_le_bytes());
}

fn append_u64(buffer: &mut Vec<u8>, value: u64) {
    buffer.extend_from_slice(&value.to_le_bytes());
}

fn write_u32(writer: &mut impl Write, value: u32) -> std::io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn write_u64(writer: &mut impl Write, value: u64) -> std::io::Result<()> {
    writer.write_all(&value.to_le_bytes())
}

fn read_u32(reader: &mut impl Read) -> std::io::Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64(reader: &mut impl Read) -> std::io::Result<u64> {
    let mut bytes = [0u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_u32_or_eof(reader: &mut impl Read) -> std::io::Result<Option<u32>> {
    let mut bytes = [0u8; 4];
    let read = reader.read(&mut bytes[..1])?;
    if read == 0 {
        return Ok(None);
    }
    reader.read_exact(&mut bytes[1..])?;
    Ok(Some(u32::from_le_bytes(bytes)))
}

fn take_u32(bytes: &[u8], offset: &mut usize) -> anyResult<u32> {
    if bytes.len().saturating_sub(*offset) < 4 {
        bail!("truncated u32 in CLMB payload");
    }
    let value = u32::from_le_bytes(bytes[*offset..*offset + 4].try_into().unwrap());
    *offset += 4;
    Ok(value)
}

fn take_u64(bytes: &[u8], offset: &mut usize) -> anyResult<u64> {
    if bytes.len().saturating_sub(*offset) < 8 {
        bail!("truncated u64 in CLMB payload");
    }
    let value = u64::from_le_bytes(bytes[*offset..*offset + 8].try_into().unwrap());
    *offset += 8;
    Ok(value)
}

#[cfg(test)]
mod clmb_tests {
    use super::*;
    use std::fs;

    #[test]
    fn clmb_round_trip_across_small_blocks() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("test.clmb");
        let contigs = vec!["ctgA".to_string(), "ctgB".to_string()];
        let expected = vec![
            ClmbRecord {
                endpoint1: encode_endpoint(0, 0).unwrap(),
                endpoint2: encode_endpoint(1, 1).unwrap(),
                distances: vec![2, 10, 100],
            },
            ClmbRecord {
                endpoint1: encode_endpoint(0, 1).unwrap(),
                endpoint2: encode_endpoint(1, 0).unwrap(),
                distances: vec![7, u32::MAX as u64 + 11],
            },
        ];
        let mut writer = ClmbWriter::create(&path, &contigs, 20, Some(2), Some(5)).unwrap();
        for record in &expected {
            writer
                .write_record(record.endpoint1, record.endpoint2, &record.distances)
                .unwrap();
        }
        writer.finish().unwrap();
        assert!(is_clmb_file(&path).unwrap());

        let mut reader = ClmbReader::open(&path).unwrap();
        assert_eq!(reader.header.contigs, contigs);
        let mut observed = Vec::new();
        while let Some(block) = reader.next_block().unwrap() {
            observed.extend(block);
        }
        assert_eq!(observed, expected);
    }

    #[test]
    fn clmb_rejects_truncated_block() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("test.clmb");
        let contigs = vec!["ctgA".to_string(), "ctgB".to_string()];
        let mut writer = ClmbWriter::create(&path, &contigs, 1024, None, None).unwrap();
        writer.write_record(0, 2, &[2, 3, 4]).unwrap();
        writer.finish().unwrap();
        let length = fs::metadata(&path).unwrap().len();
        fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .unwrap()
            .set_len(length - 1)
            .unwrap();
        let mut reader = ClmbReader::open(&path).unwrap();
        assert!(reader.next_block().is_err());
    }
}
