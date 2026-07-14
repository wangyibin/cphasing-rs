#![allow(unused)]
#![allow(dead_code)]
use anyhow::Result as anyResult;

use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use log::LevelFilter;
use rand::prelude::*;
use std::collections::{ BTreeMap, HashMap, HashSet };
use std::hash::{ BuildHasherDefault, Hasher, Hash };
use std::borrow::Cow;
use std::path::{ Path, PathBuf };
use walkdir::WalkDir;
use smallvec::{ smallvec, SmallVec };
use std::thread;
use std::io::{ BufReader, BufRead, BufWriter, Write };
use std::sync::{ Arc, Mutex};
use std::sync::atomic::{AtomicUsize, AtomicU64, Ordering};
use std::fs::{ File, OpenOptions };
use twox_hash::XxHash64;
use rand::Rng;
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};
use polars::prelude::*;
use polars::enable_string_cache;

use crate::bed::{ Bed3, Bed4 };
use crate::contacts::{ Contacts, ContactRecord };
use crate::core::{ common_reader, common_writer };
use crate::core::{ 
    BaseTable, 
    ChromSize,
    ChromSizeRecord, 
    ContigPair, ContigPair2,
    binify
};


#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Debug)]
struct Contact {
    c1: u32,
    c2: u32,
    p1: u64,
    p2: u64,
}


struct Components {
    pub comps: Vec<Vec<u32>>,
    pub found: usize,
    pub skipped_large: usize,
}

fn connected_components(adj: &HashMap<u32, HashSet<u32>>, max_component_size: usize) -> Components {
    let mut visited: HashSet<u32> = HashSet::new();
    let mut comps = Vec::new();
    let mut skipped = 0;

    for &node in adj.keys() {
        if !visited.contains(&node) {
            let mut comp = Vec::new();
            let mut stack = vec![node];
            visited.insert(node);

            while let Some(curr) = stack.pop() {
                comp.push(curr);
                if let Some(neighbors) = adj.get(&curr) {
                    for &neighbor in neighbors {
                        if visited.insert(neighbor) {
                            stack.push(neighbor);
                        }
                    }
                }
            }

            if comp.len() > max_component_size {
                skipped += 1;
            } else {
                comps.push(comp);
            }
        }
    }

    Components {
        found: comps.len() + skipped,
        comps,
        skipped_large: skipped,
    }
}


fn load_enzyme_bed_to_lapper(bed: &str) -> anyResult<(HashMap<String, Lapper<u32, u32>>, Vec<String>)> {
    use std::io::BufRead;
    let f = common_reader(bed);
    let rdr = std::io::BufReader::new(f);

    let mut chrom_to_intervals: HashMap<String, Vec<Interval<u32, u32>>> = HashMap::new();
    let mut frag_id_strs = Vec::new();
    let mut current_id = 0u32;

    for line in rdr.lines().flatten() {
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') { continue; }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 3 { continue; }
        
        let chrom = fields[0].to_string();
        let start: u32 = fields[1].parse().unwrap_or(0);
        let end: u32 = fields[2].parse().unwrap_or(0);
        let name = if fields.len() > 3 { fields[3].to_string() } else { format!("{}:{}-{}", chrom, start, end) };

        frag_id_strs.push(name);
        chrom_to_intervals.entry(chrom)
            .or_default()
            .push(Interval { start, stop: end, val: current_id });
        
        current_id += 1;
    }

    let mut chrom_lappers = HashMap::new();
    for (chrom, mut ivs) in chrom_to_intervals {
        ivs.sort_by_key(|iv| iv.start);
        chrom_lappers.insert(chrom, Lapper::new(ivs));
    }

    Ok((chrom_lappers, frag_id_strs))
}


#[inline]
fn map_to_frag_id(chrom: &str, pos: u32, lappers: &HashMap<String, Lapper<u32, u32>>) -> Option<u32> {
    if let Some(lapper) = lappers.get(chrom) {
        if let Some(iv) = lapper.find(pos, pos + 1).next() {
            return Some(iv.val);
        }
    }
    None
}

fn bron_kerbosch_pivot<F>(
    r: &mut Vec<u32>,
    p: &HashSet<u32>,
    x: &HashSet<u32>,
    neigh: &F,
    cliques: &mut Vec<Vec<u32>>,
    min_size: usize,
) where
    F: Fn(u32) -> HashSet<u32>,
{
    if p.is_empty() && x.is_empty() {
        if r.len() >= min_size {
            cliques.push(r.clone());
        }
        return;
    }

    let p_union_x: HashSet<u32> = p.union(x).copied().collect();
    if let Some(&pivot) = p_union_x.iter().next() {
        let pivot_neigh = neigh(pivot);
        let p_without_n: HashSet<u32> = p.difference(&pivot_neigh).copied().collect();

        let mut p_copy = p.clone();
        let mut x_copy = x.clone();

        for v in p_without_n {
            let v_neigh = neigh(v);
            let mut new_r = r.clone();
            new_r.push(v);

            let new_p: HashSet<u32> = p_copy.intersection(&v_neigh).copied().collect();
            let new_x: HashSet<u32> = x_copy.intersection(&v_neigh).copied().collect();

            bron_kerbosch_pivot(&mut new_r, &new_p, &new_x, neigh, cliques, min_size);

            p_copy.remove(&v);
            x_copy.insert(v);
        }
    }
}

type SmallIntVec = SmallVec<[u32; 2]>;
type SmallIntVec4 = SmallVec<[u32; 4]>;


pub const _README: &str = r#"
# .pqs file format
The _contigsizes file contains the size of each contig in the .pqs file.
The _metadata file contains the metadata of the .pqs file.
The q0, q1 directory contains the parquet files of the data.

q0 mean the mapping quality of data >= 0.
q1 mean the mapping quality of data >= 1.

.pqs/  
|-- _contigsizes 
|-- _metadata
|-- _readme 
|-- q0/
|   |-- 0.parquet
|   |-- 1.parquet
|   |-- ...
|-- q1/
|   |-- 0.parquet
|   |-- 1.parquet
|   |-- ...
    "#;


pub const _METADATA: &str = r#"
{'format-version': '0.1.0',
 'format-url': 'https://github.com/wangyibin/CPhasing',
 'creation_time': 'REPLACE',
 'is_pqs': True,
 'format': 'pairs',
 'chunksize': CHUNKSIZE,
 'is_with_mapq': True,
 'columns': ['read_idx',
             'chrom1',
             'pos1',
             'chrom2',
             'pos2',
             'strand1',
             'strand2',
             'mapq'],
 'dtypes': {'read_idx': 'U',
            'chrom1': 'U',
            'pos1': 'pos_type_lower',
            'chrom2': 'U',
            'pos2': 'pos_type_lower',
            'strand1': 'U',
            'strand2': 'U',
            'mapq': 'uint8'},
 'schema': {'read_idx': String,
            'chrom1': Categorical,
            'pos1': pos_type,
            'chrom2': Categorical,
            'pos2': pos_type,
            'strand1': Categorical,
            'strand2': Categorical,
            'mapq': UInt8}} 
"#;



fn read_metadata_counts(pqs_dir: &str) -> Option<(u64, u64)> {
    // returns (q1_records, q0_records) if valid and >0
    let path = format!("{}/_metadata_counts", pqs_dir);

    if !std::path::Path::new(&path).exists() {
        return None;
    }


    let f = common_reader(&path);
    let rdr = BufReader::new(f);

    let mut q0: Option<u64> = None;
    let mut q1: Option<u64> = None;

    for line in rdr.lines().flatten() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // expected: "q0_records\t123"  or "q1_records  456"
        let mut it = line.split_whitespace();
        let k = it.next().unwrap_or("");
        let v = it.next().unwrap_or("");
        if v.is_empty() {
            continue;
        }
        if k == "q0_records" {
            q0 = v.parse::<u64>().ok();
        } else if k == "q1_records" {
            q1 = v.parse::<u64>().ok();
        }
    }

    match (q1, q0) {
        (Some(q1), Some(q0)) => Some((q1, q0)),
        _ => None,
    }
}


fn write_metadata_counts(pqs_dir: &str, q0_records: u64, q1_records: u64) -> anyResult<()> {
    let mut wtr = common_writer(format!("{}/_metadata_counts", pqs_dir).as_str());
    writeln!(wtr, "q0_records\t{}", q0_records)?;
    writeln!(wtr, "q1_records\t{}", q1_records)?;
    wtr.flush()?;
    Ok(())
}

fn copy_metadata_counts(src_dir: &str, dst_dir: &str) -> anyResult<()> {
    let src = format!("{}/_metadata_counts", src_dir);
    let dst = format!("{}/_metadata_counts", dst_dir);
    if Path::new(&src).exists() {
        std::fs::copy(src, dst)?;
    }
    Ok(())
}

fn sample_df_by_prob(mut df: DataFrame, prob: f64, seed: u64) -> PolarsResult<DataFrame> {
    if prob >= 1.0 || df.height() == 0 {
        return Ok(df);
    }
    if prob <= 0.0 {
        return Ok(df.head(Some(0)));
    }

    let p = prob.max(0.0).min(1.0);
    let mut rng = rand::rngs::SmallRng::seed_from_u64(seed);

    let mut mask: Vec<bool> = Vec::with_capacity(df.height());
    for _ in 0..df.height() {
        mask.push(rng.gen_bool(p));
    }

    let mask = BooleanChunked::from_slice("mask".into(), &mask);
    df.filter(&mask)
}

#[derive(Debug, Clone)]
pub struct PQS {
    pub file: String,
}

impl BaseTable for PQS {
    fn new(name: &String) -> PQS {
        PQS { file: name.clone() }
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

impl PQS {

    pub fn is_pqs(&self) -> bool {
        let path = Path::new(&self.file);
        let mut is_pqs = false;
        if path.exists() {
            if path.is_dir() {
                if Path::new(format!("{}/q0", self.file).as_str()).exists()
                    && Path::new(format!("{}/q1", self.file).as_str()).exists()
                    && Path::new(format!("{}/_contigsizes", self.file).as_str()).exists() {
                    is_pqs = true;
                }
            }
        }

        return is_pqs;

    }

    pub fn to_clm_v_0_3_0(&self, 
        min_contacts: u32,
        min_quality: u8,
        output: &String,
        output_split_contacts: bool,
        output_depth: bool,
        binsize: u32,
        threads: usize,
        disable_filter: bool,
        max_depth_ratio: f64,
        max_q0_ratio: f64
    ) -> anyResult<()> {
        use hashbrown::HashMap;
        use twox_hash::XxHash64;
        use std::fmt::Write;
        
        // polars::enable_string_cache(); 
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", "1"));
        }
        
        let min_mapq = min_quality as u32;
    
        let output_prefix = if output.ends_with(".gz") {
            output.trim_end_matches(".gz").trim_end_matches(".clm").to_string()
        } else {
            output.trim_end_matches(".clm").to_string()
        };

        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        let mut q1_count = 0;
        let mut q0_only_count = 0;
        if min_mapq == 0 {
            if let Some((q1, q0_total)) = read_metadata_counts(self.file.as_str()) {
                q1_count = q1;
                q0_only_count = q0_total.saturating_sub(q1);
                log::info!(
                    "Loaded counts from _metadata_counts: q1_records={}, q0_total_records={}, q0_only_records={}",
                    q1_count, q0_total, q0_only_count
                );
            } else {
                // fallback: scan parquet
                log::info!("Scanning q0 files to calculate sampling ratio...");
                let counts: Vec<(u64, u64)> = files.par_iter().map(|file| {
                    let mut lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                        Ok(lf) => lf,
                        Err(_) => return (0, 0),
                    };
                    let has_mq = lf.collect_schema().map(|s| s.contains("mapq")).unwrap_or(false);
                    if !has_mq { return (0, 0); }

                    let df = lf.select([col("mapq")]).collect().unwrap_or_else(|_| DataFrame::empty());
                    let mq = df.column("mapq").ok().and_then(|s| s.u8().ok());
                    if mq.is_none() { return (0, 0); }
                    let mq = mq.unwrap();

                    let mut c1: u64 = 0;
                    let mut c0: u64 = 0;
                    for v in mq.into_no_null_iter() {
                        if v >= 1 { c1 += 1; } else { c0 += 1; }
                    }
                    (c1, c0)
                }).collect();
                for (q1, q0) in counts { q1_count += q1; q0_only_count += q0; }
            }
        }

        let q0_sample_prob = if min_mapq == 0 && q0_only_count > 0 {
            if max_q0_ratio == 0.0 {
                1.0
            } else {
                let target = (q1_count as f64 * max_q0_ratio).min(q0_only_count as f64);
                (target / q0_only_count as f64).min(1.0)
            }
        } else {
            1.0
        };

        if min_mapq == 0 {
                log::info!("Q1 count: {}, Q0-only count: {}, Sampling prob: {:.4}", q1_count, q0_only_count, q0_sample_prob);
        }




        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u32>().unwrap();
                contigsizes.insert(contig, size);
        }

        let mut sorted_names: Vec<_> = contigsizes.keys().cloned().collect();
        sorted_names.sort();

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = 
            sorted_names.iter().enumerate().map(|(i, k)| (k.clone(), i as u32)).collect();
        let idx_contig: HashMap<u32, String, BuildHasherDefault<XxHash64>> = 
            sorted_names.iter().enumerate().map(|(i, k)| (i as u32, k.clone())).collect();
        let idx_sizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = 
            contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), v.clone())).collect();
        let idx_contig_sizes: HashMap<u32, (String, u32), BuildHasherDefault<XxHash64>> = 
            contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), (k.clone(), v.clone()))).collect();

        let (sender, receiver) = bounded::<LazyFrame>(100);
        type FastHasher = std::hash::BuildHasherDefault<XxHash64>;
        type ShardMap = HashMap<(u32, u32), Vec<SmallIntVec>, FastHasher>;
        let num_shards = 128;
        let shards: Arc<Vec<Mutex<ShardMap>>> = Arc::new(
            (0..num_shards).map(|_| Mutex::new(ShardMap::with_capacity_and_hasher(2048, FastHasher::default()))).collect()
        );
        let contig_idx = Arc::new(contig_idx);

        files.par_iter().for_each(|file| {
            let mut lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                Ok(lf) => lf,
                Err(_) => return,
            };

            let has_mapq = lf.collect_schema().map(|s| s.contains("mapq")).unwrap_or(false);
            let select_cols = if has_mapq {
                vec![col("chrom1"), col("chrom2"), col("pos1"), col("pos2"), col("mapq")]
            } else {
                vec![col("chrom1"), col("chrom2"), col("pos1"), col("pos2")]
            };

            let df = lf.select(select_cols).collect().expect("Polars collect failed");
            if df.height() == 0 { return; }

            let c1_col = df.column("chrom1").unwrap().categorical().unwrap();
            let c2_col = df.column("chrom2").unwrap().categorical().unwrap();
            let p1_col = df.column("pos1").unwrap().u32().unwrap();
            let p2_col = df.column("pos2").unwrap().u32().unwrap();

            let mapq_col = if has_mapq { Some(df.column("mapq").unwrap().u8().unwrap()) } else { None };
            
            let rev1 = c1_col.get_rev_map();
            let phys1_ca = c1_col.physical();
            let max_id1 = phys1_ca.max().unwrap_or(0) as usize;
            let mut lookup1 = vec![u32::MAX; max_id1 + 1];
            for i in 0..=max_id1 {
                if let Some(name) = rev1.get_optional(i as u32) {
                    if let Some(&idx) = contig_idx.get(name) { lookup1[i] = idx; }
                }
            }
            
            let rev2 = c2_col.get_rev_map();
            let phys2_ca = c2_col.physical();
            let max_id2 = phys2_ca.max().unwrap_or(0) as usize;
            let mut lookup2 = vec![u32::MAX; max_id2 + 1];
            for i in 0..=max_id2 {
                if let Some(name) = rev2.get_optional(i as u32) {
                    if let Some(&idx) = contig_idx.get(name) { lookup2[i] = idx; }
                }
            }


            let mut local_map: HashMap<(u32, u32), Vec<SmallIntVec>, FastHasher> = 
                    HashMap::with_capacity_and_hasher(1024, FastHasher::default());
            let seed = 42; 
            let mut rng = rand::rngs::SmallRng::seed_from_u64(seed);
            phys1_ca.downcast_iter()
                .zip(phys2_ca.downcast_iter())
                .zip(p1_col.downcast_iter())
                .zip(p2_col.downcast_iter())
                .enumerate().for_each(|(chunk_idx, (((c1, c2), p1), p2))| {
                    let mq_chunk = mapq_col.as_ref().map(|mq| mq.downcast_iter().nth(chunk_idx).unwrap());
                    for i in 0..c1.len() {
                        if let Some(mq) = mq_chunk {
                            let m = mq.value(i);
                            if m == 0 && min_mapq == 0 {
                                if q0_sample_prob < 1.0 && !rng.gen_bool(q0_sample_prob) {
                                    continue;
                                }
                            } else if m < min_quality {
                                continue;
                            }
                        }

                        let idx1 = lookup1[c1.value(i) as usize];
                        let idx2 = lookup2[c2.value(i) as usize];
                        if idx1 != u32::MAX && idx2 != u32::MAX {
                            let key = if idx1 <= idx2 { (idx1, idx2) } else { (idx2, idx1) };
                            let val = if idx1 <= idx2 {
                                smallvec![p1.value(i), p2.value(i)]
                            } else {
                                smallvec![p2.value(i), p1.value(i)]
                            };
                            local_map.entry(key).or_default().push(val);
                        }
                    }
            });

            let mut shard_buckets: Vec<Vec<((u32, u32), Vec<SmallIntVec>)>> = (0..num_shards).map(|_| Vec::new()).collect();
            for (key, val) in local_map {
                let mut s = XxHash64::default();
                key.hash(&mut s);
                let shard_idx = (s.finish() % num_shards as u64) as usize;
                shard_buckets[shard_idx].push((key, val));
            }

            for (idx, bucket) in shard_buckets.into_iter().enumerate() {
                if bucket.is_empty() { continue; }
                let mut lock = shards[idx].lock().unwrap();
                for (key, val) in bucket {
                    lock.entry(key).or_default().extend(val);
                }
            }
        });

        let mut data = HashMap::with_capacity_and_hasher(contig_idx.len() * 2, FastHasher::default());
        let shards = Arc::try_unwrap(shards).expect("Arc unwrap failed");
        for shard in shards {
            let inner = shard.into_inner().unwrap();
            data.extend(inner); 
        }

        if output_split_contacts{
            log::info!("Calculating the distance between split contigs");
            let contacts = Contacts::new(&format!(".{}.pixels", output_prefix));
            let split_contigsizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = idx_sizes
                .par_iter()
                .map(|(k, v)| (*k, *v / 2))
                .collect(); 
           
            let writer = common_writer(format!("{}.split.contacts.gz", output_prefix.to_string()).as_str());
            let writer = Arc::new(Mutex::new(writer));
            data.par_iter().for_each(|(cp, vec) | {
                if vec.len() < min_contacts as usize {
                    return;
                }
                let (contig1, length1) = idx_contig_sizes.get(&cp.0).unwrap();
                let (contig2, length2) = idx_contig_sizes.get(&cp.1).unwrap();
                let res = vec.iter().map(
                    |x| {
                        let pos1 = x[0];
                        let pos2 = x[1];
                        let split_index1 = (pos1 / (length1 / 2)) as u8;
                        let split_index2 = (pos2 / (length2 / 2)) as u8;
                    
                        (split_index1, split_index2)
                    }
                ).collect::<Vec<_>>();

                let mut contact_hash = HashMap::with_capacity(4);
                res.iter().for_each(|(split_idx1, split_idx2)| {
                    *contact_hash.entry((split_idx1, split_idx2)).or_insert(0) += 1;
                });
                
                let mut buffer = Vec::with_capacity(4);
                contact_hash.iter().for_each(|(cp, count)| {
                    if count >= &min_contacts {
                        buffer.push(format!("{}_{}\t{}_{}\t{}\n", contig1, cp.0, contig2, cp.1, count));
                    }
                    
                });
                let buffer = buffer.join("");
                let mut writer = writer.lock().unwrap();
                writer.write_all(buffer.as_bytes()).unwrap();

            }); 

            log::info!("Successful output split contacts file `{}`", 
                                &format!("{}.split.contacts.gz", output_prefix.to_string()));

            drop(split_contigsizes);
        }

        if output_depth {
            log::info!("Calculating the depth of each contig");
        
            let mut depth: HashMap<u32, Vec<u32>> = idx_sizes
                                                        .clone()
                                                        .into_iter()
                                                        .map(|(chrom, size)| {
                                                            let num_bins = (size / binsize as u32 + 1) as usize;
                                                            (chrom, vec![0; num_bins])
                                                        }).collect();
            
            data.iter().for_each(|(cp, vec)| {
                let res = vec.iter().map(
                    |x| {
                        let pos1 = x[0];
                        let pos2 = x[1];
                        let split_index1 = (pos1 / binsize) as u32;
                        let split_index2 = (pos2 / binsize) as u32;
                        
                        (split_index1, split_index2)
                    }
                ).collect::<Vec<_>>();
                
                res.iter().for_each(|(split_index1, split_index2)| {
                    if let Some(v) = depth.get_mut(&cp.0) {
                        if let Some(v) = v.get_mut(*split_index1 as usize) {
                            *v += 1;
                        }
                
                    }

                    if let Some(v) = depth.get_mut(&cp.1) {
                        if let Some(v) = v.get_mut(*split_index2 as usize) {
                            *v += 1;
                        }
                    }

                    // *depth.get_mut(&cp.0).unwrap().get_mut(*split_index1 as usize).unwrap() += 1;
                    // *depth.get_mut(&cp.1).unwrap().get_mut(*split_index2 as usize).unwrap() += 1;
                });
                
            });

            let depth: BTreeMap<_, _> = depth.into_iter().collect();
            let writer = common_writer(format!("{}.depth", output_prefix.to_string()).as_str());
            let wtr = Arc::new(Mutex::new(writer));

            depth.par_iter().for_each(|(contig, bins)| {
                let (contig, size) = idx_contig_sizes.get(contig).unwrap();
                
                let mut buffer = Vec::with_capacity(bins.len() * 50);
                for (bin, count) in bins.iter().enumerate() {
                    let bin_start = bin * binsize as usize;
                    let mut bin_end = bin_start + binsize as usize;
    
                    if bin_end > (*size).try_into().unwrap() {
                        bin_end = *size as usize;
                    }
                    if bin_start == bin_end {
                        continue;
                    }
                    buffer.extend_from_slice(format!("{}\t{}\t{}\t{}\n", contig, bin_start, bin_end, count).as_bytes());
                }
                {
                    let mut wtr = wtr.lock().unwrap();
                    wtr.write_all(&buffer).unwrap();
                }
                
            });

            log::info!("Successful output depth file `{}`", 
                                &format!("{}.depth", output_prefix.to_string()));
      
        }

        let mut blacklist: HashSet<(u32, u32)> = HashSet::new();

        if !disable_filter && max_depth_ratio > 0.0 {
            log::info!("Identifying high depth regions...");
            
            let bin_depths: HashMap<(u32, u32), u32> = data.par_iter()
                .fold(
                    || HashMap::new(),
                    |mut acc: HashMap<(u32, u32), u32>, ((c1, c2), vec)| {
                        for pair in vec {
                            let pos1 = pair[0];
                            let pos2 = pair[1];
                            let bin1 = pos1 / binsize;
                            let bin2 = pos2 / binsize;

                            *acc.entry((*c1, bin1)).or_insert(0) += 1;
                            *acc.entry((*c2, bin2)).or_insert(0) += 1;
                        }
                        acc
                    }
                )
                .reduce(
                    || HashMap::new(),
                    |mut acc, part| {
                        for (k, v) in part {
                            *acc.entry(k).or_insert(0) += v;
                        }
                        acc
                    }
                );

            if !bin_depths.is_empty() {
                let mut depths: Vec<u32> = bin_depths.values().cloned().collect();
                depths.sort_unstable();
                let n = depths.len();
                let mid = n / 2;
                let median_depth = depths[mid] as f64;
                
                // let sum_depth: u64 = depths.iter().map(|&x| x as u64).sum();
                // let mean_depth = sum_depth as f64 / depths.len() as f64;

                // let threshold = median_depth * max_depth_ratio;
                
                let mut deviations: Vec<u32> = depths.iter()
                    .map(|&x| (x as i64 - median_depth as i64).abs() as u32)
                    .collect();
                deviations.sort_unstable();
                let mad = deviations[n / 2] as f64;
                let min_absolute_threshold = 10.0;
                let threshold = median_depth + (max_depth_ratio * mad).max(min_absolute_threshold);
                log::info!("Median bin depth: {}, Threshold: {:.2}", median_depth, threshold);

                for (key, &depth) in bin_depths.iter() {
                    if depth as f64 > threshold {
                        blacklist.insert(*key);
                    }
                }
                log::info!("Identified {} high depth bins out of {}.", blacklist.len(), bin_depths.len());
            }
        }

        let wtr = common_writer(output.as_str());
        let wtr = Arc::new(Mutex::new(wtr));
        let blacklist = Arc::new(blacklist);

        data.par_iter().for_each(|(cp, vec)| {
            if cp.0 == cp.1 { return; }

            let (contig1, length1) = idx_contig_sizes.get(&cp.0).unwrap();
            let (contig2, length2) = idx_contig_sizes.get(&cp.1).unwrap();

            let filtered_vec: Vec<&SmallIntVec> = if !disable_filter && !blacklist.is_empty() {
                vec.iter().filter(|pair| {
                    let bin1 = pair[0] / binsize;
                    let bin2 = pair[1] / binsize;
                    !blacklist.contains(&(cp.0, bin1)) && !blacklist.contains(&(cp.1, bin2))
                }).collect()
            } else {
                vec.iter().collect()
            };

            let count = filtered_vec.len();
            if count < min_contacts as usize { return; }

            let p1: Vec<u32> = filtered_vec.iter().map(|v| v[0]).collect();
            let p2: Vec<u32> = filtered_vec.iter().map(|v| v[1]).collect();

            
            let mut output_buffer = String::with_capacity(count * 50);
            let mut itoa_buf = itoa::Buffer::new();

            for i in 0..4 {
                let header = match i {
                    0 => format!("{}+ {}+\t{}\t", contig1, contig2, count),
                    1 => format!("{}+ {}-\t{}\t", contig1, contig2, count),
                    2 => format!("{}- {}+\t{}\t", contig1, contig2, count),
                    _ => format!("{}- {}-\t{}\t", contig1, contig2, count),
                };
                output_buffer.push_str(&header);

                let mut results = vec![0u32; count];
                let mut chunks_iter = results.chunks_exact_mut(8);
                let p1_chunks = p1.chunks_exact(8);
                let p2_chunks = p2.chunks_exact(8);

                for ((res_c, p1_c), p2_c) in chunks_iter.by_ref().zip(p1_chunks).zip(p2_chunks) {
                    for j in 0..8 {
                        res_c[j] = match i {
                            0 => length1.wrapping_sub(p1_c[j]).wrapping_add(p2_c[j]),
                            1 => length1.wrapping_sub(p1_c[j]).wrapping_add(*length2).wrapping_sub(p2_c[j]),
                            2 => p1_c[j].wrapping_add(p2_c[j]),
                            _ => p1_c[j].wrapping_add(*length2).wrapping_sub(p2_c[j]),
                        };
                    }
                }
                
                let rem_res = chunks_iter.into_remainder();
                let rem_p1 = p1.chunks_exact(8).remainder();
                let rem_p2 = p2.chunks_exact(8).remainder();

                for k in 0..rem_res.len() {
                    rem_res[k] = match i {
                        0 => length1.wrapping_sub(rem_p1[k]).wrapping_add(rem_p2[k]),
                        1 => length1.wrapping_sub(rem_p1[k]).wrapping_add(*length2).wrapping_sub(rem_p2[k]),
                        2 => rem_p1[k].wrapping_add(rem_p2[k]),
                        _ => rem_p1[k].wrapping_add(*length2).wrapping_sub(rem_p2[k]),
                    };
                }

                for (idx, val) in results.iter().enumerate() {
                    output_buffer.push_str(itoa_buf.format(*val));
                    if idx < count - 1 {
                        output_buffer.push(' ');
                    }
                }
                output_buffer.push('\n');
            }
            let mut writer_lock = wtr.lock().unwrap();
            writer_lock.write_all(output_buffer.as_bytes()).unwrap();

        });
      
        drop(data);
        drop(blacklist);

        {
            let mut wtr = wtr.lock().unwrap();
            wtr.flush().unwrap();
        }

        log::info!("Successful output clm file `{}`", output);

        Ok(())
    }

    pub fn to_clm(&self, 
        min_contacts: u32,
        min_quality: u8,
        output: &String,
        output_split_contacts: bool,
        output_depth: bool,
        binsize: u32,
        threads: usize,
        disable_filter: bool,
        max_depth_ratio: f64,
        max_q0_ratio: f64
    ) -> anyResult<()> {
        use hashbrown::HashMap;
        use twox_hash::XxHash64;
        use std::fmt::Write;
        use std::sync::atomic::{AtomicU32, Ordering};
        
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", "1"));
        }
        let binsize = binsize as u64;
        let min_mapq = min_quality as u32;
    
        let output_prefix = if output.ends_with(".gz") {
            output.trim_end_matches(".gz").trim_end_matches(".clm").to_string()
        } else {
            output.trim_end_matches(".clm").to_string()
        };

        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        let mut q1_count = 0;
        let mut q0_only_count = 0;
        if min_mapq == 0 {
            if let Some((q1, q0_total)) = read_metadata_counts(self.file.as_str()) {
                q1_count = q1;
                q0_only_count = q0_total.saturating_sub(q1);
                log::info!(
                    "Loaded counts from _metadata_counts: q1_records={}, q0_total_records={}, q0_only_records={}",
                    q1_count, q0_total, q0_only_count
                );
            } else {
                log::info!("Scanning q0 files to calculate sampling ratio...");
                let counts: Vec<(u64, u64)> = files.par_iter().map(|file| {
                    let mut lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                        Ok(lf) => lf,
                        Err(_) => return (0, 0),
                    };
                    let has_mq = lf.collect_schema().map(|s| s.contains("mapq")).unwrap_or(false);
                    if !has_mq { return (0, 0); }

                    let df = lf.select([col("mapq")]).collect().unwrap_or_else(|_| DataFrame::empty());
                    let mq = df.column("mapq").ok().and_then(|s| s.u8().ok());
                    if mq.is_none() { return (0, 0); }
                    let mq = mq.unwrap();

                    let mut c1: u64 = 0;
                    let mut c0: u64 = 0;
                    for v in mq.into_no_null_iter() {
                        if v >= 1 { c1 += 1; } else { c0 += 1; }
                    }
                    (c1, c0)
                }).collect();
                for (q1, q0) in counts { q1_count += q1; q0_only_count += q0; }
            }
        }

        let q0_sample_prob = if min_mapq == 0 && q0_only_count > 0 {
            if max_q0_ratio == 0.0 {
                1.0
            } else {
                let target = (q1_count as f64 * max_q0_ratio).min(q0_only_count as f64);
                (target / q0_only_count as f64).min(1.0)
            }
        } else {
            1.0
        };

        if min_mapq == 0 {
            log::info!("Q1 count: {}, Q0-only count: {}, Sampling prob: {:.4}", q1_count, q0_only_count, q0_sample_prob);
        }

        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u64>().unwrap();
            contigsizes.insert(contig, size);
        }

        let mut sorted_names: Vec<_> = contigsizes.keys().cloned().collect();
        sorted_names.sort();

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = 
            sorted_names.iter().enumerate().map(|(i, k)| (k.clone(), i as u32)).collect();
        let idx_contig: HashMap<u32, String, BuildHasherDefault<XxHash64>> = 
            sorted_names.iter().enumerate().map(|(i, k)| (i as u32, k.clone())).collect();
        let idx_sizes: HashMap<u32, u64, BuildHasherDefault<XxHash64>> = 
            contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), v.clone())).collect();
        let idx_contig_sizes: HashMap<u32, (String, u64), BuildHasherDefault<XxHash64>> = 
            contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), (k.clone(), v.clone()))).collect();

        type FastHasher = std::hash::BuildHasherDefault<XxHash64>;
        let num_shards = 128;
        let shards: Arc<Vec<Mutex<Vec<Contact>>>> = Arc::new(
            (0..num_shards).map(|_| Mutex::new(Vec::with_capacity(65536))).collect()
        );
        let contig_idx = Arc::new(contig_idx);

        files.par_iter().for_each(|file| {
            let mut lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                Ok(lf) => lf,
                Err(_) => return,
            };

            let has_mapq = lf.collect_schema().map(|s| s.contains("mapq")).unwrap_or(false);
            let select_cols = if has_mapq {
                vec![col("chrom1"), col("chrom2"), col("pos1"), col("pos2"), col("mapq")]
            } else {
                vec![col("chrom1"), col("chrom2"), col("pos1"), col("pos2")]
            };

            let df = lf.select(select_cols).collect().expect("Polars collect failed");
            if df.height() == 0 { return; }

            let c1_col = df.column("chrom1").unwrap().categorical().unwrap();
            let c2_col = df.column("chrom2").unwrap().categorical().unwrap();
            let binding = df.column("pos1").unwrap().cast(&DataType::UInt64).unwrap();
            let p1_col = binding.u64().unwrap();
            let binding = df.column("pos2").unwrap().cast(&DataType::UInt64).unwrap();
            let p2_col = binding.u64().unwrap();

            let mapq_col = if has_mapq { Some(df.column("mapq").unwrap().u8().unwrap()) } else { None };
            
            let rev1 = c1_col.get_rev_map();
            let phys1_ca = c1_col.physical();
            let max_id1 = phys1_ca.max().unwrap_or(0) as usize;
            let mut lookup1 = vec![u32::MAX; max_id1 + 1];
            for i in 0..=max_id1 {
                if let Some(name) = rev1.get_optional(i as u32) {
                    if let Some(&idx) = contig_idx.get(name) { lookup1[i] = idx; }
                }
            }
            
            let rev2 = c2_col.get_rev_map();
            let phys2_ca = c2_col.physical();
            let max_id2 = phys2_ca.max().unwrap_or(0) as usize;
            let mut lookup2 = vec![u32::MAX; max_id2 + 1];
            for i in 0..=max_id2 {
                if let Some(name) = rev2.get_optional(i as u32) {
                    if let Some(&idx) = contig_idx.get(name) { lookup2[i] = idx; }
                }
            }

            let mut local_shards: Vec<Vec<Contact>> = (0..num_shards)
                .map(|_| Vec::with_capacity(1024))
                .collect();

            let seed = 42; 
            let mut rng = rand::rngs::SmallRng::seed_from_u64(seed);
            phys1_ca.downcast_iter()
                .zip(phys2_ca.downcast_iter())
                .zip(p1_col.downcast_iter())
                .zip(p2_col.downcast_iter())
                .enumerate().for_each(|(chunk_idx, (((c1, c2), p1), p2))| {
                    let mq_chunk = mapq_col.as_ref().map(|mq| mq.downcast_iter().nth(chunk_idx).unwrap());
                    for i in 0..c1.len() {
                        if let Some(mq) = mq_chunk {
                            let m = mq.value(i);
                            if m == 0 && min_mapq == 0 {
                                if q0_sample_prob < 1.0 && !rng.gen_bool(q0_sample_prob) {
                                    continue;
                                }
                            } else if m < min_quality {
                                continue;
                            }
                        }

                        let idx1 = lookup1[c1.value(i) as usize];
                        let idx2 = lookup2[c2.value(i) as usize];
                        if idx1 != u32::MAX && idx2 != u32::MAX {
                            let (fc1, fp1, fc2, fp2) = if idx1 <= idx2 {
                                (idx1, p1.value(i), idx2, p2.value(i))
                            } else {
                                (idx2, p2.value(i), idx1, p1.value(i))
                            };

                            let contact = Contact { c1: fc1, c2: fc2, p1: fp1, p2: fp2 };
                            let mut s = XxHash64::default();
                            (fc1, fc2).hash(&mut s);
                            let shard_idx = (s.finish() % num_shards as u128 as u64) as usize;
                            local_shards[shard_idx].push(contact);
                        }
                    }
            });

            for (idx, local_vec) in local_shards.into_iter().enumerate() {
                if local_vec.is_empty() { continue; }
                let mut lock = shards[idx].lock().unwrap();
                lock.extend(local_vec);
            }
        });

        let shards = Arc::try_unwrap(shards).expect("Arc unwrap failed");
        let mut shards: Vec<Vec<Contact>> = shards.into_iter().map(|m| m.into_inner().unwrap()).collect();

        log::info!("Sorting contacts per shard...");
        shards.par_iter_mut().for_each(|shard_contacts| {
            shard_contacts.sort_unstable();
        });

        let mut depth: Vec<Vec<AtomicU32>> = Vec::new();
        let compute_depth = output_depth || (!disable_filter && max_depth_ratio > 0.0);

        if compute_depth {
            log::info!("Pre-calculating depth of each contig concurrently...");
            depth = (0..sorted_names.len())
                .map(|i| {
                    let size = idx_sizes.get(&(i as u32)).cloned().unwrap_or(0);
                    let num_bins = (size / binsize + 1) as usize;
                    let mut v = Vec::with_capacity(num_bins);
                    for _ in 0..num_bins {
                        v.push(AtomicU32::new(0));
                    }
                    v
                })
                .collect();

            shards.par_iter().for_each(|shard_contacts| {
                for contact in shard_contacts {
                    let bin1 = contact.p1 / binsize;
                    let bin2 = contact.p2 / binsize;

                    if let Some(bins) = depth.get(contact.c1 as usize) {
                        if let Some(cell) = bins.get(bin1 as usize) {
                            cell.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    if let Some(bins) = depth.get(contact.c2 as usize) {
                        if let Some(cell) = bins.get(bin2 as usize) {
                            cell.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
            });
        }

        if output_depth {
            let writer = common_writer(format!("{}.depth", output_prefix).as_str());
            let wtr = Arc::new(Mutex::new(writer));

            depth.par_iter().enumerate().for_each(|(contig_idx_val, bins)| {
                let &(ref contig, size) = idx_contig_sizes.get(&(contig_idx_val as u32)).unwrap();
                let mut buffer = Vec::with_capacity(bins.len() * 50);
                for (bin, cell) in bins.iter().enumerate() {
                    let count = cell.load(Ordering::Relaxed);
                    let bin_start = bin * binsize as usize;
                    let mut bin_end = bin_start + binsize as usize;
    
                    if bin_end > size as usize {
                        bin_end = size as usize;
                    }
                    if bin_start == bin_end {
                        continue;
                    }
                    buffer.extend_from_slice(format!("{}\t{}\t{}\t{}\n", contig, bin_start, bin_end, count).as_bytes());
                }
                if !buffer.is_empty() {
                    let mut wtr = wtr.lock().unwrap();
                    wtr.write_all(&buffer).unwrap();
                }
            });

            log::info!("Successful output depth file `{}`", &format!("{}.depth", output_prefix));
        }

        let mut blacklist: HashSet<(u32, u32)> = HashSet::new();

        if !disable_filter && max_depth_ratio > 0.0 {
            log::info!("Identifying high depth regions...");
            let mut bin_depths_flat = Vec::new();
            for (c, bins) in depth.iter().enumerate() {
                for (b, cell) in bins.iter().enumerate() {
                    let val = cell.load(Ordering::Relaxed);
                    if val > 0 {
                        bin_depths_flat.push(((c as u32, b as u32), val));
                    }
                }
            }

            if !bin_depths_flat.is_empty() {
                let mut depths: Vec<u32> = bin_depths_flat.iter().map(|&(_, val)| val).collect();
                depths.sort_unstable();
                let n = depths.len();
                let mid = n / 2;
                let median_depth = depths[mid] as f64;
                
                let mut deviations: Vec<u32> = depths.iter()
                    .map(|&x| (x as i64 - median_depth as i64).abs() as u32)
                    .collect();
                deviations.sort_unstable();
                let mad = deviations[n / 2] as f64;
                let min_absolute_threshold = 10.0;
                let threshold = median_depth + (max_depth_ratio * mad).max(min_absolute_threshold);
                log::info!("Median bin depth: {}, Threshold: {:.2}", median_depth, threshold);

                for ((c, b), depth_val) in bin_depths_flat.iter() {
                    if *depth_val as f64 > threshold {
                        blacklist.insert((*c, *b));
                    }
                }
                log::info!("Identified {} high depth bins out of {}.", blacklist.len(), bin_depths_flat.len());
            }
        }

        let split_writer = if output_split_contacts {
            let w = common_writer(format!("{}.split.contacts.gz", output_prefix).as_str());
            Some(Arc::new(Mutex::new(w)))
        } else {
            None
        };

        let wtr = common_writer(output.as_str());
        let wtr = Arc::new(Mutex::new(wtr));
        let blacklist = Arc::new(blacklist);

        log::info!("Generating final .clm output from sorted shards...");
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();

        pool.install(|| {
            shards.into_par_iter().for_each(|shard_contacts| {
                if shard_contacts.is_empty() { return; }

                let mut local_split_buf = Vec::with_capacity(256 * 1024);
                let mut local_clm_buf = Vec::with_capacity(1024 * 1024);

                let mut i = 0;
                while i < shard_contacts.len() {
                    let start_idx = i;
                    let c1 = shard_contacts[start_idx].c1;
                    let c2 = shard_contacts[start_idx].c2;

                    while i < shard_contacts.len() && shard_contacts[i].c1 == c1 && shard_contacts[i].c2 == c2 {
                        i += 1;
                    }
                    let end_idx = i;
                    let group = &shard_contacts[start_idx..end_idx];

                    let &(ref contig1, length1) = idx_contig_sizes.get(&c1).unwrap();
                    let &(ref contig2, length2) = idx_contig_sizes.get(&c2).unwrap();

                    if output_split_contacts {
                        if group.len() >= min_contacts as usize {
                            let mut contact_hash = HashMap::new();
                            for contact in group {
                                let split_index1 = (contact.p1 / (length1 / 2)) as u8;
                                let split_index2 = (contact.p2 / (length2 / 2)) as u8;
                                *contact_hash.entry((split_index1, split_index2)).or_insert(0) += 1;
                            }

                            let mut buffer = Vec::with_capacity(4);
                            for ((s_idx1, s_idx2), count) in contact_hash {
                                if count >= min_contacts {
                                    buffer.push(format!("{}_{}\t{}_{}\t{}\n", contig1, s_idx1, contig2, s_idx2, count));
                                }
                            }
                            if !buffer.is_empty() {
                                let buffer_str = buffer.join("");
                                local_split_buf.extend_from_slice(buffer_str.as_bytes());
                            }
                        }
                    }
                    if c1 == c2 { continue; }

                    let filtered_vec: Vec<&Contact> = if !disable_filter && !blacklist.is_empty() {
                        group.iter().filter(|contact| {
                            let bin1 = contact.p1 / binsize;
                            let bin2 = contact.p2 / binsize;
                            !blacklist.contains(&(c1, bin1.try_into().unwrap())) && !blacklist.contains(&(c2, bin2.try_into().unwrap()))
                        }).collect()
                    } else {
                        group.iter().collect()
                    };

                    let count = filtered_vec.len();
                    if count < min_contacts as usize { continue; }

                    let mut output_buffer = String::with_capacity(count * 50);
                    let mut itoa_buf = itoa::Buffer::new();

                    for ori in 0..4 {
                        let header = match ori {
                            0 => format!("{}+ {}+\t{}\t", contig1, contig2, count),
                            1 => format!("{}+ {}-\t{}\t", contig1, contig2, count),
                            2 => format!("{}- {}+\t{}\t", contig1, contig2, count),
                            _ => format!("{}- {}-\t{}\t", contig1, contig2, count),
                        };
                        output_buffer.push_str(&header);
                        let len1 = length1 as u64;
                        let len2 = length2 as u64;

                        for (j, contact) in filtered_vec.iter().enumerate() {
                            let cp1 = contact.p1 as u64;
                            let cp2 = contact.p2 as u64;
                            let val = match ori {
                                0 => len1.saturating_sub(cp1).saturating_add(cp2),
                                1 => len1.saturating_sub(cp1).saturating_add(len2).saturating_sub(cp2),
                                2 => cp1.saturating_add(cp2),
                                _ => cp1.saturating_add(len2).saturating_sub(cp2),
                                // 0 => length1.wrapping_sub(contact.p1).wrapping_add(contact.p2),
                                // 1 => length1.wrapping_sub(contact.p1).wrapping_add(length2).wrapping_sub(contact.p2),
                                // 2 => contact.p1.wrapping_add(contact.p2),
                                // _ => contact.p1.wrapping_add(length2).wrapping_sub(contact.p2),
                            };

                            output_buffer.push_str(itoa_buf.format(val));
                            if j < count - 1 {
                                output_buffer.push(' ');
                            }
                        }
                        output_buffer.push('\n');
                    }
                    local_clm_buf.extend_from_slice(output_buffer.as_bytes());
                }

                if !local_split_buf.is_empty() {
                    let mut writer = split_writer.as_ref().unwrap().lock().unwrap();
                    writer.write_all(&local_split_buf).unwrap();
                }
                if !local_clm_buf.is_empty() {
                    let mut writer_lock = wtr.lock().unwrap();
                    writer_lock.write_all(&local_clm_buf).unwrap();
                }
            });
        });

        if let Some(w) = split_writer {
            let mut writer = w.lock().unwrap();
            writer.flush().unwrap();
        }

        {
            let mut writer = wtr.lock().unwrap();
            writer.flush().unwrap();
        }

        log::info!("Successful output clm file `{}`", output);

        Ok(())
    }

    pub fn to_mnd(&self, min_quality: u8, output: &String) -> anyResult<()> {
        enable_string_cache();
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", 10));
        }
        
        let min_mapq = min_quality as u32;
        // get prefix of parquet_dir
        let output_prefix = if output.ends_with(".gz") {
            output.trim_end_matches(".gz").to_string()
        } else {
            output.to_string()
        };

        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };


        let results = files.into_par_iter().map(|file| {
            let mut df = LazyFrame::scan_parquet(file,  ScanArgsParquet::default()).unwrap();

            if min_mapq > 1 {
                df = df.clone().filter(col("mapq").gt_eq(min_mapq));
            }
        
            let result = df.select(
                &[
                    lit(0i32).alias("strand1"),
                    col("chrom1"),
                    col("pos1"),
                    lit(0u32).alias("frag1"),
                    lit(0i32).alias("strand2"),
                    col("chrom2"),
                    col("pos2"),
                    lit(1u32).alias("frag2"),
                    col("mapq").alias("mapq1"),
                    lit("-").alias("cigar1"),
                    lit("-").alias("sequence1"),
                    col("mapq").alias("mapq2"),
                    lit("-").alias("cigar2"),
                    lit("-").alias("sequence2"),
                    lit("-").alias("readname1"),
                    lit("-").alias("readname2"),
                ]
            );
            result
        }).collect::<Vec<_>>();

        let mut file = File::create(output.as_str()).unwrap();

        if !results.is_empty() {
            let mut df = results[0].clone().collect().unwrap();
            CsvWriter::new(&mut file)
                .include_header(false)
                .with_separator(b' ')
                .finish(&mut df)
                .unwrap();
        }

        for result in results.iter().skip(1) {
            let mut file = OpenOptions::new().append(true).open(output.as_str()).unwrap();
            let mut df = result.clone().collect().unwrap();
            CsvWriter::new(&mut file)
                .include_header(false)
                .with_separator(b' ')
                .finish(&mut df)
                .unwrap();
        }

        log::info!("Successful output mnd file `{}`", output);

        Ok(())
    }

    pub fn to_contacts(&self, min_contacts: u32, min_quality: u8, output: &String) -> anyResult<()> {
        use hashbrown::HashMap;
        polars::enable_string_cache(); 
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", 4));
        }
        
        let min_mapq = min_quality as u32;
       
        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u32>().unwrap();
                contigsizes.insert(contig, size);
        }

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (k.clone(), i as u32)).collect();
        let idx_contig: HashMap<u32, String, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (i as u32, k.clone())).collect();
        let idx_sizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), v.clone())).collect();

        let (sender, receiver) = bounded::<LazyFrame>(100);
        let data = Arc::new(Mutex::new(HashMap::new()));

       

        let consumer_handles: Vec<_> = (0..8).map(|_| {
            let receiver = receiver.clone();
            let data = Arc::clone(&data);
            let contig_idx = contig_idx.clone();

            thread::spawn(move || {
                
                while let Ok(df) = receiver.recv() {
                    let df = df.collect().unwrap();
                    let mut local_data: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
                    let cat_col1 = df.column("chrom1").unwrap().categorical().unwrap();
                    let rev_map1 = cat_col1.get_rev_map();
    
                    let cat_col2 = df.column("chrom2").unwrap().categorical().unwrap();
                    let rev_map2 = cat_col2.get_rev_map();
    
                    let nrows = df.height();
    
                    for idx in 0..nrows {
                        let row = df.get(idx).unwrap();
                        let chrom1 = match row.get(0) {
                            Some(AnyValue::Categorical(v, _, _)) => Some(v),
                            _ => None
                        };
    
                        let chrom2 = match row.get(1) {
                            Some(AnyValue::Categorical(v, _, _)) => Some(v),
                            _ => None
                        };
    
                        let count = match row.get(2) {
                            Some(AnyValue::UInt32(v)) => Some(v),
                            _ => None
                        };
    
                      
                        if let (Some(chrom1), Some(chrom2), Some(count)) = (chrom1, chrom2, count) {
                            let chrom1 = rev_map1.get(*chrom1);
                            let chrom2 = rev_map2.get(*chrom2);
                            let chrom1 = contig_idx.get(chrom1).unwrap();
                            let chrom2 = contig_idx.get(chrom2).unwrap();
                            
                            
                            local_data.entry((*chrom1, *chrom2)).or_insert(Vec::new()).push(*count);
                        }
                    }

                    let mut data = data.lock().unwrap();
                    for (key, value) in local_data {
                        data.entry(key).or_insert(Vec::new()).extend(value);
                    }
             
                }
            })
        }).collect();


        // for handle in producer_handles {
        //     handle.join().unwrap();
        // }
        
       
        for file in files.into_iter() {
            let df = match LazyFrame::scan_parquet(&file, ScanArgsParquet::default()) {
                Ok(df) => df,
                Err(e) => {
                    log::warn!("Empty file: {:?}", file);
                    continue
                }
            };

            let df = if min_mapq > 1 {
                df.filter(col("mapq").gt_eq(min_mapq))
            } else {
                df
            };

            let result = df.group_by(["chrom1", "chrom2"])
                .agg(&[col("pos1"), col("pos2")])
                .with_column(
                    col("pos1").arr().0.apply(
                        |s| {
                            let ca = s.list().unwrap();
                            let mut vec = Vec::with_capacity(ca.len());
                            for i in 0..ca.len() {
                                let val = ca.get(i).unwrap();
                                
                                vec.push(val.len() as u32);
                            }
                            

                            Ok(Some(Series::new("count1".into(), vec).into()))
                        },
                        GetOutput::from_type(DataType::UInt32)
                    ).alias("count")
                ).select(
                    &[
                        col("chrom1"),
                        col("chrom2"),
                        col("count")
                    ]
                );

            sender.send(result).unwrap();

        }

        drop(sender);

        for handle in consumer_handles {
            handle.join().unwrap();
        }
       
        let data = Arc::try_unwrap(data).unwrap().into_inner().unwrap();


        let mut wtr = common_writer(output.as_str());
        for (cp, vec) in data {
            let count = vec.iter().sum::<u32>();    
            if count < min_contacts {
                continue;
            }
            let contig1 = idx_contig.get(&cp.0).unwrap();
            let contig2 = idx_contig.get(&cp.1).unwrap();
            let buffer = format!("{}\t{}\t{}\n", contig1, contig2, count);
            
            wtr.write_all(buffer.as_bytes()).unwrap();
        }
       

        Ok(())
       
    }

    pub fn to_split_contacts(&self, min_contacts: u32, split_num: u32, min_quality: u8, output: &String) {
        polars::enable_string_cache(); 
        use hashbrown::HashMap;
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", 4));
        }
        
        let min_mapq = min_quality as u32;
       
        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u32>().unwrap();
                contigsizes.insert(contig, size);
        }

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (k.clone(), i as u32)).collect();
        let idx_contig: HashMap<u32, String, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (i as u32, k.clone())).collect();
        let idx_sizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), v.clone())).collect();

        let (sender, receiver) = bounded::<LazyFrame>(100);
        let data = Arc::new(Mutex::new(HashMap::new()));

        let consumer_handles: Vec<_> = (0..8).map(|_| {
            let receiver = receiver.clone();
            let data = Arc::clone(&data);
            let contig_idx = contig_idx.clone();

            thread::spawn(move || {
                
                while let Ok(df) = receiver.recv() {
                    let df = df.collect().unwrap();
                    let mut local_data: HashMap<(u32, u32), Vec<SmallIntVec>> = HashMap::new();
                    let cat_col1 = df.column("chrom1").unwrap().categorical().unwrap();
                    let rev_map1 = cat_col1.get_rev_map();
    
                    let cat_col2 = df.column("chrom2").unwrap().categorical().unwrap();
                    let rev_map2 = cat_col2.get_rev_map();
    
                    let nrows = df.height();
    
                    for idx in 0..nrows {
                        let row = df.get(idx).unwrap();
                        let chrom1 = match row.get(0) {
                            Some(AnyValue::Categorical(v, _, _)) => Some(v),
                            _ => None
                        };
    
                        let chrom2 = match row.get(1) {
                            Some(AnyValue::Categorical(v, _, _)) => Some(v),
                            _ => None
                        };
    

                        let pos1 = match row.get(2) {
                            Some(AnyValue::List(v)) => Some(v),
                            _ => None
                        };
    
                        let pos2 = match row.get(3) {
                            Some(AnyValue::List(v)) => Some(v),
                            _ => None
                        };
    
                        if let (Some(chrom1), Some(chrom2), Some(pos1), Some(pos2)) = (chrom1, chrom2, pos1, pos2) {
                            let chrom1 = rev_map1.get(*chrom1);
                            let chrom2 = rev_map2.get(*chrom2);
                            let chrom1 = contig_idx.get(chrom1).unwrap();
                            let chrom2 = contig_idx.get(chrom2).unwrap();
                            
                            let mut vec: Vec<SmallIntVec> = Vec::new();

                            for (p1, p2) in pos1.u32().unwrap().iter().zip(pos2.u32().unwrap().iter()) {
                                vec.push(smallvec![p1.unwrap(), p2.unwrap()]);
                            }
                            
                            local_data.entry((*chrom1, *chrom2)).or_insert(Vec::new()).extend(vec);
                        }
                    }

                    let mut data = data.lock().unwrap();
                    for (key, value) in local_data {
                        data.entry(key).or_insert(Vec::new()).extend(value);
                    }
             
                }
            })
        }).collect();

        for file in files.into_iter() {
            let df = match LazyFrame::scan_parquet(&file, ScanArgsParquet::default()) {
                Ok(df) => df,
                Err(e) => {
                    log::warn!("Empty file: {:?}", file);
                    continue
                }
            };

            let df = if min_mapq > 1 {
                df.filter(col("mapq").gt_eq(min_mapq))
            } else {
                df
            };

            let result = df.group_by(["chrom1", "chrom2"])
                .agg(&[col("pos1"), col("pos2")]);

            sender.send(result).unwrap();

        }

        drop(sender);
        for handle in consumer_handles {
            handle.join().unwrap();
        }
        let data = Arc::try_unwrap(data).unwrap().into_inner().unwrap();
        log::info!("Calculating the split contacts");

        let writer = common_writer(output.as_str());
        let writer = Arc::new(Mutex::new(writer));
        
        data.par_iter().for_each(|(cp, vec)| {
            let (c1, c2) = cp;
            let size1 = idx_sizes.get(c1).unwrap();
            let size2 = idx_sizes.get(c2).unwrap();
            let name1 = idx_contig.get(c1).unwrap();
            let name2 = idx_contig.get(c2).unwrap();

            let split_size1 = if *size1 > split_num { *size1 / split_num } else { 1 };
            let split_size2 = if *size2 > split_num { *size2 / split_num } else { 1 };

            let mut contact_hash = HashMap::new();
            
            for pair in vec {
                let p1 = pair[0];
                let p2 = pair[1];
                

                let idx1 = (p1 / split_size1).min(split_num - 1);
                let idx2 = (p2 / split_size2).min(split_num - 1);
                
                *contact_hash.entry((idx1, idx2)).or_insert(0) += 1;
            }
            
            let mut buffer = Vec::with_capacity(contact_hash.len());
            contact_hash.iter().for_each(|((idx1, idx2), count)| {
                if count >= &min_contacts {
                    buffer.push(format!("{}_{}\t{}_{}\t{}\n", name1, idx1, name2, idx2, count));
                }
            });
            
            let buffer = buffer.join("");
            let mut writer = writer.lock().unwrap();
            writer.write_all(buffer.as_bytes()).unwrap();

        });

        log::info!("Successful output split contacts file `{}`", output);
    }

    pub fn to_depth(&self, binsize: u32, min_quality: u8, output: &String) {
        use hashbrown::HashMap;
        polars::enable_string_cache(); 
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", 2));
        }
        
        let min_mapq = min_quality as u32;
        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u32>().unwrap();
                contigsizes.insert(contig, size);
        }

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (k.clone(), i as u32)).collect();
        let idx_contig: HashMap<u32, String, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (i as u32, k.clone())).collect();
        let idx_sizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), v.clone())).collect();


        let (sender, receiver) = bounded::<LazyFrame>(100);
        let data = Arc::new(Mutex::new(HashMap::new()));


        let consumer_handles: Vec<_> = (0..8).map(|_| {
            let receiver = receiver.clone();
            let data = Arc::clone(&data);
            let contig_idx = contig_idx.clone();

            thread::spawn(move || {
                
                while let Ok(df) = receiver.recv() {
                    let df = df.collect().unwrap();
                    let mut local_data: HashMap<u32, Vec<u32>> = HashMap::new();
                    let cat_col1 = df.column("chrom1").unwrap().categorical().unwrap();
                    let rev_map1 = cat_col1.get_rev_map();
    
                    let cat_col2 = df.column("chrom2").unwrap().categorical().unwrap();
                    let rev_map2 = cat_col2.get_rev_map();
    
                    let nrows = df.height();
    
                    for idx in 0..nrows {
                        let row = df.get(idx).unwrap();
                        let chrom1 = match row.get(0) {
                            Some(AnyValue::Categorical(v, _, _)) => Some(v),
                            _ => None
                        };
    
                        let chrom2 = match row.get(1) {
                            Some(AnyValue::Categorical(v, _, _)) => Some(v),
                            _ => None
                        };
    
                        let pos1 = match row.get(2) {
                            Some(AnyValue::List(v)) => Some(v),
                            _ => None
                        };
    
                        let pos2 = match row.get(3) {
                            Some(AnyValue::List(v)) => Some(v),
                            _ => None
                        };
    
                        if let (Some(chrom1), Some(chrom2), Some(pos1), Some(pos2)) = (chrom1, chrom2, pos1, pos2) {
                            let chrom1 = rev_map1.get(*chrom1);
                            let chrom2 = rev_map2.get(*chrom2);
                            let chrom1 = match contig_idx.get(chrom1) {
                                Some(v) => v,
                                None => continue,
                            };
                            let chrom2 = match contig_idx.get(chrom2) {
                                Some(v) => v,
                                None => continue,
                            };
                           
                            let mut vec1: Vec<u32> = Vec::new();
                            for p1 in pos1.u32().unwrap().iter() {
                                vec1.push(p1.unwrap());
                            }
                            
                            let mut vec2: Vec<u32> = Vec::new();
                            for p2 in pos2.u32().unwrap().iter() {
                                vec2.push(p2.unwrap());
                            }

                            local_data.entry(*chrom1).or_insert(Vec::new()).extend(vec1);
                            local_data.entry(*chrom2).or_insert(Vec::new()).extend(vec2);
                            
                        }
                    }

                    let mut data = data.lock().unwrap();
                    for (key, value) in local_data {
                        data.entry(key).or_insert(Vec::new()).extend(value);
                    }
             
                }
            })
        }).collect();


        // for handle in producer_handles {

        //     handle.join().unwrap();
            
        // }

        for file in files.into_iter() {
            let df = match LazyFrame::scan_parquet(&file, ScanArgsParquet::default()) {
                Ok(df) => df,
                Err(e) => {
                    log::warn!("Empty file: {:?}", file);
                    continue
                }
            };

            let df = if min_mapq > 1 {
                df.filter(col("mapq").gt_eq(min_mapq))
            } else {
                df
            };

            let result = df.with_column(
                    col("pos1") / binsize.into()
                )
                .with_column(
                    col("pos2") / binsize.into()
                )
                .group_by(["chrom1", "chrom2"])
                .agg(&[col("pos1"), col("pos2")]);

            sender.send(result).unwrap();
     
        }
    
        drop(sender);

        for handle in consumer_handles {
            handle.join().unwrap();
        }

        let data = Arc::try_unwrap(data).unwrap().into_inner().unwrap();

        log::info!("Calculating the depth of each contig");
        
        let mut depth: HashMap<u32, Vec<u32>> = idx_sizes
                                                    .clone()
                                                    .into_iter()
                                                    .map(|(chrom, size)| {
                                                        let num_bins = (size / binsize as u32 + 1) as usize;
                                                        (chrom, vec![0; num_bins])
                                                    }).collect();
        
        data.iter().for_each(|(cp, vec)| {
            vec.iter().for_each(|pos| {
                *depth.get_mut(cp).unwrap().get_mut(*pos as usize).unwrap() += 1;
               
            });

            
            
        });

        let depth: BTreeMap<_, _> = depth.into_iter().collect();
        let writer = common_writer(output);
        let wtr = Arc::new(Mutex::new(writer));

        depth.par_iter().for_each(|(contig, bins)| {
            let size = idx_sizes.get(contig).unwrap_or(&0);
            let contig = idx_contig.get(contig).unwrap();
            
            let mut buffer = Vec::with_capacity(bins.len() * 50);
            for (bin, count) in bins.iter().enumerate() {
                let bin_start = bin * binsize as usize;
                let mut bin_end = bin_start + binsize as usize;

                if bin_end > (*size).try_into().unwrap() {
                    bin_end = *size as usize;
                }
                if bin_start == bin_end {
                    continue;
                }
                buffer.extend_from_slice(format!("{}\t{}\t{}\t{}\n", contig, bin_start, bin_end, count).as_bytes());
            }
            {
                let mut wtr = wtr.lock().unwrap();
                wtr.write_all(&buffer).unwrap();
            }
            
        });

        log::info!("Successful output depth file `{}`", output);

    }

    pub fn break_contigs(&mut self, break_bed: &String, output: &String) {
        // enable_string_cache();
        polars::enable_string_cache(); 
        type IvString = Interval<usize, String>;

        let bed = Bed4::new(break_bed);
        let interval_hash = bed.to_interval_hash();
        let breaked_contigs = interval_hash.keys().map(|x| x.clone()).collect::<Vec<String>>();
        let breaked_series = Series::new("contig".into(), breaked_contigs);


        let files = collect_parquet_files(format!("{}/q0", self.file).as_str());

        let _ = std::fs::create_dir_all(output);
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));

        
        let _ = std::fs::copy(format!("{}/_metadata", self.file), format!("{}/_metadata", output));
        let _ = std::fs::copy(format!("{}/_readme", self.file), format!("{}/_readme", output));


        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u32>().unwrap();
                contigsizes.insert(contig, size);
        }

        let contigsizes_data: HashMap<String, u32> = contigsizes.clone();
        let mut new_contigsizes_data = contigsizes_data.clone();
        for contig in contigsizes_data.keys() {
            if interval_hash.contains_key(contig) {
                new_contigsizes_data.remove(contig);
                let interval = interval_hash.get(contig).unwrap();
                let res = interval.iter().collect::<Vec<_>>();
                for sub_res in res.iter() {
                    let new_contigs = &sub_res.val;
                    let size = sub_res.stop - sub_res.start + 1;
                    new_contigsizes_data.insert(new_contigs.clone(), size.try_into().unwrap());
                }
            }
        }


        let mut writer = common_writer(format!("{}/_contigsizes", output).as_str());
        for (contig, size) in new_contigsizes_data {
            let buffer = format!("{}\t{}\n", contig, size);
            writer.write_all(buffer.as_bytes()).unwrap();
        }

        let _ = writer.flush();

        let output_q_dir = format!("{}/q0", output);
        let copy_q_dir = format!("{}/q1", output);


        let results = files.chunks(500).for_each(|file_chunk| { file_chunk.into_par_iter().for_each(|file| {
            let df = LazyFrame::scan_parquet(file.clone(),  
                        ScanArgsParquet::default()).unwrap();
            
            let mut df = df.collect().unwrap();
            
            let cat_col1 = df.column("chrom1").unwrap().categorical().unwrap();
            let rev_map1 = cat_col1.get_rev_map();

            let cat_col2 = df.column("chrom2").unwrap().categorical().unwrap();
            let rev_map2 = cat_col2.get_rev_map();
          
            let cat_col = df.column("strand1").unwrap().categorical().unwrap();
            let rev_map3 = cat_col.get_rev_map();

            let cat_col = df.column("strand2").unwrap().categorical().unwrap();
            let rev_map4 = cat_col.get_rev_map();

            let phys_map1 = cat_col1.physical();
            let phys_map2 = cat_col2.physical();
            
            let mask1 = phys_map1.iter().map(|x| {
                let chrom = rev_map1.get(x.unwrap());
                interval_hash.contains_key(chrom)
            }).collect::<Vec<_>>();
            
            let mask2 = phys_map2.iter().map(|x| {
                let chrom = rev_map2.get(x.unwrap());
                interval_hash.contains_key(chrom)
            }).collect::<Vec<_>>();

            let mask_any = mask1.iter().zip(mask2.iter()).map(|(x, y)| *x || *y).collect::<Vec<_>>();
            // let is_filter = mask.iter().all(|x| !x);
            if mask_any.iter().all(|x| !x) {
                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let new_file = format!("{}/{}", output_q_dir, file_name);
                let mut new_file = File::create(new_file).unwrap();
                ParquetWriter::new(&mut new_file)
                    .finish(&mut df)
                    .unwrap();
            
                let mut df = df.lazy().filter(
                    col("mapq").gt_eq(1)
                ).collect().unwrap();
            
                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let new_file = format!("{}/{}", copy_q_dir, file_name);
                let mut new_file = File::create(new_file).unwrap();
                ParquetWriter::new(&mut new_file)
                    .finish(&mut df)
                    .unwrap();
            } else {
                let nrows = df.height();

                let mut data: Vec<bool> = Vec::new();
            
                let mut chrom1_vec = Vec::new();
                let mut chrom2_vec = Vec::new();
                let mut pos1_vec = Vec::new();
                let mut pos2_vec = Vec::new();
                let mut read_id_vec = Vec::new();
                let mut strand1_vec = Vec::new();
                let mut strand2_vec = Vec::new();
                let mut mapq_vec = Vec::new();


                for idx in 0..nrows {
                    let row = df.get(idx).unwrap();
                    let mut keep = true;
                    let chrom1 = match row.get(1) {
                        Some(AnyValue::Categorical(v, _, _)) => Some(v),
                        _ => None
                    };

                    let chrom2 = match row.get(3) {
                        Some(AnyValue::Categorical(v, _, _)) => Some(v),
                        _ => None
                    };

                    let pos1 = match row.get(2) {
                        Some(AnyValue::UInt32(v)) => Some(v),
                        _ => None
                    };
                    
                    let pos2 = match row.get(4) {
                        Some(AnyValue::UInt32(v)) => Some(v),
                        _ => None
                    };

                    if let (Some(chrom1), Some(chrom2), Some(pos1), Some(pos2)) = (chrom1, chrom2, pos1, pos2) {
                        let chrom1 = rev_map1.get(*chrom1);
                        let chrom2 = rev_map2.get(*chrom2);
                        let pos1 = *pos1 as usize;
                        let pos2 = *pos2 as usize;

                        let is_break_contig1 = interval_hash.contains_key(chrom1);
                        let is_break_contig2 = interval_hash.contains_key(chrom2);

                        if is_break_contig1 || is_break_contig2 {
                            
                            let strand1 = match row.get(5) {
                                Some(AnyValue::Categorical(v, _, _)) => Some(v),
                                _ => None
                            };
            
                            let strand2 = match row.get(6) {
                                Some(AnyValue::Categorical(v, _, _)) => Some(v),
                                _ => None
                            };
            
                            let mapq = match row.get(7) {
                                Some(AnyValue::UInt8(v)) => Some(v),
                                _ => None
                            };

                            let read_idx = match row.get(0) {
                                Some(AnyValue::String(v)) => Some(v),
                                _ => None
                            };
                            
                            let pos1 = pos1 as usize;
                            let pos2 = pos2 as usize;

                            let new1 = match is_break_contig1 {
                                true => {
                                    let interval = interval_hash.get(chrom1).unwrap();
                                    let res = interval.find(pos1, pos1 + 1).collect::<Vec<_>>();
                                    if res.len() > 0 {
                                        let new_pos = pos1 - res[0].start + 1;
                                        let new_chrom = res[0].val.clone();
                                        Some((new_chrom, new_pos as u32))
                                    } else {
                                        None
                                    }
                                },
                                false => Some((chrom1.to_owned(), pos1 as u32))
                            };

                            let new2 = match is_break_contig2 {
                                true => {
                                    let interval = interval_hash.get(chrom2).unwrap();
                                    let res = interval.find(pos2, pos2 + 1).collect::<Vec<_>>();
                                    if res.len() > 0 {
                                        let new_pos = pos2 - res[0].start + 1;
                                        let new_chrom = res[0].val.clone();
                                        Some((new_chrom, new_pos as u32 ))
                                    } else {
                                        None
                                    }
                                
                                },
                                false => Some((chrom2.to_owned(), pos2 as u32))
                            };
                        
                            if let (Some((new_chrom1, new_pos1)), Some((new_chrom2, new_pos2))) = (new1, new2) {

                                if let Some(read_idx) = read_idx {
                                    read_id_vec.push(read_idx.to_string());
                                }
                
                                if let (Some(strand1), Some(strand2), Some(mapq)) = (strand1, strand2, mapq) {
                                    strand1_vec.push(rev_map3.get(*strand1));
                                    strand2_vec.push(rev_map4.get(*strand2));
                                    mapq_vec.push(mapq.clone());
                                }
                
                                if new_chrom1 <= new_chrom2 {
                                    chrom1_vec.push(new_chrom1);
                                    pos1_vec.push(new_pos1);
                                    chrom2_vec.push(new_chrom2);
                                    pos2_vec.push(new_pos2);
                                } else {
                                    chrom1_vec.push(new_chrom2);
                                    pos1_vec.push(new_pos2);
                                    chrom2_vec.push(new_chrom1);
                                    pos2_vec.push(new_pos1);
                                }
                
                                keep = false;
                                
                            } else {
                                keep = true;
                            }
                            
                            
                        } else {
                            keep = true;
                        }
                    } else {
                        keep = true;
                    }
                    
                    data.push(keep);
                } 

                let df2 = df![
                    "read_idx" => read_id_vec,
                    "chrom1" => chrom1_vec,
                    "pos1" => pos1_vec,
                    "chrom2" => chrom2_vec,
                    "pos2" => pos2_vec,
                    "strand1" => strand1_vec,
                    "strand2" => strand2_vec,
                    "mapq" => mapq_vec
                ].unwrap();

                let data = Series::new("break".into(), data);
                let df = df.filter(
                    data.bool().unwrap()
                ).unwrap()
                .lazy()
                .with_column(col("chrom1").cast(DataType::String))
                .with_column(col("chrom2").cast(DataType::String))
                .with_column(col("strand1").cast(DataType::String))
                .with_column(col("strand2").cast(DataType::String))
                .collect().unwrap();

                let n_rows = df2.height();
                let df = if n_rows > 0 {
                    df.vstack(&df2).unwrap()
                } else {
                    df
                };

                let mut df = df.lazy().with_column(
                    col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                ).with_column(
                    col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                ).with_column(
                    col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                ).with_column(
                    col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                )
                .collect().unwrap();

                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let new_file = format!("{}/{}", output_q_dir, file_name);
                let mut new_file = File::create(new_file).unwrap();
                ParquetWriter::new(&mut new_file)
                    .finish(&mut df)
                    .unwrap();

                let mut df = df.lazy().filter(
                    col("mapq").gt_eq(1)
                ).collect().unwrap();

                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let new_file = format!("{}/{}", copy_q_dir, file_name);
                let mut new_file = File::create(new_file).unwrap();
                ParquetWriter::new(&mut new_file)
                    .finish(&mut df)
                    .unwrap();
            }


            });
                
        });

        let _ = copy_metadata_counts(self.file.as_str(), output.as_str());


    }

    pub fn intersect(&self, hcr_bed: &String, invert: bool,
                    min_mapq: u8, threads: usize,
                    max_q0_ratio: f64,
                    output: &String) -> anyResult<()> {
        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", threads));
        }
        
        
        enable_string_cache();  
        let bed = Bed3::new(hcr_bed);
        let interval_hash = Arc::new(bed.to_interval_hash());
        
        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        let mut q1_count: u64 = 0;
        let mut q0_only_count: u64 = 0;
        if min_mapq == 0 {
            if let Some((q1, q0_total)) = read_metadata_counts(self.file.as_str()) {
                q1_count = q1;
                q0_only_count = q0_total.saturating_sub(q1);
                log::info!(
                    "[intersect] Loaded counts: q1_records={}, q0_total_records={}, q0_only_records={}",
                    q1_count, q0_total, q0_only_count
                );
            } else {
                log::warn!("[intersect] Missing _metadata_counts; fallback to keep all q0-only (no sampling).");
            }
        }

        let q0_sample_prob: f64 = if min_mapq == 0 && q0_only_count > 0 {
            if max_q0_ratio <= 0.0 {
                1.0
            } else {
                let target = (q1_count as f64 * max_q0_ratio).min(q0_only_count as f64);
                (target / q0_only_count as f64).min(1.0)
            }
        } else {
            1.0
        };

        if min_mapq == 0 {
            log::info!(
                "[intersect] max_q0_ratio={}, q0_sample_prob={:.4}",
                max_q0_ratio, q0_sample_prob
            );
        }
        log::info!("Calculating the intersection of contacts with regions");

        let _ = std::fs::create_dir_all(output);
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));

        let _ = std::fs::copy(format!("{}/_contigsizes", self.file), format!("{}/_contigsizes", output));
        let _ = std::fs::copy(format!("{}/_metadata", self.file), format!("{}/_metadata", output));
        let _ = std::fs::copy(format!("{}/_readme", self.file), format!("{}/_readme", output));

        let q0_total = Arc::new(std::sync::atomic::AtomicU64::new(0));
        let q1_total = Arc::new(std::sync::atomic::AtomicU64::new(0));

        let output_q_dir = if min_mapq == 0 { format!("{}/q0", output) } else { format!("{}/q1", output) };
        let copy_q_dir = if min_mapq == 0 { format!("{}/q1", output) } else { format!("{}/q0", output) };

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .stack_size(16 * 1024 * 1024)
            .build()
            .unwrap();
        pool.install(|| {
            files.par_iter().for_each(|file| {

            let mut df = LazyFrame::scan_parquet(file.clone(), ScanArgsParquet::default()).unwrap();

            if min_mapq > 1 {
                df = df.filter(col("mapq").gt_eq(min_mapq));
            }

            
            let df_coords = df.clone()
                    .select([col("chrom1"), col("pos1"), col("chrom2"), col("pos2")])
                    .collect()
                    .unwrap();
                
            let nrows = df_coords.height();
            if nrows == 0 { return; }

            let chrom1_col = df_coords.column("chrom1").unwrap().categorical().unwrap();
            let chrom2_col = df_coords.column("chrom2").unwrap().categorical().unwrap();
            let pos1_col = df_coords.column("pos1").unwrap().u32().unwrap();
            let pos2_col = df_coords.column("pos2").unwrap().u32().unwrap();


            let rev_map1 = chrom1_col.get_rev_map();
            let rev_map2 = chrom2_col.get_rev_map();

            let mut lapper_cache1 = HashMap::new();
            if let Ok(unique_phys1) = chrom1_col.physical().unique() {
                for phys in unique_phys1.into_iter().flatten() {
                    lapper_cache1.insert(phys, interval_hash.get(rev_map1.get(phys)));
                }
            }

            let mut lapper_cache2 = HashMap::new();
            if let Ok(unique_phys2) = chrom2_col.physical().unique() {
                for phys in unique_phys2.into_iter().flatten() {
                    lapper_cache2.insert(phys, interval_hash.get(rev_map2.get(phys)));
                }
            }

            let mask: BooleanChunked = (0..nrows)
            .map(|i| {
                let c1 = chrom1_col.physical().get(i);
                let p1 = pos1_col.get(i);
                let c2 = chrom2_col.physical().get(i);
                let p2 = pos2_col.get(i);

                if let (Some(c1_idx), Some(p1_val), Some(c2_idx), Some(p2_val)) = (c1, p1, c2, p2) {
                    let p1_val = p1_val as usize;
                    let p2_val = p2_val as usize;

                    let is_in = lapper_cache1.get(&c1_idx).and_then(|&l| l).map_or(false, |l1| {
                        l1.find(p1_val, p1_val + 1).next().is_some() && 
                        lapper_cache2.get(&c2_idx).and_then(|&l| l).map_or(false, |l2| {
                            l2.find(p2_val, p2_val + 1).next().is_some()
                        })
                    });
                    is_in ^ invert
                } else {
                    false
                }
            })
            .collect();
            let mask = mask.with_name("mask".into());
                
            
            drop(df_coords);
            drop(lapper_cache1);
            drop(lapper_cache2);

            let full_df = df.collect().unwrap();
            let mut filtered_df = full_df.filter(&mask).unwrap();
            drop(full_df);

            // let mask_s = mask.into_series();
            // let mut lf = LazyFrame::scan_parquet(file.clone(), ScanArgsParquet::default()).unwrap();
            // let original_columns: Vec<Expr> = lf.collect_schema()
            //     .unwrap()
            //     .iter_names()
            //     .map(|n| col(n.as_str()))
            //     .collect();

            // if min_mapq > 1 {
            //     lf = lf.filter(col("mapq").gt_eq(min_mapq));
            // }

            // let mut filtered_df = lf
            //     .with_column(lit(mask_s).alias("mask"))
            //     .filter(col("mask"))
            //     .select(original_columns)
                
            //     .collect()
            //     .unwrap();

            // let mut df_full = LazyFrame::scan_parquet(file.clone(), ScanArgsParquet::default())
            //         .unwrap();
                
            // if min_mapq > 1 {
            //     df_full = df_full.filter(col("mapq").gt_eq(min_mapq));
            // }

            // let df_full = df_full.collect().unwrap();
            // let mut filtered_df = df_full.filter(&mask).unwrap();

 
            // drop(df_full);

            let q0_total_cl = Arc::clone(&q0_total);
            let q1_total_cl = Arc::clone(&q1_total);

                
            let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
            let main_output_path = format!("{}/{}", output_q_dir, file_name);
            let mut main_file = File::create(&main_output_path).unwrap();
            ParquetWriter::new(&mut main_file).finish(&mut filtered_df).unwrap();

            if min_mapq == 0 {

                let mut df_q1 = filtered_df
                    .filter(&filtered_df.column("mapq").unwrap().u8().unwrap().gt_eq(1))  // keep mapq>=1
                    .unwrap();

                let mut df_q0only = filtered_df
                    .filter(&filtered_df.column("mapq").unwrap().u8().unwrap().equal(0u8)) // keep mapq==0
                    .unwrap();

                if q0_sample_prob < 1.0 && df_q0only.height() > 0 {
                    df_q0only = sample_df_by_prob(df_q0only, q0_sample_prob, 0).unwrap();
                }

                let mut df_q0_out = if df_q0only.height() > 0 {
                    df_q1.vstack(&df_q0only).unwrap()
                } else {
                    df_q1.clone()
                };

                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let q0_out_path = format!("{}/{}", output_q_dir, file_name);
                let mut q0_file = File::create(&q0_out_path).unwrap();
                ParquetWriter::new(&mut q0_file).finish(&mut df_q0_out).unwrap();

                let q1_out_path = format!("{}/{}", copy_q_dir, file_name);
                let mut q1_file = File::create(&q1_out_path).unwrap();
                ParquetWriter::new(&mut q1_file).finish(&mut df_q1).unwrap();

                q0_total_cl.fetch_add(df_q0_out.height() as u64, std::sync::atomic::Ordering::Relaxed);
                q1_total_cl.fetch_add(df_q1.height() as u64, std::sync::atomic::Ordering::Relaxed);

              
            }  else {
             

                q1_total_cl.fetch_add(filtered_df.height() as u64, std::sync::atomic::Ordering::Relaxed);
                q0_total_cl.fetch_add(filtered_df.height() as u64, std::sync::atomic::Ordering::Relaxed);

                let q0_output_path = format!("{}/{}", copy_q_dir, file_name);
                let mut q0_file = File::create(q0_output_path).unwrap();
                ParquetWriter::new(&mut q0_file).finish(&mut filtered_df).unwrap();
            }
            // let new_file_path = format!("{}/{}", output_q_dir, file_name);
            // let mut new_file = File::create(new_file_path).unwrap();
            // ParquetWriter::new(&mut new_file).finish(&mut filtered_df).unwrap();
        });
        });
    
        
        // let files = collect_parquet_files(format!("{}", output_q_dir).as_str());
        // files.par_iter().for_each(|file| {
        //     let file_name = Path::new(file).file_name().unwrap().to_str().unwrap();
        //     let new_file = format!("{}/{}", copy_q_dir, file_name);
        //     let _ = std::fs::copy(file, new_file);
        // });

        let q0_n = q0_total.load(std::sync::atomic::Ordering::Relaxed);
        let q1_n = q1_total.load(std::sync::atomic::Ordering::Relaxed);
        write_metadata_counts(output.as_str(), q0_n, q1_n)?;

        log::info!("Successful output intersect file `{}`", output);
        Ok(())
    }

    pub fn dup_v0_3_0(&self, collapsed_list: &String, seed: usize, output: &String) {
        enable_string_cache();

        let output = match output.as_str() {
            "-" => {
                let mut output = self.file.clone();
                output.push_str("_dup");
                output
            },
            _ => output.clone()
        };

        let seed_bytes = seed.to_ne_bytes();
        let mut seed_array = [0u8; 32];
        for (i, byte) in seed_bytes.iter().enumerate() {
            seed_array[i] = *byte;
        }
        
        let seed_bytes2 = (seed * 2).to_ne_bytes();
        let mut seed_array2 = [0u8; 32];
        for (i, byte) in seed_bytes2.iter().enumerate() {
            seed_array2[i] = *byte;
        }

        let files = collect_parquet_files(format!("{}/q0", self.file).as_str());
        let _ = std::fs::create_dir_all(output.clone());
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));

        let _ = std::fs::copy(format!("{}/_contigsizes", self.file), format!("{}/_contigsizes", output));
        let _ = std::fs::copy(format!("{}/_metadata", self.file), format!("{}/_metadata", output));
        let _ = std::fs::copy(format!("{}/_readme", self.file), format!("{}/_readme", output));

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
            collapsed_contigs.entry(contig1.clone()).or_insert(vec![contig1]).push(contig2);
        }    

        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u32>().unwrap();
                contigsizes.insert(contig, size);
        }

       
        let mut wtr = common_writer(format!("{}/_contigsizes", output).as_str());
        for (contig, size) in contigsizes {
            if let Some(collapsed_contigs) = collapsed_contigs.get(&contig) {
                for collapsed_contig in collapsed_contigs {
                    let buffer = format!("{}\t{}\n", collapsed_contig, size);
                    wtr.write_all(buffer.as_bytes()).unwrap();
                }
            } else {
                let buffer = format!("{}\t{}\n", contig, size);
                wtr.write_all(buffer.as_bytes()).unwrap();
            }
            
        }
        
        let output_q_dir = format!("{}/q0", output);
        let copy_q_dir = format!("{}/q1", output);
        
        let chunksize = (files.len() / 10).max(1);
        files.par_chunks(chunksize).enumerate().for_each(|(chunk_idx, file_chunk)| {
            // log::info!("Processing chunk {}: {}", chunk_idx, file_chunk.len());
            let mut rng = StdRng::from_seed(seed_array);
            let mut rng2 = StdRng::from_seed(seed_array2);
            file_chunk.iter().for_each(|file| {
                
                let df = LazyFrame::scan_parquet(file.clone(),  ScanArgsParquet::default()).unwrap();
                let mut df = df.collect().unwrap();

                let c1_cat = df.column("chrom1").unwrap().categorical().unwrap();
                let c2_cat = df.column("chrom2").unwrap().categorical().unwrap();
                
                let rev1 = c1_cat.get_rev_map();
                let rev2 = c2_cat.get_rev_map();
                
                let phys1 = c1_cat.physical();
                let phys2 = c2_cat.physical();
                
                let max_id1 = phys1.max().unwrap_or(0) as usize;
                let max_id2 = phys2.max().unwrap_or(0) as usize;

                let mut lookup1: Vec<Option<&Vec<String>>> = vec![None; max_id1 + 1];
                for i in 0..=max_id1 {
                    if let Some(s) = rev1.get_optional(i as u32) {
                        lookup1[i] = collapsed_contigs.get(s);
                    }
                }
                
                let mut lookup2: Vec<Option<&Vec<String>>> = vec![None; max_id2 + 1];
                for i in 0..=max_id2 {
                    if let Some(s) = rev2.get_optional(i as u32) {
                        lookup2[i] = collapsed_contigs.get(s);
                    }
                }

                let nrows = df.height();
                let mut new_c1 = Vec::with_capacity(nrows);
                let mut new_c2 = Vec::with_capacity(nrows);

                for p in phys1.into_no_null_iter() {
                    let p_usize = p as usize;
                    let s = match lookup1[p_usize] {
                        Some(v) => v[rng.gen_range(0..v.len())].as_str(),
                        None => rev1.get(p),
                    };
                    new_c1.push(s);
                }

                for p in phys2.into_no_null_iter() {
                    let p_usize = p as usize;
                    let s = match lookup2[p_usize] {
                        Some(v) => v[rng2.gen_range(0..v.len())].as_str(),
                        None => rev2.get(p),
                    };
                    new_c2.push(s);
                }

                let chrom1_series = Series::new("chrom1".into(), new_c1);
                let chrom2_series = Series::new("chrom2".into(), new_c2);
                
                df.with_column(chrom1_series).unwrap();
                df.with_column(chrom2_series).unwrap();

                df = df.lazy().with_column(
                    col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                ).with_column(
                    col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                ).collect().unwrap();

                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let new_file = format!("{}/{}", output_q_dir, file_name);
                let mut new_file = File::create(new_file).unwrap();
                ParquetWriter::new(&mut new_file)
                    .finish(&mut df)
                    .unwrap();

                let mut df = df.lazy().filter(
                    col("mapq").gt_eq(1)
                ).collect().unwrap();

                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                
                let new_file = format!("{}/{}", copy_q_dir, file_name);
                let mut new_file = File::create(new_file).unwrap();
                ParquetWriter::new(&mut new_file)
                    .finish(&mut df)
                    .unwrap();
            
            })

        });

        let _ = copy_metadata_counts(self.file.as_str(), output.as_str());
        
        log::info!("Successful output a contig duplicated pairs file into {}", output);
    }


    pub fn dup(&self, collapsed_list: &String, seed: usize, output: &String) {
        enable_string_cache();

        let output = match output.as_str() {
            "-" => {
                let mut output = self.file.clone();
                output.push_str("_dup");
                output
            },
            _ => output.clone()
        };

        let seed_bytes = seed.to_ne_bytes();
        let mut seed_array = [0u8; 32];
        for (i, byte) in seed_bytes.iter().enumerate() {
            seed_array[i] = *byte;
        }
        
        let seed_bytes2 = (seed * 2).to_ne_bytes();
        let mut seed_array2 = [0u8; 32];
        for (i, byte) in seed_bytes2.iter().enumerate() {
            seed_array2[i] = *byte;
        }

        let files = collect_parquet_files(format!("{}/q0", self.file).as_str());
        let _ = std::fs::create_dir_all(output.clone());
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));

        let _ = std::fs::copy(format!("{}/_contigsizes", self.file), format!("{}/_contigsizes", output));
        let _ = std::fs::copy(format!("{}/_metadata", self.file), format!("{}/_metadata", output));
        let _ = std::fs::copy(format!("{}/_readme", self.file), format!("{}/_readme", output));

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
            collapsed_contigs.entry(contig1.clone()).or_insert(vec![contig1]).push(contig2);
        }    

        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = common_reader(&contigsize_file);
        let mut contigsizes = HashMap::new();
        for record in reader.lines() {
            let record = record.unwrap();
            let record = record.split("\t").collect::<Vec<&str>>();
            let contig = record.get(0).unwrap().to_string();
            let size = record.get(1).unwrap().parse::<u32>().unwrap();
                contigsizes.insert(contig, size);
        }

        
        let mut wtr = common_writer(format!("{}/_contigsizes", output).as_str());
        for (contig, size) in contigsizes {
            if let Some(collapsed_contigs) = collapsed_contigs.get(&contig) {
                for collapsed_contig in collapsed_contigs {
                    let buffer = format!("{}\t{}\n", collapsed_contig, size);
                    wtr.write_all(buffer.as_bytes()).unwrap();
                }
            } else {
                let buffer = format!("{}\t{}\n", contig, size);
                wtr.write_all(buffer.as_bytes()).unwrap();
            }
            
        }
        
        let output_q_dir = format!("{}/q0", output);
        let copy_q_dir = format!("{}/q1", output);
        
        // files.par_iter().enumerate().for_each(|(file_idx, file)| {
        //     let mut rng = StdRng::from_seed(seed_array);
        //     let mut rng2 = StdRng::from_seed(seed_array2);
            
        //     let df = LazyFrame::scan_parquet(file.clone(),  ScanArgsParquet::default()).unwrap();
        //     let mut df = df.collect().unwrap();

        //     let c1_cat = df.column("chrom1").unwrap().categorical().unwrap();
        //     let c2_cat = df.column("chrom2").unwrap().categorical().unwrap();
            
        //     let rev1 = c1_cat.get_rev_map();
        //     let rev2 = c2_cat.get_rev_map();
            
        //     let phys1 = c1_cat.physical();
        //     let phys2 = c2_cat.physical();
            
        //     let max_id1 = phys1.max().unwrap_or(0) as usize;
        //     let max_id2 = phys2.max().unwrap_or(0) as usize;

        //     let mut has_collapsed = false;
        //     for i in 0..=max_id1 {
        //         if let Some(s) = rev1.get_optional(i as u32) {
        //             if collapsed_contigs.contains_key(s) {
        //                 has_collapsed = true;
        //                 break;
        //             }
        //         }
        //     }
        //     if !has_collapsed {
        //         for i in 0..=max_id2 {
        //             if let Some(s) = rev2.get_optional(i as u32) {
        //                 if collapsed_contigs.contains_key(s) {
        //                     has_collapsed = true;
        //                     break;
        //                 }
        //             }
        //         }
        //     }

        //     let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
        //     let new_file = format!("{}/{}", output_q_dir, file_name);

        //     if !has_collapsed {
        //         let mut out_f = File::create(new_file).unwrap();
        //         ParquetWriter::new(&mut out_f).finish(&mut df).unwrap();

        //         let mut df_q1 = df.lazy().filter(col("mapq").gt_eq(1)).collect().unwrap();
        //         let copy_file = format!("{}/{}", copy_q_dir, file_name);
        //         let mut out_copy_f = File::create(copy_file).unwrap();
        //         ParquetWriter::new(&mut out_copy_f).finish(&mut df_q1).unwrap();
        //         return;
        //     }

        //     let mut lookup1: Vec<Option<&Vec<String>>> = vec![None; max_id1 + 1];
        //     for i in 0..=max_id1 {
        //         if let Some(s) = rev1.get_optional(i as u32) {
        //             lookup1[i] = collapsed_contigs.get(s);
        //         }
        //     }
            
        //     let mut lookup2: Vec<Option<&Vec<String>>> = vec![None; max_id2 + 1];
        //     for i in 0..=max_id2 {
        //         if let Some(s) = rev2.get_optional(i as u32) {
        //             lookup2[i] = collapsed_contigs.get(s);
        //         }
        //     }

        //     let lookup_fast1: Vec<&str> = (0..=max_id1).map(|i| rev1.get_optional(i as u32).unwrap_or("")).collect();
        //     let lookup_fast2: Vec<&str> = (0..=max_id2).map(|i| rev2.get_optional(i as u32).unwrap_or("")).collect();
        //     let nrows = df.height();
        //     let mut new_c1 = Vec::with_capacity(nrows);
        //     let mut new_c2 = Vec::with_capacity(nrows);

        //     for p in phys1.into_no_null_iter() {
        //         let p_usize = p as usize;
        //         let s = match lookup1[p_usize] {
        //             Some(v) => v[rng.gen_range(0..v.len())].as_str(),
        //             None => lookup_fast1[p_usize],
        //         };
        //         new_c1.push(s);
        //     }

        //     for p in phys2.into_no_null_iter() {
        //         let p_usize = p as usize;
        //         let s = match lookup2[p_usize] {
        //             Some(v) => v[rng2.gen_range(0..v.len())].as_str(),
        //             None => lookup_fast2[p_usize],
        //         };
        //         new_c2.push(s);
        //     }

        //     let chrom1_series = Series::new("chrom1".into(), new_c1);
        //     let chrom2_series = Series::new("chrom2".into(), new_c2);
            
        //     df.with_column(chrom1_series).unwrap();
        //     df.with_column(chrom2_series).unwrap();

        //     df = df.lazy().with_column(
        //         col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
        //     ).with_column(
        //         col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
        //     ).collect().unwrap();

        //     let mut out_f = File::create(new_file).unwrap();
        //     ParquetWriter::new(&mut out_f)
        //         .finish(&mut df)
        //         .unwrap();

        //     let mut df = df.lazy().filter(
        //         col("mapq").gt_eq(1)
        //     ).collect().unwrap();
            
        //     let copy_file = format!("{}/{}", copy_q_dir, file_name);
        //     let mut out_copy_f = File::create(copy_file).unwrap();
        //     ParquetWriter::new(&mut out_copy_f)
        //         .finish(&mut df)
        //         .unwrap();
        // });
        let chunksize = (files.len() / 8).max(1);
        files.par_chunks(chunksize).enumerate().for_each(|(chunk_idx, file_chunk)| {
            file_chunk.iter().for_each(|file| {
                let mut rng = StdRng::from_seed(seed_array);
                let mut rng2 = StdRng::from_seed(seed_array2);
                
                let df = LazyFrame::scan_parquet(file.clone(),  ScanArgsParquet::default()).unwrap();
                let mut df = df.collect().unwrap();

                let c1_cat = df.column("chrom1").unwrap().categorical().unwrap();
                let c2_cat = df.column("chrom2").unwrap().categorical().unwrap();
                
                let rev1 = c1_cat.get_rev_map();
                let rev2 = c2_cat.get_rev_map();
                
                let phys1 = c1_cat.physical();
                let phys2 = c2_cat.physical();
                
                let max_id1 = phys1.max().unwrap_or(0) as usize;
                let max_id2 = phys2.max().unwrap_or(0) as usize;

                let mut has_collapsed = false;
                for i in 0..=max_id1 {
                    if let Some(s) = rev1.get_optional(i as u32) {
                        if collapsed_contigs.contains_key(s) {
                            has_collapsed = true;
                            break;
                        }
                    }
                }
                if !has_collapsed {
                    for i in 0..=max_id2 {
                        if let Some(s) = rev2.get_optional(i as u32) {
                            if collapsed_contigs.contains_key(s) {
                                has_collapsed = true;
                                break;
                            }
                        }
                    }
                }

                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let new_file = format!("{}/{}", output_q_dir, file_name);

                if !has_collapsed {
                    let mut out_f = File::create(&new_file).unwrap();
                    ParquetWriter::new(&mut out_f)
                        .finish(&mut df)
                        .unwrap();

                    let mut df_q1 = df.lazy().filter(
                        col("mapq").gt_eq(1)
                    ).collect().unwrap();
                    
                    let copy_file = format!("{}/{}", copy_q_dir, file_name);
                    let mut out_copy_f = File::create(copy_file).unwrap();
                    ParquetWriter::new(&mut out_copy_f)
                        .finish(&mut df_q1)
                        .unwrap();
                    
                    return;
                }

                let mut lookup1: Vec<Option<&Vec<String>>> = vec![None; max_id1 + 1];
                for i in 0..=max_id1 {
                    if let Some(s) = rev1.get_optional(i as u32) {
                        lookup1[i] = collapsed_contigs.get(s);
                    }
                }
                
                let mut lookup2: Vec<Option<&Vec<String>>> = vec![None; max_id2 + 1];
                for i in 0..=max_id2 {
                    if let Some(s) = rev2.get_optional(i as u32) {
                        lookup2[i] = collapsed_contigs.get(s);
                    }
                }

                let lookup_fast1: Vec<&str> = (0..=max_id1).map(|i| rev1.get_optional(i as u32).unwrap_or("")).collect();
                let lookup_fast2: Vec<&str> = (0..=max_id2).map(|i| rev2.get_optional(i as u32).unwrap_or("")).collect();
                let nrows = df.height();
                let mut new_c1 = Vec::with_capacity(nrows);
                let mut new_c2 = Vec::with_capacity(nrows);

                for p in phys1.into_no_null_iter() {
                    let p_usize = p as usize;
                    let s = match lookup1[p_usize] {
                        Some(v) => v[rng.gen_range(0..v.len())].as_str(),
                        None => lookup_fast1[p_usize],
                    };
                    new_c1.push(s);
                }

                for p in phys2.into_no_null_iter() {
                    let p_usize = p as usize;
                    let s = match lookup2[p_usize] {
                        Some(v) => v[rng2.gen_range(0..v.len())].as_str(),
                        None => lookup_fast2[p_usize],
                    };
                    new_c2.push(s);
                }

                let chrom1_series = Series::new("chrom1".into(), new_c1);
                let chrom2_series = Series::new("chrom2".into(), new_c2);
                
                df.with_column(chrom1_series).unwrap();
                df.with_column(chrom2_series).unwrap();

                df = df.lazy().with_column(
                    col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                ).with_column(
                    col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical))
                ).collect().unwrap();

                let mut out_f = File::create(new_file).unwrap();
                ParquetWriter::new(&mut out_f)
                    .finish(&mut df)
                    .unwrap();

                let mut df = df.lazy().filter(
                    col("mapq").gt_eq(1)
                ).collect().unwrap();
                
                let copy_file = format!("{}/{}", copy_q_dir, file_name);
                let mut out_copy_f = File::create(copy_file).unwrap();
                ParquetWriter::new(&mut out_copy_f)
                    .finish(&mut df)
                    .unwrap();
            });
        });
        let _ = copy_metadata_counts(self.file.as_str(), output.as_str());
        
        log::info!("Successful output a contig duplicated pairs file into {}", output);
    }
    
    pub fn to_bam(&self, min_quality: u8, output: &String, threads: usize) -> anyResult<()> {
        use rust_htslib::bam::{self, Header, HeaderView, Record, Writer, header::HeaderRecord};
        use rust_htslib::bam::record::{Cigar, CigarString};
        use polars::prelude::*;
        use hashbrown::HashMap;

        polars::enable_string_cache(); 
        
        let min_mapq = min_quality as u32;
        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        log::info!("Parsing _contigsizes for BAM header ...");
        let contigsize_file = format!("{}/_contigsizes", self.file);
        let reader = crate::core::common_reader(&contigsize_file);
        
        let mut bam_header = Header::new();
        let mut hd_record = HeaderRecord::new(b"HD");
        hd_record.push_tag(b"VN", &"1.6");
        hd_record.push_tag(b"SO", &"unsorted");
        bam_header.push_record(&hd_record);

        let mut contig_idx: HashMap<String, i32> = HashMap::new();
        for (i, record) in reader.lines().enumerate() {
            let record = record.unwrap();
            let mut s = record.split('\t');
            if let (Some(contig), Some(size_str)) = (s.next(), s.next()) {
                let size: u64 = size_str.parse().unwrap_or(0);
                let mut sq_record = HeaderRecord::new(b"SQ");
                sq_record.push_tag(b"SN", contig);
                sq_record.push_tag(b"LN", &size);
                bam_header.push_record(&sq_record);
                contig_idx.insert(contig.to_string(), i as i32);
            }
        }

        let mut pg_record = HeaderRecord::new(b"PG");
        pg_record.push_tag(b"ID", &"pqs2bam");
        pg_record.push_tag(b"PN", &"cphasing-rs");
        pg_record.push_tag(b"VN", env!("CARGO_PKG_VERSION"));
        bam_header.push_record(&pg_record);

        let bam_header_view = HeaderView::from_header(&bam_header);
        let mut wtr = Writer::from_path(output, &bam_header, bam::Format::Bam).unwrap();
        let _ = wtr.set_threads(threads);

        let dummy_cigar = CigarString(vec![Cigar::Match(1)]);

        log::info!("Writing mapping records to {} ...", output);
        for file in files.into_iter() {
            let df = match LazyFrame::scan_parquet(&file, ScanArgsParquet::default()) {
                Ok(df) => df,
                Err(_) => {
                    log::warn!("Skip empty or unparsable file: {:?}", file);
                    continue;
                }
            };
            
            let df = if min_mapq > 1 {
                df.filter(col("mapq").gt_eq(min_mapq))
            } else {
                df
            }.collect().unwrap();

            let nrows = df.height();
            if nrows == 0 { continue; }

            let cat_col1 = df.column("chrom1").unwrap().categorical().unwrap();
            let rev_map1 = cat_col1.get_rev_map();
            
            let cat_col2 = df.column("chrom2").unwrap().categorical().unwrap();
            let rev_map2 = cat_col2.get_rev_map();
            
            let strand_col1 = df.column("strand1").unwrap().categorical().unwrap();
            let rev_strand1 = strand_col1.get_rev_map();
            
            let strand_col2 = df.column("strand2").unwrap().categorical().unwrap();
            let rev_strand2 = strand_col2.get_rev_map();

            for i in 0..nrows {
                let row = df.get(i).unwrap();
                let read_idx = match row.get(0) {
                    Some(AnyValue::String(v)) => v,
                    _ => continue
                };
                let chrom1 = match row.get(1) {
                    Some(AnyValue::Categorical(v, _, _)) => rev_map1.get(*v),
                    _ => continue
                };
                let pos1: i64 = match row.get(2) {
                    Some(AnyValue::UInt32(v)) => *v as i64,
                    Some(AnyValue::UInt64(v)) => *v as i64,
                    Some(AnyValue::Int32(v))  => *v as i64,
                    Some(AnyValue::Int64(v))  => *v as i64,
                    _ => continue
                };
                let chrom2 = match row.get(3) {
                    Some(AnyValue::Categorical(v, _, _)) => rev_map2.get(*v),
                    _ => continue
                };
                let pos2: i64 = match row.get(4) {
                    Some(AnyValue::UInt32(v)) => *v as i64,
                    Some(AnyValue::UInt64(v)) => *v as i64,
                    Some(AnyValue::Int32(v))  => *v as i64,
                    Some(AnyValue::Int64(v))  => *v as i64,
                    _ => continue
                };
                let strand1 = match row.get(5) {
                    Some(AnyValue::Categorical(v, _, _)) => rev_strand1.get(*v),
                    _ => continue
                };
                let strand2 = match row.get(6) {
                    Some(AnyValue::Categorical(v, _, _)) => rev_strand2.get(*v),
                    _ => continue
                };
                let mapq = match row.get(7) {
                    Some(AnyValue::UInt8(v)) => *v,
                    _ => 60
                };

                let tid1 = *contig_idx.get(chrom1).unwrap_or(&-1);
                let tid2 = *contig_idx.get(chrom2).unwrap_or(&-1);

                if tid1 == -1 || tid2 == -1 { continue; }

                let rev1 = strand1 == "-";
                let rev2 = strand2 == "-";
                
                // Pair-format is 1-based, however BAM pos is 0-based. Subtract 1 and clamp to at least 0.
                let p1_0based = (pos1 - 1).max(0);
                let p2_0based = (pos2 - 1).max(0);

                let mut rec1 = Record::new();
                let mut rec2 = Record::new();

                // Read 1 Setup
                rec1.set(read_idx.as_bytes(), Some(&dummy_cigar), b"N", b"!");
                rec1.set_tid(tid1);
                rec1.set_pos(p1_0based); 
                rec1.set_mtid(tid2);
                rec1.set_mpos(p2_0based);
                rec1.set_mapq(mapq);
                let flag1 = 0x01 | 0x40 | (if rev1 { 0x10 } else { 0 }) | (if rev2 { 0x20 } else { 0 });
                rec1.set_flags(flag1);

                // Read 2 Setup
                rec2.set(read_idx.as_bytes(), Some(&dummy_cigar), b"N", b"!");
                rec2.set_tid(tid2);
                rec2.set_pos(p2_0based); 
                rec2.set_mtid(tid1);
                rec2.set_mpos(p1_0based);
                rec2.set_mapq(mapq);
                let flag2 = 0x01 | 0x80 | (if rev2 { 0x10 } else { 0 }) | (if rev1 { 0x20 } else { 0 });
                rec2.set_flags(flag2);

                wtr.write(&rec1).unwrap();
                wtr.write(&rec2).unwrap();
            }
        }

        log::info!("Successful output bam file `{}`", output);
        Ok(())
    }
 
    pub fn to_porec(
        &self,
        min_quality: u8,
        output: &String,
        enzyme_bed: &String,
        min_edge_support: u32,
        min_clique_size: usize,
        max_component_size: usize,
        threads: usize,
    ) -> anyResult<()> {
        use polars::prelude::*;
        use rayon::prelude::*;
        use rust_lapper::{Interval, Lapper};

        enable_string_cache();

        unsafe {
            std::env::set_var("POLARS_MAX_THREADS", format!("{}", threads));
        }

        let (chrom_lappers, frag_id_strs) = load_enzyme_bed_to_lapper(enzyme_bed)?;

        let min_mapq = min_quality as u32;
        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        log::info!(
            "[to_porec] start: files={}, enzyme_bed={}, min_mapq={}, min_edge_support={}, min_clique_size={}, max_component_size={}",
            files.len(),
            enzyme_bed,
            min_mapq,
            min_edge_support,
            min_clique_size,
            max_component_size
        );

        type Node = u32;
        type EdgeKey = (Node, Node);

        let n_shards = 128usize;
        let shards_adj: Arc<Vec<Mutex<HashMap<Node, HashSet<Node>>>>> = Arc::new(
            (0..n_shards).map(|_| Mutex::new(HashMap::new())).collect()
        );
        let shards_sup: Arc<Vec<Mutex<HashMap<EdgeKey, u32>>>> = Arc::new(
            (0..n_shards).map(|_| Mutex::new(HashMap::new())).collect()
        );

        let shard_of = |a: Node, b: Node| -> usize {
            let mut h = XxHash64::default();
            a.hash(&mut h);
            b.hash(&mut h);
            (h.finish() as usize) % n_shards
        };

        files.par_iter().for_each(|file| {
            let mut lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                Ok(lf) => lf,
                Err(_) => return,
            };

            let mut has_mapq = lf.collect_schema().map(|s| s.contains("mapq")).unwrap_or(false);
            let mut lf = lf.select([
                col("chrom1"),
                col("pos1"),
                col("chrom2"),
                col("pos2"),
                if has_mapq { col("mapq") } else { lit(255u8).alias("mapq") },
            ]);

            if min_mapq > 0 {
                lf = lf.filter(col("mapq").cast(DataType::UInt32).gt_eq(lit(min_mapq)));
            }

            let df = match lf.collect() {
                Ok(df) => df,
                Err(_) => return,
            };
            if df.height() == 0 {
                return;
            }

            let chrom1 = match df.column("chrom1").and_then(|s| s.categorical()) {
                Ok(x) => x,
                Err(_) => return,
            };
            let chrom2 = match df.column("chrom2").and_then(|s| s.categorical()) {
                Ok(x) => x,
                Err(_) => return,
            };
            let pos1 = match df.column("pos1").and_then(|s| s.u32()) {
                Ok(x) => x,
                Err(_) => return,
            };
            let pos2 = match df.column("pos2").and_then(|s| s.u32()) {
                Ok(x) => x,
                Err(_) => return,
            };

            let rev1 = chrom1.get_rev_map();
            let rev2 = chrom2.get_rev_map();

            let pchrom1 = chrom1.physical();
            let pchrom2 = chrom2.physical();

            let mut local_edges: HashMap<EdgeKey, u32> = HashMap::new();
            let n = df.height();

            for i in 0..n {
                let c1p = pchrom1.get(i);
                let c2p = pchrom2.get(i);
                let p1 = pos1.get(i);
                let p2 = pos2.get(i);
                if c1p.is_none() || c2p.is_none() || p1.is_none() || p2.is_none() {
                    continue;
                }

                let c1 = rev1.get(c1p.unwrap());
                let c2 = rev2.get(c2p.unwrap());

                let p1_1 = p1.unwrap();
                let p2_1 = p2.unwrap();
                if p1_1 == 0 || p2_1 == 0 {
                    continue;
                }
                let p1_0 = p1_1 - 1;
                let p2_0 = p2_1 - 1;

                let f1 = map_to_frag_id(c1, p1_0, &chrom_lappers);
                let f2 = map_to_frag_id(c2, p2_0, &chrom_lappers);
                let (Some(f1), Some(f2)) = (f1, f2) else { continue; };
                if f1 == f2 {
                    continue;
                }

                let (a, b) = if f1 < f2 { (f1, f2) } else { (f2, f1) };
                *local_edges.entry((a, b)).or_insert(0) += 1;
            }

            for ((a, b), c) in local_edges {
                let sid = shard_of(a, b);
                {
                    let mut sup = shards_sup[sid].lock().unwrap();
                    *sup.entry((a, b)).or_insert(0) += c;
                }
            }
        });

        (0..n_shards).into_par_iter().for_each(|sid| {
            let sup = {
                let mut lock = shards_sup[sid].lock().unwrap();
                std::mem::take(&mut *lock)
            };
            if sup.is_empty() { return; }

            let mut adj_local: HashMap<Node, HashSet<Node>> = HashMap::new();
            for ((a, b), c) in sup {
                if c < min_edge_support { continue; }
                adj_local.entry(a).or_default().insert(b);
                adj_local.entry(b).or_default().insert(a);
            }

            for (node, neighs) in adj_local {
                let nid = (node as usize) % n_shards;
                let mut lock = shards_adj[nid].lock().unwrap();
                lock.entry(node).or_default().extend(neighs);
            }
        });

        let mut adj: HashMap<Node, HashSet<Node>> = HashMap::new();
        let shards_adj = Arc::try_unwrap(shards_adj).expect("adj Arc unwrap failed");
        for m in shards_adj {
            let inner = m.into_inner().unwrap();
            for (k, v) in inner {
                adj.entry(k).or_default().extend(v);
            }
        }

        log::info!("[to_porec] graph nodes={}, output={}", adj.len(), output);

        let components = connected_components(&adj, max_component_size);
        log::info!(
            "[to_porec] components={}, skipped_large={}",
            components.found,
            components.skipped_large
        );

        let mut wtr = common_writer(output.as_str());
        let wtr = Arc::new(Mutex::new(&mut wtr));

        components.comps.into_par_iter().enumerate().for_each(|(ci, comp)| {
            let mut comp_nodes = comp;
            comp_nodes.sort_unstable();

            let comp_set: HashSet<Node> = comp_nodes.iter().copied().collect();
            if comp_nodes.len() < min_clique_size { return; }

            let neigh = |v: Node| -> HashSet<Node> {
                adj.get(&v)
                    .map(|hs| hs.iter().filter(|&&u| comp_set.contains(&u)).copied().collect())
                    .unwrap_or_default()
            };

            let all_vertices: HashSet<Node> = comp_set.clone();
            let mut cliques: Vec<Vec<Node>> = Vec::new();

            bron_kerbosch_pivot(
                &mut Vec::new(),
                &all_vertices,
                &HashSet::new(),
                &neigh,
                &mut cliques,
                min_clique_size,
            );

            if cliques.is_empty() { return; }

            let mut out = String::new();
            for (ki, mut clique) in cliques.into_iter().enumerate() {
                clique.sort_unstable();
                if clique.len() < min_clique_size { continue; }
                
                let read_id = format!("v_read_comp{}_cl{}", ci, ki);
                out.push_str(&read_id);
                for &fid in &clique {
                    let name = frag_id_strs.get(fid as usize).map(|s| s.as_str()).unwrap_or("NA:0");
                    out.push_str(&format!("\t{}", name));
                }
                out.push('\n');
            }

            if !out.is_empty() {
                let mut lock = wtr.lock().unwrap();
                let _ = lock.write_all(out.as_bytes());
            }
        });
        {
            let mut w = wtr.lock().unwrap();
            w.flush()?;
        }

        log::info!("[to_porec] done: {}", output);
        Ok(())
    }
}

fn collect_parquet_files(parquets_dir: &str) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for entry in WalkDir::new(parquets_dir) {
        match entry {
            Ok(entry) => {
                if entry.path().is_file() {
                    let path = entry.path();
                    if let Some(ext) = path.extension() {

                        if ext == "parquet" || ext == "pq" {
                            files.push(path.to_path_buf());
                        }
                    }
                }
            }
            Err(e) => {
                eprintln!("Error reading directory entry: {}", e);
            }
        }
    }
    files
}


pub fn merge_pqs(input: Vec<&String>, output: &String) {

    let _ = std::fs::create_dir_all(output);
    let _ = std::fs::create_dir_all(format!("{}/q0", output));
    let _ = std::fs::create_dir_all(format!("{}/q1", output));

    let _ = std::fs::copy(format!("{}/_contigsizes", input[0]), format!("{}/_contigsizes", output));
    let _ = std::fs::copy(format!("{}/_metadata", input[0]), format!("{}/_metadata", output));
    let _ = std::fs::copy(format!("{}/_readme", input[0]), format!("{}/_readme", output));


    let mut q0_pq_files = Vec::new();
    let mut q1_pq_files = Vec::new();

    let mut q0_total: u64 = 0;
    let mut q1_total: u64 = 0;

    for file in input {
        let p = PQS::new(file);
        if !p.is_pqs() {
            log::warn!("{} is not a valid PQS directory, skipped.", file);
            continue;
        }

        if let Some((q1, q0)) = read_metadata_counts(file.as_str()) {
            q0_total = q0_total.saturating_add(q0);
            q1_total = q1_total.saturating_add(q1);
        } else {
            log::warn!("No _metadata_counts in {}, merged output will miss these counts unless regenerated.", file);
        }


        let q0 = format!("{}/q0", file);
        let q1 = format!("{}/q1", file);

        let q0_files = collect_parquet_files(q0.as_str());
        let q1_files = collect_parquet_files(q1.as_str());

        q0_pq_files.extend(q0_files);
        q1_pq_files.extend(q1_files);
    }

    rayon::join(
        || {
            q0_pq_files.par_iter().enumerate().for_each(|(idx, file)| {
                let output_file = format!("{}/q0/{}.parquet", output, idx);
                let _ = std::fs::copy(file, output_file);
            });
        },
        || {
            q1_pq_files.par_iter().enumerate().for_each(|(idx, file)| {
                let output_file = format!("{}/q1/{}.parquet", output, idx);
                let _ = std::fs::copy(file, output_file);
            });
        }
    );

    if q0_total > 0 || q1_total > 0 {
        if let Err(e) = write_metadata_counts(output.as_str(), q0_total, q1_total) {
            log::warn!("Failed to write _metadata_counts for merged PQS: {}", e);
        }
    } else {
        log::warn!("Merged PQS: no counts found from inputs; _metadata_counts not written.");
    }
    // q0_pq_files.into_par_iter().enumerate().for_each(|(idx, file)| {
    //     let output_file = format!("{}/q0/{}.parquet", output, idx);
    //     let _ = std::fs::copy(file, output_file);
    // });


    // q1_pq_files.into_par_iter().enumerate().for_each(|(idx, file)| {
    //     let output_file = format!("{}/q1/{}.parquet", output, idx);
    //     let _ = std::fs::copy(file, output_file);
    // });

}



#[derive(Debug, Clone)]
struct ContigMapRec {
    contig: String,
    start: u32,
    end: u32,
}

fn parse_contig_bed_to_lapper(bed: &str) -> anyResult<(HashMap<String, Lapper<u32, usize>>, Vec<ContigMapRec>, HashMap<String, u32>)> {
    // Return:
    // 1) per-chrom lapper of intervals -> index into `recs`
    // 2) recs: idx -> (contig, start, end)
    // 3) contig_sizes: contig -> length
    let f = common_reader(bed);
    let rdr = BufReader::new(f);

    let mut recs: Vec<ContigMapRec> = Vec::new();
    let mut chrom_to_intervals: HashMap<String, Vec<Interval<u32, usize>>> = HashMap::new();
    let mut contig_sizes: HashMap<String, u32> = HashMap::new();

    for (ln, line) in rdr.lines().enumerate() {
        let line = line?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') { continue; }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(anyhow::anyhow!("Invalid BED at line {}: expected 4 cols, got {}", ln + 1, fields.len()));
        }
        let chrom = fields[0].to_string();
        let start: u32 = fields[1].parse()?;
        let end: u32 = fields[2].parse()?;
        let contig = fields[3].to_string();
        if end <= start {
            return Err(anyhow::anyhow!("Invalid interval at line {}: start>=end ({}-{})", ln + 1, start, end));
        }

        let idx = recs.len();
        recs.push(ContigMapRec { contig: contig.clone(), start, end });

        chrom_to_intervals.entry(chrom)
            .or_default()
            .push(Interval { start, stop: end, val: idx });

        contig_sizes.insert(contig, end - start);
    }

    let mut chrom_lappers: HashMap<String, Lapper<u32, usize>> = HashMap::new();
    for (chrom, mut ivs) in chrom_to_intervals {
        ivs.sort_by_key(|iv| iv.start);
        chrom_lappers.insert(chrom, Lapper::new(ivs));
    }
    Ok((chrom_lappers, recs, contig_sizes))
}

#[inline]
fn map_chrpos_to_contig(
    chrom: &str,
    pos_1based: u32,
    chrom_lappers: &HashMap<String, Lapper<u32, usize>>,
    recs: &Vec<ContigMapRec>,
) -> Option<(String, u32)> {
    let p0 = pos_1based.checked_sub(1)?;
    let lap = chrom_lappers.get(chrom)?;
    let mut it = lap.find(p0, p0 + 1);
    let hit = it.next()?;

    // hit.val: usize (not a reference)
    let rec = &recs[hit.val];

    let contig_pos_1based = p0 - rec.start + 1;
    Some((rec.contig.clone(), contig_pos_1based))
}


#[derive(Debug, Clone)]
struct BedContigIv {
    start: u32,
    end: u32,
    contig_id: u32,
    contig_start: u32, // for converting chr pos to contig pos
}


fn parse_contig_bed_index(bed: &str) -> anyResult<(HashMap<String, Vec<BedContigIv>>, Vec<String>, HashMap<String, u32>)> {
    let f = common_reader(bed);
    let rdr = BufReader::new(f);

    let mut contig_name_to_id: HashMap<String, u32> = HashMap::new();
    let mut contigs: Vec<String> = Vec::new();
    let mut contig_sizes: HashMap<String, u32> = HashMap::new();

    let mut chrom_ivs: HashMap<String, Vec<BedContigIv>> = HashMap::new();

    for (ln, line) in rdr.lines().enumerate() {
        let line = line?;
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') { continue; }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(anyhow::anyhow!("Invalid BED at line {}: expected 4 cols, got {}", ln + 1, fields.len()));
        }
        let chrom = fields[0].to_string();
        let start: u32 = fields[1].parse()?;
        let end: u32 = fields[2].parse()?;
        let contig = fields[3].to_string();
        if end <= start {
            return Err(anyhow::anyhow!("Invalid interval at line {}: start>=end ({}-{})", ln + 1, start, end));
        }

        let contig_id = *contig_name_to_id.entry(contig.clone()).or_insert_with(|| {
            let id = contigs.len() as u32;
            contigs.push(contig.clone());
            id
        });

        chrom_ivs.entry(chrom)
            .or_default()
            .push(BedContigIv { start, end, contig_id, contig_start: start });

        contig_sizes.insert(contig, end - start);
    }

    // sort each chrom intervals by start for binary search
    for (_chrom, ivs) in chrom_ivs.iter_mut() {
        ivs.sort_by_key(|x| x.start);
    }

    Ok((chrom_ivs, contigs, contig_sizes))
}

#[inline]
fn map_chrpos_to_contig_fast(
    chrom: &str,
    pos_1based: u32,
    chrom_ivs: &HashMap<String, Vec<BedContigIv>>,
) -> Option<(u32, u32)> {
    let p0 = pos_1based.checked_sub(1)?;
    let ivs = chrom_ivs.get(chrom)?;
    // find rightmost interval with start <= p0
    let mut lo = 0usize;
    let mut hi = ivs.len();
    while lo < hi {
        let mid = (lo + hi) >> 1;
        if ivs[mid].start <= p0 { lo = mid + 1; } else { hi = mid; }
    }
    if lo == 0 { return None; }
    let iv = &ivs[lo - 1];
    if p0 < iv.end {
        let contig_pos_1based = p0 - iv.contig_start + 1;
        Some((iv.contig_id, contig_pos_1based))
    } else {
        None
    }
}

pub fn chr_pqs_to_contig_pqs(
    input_pqs: &String,
    contig_bed: &String,
    output: &String,
    threads: usize,
) -> anyResult<()> {
    enable_string_cache();

    let p = PQS::new(input_pqs);
    if !p.is_pqs() {
        return Err(anyhow::anyhow!("Input is not a valid PQS dir: {}", input_pqs));
    }

    let (chrom_lappers, recs, contig_sizes) = parse_contig_bed_to_lapper(contig_bed)?;

    // prepare output dirs
    std::fs::create_dir_all(output)?;
    std::fs::create_dir_all(format!("{}/q0", output))?;
    std::fs::create_dir_all(format!("{}/q1", output))?;
    let _ = std::fs::copy(format!("{}/_metadata", input_pqs), format!("{}/_metadata", output));
    let _ = std::fs::copy(format!("{}/_readme", input_pqs), format!("{}/_readme", output));

    // write new _contigsizes
    {
        let mut w = common_writer(format!("{}/_contigsizes", output).as_str());
        let mut names: Vec<_> = contig_sizes.keys().cloned().collect();
        names.sort();
        for name in names {
            let sz = contig_sizes.get(&name).unwrap();
            writeln!(w, "{}\t{}", name, sz)?;
        }
        w.flush()?;
    }

    // only read q0
    let q0_files = collect_parquet_files(format!("{}/q0", input_pqs).as_str());

    let q0_total = Arc::new(AtomicU64::new(0));
    let q1_total = Arc::new(AtomicU64::new(0));

    let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build()?;
    pool.install(|| -> anyResult<()> {
        q0_files.par_iter().enumerate().try_for_each(|(idx, file)| -> anyResult<()> {
            let lf = LazyFrame::scan_parquet(file, ScanArgsParquet::default())
                .map_err(|e| anyhow::anyhow!("scan_parquet failed: {:?}: {}", file, e))?;

            let df = lf.select([
                    col("read_idx"),
                    col("chrom1"),
                    col("pos1"),
                    col("chrom2"),
                    col("pos2"),
                    col("strand1"),
                    col("strand2"),
                    col("mapq"),
                ])
                .collect()
                .map_err(|e| anyhow::anyhow!("collect failed: {:?}: {}", file, e))?;

            // empty input -> write empty outputs
            if df.height() == 0 {
                let out_q0 = format!("{}/q0/{}.parquet", output, idx);
                let out_q1 = format!("{}/q1/{}.parquet", output, idx);
                let mut f0 = File::create(&out_q0)?;
                let mut f1 = File::create(&out_q1)?;
                let mut empty0 = df.clone();
                let mut empty1 = df;
                ParquetWriter::new(&mut f0).finish(&mut empty0)?;
                ParquetWriter::new(&mut f1).finish(&mut empty1)?;
                return Ok(());
            }

            // Convert to concrete types for fast loop
            let read_idx = df.column("read_idx")?.str()?;
            let chrom1 = df.column("chrom1")?.categorical()?;
            let pos1 = df.column("pos1")?.u32()?;
            let chrom2 = df.column("chrom2")?.categorical()?;
            let pos2 = df.column("pos2")?.u32()?;
            let strand1 = df.column("strand1")?.categorical()?;
            let strand2 = df.column("strand2")?.categorical()?;
            let mapq = df.column("mapq")?.u8()?;

            let rev_c1 = chrom1.get_rev_map();
            let rev_c2 = chrom2.get_rev_map();
            let rev_s1 = strand1.get_rev_map();
            let rev_s2 = strand2.get_rev_map();

            let n = df.height();
            let mut out_read: Vec<String> = Vec::with_capacity(n);
            let mut out_c1: Vec<String> = Vec::with_capacity(n);
            let mut out_p1: Vec<u32> = Vec::with_capacity(n);
            let mut out_c2: Vec<String> = Vec::with_capacity(n);
            let mut out_p2: Vec<u32> = Vec::with_capacity(n);
            let mut out_s1: Vec<String> = Vec::with_capacity(n);
            let mut out_s2: Vec<String> = Vec::with_capacity(n);
            let mut out_mq: Vec<u8> = Vec::with_capacity(n);

            for i in 0..n {
                let rid = read_idx.get(i);
                let c1_phys = chrom1.physical().get(i);
                let c2_phys = chrom2.physical().get(i);
                let p1 = pos1.get(i);
                let p2 = pos2.get(i);
                let s1_phys = strand1.physical().get(i);
                let s2_phys = strand2.physical().get(i);
                let mq = mapq.get(i);

                if rid.is_none() || c1_phys.is_none() || c2_phys.is_none() || p1.is_none() || p2.is_none()
                    || s1_phys.is_none() || s2_phys.is_none() || mq.is_none() {
                    continue;
                }

                let chrom1_name = rev_c1.get(c1_phys.unwrap());
                let chrom2_name = rev_c2.get(c2_phys.unwrap());
                let p1 = p1.unwrap();
                let p2 = p2.unwrap();

                let m1 = map_chrpos_to_contig(chrom1_name, p1, &chrom_lappers, &recs);
                let m2 = map_chrpos_to_contig(chrom2_name, p2, &chrom_lappers, &recs);
                if m1.is_none() || m2.is_none() {
                    continue;
                }
                let (nc1, np1) = m1.unwrap();
                let (nc2, np2) = m2.unwrap();

                let (fc1, fp1, fc2, fp2) = if nc1 <= nc2 {
                    (nc1, np1, nc2, np2)
                } else {
                    (nc2, np2, nc1, np1)
                };

                out_read.push(rid.unwrap().to_string());
                out_c1.push(fc1);
                out_p1.push(fp1);
                out_c2.push(fc2);
                out_p2.push(fp2);
                out_s1.push(rev_s1.get(s1_phys.unwrap()).to_string());
                out_s2.push(rev_s2.get(s2_phys.unwrap()).to_string());
                out_mq.push(mq.unwrap());
            }

            let mut out_df = df![
                "read_idx" => out_read,
                "chrom1" => out_c1,
                "pos1" => out_p1,
                "chrom2" => out_c2,
                "pos2" => out_p2,
                "strand1" => out_s1,
                "strand2" => out_s2,
                "mapq" => out_mq
            ]?;

            out_df = out_df.lazy()
                .with_column(col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .collect()?;

            // write q0
            let out_q0 = format!("{}/q0/{}.parquet", output, idx);
            let mut f0 = File::create(&out_q0)?;
            ParquetWriter::new(&mut f0).finish(&mut out_df)?;
            q0_total.fetch_add(out_df.height() as u64, Ordering::Relaxed);

            // derive and write q1 from q0 output
            let mut out_df_q1 = out_df
                .lazy()
                .filter(col("mapq").gt_eq(lit(1u8)))
                .collect()?;
            let out_q1 = format!("{}/q1/{}.parquet", output, idx);
            let mut f1 = File::create(&out_q1)?;
            ParquetWriter::new(&mut f1).finish(&mut out_df_q1)?;
            q1_total.fetch_add(out_df_q1.height() as u64, Ordering::Relaxed);

            Ok(())
        })?;
        Ok(())
    })?;

    write_metadata_counts(
        output.as_str(),
        q0_total.load(Ordering::Relaxed),
        q1_total.load(Ordering::Relaxed),
    )?;
    Ok(())
}


pub fn downsample_pqs(
    input_pqs: &String,
    output: &String,
    n: Option<u64>,
    prob: Option<f64>,
    seed: u64,
    min_mapq: u8,
    threads: usize,
) -> anyResult<()> {
    enable_string_cache();

    let p = PQS::new(input_pqs);
    if !p.is_pqs() {
        return Err(anyhow::anyhow!("Input is not a valid PQS dir: {}", input_pqs));
    }

    if n.is_none() && prob.is_none() {
        return Err(anyhow::anyhow!("downsample_pqs: either n or prob must be provided"));
    }
    if let Some(p) = prob {
        if !(0.0..=1.0).contains(&p) {
            return Err(anyhow::anyhow!("downsample_pqs: prob must be in [0,1], got {}", p));
        }
    }

    // prepare output dirs
    std::fs::create_dir_all(output)?;
    std::fs::create_dir_all(format!("{}/q0", output))?;
    std::fs::create_dir_all(format!("{}/q1", output))?;

    // copy metadata/readme/contigsizes
    let _ = std::fs::copy(format!("{}/_contigsizes", input_pqs), format!("{}/_contigsizes", output));
    let _ = std::fs::copy(format!("{}/_metadata", input_pqs), format!("{}/_metadata", output));
    let _ = std::fs::copy(format!("{}/_readme", input_pqs), format!("{}/_readme", output));

    let q0_files = collect_parquet_files(format!("{}/q0", input_pqs).as_str());
    if q0_files.is_empty() {
        return Err(anyhow::anyhow!("No parquet files found under {}/q0", input_pqs));
    }

    // ---------- Pass 1: count eligible rows ----------
    // Eligibility: (min_mapq==0) or (mapq >= min_mapq)
    let total_eligible = Arc::new(AtomicU64::new(0));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(16 * 1024 * 1024)
        .build()?;

    pool.install(|| {
        q0_files.par_iter().for_each(|file| {
            let lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                Ok(lf) => lf,
                Err(_) => return,
            };

            // Only need mapq column
            let df = match lf.select([col("mapq")]).collect() {
                Ok(df) => df,
                Err(_) => return,
            };
            if df.height() == 0 {
                return;
            }

            let mq = match df.column("mapq").ok().and_then(|s| s.u8().ok()) {
                Some(x) => x,
                None => return,
            };

            let mut c: u64 = 0;
            if min_mapq == 0 {
                c = mq.len() as u64;
            } else {
                for v in mq.into_no_null_iter() {
                    if v >= min_mapq {
                        c += 1;
                    }
                }
            }
            total_eligible.fetch_add(c, Ordering::Relaxed);
        });
    });

    let total_eligible = total_eligible.load(Ordering::Relaxed);
    if total_eligible == 0 {
        // still write empty outputs (and metadata_counts=0)
        for (i, _file) in q0_files.iter().enumerate() {
            let out_q0 = format!("{}/q0/{}.parquet", output, i);
            let out_q1 = format!("{}/q1/{}.parquet", output, i);
            let mut f0 = File::create(&out_q0)?;
            let mut f1 = File::create(&out_q1)?;
            let mut empty = DataFrame::empty();
            ParquetWriter::new(&mut f0).finish(&mut empty)?;
            ParquetWriter::new(&mut f1).finish(&mut empty)?;
        }
        write_metadata_counts(output.as_str(), 0, 0)?;
        return Ok(());
    }

    // ---------- Determine sampling parameters ----------
    // If `n` is set, compute keep_prob = min(1, n/total_eligible), then do prob sampling.
    // This is deterministic and parallel-friendly, but keeps approximately n rows (binomial variance).
    // For a strict "exactly n" sampler, we would need a reservoir across all parquet rows.
    let keep_prob: f64 = if let Some(n_keep) = n {
        if n_keep >= total_eligible {
            1.0
        } else {
            (n_keep as f64) / (total_eligible as f64)
        }
    } else {
        prob.unwrap_or(1.0)
    };

    log::info!(
        "[downsample] eligible_records={}, target_n={:?}, prob_arg={:?}, effective_keep_prob={:.6}, min_mapq={}",
        total_eligible, n, prob, keep_prob, min_mapq
    );

    // ---------- Pass 2: write downsampled PQS ----------
    let out_q0_total = Arc::new(AtomicU64::new(0));
    let out_q1_total = Arc::new(AtomicU64::new(0));

    pool.install(|| {
        q0_files.par_iter().enumerate().for_each(|(idx, file)| {
            let lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                Ok(lf) => lf,
                Err(_) => return,
            };

            // select minimal required columns (keep full schema if you prefer, but this is consistent with other converters)
            let df = match lf.select([
                col("read_idx"),
                col("chrom1"),
                col("pos1"),
                col("chrom2"),
                col("pos2"),
                col("strand1"),
                col("strand2"),
                col("mapq"),
            ]).collect() {
                Ok(df) => df,
                Err(_) => return,
            };

            if df.height() == 0 {
                // write empty
                let out_q0 = format!("{}/q0/{}.parquet", output, idx);
                let out_q1 = format!("{}/q1/{}.parquet", output, idx);
                let _ = File::create(&out_q0).and_then(|mut f| ParquetWriter::new(&mut f).finish(&mut DataFrame::empty()).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e)));
                let _ = File::create(&out_q1).and_then(|mut f| ParquetWriter::new(&mut f).finish(&mut DataFrame::empty()).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e)));
                return;
            }

            // We downsample by building a boolean mask, using per-row stable seeding:
            // Randomness is independent of row ordering and parallelism.
            //
            // key = (global_seed, file_idx, row_idx, mapq, chrom1_phys, chrom2_phys, pos1, pos2)
            // -> deterministic pseudo-random u64 -> uniform in [0,1)
            //
            // This avoids having to share RNG across threads.
            let mq = match df.column("mapq").ok().and_then(|s| s.u8().ok()) {
                Some(x) => x,
                None => return,
            };

            let chrom1_cat = df.column("chrom1").ok().and_then(|s| s.categorical().ok());
            let chrom2_cat = df.column("chrom2").ok().and_then(|s| s.categorical().ok());
            let pos1 = df.column("pos1").ok().and_then(|s| s.u32().ok());
            let pos2 = df.column("pos2").ok().and_then(|s| s.u32().ok());
            if chrom1_cat.is_none() || chrom2_cat.is_none() || pos1.is_none() || pos2.is_none() {
                return;
            }
            let chrom1_phys = chrom1_cat.unwrap().physical();
            let chrom2_phys = chrom2_cat.unwrap().physical();
            let pos1 = pos1.unwrap();
            let pos2 = pos2.unwrap();

            // Build mask
            let nrows = df.height();
            let mut mask: Vec<bool> = Vec::with_capacity(nrows);

            for i in 0..nrows {
                let mqv = mq.get(i);
                let c1 = chrom1_phys.get(i);
                let c2 = chrom2_phys.get(i);
                let p1 = pos1.get(i);
                let p2 = pos2.get(i);

                if mqv.is_none() || c1.is_none() || c2.is_none() || p1.is_none() || p2.is_none() {
                    mask.push(false);
                    continue;
                }
                let mqv = mqv.unwrap();
                if min_mapq > 0 && mqv < min_mapq {
                    mask.push(false);
                    continue;
                }
                if keep_prob >= 1.0 {
                    mask.push(true);
                    continue;
                }
                if keep_prob <= 0.0 {
                    mask.push(false);
                    continue;
                }

                // deterministic hash -> uniform
                let mut h = XxHash64::with_seed(seed);
                // mix file/row coordinates + some content to reduce accidental correlations
                h.write_u64(idx as u64);
                h.write_u64(i as u64);
                h.write_u8(mqv);
                h.write_u32(c1.unwrap());
                h.write_u32(c2.unwrap());
                h.write_u32(p1.unwrap());
                h.write_u32(p2.unwrap());
                let x = h.finish();
                // convert top 53 bits to f64 in [0,1)
                let u = ((x >> 11) as f64) / ((1u64 << 53) as f64);
                mask.push(u < keep_prob);
            }

            let mask = BooleanChunked::from_slice("mask".into(), &mask);
            let mut out_df = match df.filter(&mask) {
                Ok(df) => df,
                Err(_) => return,
            };

            // Ensure categorical back (may already be categorical)
            out_df = match out_df.lazy()
                .with_column(col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .collect()
            {
                Ok(df) => df,
                Err(_) => return,
            };

            // write q0
            let out_q0_path = format!("{}/q0/{}.parquet", output, idx);
            if let Ok(mut f0) = File::create(&out_q0_path) {
                let _ = ParquetWriter::new(&mut f0).finish(&mut out_df.clone());
            }
            let out_q0_n = out_df.height() as u64;

            // derive q1 as mapq>=1 subset of q0 output
            let mut out_df_q1 = match out_df.lazy().filter(col("mapq").gt_eq(lit(1u8))).collect() {
                Ok(df) => df,
                Err(_) => DataFrame::empty(),
            };
            let out_q1_path = format!("{}/q1/{}.parquet", output, idx);
            if let Ok(mut f1) = File::create(&out_q1_path) {
                let _ = ParquetWriter::new(&mut f1).finish(&mut out_df_q1);
            }

            out_q0_total.fetch_add(out_q0_n, Ordering::Relaxed);
           
            out_q1_total.fetch_add(out_df_q1.height() as u64, Ordering::Relaxed);
        });
    });

    write_metadata_counts(
        output.as_str(),
        out_q0_total.load(Ordering::Relaxed),
        out_q1_total.load(Ordering::Relaxed),
    )?;

    Ok(())
}


fn parse_prune_table(table_path: &str) -> anyResult<HashSet<(String, String)>> {
    let f = common_reader(table_path);
    let rdr = BufReader::new(f);
    let mut set = HashSet::new();
    for line in rdr.lines().flatten() {
        let s = line.trim();
        if s.is_empty() || s.starts_with('#') { continue; }
        let fields: Vec<&str> = s.split_whitespace().collect();
        if fields.len() < 2 { continue; }
        let c1 = fields[0].to_string();
        let c2 = fields[1].to_string();
        if c1 <= c2 {
            set.insert((c1, c2));
        } else {
            set.insert((c2, c1));
        }
    }
    Ok(set)
}

pub fn prune_pqs(
    input_pqs: &String,
    prune_table: &String,
    output: &String,
    threads: usize,
) -> anyResult<()> {
    enable_string_cache();

    let p = PQS::new(input_pqs);
    if !p.is_pqs() {
        return Err(anyhow::anyhow!("Input is not a valid PQS dir: {}", input_pqs));
    }

    let blacklist = parse_prune_table(prune_table)?;
    log::info!("Loaded {} contig pairs to prune from {}", blacklist.len(), prune_table);

    std::fs::create_dir_all(output)?;
    std::fs::create_dir_all(format!("{}/q0", output))?;
    std::fs::create_dir_all(format!("{}/q1", output))?;

    let _ = std::fs::copy(format!("{}/_contigsizes", input_pqs), format!("{}/_contigsizes", output));
    let _ = std::fs::copy(format!("{}/_metadata", input_pqs), format!("{}/_metadata", output));
    let _ = std::fs::copy(format!("{}/_readme", input_pqs), format!("{}/_readme", output));

    let q0_files = collect_parquet_files(format!("{}/q0", input_pqs).as_str());
    if q0_files.is_empty() {
        return Err(anyhow::anyhow!("No parquet files found under {}/q0", input_pqs));
    }

    let out_q0_total = Arc::new(AtomicU64::new(0));
    let out_q1_total = Arc::new(AtomicU64::new(0));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(16 * 1024 * 1024)
        .build()?;

    pool.install(|| {
        q0_files.par_iter().enumerate().for_each(|(idx, file)| {
            let lf = match LazyFrame::scan_parquet(file, ScanArgsParquet::default()) {
                Ok(lf) => lf,
                Err(_) => return,
            };

            let df = match lf.select([
                col("read_idx"),
                col("chrom1"),
                col("pos1"),
                col("chrom2"),
                col("pos2"),
                col("strand1"),
                col("strand2"),
                col("mapq"),
            ]).collect() {
                Ok(df) => df,
                Err(_) => return,
            };

            if df.height() == 0 {
                let out_q0 = format!("{}/q0/{}.parquet", output, idx);
                let out_q1 = format!("{}/q1/{}.parquet", output, idx);
                let _ = File::create(&out_q0).and_then(|mut f| ParquetWriter::new(&mut f).finish(&mut DataFrame::empty()).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e)));
                let _ = File::create(&out_q1).and_then(|mut f| ParquetWriter::new(&mut f).finish(&mut DataFrame::empty()).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e)));
                return;
            }

            let chrom1_cat = df.column("chrom1").unwrap().categorical().unwrap();
            let chrom2_cat = df.column("chrom2").unwrap().categorical().unwrap();

            let rev1 = chrom1_cat.get_rev_map();
            let rev2 = chrom2_cat.get_rev_map();

            let mut name1_to_phys = HashMap::new();
            let max_id1 = chrom1_cat.physical().max().unwrap_or(0) as usize;
            for i in 0..=max_id1 {
                if let Some(name) = rev1.get_optional(i as u32) {
                    name1_to_phys.insert(name, i as u32);
                }
            }

            let mut name2_to_phys = HashMap::new();
            let max_id2 = chrom2_cat.physical().max().unwrap_or(0) as usize;
            for i in 0..=max_id2 {
                if let Some(name) = rev2.get_optional(i as u32) {
                    name2_to_phys.insert(name, i as u32);
                }
            }

            let mut blacklist_phys = HashSet::new();
            for (c1, c2) in blacklist.iter() {
                if let (Some(&p1), Some(&p2)) = (name1_to_phys.get(c1.as_str()), name2_to_phys.get(c2.as_str())) {
                    blacklist_phys.insert((p1, p2));
                }
                if let (Some(&p1), Some(&p2)) = (name1_to_phys.get(c2.as_str()), name2_to_phys.get(c1.as_str())) {
                    blacklist_phys.insert((p1, p2));
                }
            }

            let nrows = df.height();
            let mut mask: Vec<bool> = Vec::with_capacity(nrows);

            let chrom1_phys = chrom1_cat.physical();
            let chrom2_phys = chrom2_cat.physical();

            for i in 0..nrows {
                let c1 = chrom1_phys.get(i);
                let c2 = chrom2_phys.get(i);

                if c1.is_none() || c2.is_none() {
                    mask.push(false);
                    continue;
                }
                let c1_val = c1.unwrap();
                let c2_val = c2.unwrap();

                mask.push(!blacklist_phys.contains(&(c1_val, c2_val)));
            }

            let mask = BooleanChunked::from_slice("mask".into(), &mask);
            let mut out_df = match df.filter(&mask) {
                Ok(df) => df,
                Err(_) => return,
            };

            out_df = match out_df.lazy()
                .with_column(col("chrom1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("chrom2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("strand1").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .with_column(col("strand2").cast(DataType::Categorical(None, CategoricalOrdering::Physical)))
                .collect()
            {
                Ok(df) => df,
                Err(_) => return,
            };

            let out_q0_path = format!("{}/q0/{}.parquet", output, idx);
            if let Ok(mut f0) = File::create(&out_q0_path) {
                let _ = ParquetWriter::new(&mut f0).finish(&mut out_df.clone());
            }
            let out_q0_n = out_df.height() as u64;

            let mut out_df_q1 = match out_df.lazy().filter(col("mapq").gt_eq(lit(1u8))).collect() {
                Ok(df) => df,
                Err(_) => DataFrame::empty(),
            };
            let out_q1_path = format!("{}/q1/{}.parquet", output, idx);
            if let Ok(mut f1) = File::create(&out_q1_path) {
                let _ = ParquetWriter::new(&mut f1).finish(&mut out_df_q1);
            }

            out_q0_total.fetch_add(out_q0_n, Ordering::Relaxed);
            out_q1_total.fetch_add(out_df_q1.height() as u64, Ordering::Relaxed);
        });
    });

    write_metadata_counts(
        output.as_str(),
        out_q0_total.load(Ordering::Relaxed),
        out_q1_total.load(Ordering::Relaxed),
    )?;

    log::info!("Successful output pruned PQS to `{}`", output);
    Ok(())
}