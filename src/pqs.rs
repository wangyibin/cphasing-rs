#[allow(unused)]
use anyhow::Result as anyResult;

use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use log::LevelFilter;
use std::collections::{ BTreeMap, HashMap, HashSet };
use std::hash::BuildHasherDefault;
use std::borrow::Cow;
use std::path::{ Path, PathBuf };
use walkdir::WalkDir;
use smallvec::{ smallvec, SmallVec };
use std::thread;
use std::io::{ BufReader, BufRead, BufWriter, Write };
use std::sync::{ Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::fs::{ File, OpenOptions };
use twox_hash::XxHash64;
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

    pub fn to_clm(&self, 
        min_contacts: u32,
        min_quality: u8,
        output: &String,
        output_split_contacts: bool,
        output_depth: bool,
        binsize: u32,
        threads: usize) -> anyResult<()> {
        // use arrow::array::UInt32Array;
        use hashbrown::HashMap;
        use std::fmt::Write;


        std::env::set_var("POLARS_MAX_THREADS", format!("{}", threads));
        let min_mapq = min_quality as u32;
        // get prefix of parquet_dir
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

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = 
                contigsizes.keys().enumerate().map(|(i, k)| (k.clone(), i as u32)).collect();
        let idx_contig: HashMap<u32, String, BuildHasherDefault<XxHash64>> = 
                contigsizes.keys().enumerate().map(|(i, k)| (i as u32, k.clone())).collect();
        let idx_sizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = 
                contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), v.clone())).collect();
        let idx_contig_sizes: HashMap<u32, (String, u32), BuildHasherDefault<XxHash64>> = 
                contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), (k.clone(), v.clone()))).collect();

        let (sender, receiver) = bounded(100);
        let data = Arc::new(Mutex::new(HashMap::new()));

        let producer_handles: Vec<_> = files.into_par_iter().map(|file| {
            let sender = sender.clone();
            thread::spawn(move || {
                let df = match LazyFrame::scan_parquet(&file, ScanArgsParquet::default()) {
                    Ok(df) => df,
                    Err(e) => {
                        log::warn!("Empty file: {:?}", file);
                        return;
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
            })
        }).collect();

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
                    
                    // let phys_map1 = cat_col1.physical();
                    // let phys_map2 = cat_col2.physical();
                    // let pos1_arr = df.column("pos1").unwrap().list().unwrap();
                    // let pos2_arr = df.column("pos2").unwrap().list().unwrap();

                    // for idx in 0..nrows {
                    //     let chrom1 = phys_map1.get(idx);
                    //     let chrom1 = rev_map1.get(chrom1.unwrap());

                    //     let chrom2 = phys_map2.get(idx);
                    //     let chrom2 = rev_map2.get(chrom2.unwrap());
                    //     let pos1 = pos1_arr.get(idx).unwrap();
                    //     let pos2 = pos2_arr.get(idx).unwrap();
                    
                        
                    //     let chrom1 = match contig_idx.get(chrom1) {
                    //         Some(v) => v,
                    //         None => {
                    //             log::warn!("Could not found {:?} in _contigsizes", chrom1);
                    //             continue
                    //         },
                    //     };
                    //     let chrom2 = match contig_idx.get(chrom2) {
                    //         Some(v) => v,
                    //         None => {
                    //             log::warn!("Could not found {:?} in _contigsizes", chrom2);
                    //             continue
                    //         },
                    //     };
                        
                    //     let mut vec: Vec<SmallIntVec> = Vec::new();
                        
                    //     let array1 = pos1.as_any().downcast_ref::<UInt32Array>().unwrap();
                    //     let array2= pos2.as_any().downcast_ref::<UInt32Array>().unwrap();
                    //     for (p1, p2) in array1.iter().zip(array2.iter()) {
                    //         vec.push(smallvec![p1.unwrap(), p2.unwrap()]);
                    //     }
                    //     local_data.entry((*chrom1, *chrom2)).or_insert(Vec::new()).extend(vec);

                    // }

                    // let local_data = (0..nrows).into_par_iter().filter_map(|idx| {
                    //     let row = df.get(idx).unwrap();
                    //     let chrom1 = match row.get(0) {
                    //         Some(AnyValue::Categorical(v, _, _)) => Some(v),
                    //         _ => None
                    //     };
    
                    //     let chrom2 = match row.get(1) {
                    //         Some(AnyValue::Categorical(v, _, _)) => Some(v),
                    //         _ => None
                    //     };
    
                    //     let pos1 = match row.get(2) {
                    //         Some(AnyValue::List(v)) => Some(v),
                    //         _ => None
                    //     };
    
                    //     let pos2 = match row.get(3) {
                    //         Some(AnyValue::List(v)) => Some(v),
                    //         _ => None
                    //     };
    
                    //     if let (Some(chrom1), Some(chrom2), Some(pos1), Some(pos2)) = (chrom1, chrom2, pos1, pos2) {
                    //         let chrom1 = rev_map1.get(*chrom1);
                    //         let chrom2 = rev_map2.get(*chrom2);
                    //         let chrom1 = match contig_idx.get(chrom1) {
                    //             Some(v) => v,
                    //             None => {
                    //                 log::warn!("Could not found {:?} in _contigsizes", chrom1);
                    //                 return None
                    //             },
                    //         };
                    //         let chrom2 = match contig_idx.get(chrom2) {
                    //             Some(v) => v,
                    //             None => {
                    //                 log::warn!("Could not found {:?} in _contigsizes", chrom2);
                    //                 return None
                    //             },
                    //         };
                    //         let mut vec: Vec<SmallIntVec> = Vec::new();
                    //         for (p1, p2) in pos1.u32().unwrap().iter().zip(pos2.u32().unwrap().iter()) {
                    //             vec.push(smallvec![p1.unwrap(), p2.unwrap()]);
                    //         }

                    //         Some((*chrom1, *chrom2, vec))
                    //     } else {
                    //         None
                    //     }
                    // }).collect::<Vec<_>>()
                    // let mut data = data.lock().unwrap();
                    // for (chrom1, chrom2, vec) in local_data {
                    //     data.entry((chrom1, chrom2)).or_insert(Vec::new()).extend(vec);
                    // }

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
                                None => {
                                    log::warn!("Could not found {:?} in _contigsizes", chrom1);
                                    continue
                                },
                            };
                            let chrom2 = match contig_idx.get(chrom2) {
                                Some(v) => v,
                                None => {
                                    log::warn!("Could not found {:?} in _contigsizes", chrom2);
                                    continue
                                },
                            };
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

        for handle in producer_handles {
            handle.join().unwrap();
        }
    
        drop(sender);

        for handle in consumer_handles {
            handle.join().unwrap();
        }
    
        let data = Arc::try_unwrap(data).unwrap().into_inner().unwrap();

        if output_split_contacts{
            log::info!("Calculating the distance between split contigs");
            let contacts = Contacts::new(&format!(".{}.pixels", output_prefix));
            let split_contigsizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = idx_sizes
                .par_iter()
                .map(|(k, v)| (*k, *v / 2))
                .collect(); 
    
            let writer = common_writer(format!("{}.split.contacts", output_prefix.to_string()).as_str());
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
                                &format!("{}.split.contacts", output_prefix.to_string()));

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

        let wtr = common_writer(output.as_str());
        let wtr = Arc::new(Mutex::new(wtr));
        data.par_iter().for_each(|(cp, vec)| {
            if cp.0 == cp.1 {
                return;
            }

            let (contig1, length1) = idx_contig_sizes.get(&cp.0).unwrap();
            let (contig2, length2) = idx_contig_sizes.get(&cp.1).unwrap();

            let res = vec.iter().map(
                |x| {
                    let pos1 = x[0];
                    let pos2 = x[1];
                    smallvec![
                        length1.wrapping_sub(pos1).wrapping_add(pos2),                      // ctg1+ ctg2+
                        length1.wrapping_sub(pos1).wrapping_add(*length2).wrapping_sub(pos2), // ctg1+ ctg2-
                        pos1.wrapping_add(pos2),                                             // ctg1- ctg2+
                        pos1.wrapping_add(*length2).wrapping_sub(pos2),                      // ctg1- ctg2-
                    ]
                }
            ).collect::<Vec<_>>();

            let zipped: Vec<Vec<u32>> = res[0].iter().enumerate().map(|(i, _)| {
                res.iter().map(|x: &SmallIntVec4| x[i]).collect::<Vec<_>>()
            }).collect::<Vec<_>>();

            let count = zipped[0].len();
            if count < min_contacts as usize {
                return
            }

            let mut buffer = Vec::with_capacity(4);
            for (i, res1) in zipped.iter().enumerate() {
                let res1 = res1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ");
                let line = match i {
                    0 => format!("{}+ {}+\t{}\t{}\n", contig1, contig2, count, res1),
                    1 => format!("{}+ {}-\t{}\t{}\n", contig1, contig2, count, res1),
                    2 => format!("{}- {}+\t{}\t{}\n", contig1, contig2, count, res1),
                    3 => format!("{}- {}-\t{}\t{}\n", contig1, contig2, count, res1),
                    _ => panic!("Error: Invalidindex"),
                };

                buffer.extend_from_slice(line.as_bytes());
            }
            {
                let mut wtr = wtr.lock().unwrap();
                wtr.write_all(&buffer).unwrap();
            }

        });


        log::info!("Successful output clm file `{}`", output);

        let mut wtr = wtr.lock().unwrap();
        wtr.flush().unwrap();

        Ok(())
    }

    

    pub fn to_mnd(&self, min_quality: u8, output: &String) -> anyResult<()> {
        std::env::set_var("POLARS_MAX_THREADS", format!("{}", 10));
        let min_mapq = min_quality as u32;
        // get prefix of parquet_dir
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


        let mut results = files.into_par_iter().map(|file| {
            let mut df = LazyFrame::scan_parquet(file,  ScanArgsParquet::default()).unwrap();

            if min_mapq > 1 {
                df = df.clone().filter(col("mapq").gt_eq(min_mapq));
                
            }
           
          

            let result = df.select(
                &[
                    lit(1i32).alias("strand1"),
                    col("chrom1"),
                    col("pos1"),
                    lit(0u32).alias("frag1"),
                    lit(1i32).alias("strand2"),
                    col("chrom2"),
                    col("pos2"),
                    lit(0u32).alias("frag2"),
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
        let mut df = results.remove(0).collect().unwrap();
        CsvWriter::new(&mut file)
            .include_header(false)
            .with_separator(b' ')
            .finish(&mut df)
            .unwrap();
       

        for i in 1..results.len() {
            let mut file = OpenOptions::new().append(true).open(output.as_str()).unwrap();
            let mut df = results.remove(0).collect().unwrap();
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
        std::env::set_var("POLARS_MAX_THREADS", format!("{}", 4));
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

        let (sender, receiver) = bounded(100);
        let data = Arc::new(Mutex::new(HashMap::new()));

        let producer_handles: Vec<_> = files.into_iter().map(|file| {
            let sender = sender.clone();
            thread::spawn(move || {
                let df = match LazyFrame::scan_parquet(&file, ScanArgsParquet::default()) {
                    Ok(df) => df,
                    Err(e) => {
                        log::warn!("Empty file: {:?}", file);
                        return;
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
                                

                                Ok(Some(Series::new("count1", vec)))
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
            })
        }).collect();

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


        for handle in producer_handles {
            handle.join().unwrap();
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

    pub fn to_depth(&self, binsize: u32, min_quality: u8, output: &String) {
        use hashbrown::HashMap;
        std::env::set_var("POLARS_MAX_THREADS", format!("{}", 2));
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


        let (sender, receiver) = bounded(100);
        let data = Arc::new(Mutex::new(HashMap::new()));

        let producer_handles: Vec<_> = files.into_iter().map(|file| {
            let sender = sender.clone();
            thread::spawn(move || {
                let df = match LazyFrame::scan_parquet(&file, ScanArgsParquet::default()) {
                    Ok(df) => df,
                    Err(e) => {
                        log::warn!("Empty file: {:?}", file);
                        return;
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
            })
        }).collect();

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
                            let chrom1 = contig_idx.get(chrom1).unwrap();
                            let chrom2 = contig_idx.get(chrom2).unwrap();
                           
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


        for handle in producer_handles {
            handle.join().unwrap();
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
        enable_string_cache();
        type IvString = Interval<usize, String>;

        let bed = Bed4::new(break_bed);
        let interval_hash = bed.to_interval_hash();
        

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

        let output_q_dir = format!("{}/q0", output);
        let copy_q_dir = format!("{}/q1", output);

        let results = files.into_par_iter().for_each(|file| {
            let df = LazyFrame::scan_parquet(file.clone(),  
                    ScanArgsParquet::default()).unwrap();
            
            let df = df.collect().unwrap();

            let cat_col = df.column("chrom1").unwrap().categorical().unwrap();
            let rev_map1 = cat_col.get_rev_map();

            let cat_col = df.column("chrom2").unwrap().categorical().unwrap();
            let rev_map2 = cat_col.get_rev_map();

            let cat_col = df.column("strand1").unwrap().categorical().unwrap();
            let rev_map3 = cat_col.get_rev_map();

            let cat_col = df.column("strand2").unwrap().categorical().unwrap();
            let rev_map4 = cat_col.get_rev_map();

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
                            
                           
                        } else {
                            continue
                        }
                        
                        if let Some(read_idx) = read_idx {
                                read_id_vec.push(read_idx.clone());
                        }

                        if let (Some(strand1), Some(strand2), Some(mapq)) = (strand1, strand2, mapq) {
                            strand1_vec.push(rev_map3.get(*strand1));
                            strand2_vec.push(rev_map4.get(*strand2));
                            mapq_vec.push(mapq.clone());
                        }

                        data.push(false);
                    } else {
                        data.push(true);
                    }
                }

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

            let data = Series::new("break", data);
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
        });

    }

    pub fn intersect(&self, hcr_bed: &String, invert: bool,
                    min_mapq: u8, output: &String) -> anyResult<()> {

        std::env::set_var("POLARS_MAX_THREADS", format!("{}", 10));
        enable_string_cache();   
        let bed = Bed3::new(hcr_bed);
        let interval_hash = bed.to_interval_hash();
        
        let files = if min_mapq == 0 {
            collect_parquet_files(format!("{}/q0", self.file).as_str())
        } else {
            collect_parquet_files(format!("{}/q1", self.file).as_str())
        };

        log::info!("Calculating the intersection of contacts with regions");

        let _ = std::fs::create_dir_all(output);
        let _ = std::fs::create_dir_all(format!("{}/q0", output));
        let _ = std::fs::create_dir_all(format!("{}/q1", output));

        let _ = std::fs::copy(format!("{}/_contigsizes", self.file), format!("{}/_contigsizes", output));
        let _ = std::fs::copy(format!("{}/_metadata", self.file), format!("{}/_metadata", output));
        let _ = std::fs::copy(format!("{}/_readme", self.file), format!("{}/_readme", output));

        let output_q_dir = if min_mapq == 0 {
            format!("{}/q0", output)
        } else {
            format!("{}/q1", output)
        };

        let copy_q_dir = if min_mapq == 0 {
            format!("{}/q1", output)
        } else {
            format!("{}/q0", output)
        };

        
        files.chunks(500).for_each(|file_chunk| {
            file_chunk.par_iter().for_each(|file| {
                let mut df = LazyFrame::scan_parquet(file.clone(),  ScanArgsParquet::default()).unwrap();

                if min_mapq > 1 {
                    df = df.clone().filter(col("mapq").gt_eq(min_mapq));
                }

                let df = df.collect().unwrap();

                let cat_col = df.column("chrom1").unwrap().categorical().unwrap();
                let rev_map1 = cat_col.get_rev_map();

                let cat_col = df.column("chrom2").unwrap().categorical().unwrap();
                let rev_map2 = cat_col.get_rev_map();

                let nrows = df.height();

                let mut data: Vec<bool> = Vec::with_capacity(nrows);

                for idx in 0..nrows {
                    let row = df.get(idx).unwrap();
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

                        let is_in_regions = interval_hash.get(chrom1).map_or(false, |interval1| {
                            interval_hash.get(chrom2).map_or(false, |interval2| {
                                interval1.count(pos1, pos1 + 1) > 0 && interval2.count(pos2, pos2 + 1) > 0
                            })
                        });
                        if is_in_regions ^ invert {
                            data.push(true);
                        } else {
                            data.push(false);
                        }
                    }
                }

                let data = Series::new("intersect", data);
                
                let mut df = df.filter(
                    data.bool().unwrap()
                ).unwrap();
                
                let file_name = Path::new(&file).file_name().unwrap().to_str().unwrap();
                let new_file = format!("{}/{}", output_q_dir, file_name);
        
                let mut new_file = File::create(new_file).unwrap();
                ParquetWriter::new(&mut new_file)
                    .finish(&mut df)
                    .unwrap();
            })
        });
        
        let files = collect_parquet_files(format!("{}", output_q_dir).as_str());

        files.par_iter().for_each(|file| {
            let file_name = Path::new(file).file_name().unwrap().to_str().unwrap();
            let new_file = format!("{}/{}", copy_q_dir, file_name);
            let _ = std::fs::copy(file, new_file);
        });

        log::info!("Successful output intersect file `{}`", output);
        
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

    for file in input {
        let p = PQS::new(file);
        if !p.is_pqs() {
            log::warn!("{} is not a valid PQS directory, skipped.", file);
            continue;
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
    // q0_pq_files.into_par_iter().enumerate().for_each(|(idx, file)| {
    //     let output_file = format!("{}/q0/{}.parquet", output, idx);
    //     let _ = std::fs::copy(file, output_file);
    // });


    // q1_pq_files.into_par_iter().enumerate().for_each(|(idx, file)| {
    //     let output_file = format!("{}/q1/{}.parquet", output, idx);
    //     let _ = std::fs::copy(file, output_file);
    // });

}