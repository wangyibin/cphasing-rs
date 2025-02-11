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
        use hashbrown::HashMap;
        enable_string_cache();

        std::env::set_var("POLARS_MAX_THREADS", format!("{}", 1));
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

        let contig_idx: HashMap<String, u32, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (k.clone(), i as u32)).collect();
        let idx_contig: HashMap<u32, String, BuildHasherDefault<XxHash64>> = contigsizes.keys().enumerate().map(|(i, k)| (i as u32, k.clone())).collect();
        let idx_sizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = contigsizes.iter().map(|(k, v)| (contig_idx.get(k).unwrap().clone(), v.clone())).collect();

        // let mut hashmap = HashMap::new();
        // let mut results = Vec::new();
        // for file in files {
        //     let mut df = LazyFrame::scan_parquet(file,  ScanArgsParquet::default()).unwrap();

        //     if min_mapq > 1 {
        //         df = df.clone().filter(col("mapq").gt_eq(min_mapq));
            
        //     }

        //     let mut result = df.group_by(["chrom1", "chrom2"])
        //                         .agg(vec![col("pos1"), col("pos2")])
        //                         .collect().unwrap();

        //     results.push(result);

        // }

        
        let mut results = files.into_par_iter().filter_map(|file| {
            let mut df = match LazyFrame::scan_parquet(file.clone(),  ScanArgsParquet::default()) {
                Ok(df) => df,
                Err(e) => {
                    log::warn!("Empty file: {:?}", file);
                    return None;
                }
            };

            if min_mapq > 1 {
                df = df.clone().filter(col("mapq").gt_eq(min_mapq));
                
            }

            let result = df.group_by(["chrom1", "chrom2"])
                        .agg(vec![col("pos1"), col("pos2")])
                        .collect().unwrap();

            // let result = df.collect().unwrap();
            Some(result)
        }).collect::<Vec<_>>();

        let mut data: HashMap<(u32, u32), Vec<SmallIntVec>> = HashMap::new();

        for df in &results {
           
            let cat_col = df.column("chrom1").unwrap().categorical().unwrap();
            let rev_map1 = cat_col.get_rev_map();

            let cat_col = df.column("chrom2").unwrap().categorical().unwrap();
            let rev_map2 = cat_col.get_rev_map();

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
                    // let vec: SmallIntVec = smallvec![*pos1, *pos2];
                    let mut vec: Vec<SmallIntVec> = Vec::new();
                    for (p1, p2) in pos1.u32().unwrap().iter().zip(pos2.u32().unwrap().iter()) {
                        vec.push(smallvec![p1.unwrap(), p2.unwrap()]);
                    }

                    data.entry((*chrom1, *chrom2)).or_insert(Vec::new()).extend(vec);
                }

            }
        }

        if output_split_contacts{
            log::info!("Calculating the distance between split contigs");
            let mut contacts = Contacts::new(&format!(".{}.pixels", output_prefix));
            let split_contigsizes: HashMap<u32, u32, BuildHasherDefault<XxHash64>> = idx_sizes
                .par_iter()
                .map(|(k, v)| (*k, *v / 2))
                .collect(); 
    
            let writer = common_writer(format!("{}.split.contacts", output_prefix.to_string()).as_str());
            let mut writer = Arc::new(Mutex::new(writer));
            data.par_iter().for_each(|(cp, vec) | {
                if vec.len() < min_contacts as usize {
                    return;
                }
                let length1 = idx_sizes.get(&cp.0).unwrap();
                let length2 = idx_sizes.get(&cp.1).unwrap();
                let contig1 = idx_contig.get(&cp.0).unwrap();
                let contig2 = idx_contig.get(&cp.1).unwrap();
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
                    *depth.get_mut(&cp.0).unwrap().get_mut(*split_index1 as usize).unwrap() += 1;
                    *depth.get_mut(&cp.1).unwrap().get_mut(*split_index2 as usize).unwrap() += 1;
                });
                
            });

            let depth: BTreeMap<_, _> = depth.into_iter().collect();
            let writer = common_writer(format!("{}.depth", output_prefix.to_string()).as_str());
            let mut wtr = Arc::new(Mutex::new(writer));

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

            log::info!("Successful output depth file `{}`", 
                                &format!("{}.depth", output_prefix.to_string()));
      
        }

        let mut wtr = common_writer(output.as_str());
        let wtr = Arc::new(Mutex::new(wtr));
        data.par_iter().for_each(|(cp, vec)| {
            if cp.0 == cp.1 {
                return;
            }

            let length1 = idx_sizes.get(&cp.0).unwrap();
            let length2 = idx_sizes.get(&cp.1).unwrap();

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
            let contig1 = idx_contig.get(&cp.0).unwrap();
            let contig2 = idx_contig.get(&cp.1).unwrap();

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


        // results.into_par_iter()
        //     .map(|df| {
        //         let cat_col1 = df.column("chrom1").unwrap().categorical().unwrap();
        //         let rev_map1 = cat_col1.get_rev_map();

        //         let cat_col2 = df.column("chrom2").unwrap().categorical().unwrap();
        //         let rev_map2 = cat_col2.get_rev_map();

        //         let nrows = df.height();
        //         let mut local_data: HashMap<(u32, u32), Vec<SmallVec<[u32; 2]>>> = HashMap::new();

        //         for idx in 0..nrows {
        //             let row = df.get(idx).unwrap();
        //             let chrom1 = match row.get(0) {
        //                 Some(AnyValue::Categorical(v, _, _)) => Some(v),
        //                 _ => None
        //             };

        //             let chrom2 = match row.get(1) {
        //                 Some(AnyValue::Categorical(v, _, _)) => Some(v),
        //                 _ => None
        //             };

        //             let pos1 = match row.get(2) {
        //                 Some(AnyValue::List(v)) => Some(v),
        //                 _ => None
        //             };

        //             let pos2 = match row.get(3) {
        //                 Some(AnyValue::List(v)) => Some(v),
        //                 _ => None
        //             };

        //             if let (Some(chrom1), Some(chrom2), Some(pos1), Some(pos2)) = (chrom1, chrom2, pos1, pos2) {
        //                 let mut vec: Vec<SmallIntVec> = Vec::new();
        //                 for (p1, p2) in pos1.u32().unwrap().iter().zip(pos2.u32().unwrap().iter()) {
        //                     vec.push(smallvec![p1.unwrap(), p2.unwrap()]);
        //                 }

        //                 local_data.entry((*chrom1, *chrom2)).or_insert(vec);
        //             }
        //         }

        //         local_data
        //     })
        //     .reduce(
        //         || HashMap::new(),
        //         |mut acc, local_data| {
        //             for (key, value) in local_data {
        //                 acc.entry(key).or_insert_with(Vec::new).extend(value);
        //             }
        //             acc
        //         },
        //     );





        // if !low_memory {
        //     // let mut concatenated = concat(
        //     //     results.iter().map(|x| x.lazy()).collect::<Vec<_>>(),
        //     //     UnionArgs::default()
        //     // ).unwrap();
        //     let mut concatenated = results.remove(0);
        //     for result in results {
        //         concatenated = concatenated.vstack(&result).unwrap();
        //     }

        //     let mut results = concatenated
        //                         .lazy()
        //                         .group_by(["chrom1", "chrom2"])
        //                         .agg(vec![col("pos1"), col("pos2")]);
        //                         // .collect().unwrap();



        //     let mut results = results.with_column(
        //         col("pos1")
        //             .arr().0
        //             .apply(
        //                 |nested_list_series| {
        //                     let ca = nested_list_series.list().unwrap();

        //                     let flattened: Vec<Option<Series>> = ca
        //                     .into_iter() 
        //                     .map(|opt_sublist| {
        //                         if let Some(sublist) = opt_sublist {
        //                             let merged: Vec<u32> = sublist
        //                                 .list()
        //                                 .unwrap()
        //                                 .into_no_null_iter()
        //                                 .flat_map(|inner_s| {
        //                                     let s = inner_s.u32().unwrap().clone();
        //                                     s.into_no_null_iter().collect::<Vec<u32>>()
        //                                 })
        //                                 .collect();

        //                             Some(Series::new("", merged))
        //                         } else {
        //                             None
        //                         }
        //                     })
        //                     .collect();
                    
        //                     let lc = ListChunked::from_iter(flattened);

        //                     Ok(Some(lc.into_series()))
        //                 },
        //                 GetOutput::from_type(DataType::List(Box::new(DataType::UInt32)))
        //             )
        //             .alias("pos1"))
        //     .with_column(
        //         col("pos2")
        //             .arr().0
        //             .apply(
        //                 |nested_list_series| {

        //                     let ca = nested_list_series.list().unwrap();

        //                     let flattened: Vec<Option<Series>> = ca
        //                     .into_iter() 
        //                     .map(|opt_sublist| {
        //                         if let Some(sublist) = opt_sublist {
        //                             let merged: Vec<u32> = sublist
        //                                 .list()
        //                                 .unwrap()
        //                                 .into_no_null_iter()
        //                                 .flat_map(|inner_s| {
        //                                     let s = inner_s.u32().unwrap().clone();
        //                                     s.into_no_null_iter().collect::<Vec<u32>>()
        //                                 })
        //                                 .collect();

        //                             Some(Series::new("", merged))
        //                         } else {
        //                             None
        //                         }
        //                     })
        //                     .collect();

                    
        //                     let lc = ListChunked::from_iter(flattened);

        //                     Ok(Some(lc.into_series()))
        //                 },
        //                 GetOutput::from_type(DataType::List(Box::new(DataType::UInt32)))
        //             )
        //             .alias("pos2")
        //     )
        //     // .with_column(
        //     //     col("pos1").arr().0
        //     //         .apply(
        //     //             |s| {
        //     //                 let ca = s.list().unwrap();
        //     //                 // let mut vec = Vec::with_capacity(ca.len());
        //     //                 // for opt_sublis in ca {
        //     //                 //     if let Some(sublis) = opt_sublis {
        //     //                 //         vec.push(sublis.len() as u32);
        //     //                 //     } else {
        //     //                 //         vec.push(0);
        //     //                 //     }
        //     //                 // }
        //     //                 let vec: Vec<u32> = ca
        //     //                         .into_iter()
        //     //                         .map(|opt_sublist| {
        //     //                             if let Some(sublist) = opt_sublist {
        //     //                                 sublist.len() as u32
        //     //                             } else {
        //     //                                 0
        //     //                             }
        //     //                         })
        //     //                         .collect::<Vec<_>>()
        //     //                         .into_par_iter()
        //     //                         .collect();
        //     //                 Ok(Some(Series::new("", vec)))
        //     //             },
        //     //             GetOutput::from_type(DataType::UInt32)
        //     //         ).alias("count")
        //     // )
        //     .collect().unwrap();


        //     let cat_col = results.column("chrom1").unwrap().categorical().unwrap();
        //     let rev_map1 = cat_col.get_rev_map();

        //     // for idx in &cat_col.physical().unique().unwrap() {
        //     //     if let Some(idx) = idx {
        //     //         println!("{:?}, {:?}", idx, rev_map.get(idx));
        //     //     };
        //     // }

        //     let cat_col = results.column("chrom2").unwrap().categorical().unwrap();
        //     let rev_map2 = cat_col.get_rev_map();

        //     // for idx in &cat_col.physical().unique().unwrap() {
        //     //     if let Some(idx) = idx {
        //     //         println!("{:?}, {:?}", idx,  rev_map.get(idx));
        //     //     };
        //     // }


        //     // let nrows = results.height();
        //     // for idx in 0..nrows {
        //     //     let row = results.get(idx).unwrap();
        //     //     let chrom1 = match row.get(0) {
        //     //         Some(AnyValue::Categorical(v, _, _)) => Some(v),
        //     //         _ => None
        //     //     };
            
        //     //     let chrom2 = match row.get(1) {
        //     //         Some(AnyValue::Categorical(v, _, _)) => Some(v),
        //     //         _ => None
        //     //     };

        //     //     let pos1 = match row.get(2) {
        //     //         Some(AnyValue::List(v)) => Some(v),
        //     //         _ => None
        //     //     };
            
        //     //     let pos2 = match row.get(3) {
        //     //         Some(AnyValue::List(v)) => Some(v),
        //     //         _ => None
        //     //     };
            
        //     //     if let (Some(chrom1), Some(chrom2)) = (chrom1, chrom2) {
        //     //         // println!("{:?} {:?}", rev_map1.get(*chrom1), rev_map2.get(*chrom2));
        //     //     }
            
            
        //     //     for (p1, p2) in pos1.unwrap().iter().zip(pos2.unwrap().iter()) {
        //     //         // println!("{:?} {:?}", p1, p2);
        //     //     }
            
        //     // }
        //     println!("{:?}", results);
        // }



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

        Ok(())
    }

    pub fn intersect(&self, hcr_bed: &String, invert: bool,
                    min_mapq: u8, output: &String) -> anyResult<()> {
                        
        type IvU8 = Interval<usize, u8>;
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

        let mut results = files.into_par_iter().map(|file| {
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

            let mut data: Vec<bool> = Vec::new();

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
        
        }
        ).collect::<Vec<_>>();
        

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

    // mkdir -p output
    let _ = std::fs::create_dir_all(output);
    let _ = std::fs::create_dir_all(format!("{}/q0", output));
    let _ = std::fs::create_dir_all(format!("{}/q1", output));

    // mv input/_contigsizes output/_contigsizes
    let _ = std::fs::copy(format!("{}/_contigsizes", input[0]), format!("{}/_contigsizes", output));
    let _ = std::fs::copy(format!("{}/_metadata", input[0]), format!("{}/_metadata", output));
    let _ = std::fs::copy(format!("{}/_readme", input[0]), format!("{}/_metadata", output));

    // let idx = AtomicUsize::new(0);
    // let idx2 = AtomicUsize::new(0);

    // input.par_iter().for_each(|file| {
    //     let p = PQS::new(file);
    //     if !p.is_pqs() {
    //         log::warn!("{} is not a valid PQS directory, skipped.", file);
    //         return;
    //     }

    //     let q0 = format!("{}/q0", file);
    //     let q1 = format!("{}/q1", file);

    //     let q0_files = collect_parquet_files(q0.as_str());
    //     let q1_files = collect_parquet_files(q1.as_str());

    //     q0_files.par_iter().for_each(|q0_file| {
    //         let current = idx.fetch_add(1, Ordering::SeqCst);
    //         let output_file = format!("{}/q0/{}.q0.parquet", output, current);
    //         let _ = std::fs::copy(q0_file, output_file);
    //     });

    //     q1_files.par_iter().for_each(|q1_file| {
    //         let current = idx2.fetch_add(1, Ordering::SeqCst);
    //         let output_file = format!("{}/q1/{}.q1.parquet", output, current);
    //         let _ = std::fs::copy(q1_file, output_file);
    //     });
    // });

    let mut idx = 0;
    let mut idx2 = 0;
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

        for q0_file in q0_files {
            
            let output_file = format!("{}/q0/{}.q0.parquet", output, idx);
            let _ = std::fs::copy(q0_file, output_file);
            idx += 1;
        }

        for q1_file in q1_files {
            let output_file = format!("{}/q1/{}.q1.parquet", output, idx2);
            let _ = std::fs::copy(q1_file, output_file);
            idx2 += 1;
        }
    }
}