#[allow(unused)]
use anyhow::Result as anyResult;

use crossbeam_channel::{unbounded, bounded, Receiver, Sender};
use log::LevelFilter;
use std::collections::{ BTreeMap, HashMap, HashSet };
use std::borrow::Cow;
use std::path::{ Path, PathBuf };
use walkdir::WalkDir;

use std::fs::File;
use rayon::prelude::*;


use crate::core::{ common_reader, common_writer };
use crate::core::{ 
    BaseTable, 
    ChromSize,
    ChromSizeRecord, 
    ContigPair, ContigPair2,
    binify
};




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


pub fn parquets2clm(parquets_dir: &String) -> anyResult<()> {
    use polars::prelude::*;
    use polars::enable_string_cache;
    enable_string_cache();

    std::env::set_var("POLARS_MAX_THREADS", "10");

    let min_mapq = 0;
    let files = collect_parquet_files(parquets_dir);
    // let mut hashmap = HashMap::new();
    let mut results = Vec::new();
    for file in files {
        // let mut file = File::open(file).unwrap();
        let mut df = LazyFrame::scan_parquet(file,  ScanArgsParquet::default()).unwrap();
        
        if min_mapq > 1 {
            df = df.clone().filter(col("mapq").gt_eq(min_mapq));
            
        }
        
        let mut result = df.group_by(["chrom1", "chrom2"]).agg(vec![col("pos1"), col("pos2")]);
      
        results.push(result);


    }

    let mut concatenated = concat(
        results.iter().map(|x| x.clone()).collect::<Vec<_>>(),
        UnionArgs::default()
    ).unwrap();
    // let mut concatenated = results.remove(0);
    // for result in results {
    //     concatenated = concatenated.vstack(&result).unwrap();
    // }
   

    let mut results = concatenated
                        // .lazy()
                        .group_by(["chrom1", "chrom2"])
                        .agg(vec![col("pos1"), col("pos2")]);
                        // .collect().unwrap();

    

    let mut results = results.with_column(
        col("pos1")
            .arr().0
            .apply(
                |nested_list_series| {

                    let ca = nested_list_series.list().unwrap();
                    // let mut flattened = Vec::with_capacity(ca.len());

                    // for opt_sublist in ca.into_iter() {
                        
                    //     if let Some(sublist) = opt_sublist {

                    //         let merged: Vec<u32> = sublist
                    //             .list()
                    //             .unwrap()
                    //             .into_no_null_iter()
                    //             .flat_map(|inner_s| {
                    //                 let s = inner_s.u32().unwrap().clone();
                    //                 s.into_no_null_iter().collect::<Vec<u32>>()
                    //             })
                    //             .collect();

                    //         flattened.push(Some(Series::new("", merged)));
                            
                    //     } else {
                    //         flattened.push(None);
                    //     }
                    // }

                    let flattened: Vec<Option<Series>> = ca
                    .into_iter() 
                    .map(|opt_sublist| {
                        if let Some(sublist) = opt_sublist {
                            let merged: Vec<u32> = sublist
                                .list()
                                .unwrap()
                                .into_no_null_iter()
                                .flat_map(|inner_s| {
                                    let s = inner_s.u32().unwrap().clone();
                                    s.into_no_null_iter().collect::<Vec<u32>>()
                                })
                                .collect();

                            Some(Series::new("", merged))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .collect();


                   
                    let lc = ListChunked::from_iter(flattened);
    
                    Ok(Some(lc.into_series()))
                },
                GetOutput::from_type(DataType::List(Box::new(DataType::UInt32)))
            )
            .alias("pos1"))
    .with_column(
        col("pos2")
            .arr().0
            .apply(
                |nested_list_series| {

                    let ca = nested_list_series.list().unwrap();
                    // let mut flattened = Vec::with_capacity(ca.len());
                    // for opt_sublist in ca.into_iter() { 
                    //     if let Some(sublist) = opt_sublist {
                    //         let merged: Vec<u32> = sublist
                    //             .list()
                    //             .unwrap()
                    //             .into_no_null_iter()
                    //             .flat_map(|inner_s| {
                    //                 let s = inner_s.u32().unwrap().clone();
                    //                 s.into_no_null_iter().collect::<Vec<u32>>()
                    //             })
                    //             .collect();

                    //         flattened.push(Some(Series::new("", merged)));
                            
                    //     } else {
                    //         flattened.push(None);
                    //     }
                    // }

                    let flattened: Vec<Option<Series>> = ca
                    .into_iter() 
                    .map(|opt_sublist| {
                        if let Some(sublist) = opt_sublist {
                            let merged: Vec<u32> = sublist
                                .list()
                                .unwrap()
                                .into_no_null_iter()
                                .flat_map(|inner_s| {
                                    let s = inner_s.u32().unwrap().clone();
                                    s.into_no_null_iter().collect::<Vec<u32>>()
                                })
                                .collect();

                            Some(Series::new("", merged))
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>()
                    .into_par_iter()
                    .collect();

                   
                    let lc = ListChunked::from_iter(flattened);
    
                    Ok(Some(lc.into_series()))
                },
                GetOutput::from_type(DataType::List(Box::new(DataType::UInt32)))
            )
            .alias("pos2")
    )
    // .with_column(
    //     col("pos1").arr().0
    //         .apply(
    //             |s| {
    //                 let ca = s.list().unwrap();
    //                 // let mut vec = Vec::with_capacity(ca.len());
    //                 // for opt_sublis in ca {
    //                 //     if let Some(sublis) = opt_sublis {
    //                 //         vec.push(sublis.len() as u32);
    //                 //     } else {
    //                 //         vec.push(0);
    //                 //     }
    //                 // }
    //                 let vec: Vec<u32> = ca
    //                         .into_iter()
    //                         .map(|opt_sublist| {
    //                             if let Some(sublist) = opt_sublist {
    //                                 sublist.len() as u32
    //                             } else {
    //                                 0
    //                             }
    //                         })
    //                         .collect::<Vec<_>>()
    //                         .into_par_iter()
    //                         .collect();
    //                 Ok(Some(Series::new("", vec)))
    //             },
    //             GetOutput::from_type(DataType::UInt32)
    //         ).alias("count")
    // )
    .collect().unwrap();
   
    
    let cat_col = results.column("chrom1").unwrap().categorical().unwrap();
    let rev_map1 = cat_col.get_rev_map();

    // for idx in &cat_col.physical().unique().unwrap() {
    //     if let Some(idx) = idx {
    //         println!("{:?}, {:?}", idx, rev_map.get(idx));
    //     };
    // }

    let cat_col = results.column("chrom2").unwrap().categorical().unwrap();
    let rev_map2 = cat_col.get_rev_map();

    // for idx in &cat_col.physical().unique().unwrap() {
    //     if let Some(idx) = idx {
    //         println!("{:?}, {:?}", idx,  rev_map.get(idx));
    //     };
    // }


    // let nrows = results.height();
    // for idx in 0..nrows {
    //     let row = results.get(idx).unwrap();
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
        
    //     if let (Some(chrom1), Some(chrom2)) = (chrom1, chrom2) {
    //         // println!("{:?} {:?}", rev_map1.get(*chrom1), rev_map2.get(*chrom2));
    //     }
        
        
    //     for (p1, p2) in pos1.unwrap().iter().zip(pos2.unwrap().iter()) {
    //         // println!("{:?} {:?}", p1, p2);
    //     }
        
    // }
    

    println!("{:?}", results);


    Ok(())
}


