#![allow(dead_code)]
#![allow(unused_variables, unused_assignments)]
use std::io::BufReader;
use std::io::BufRead;
use std::path::Path;
use indexmap::IndexMap;
use rayon::ThreadPoolBuilder;
use cphasing::alleles::AllelesFasta;
use cphasing::aligner::read_bam;
use cphasing::bam::*;
use cphasing::cli::cli;
use cphasing::clm::{Clm, merge_clm};
use cphasing::core::{ 
    BaseTable,  common_reader, 
    common_writer, ContigPair,
    check_program};
use cphasing::count_re::CountRE;
// use cphasing::contacts::Contacts;
use cphasing::cutsite::cut_site;
use cphasing::fastx::{ Fastx, split_fastq };
use cphasing::methy::{ modbam2fastq, modify_fasta };
use cphasing::paf::PAFTable;
use cphasing::pairs::*;
use cphasing::porec::{
        PoreCTable, merge_porec_tables };
use cphasing::kprune::{ PruneTable, KPruner };
use cphasing::prune::{ Pruner };
use cphasing::pqs::*;
use cphasing::simulation::{ 
        simulation_from_split_read, simulate_porec,
        simulate_hic };
use cphasing::order::*;
use cphasing::orientation::*;
use cphasing::splitcontacts::*;
use std::collections::{ HashMap, HashSet };
use std::io::Write; 
use chrono::Local;
use env_logger::Builder;
use log::LevelFilter;

// #[cfg(not(all(target_os = "linux", target_arch = "aarch64")))]
// use jemallocator::Jemalloc;

// #[cfg(not(all(target_os = "linux", target_arch = "aarch64")))]
// #[global_allocator]
// static GLOBAL: Jemalloc = Jemalloc;

use tikv_jemallocator::Jemalloc;
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

fn main() {

    Builder::new()
        .format(|buf, record| {
            writeln!(buf, "{} [{}] - {}", 
                Local::now().format("%Y-%m-%dT%H:%M:%S"),
                record.level(),
                record.args()
            )
        })
        .filter(None, LevelFilter::Info)
        .init();

    let matches = cli().get_matches();
    match matches.subcommand() {
        Some(("aligner", sub_matches)) => {
            let fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let input_bam = sub_matches.get_one::<String>("BAM").expect("required");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let min_prob = sub_matches.get_one::<f32>("MIN_PROB").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            check_program("minimap2");

            let fa = Fastx::new(&fasta);
            let seqs = fa.get_chrom_seqs().unwrap();
            read_bam(&input_bam, &seqs, *min_quality, *min_prob, &output);
        }
        Some(("realign", sub_matches)) => {
            use cphasing::realign::read_bam;
            use cphasing::realign::read_paf;
            let input = sub_matches.get_one::<String>("PAF").expect("required");
            let import_format = sub_matches.get_one::<String>("FORMAT").expect("error");
            // let contacts = sub_matches.get_one::<String>("CONTACTS").expect("required");
            let contacts = if let Some(c) = sub_matches.get_one::<String>("CONTACTS") {
                Some(c.to_string())
            } else {
                None
            };
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let min_mapq = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();


            match import_format.as_str() {
                "bam" => {
                    read_bam(&input, *min_mapq, &output, *threads, contacts);
                }
                "paf" => {
                    read_paf(&input, *min_mapq, &output);
                }
                _ => {
                    eprintln!("No such format.");
                }
            }

        }
        Some(("methalign", submathes)) => {
            use cphasing::methalign::parse_bam;
            let input_bam = submathes.get_one::<String>("BAM").expect("required");
            let fasta_opt: Option<String> = submathes
                .get_one::<String>("FASTA")
                .and_then(|s| if s == "none" || s == "-" { None } else { Some(s.to_string()) });
            let bedgraph_opt: Option<String> = submathes
                .get_one::<String>("BEDGRAPH")
                .and_then(|s| if s == "none" || s == "-" { None } else { Some(s.to_string()) });
            let match_score = submathes.get_one::<i32>("MATCH_SCORE").expect("error");
            let ref_penalty = submathes.get_one::<i32>("REF_PENALTY").expect("error");
            let read_penalty = submathes.get_one::<i32>("READ_PENALTY").expect("error");
            let ref_prob_cutoff = submathes.get_one::<f64>("REF_PROB_CUTOFF").expect("error");
            let prob_cutoff = submathes.get_one::<u8>("PROB_CUTOFF").expect("error");
            let designate_mapq = submathes.get_one::<u8>("DESIGNATE_MAPQ").expect("error");
            let is_set_y = submathes.get_one::<bool>("IS_SET_Y").expect("error");
            let cpg = submathes.get_one::<bool>("CPG").expect("error");
            let output_secondary = submathes.get_one::<bool>("OUTPUT_SECONDARY").expect("error");
            let output_bam = submathes.get_one::<String>("OUTPUT").expect("error");
            let threads = submathes.get_one::<usize>("THREADS").expect("error");
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();
            
            parse_bam(&input_bam, 
                      &fasta_opt,
                      &bedgraph_opt,
                      *match_score,
                      *ref_penalty,
                      *read_penalty,
                      *ref_prob_cutoff,
                      *prob_cutoff,
                      *designate_mapq,
                      *is_set_y,
                      *cpg,
                      *output_secondary,
                      &output_bam,
                      *threads);
        }
        Some(("alleles", sub_matches)) => {

            let fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let kmer_size = sub_matches.get_one::<usize>("K").expect("error");
            let window_size = sub_matches.get_one::<usize>("W").expect("error");
            let minimum_similarity = sub_matches.get_one::<f64>("M").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let mut alleles = AllelesFasta::new(&fasta);
            alleles.run(*kmer_size, *window_size, *minimum_similarity, output);

        }

        Some(("kprune", sub_matches)) => {
            use rayon::prelude::*;
            
            let alleletable = sub_matches.get_one::<String>("ALLELETABLE").expect("required");
            let contacts = sub_matches.get_one::<String>("CONTACTS").expect("required");
            
            let prunetable = sub_matches.get_one::<String>("PRUNETABLE").expect("required");
            let count_re_opt: Option<String> = sub_matches
                .get_one::<String>("COUNTRE")
                .and_then(|s| if s == "none" || s == "-" { None } else { Some(s.to_string()) });
            let method = sub_matches.get_one::<String>("METHOD").expect("error");
            let normalization_method = sub_matches.get_one::<String>("NORMALIZATION_METHOD").expect("error");
            let whitelist = sub_matches.get_one::<String>("WHITELIST").expect("error");
            let partial_whitelist = sub_matches.get_one::<bool>("PARTIAL_WHITELIST").expect("error");
            let first_cluster = sub_matches.get_one::<String>("FIRST_CLUSTER").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");

            assert!(method == "fast" || method == "precise" || method == "greedy", 
                     "method must be in ['fast', 'precise', 'greedy']");
     
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let mut whitehash: HashSet<String> = HashSet::new();

            if whitelist != "none" {
                let reader = common_reader(&whitelist);
                let reader = BufReader::new(reader);
                for line in reader.lines() {
                    let line = line.unwrap();
                    whitehash.insert(line);
                }
            }


            let mut first_cluster_hashmap: HashMap<String, HashSet<String>> = HashMap::new();
            if first_cluster != "none" {
                let reader = common_reader(&first_cluster);
                let reader = BufReader::new(reader);
              
                let pairs: Vec<(String, HashSet<String>)> = reader
                    .lines()
                    .par_bridge()
                    .filter_map(|r| r.ok())
                    .map(|line| {
                       
                        let mut parts = line.splitn(3, '\t');
                        
                        let group = parts.next().unwrap_or("").to_string();
                        let _count = parts.next().unwrap_or("");
                        let contig_field = parts.next().unwrap_or("");

                       
                        let mut set = HashSet::new();
                        if whitehash.is_empty() {
                            for s in contig_field.split_ascii_whitespace() {
                                if !s.is_empty() {
                                    set.insert(s.to_string());
                                }
                            }
                        } else {
   
                            for s in contig_field.split_ascii_whitespace() {
                                if !s.is_empty() && whitehash.contains(s) {
                                    set.insert(s.to_string());
                                }
                            }
                        }
                        (group, set)
                    })
                    .filter(|(_, set)| !set.is_empty()) 
                    .collect();

      
                first_cluster_hashmap = pairs.into_iter().collect();

            }

            if first_cluster != "none" {
                let mut writer = common_writer(&prunetable);
                log::set_max_level(log::LevelFilter::Off);
                let results: Vec<(String, String, u32, u32)> = first_cluster_hashmap
                    .into_par_iter()
                    .map(|(k, v)| {
                        let tmp_output = format!("{}.kprune.tmp", k);
                   
                        let mut tmp_writer = common_writer(&tmp_output);

                        log::info!("Pruning cluster `{}`", k);

                        let mut kpruner = KPruner::new(
                            &alleletable, &contacts, &prunetable, &count_re_opt,
                            normalization_method);

                        let white_refs: HashSet<&String> = v.iter().collect();


                        kpruner.prune(&method.as_str(), &white_refs, &mut tmp_writer, *partial_whitelist);
                        (k, tmp_output, kpruner.allelic_counts, kpruner.cross_allelic_counts)
                    })
                    .collect();

                log::set_max_level(log::LevelFilter::Info);

                let mut allelic_counts: u32 = 0;
                let mut cross_allelic_counts: u32 = 0;

                for (k, tmp_output, ac, cac) in results {
                    let reader = common_reader(&tmp_output);
                    let reader = BufReader::new(reader);
                    for line in reader.lines() {
                        let line = line.unwrap();
                        writer.write(format!("{}\n", line).as_bytes()).unwrap();
                    }
                    std::fs::remove_file(&tmp_output).unwrap();

                    allelic_counts += ac;
                    cross_allelic_counts += cac;
                }

                log::info!("Allelic counts: {}", allelic_counts);
                log::info!("Cross-allelic counts: {}", cross_allelic_counts);
            } else {
                let mut writer = common_writer(&prunetable);
                let mut whitehash2: HashSet<&String> = HashSet::new();
                for x in whitehash.iter() {
                    whitehash2.insert(x);
                }

                let mut kpruner = KPruner::new(&alleletable, &contacts, &prunetable, &count_re_opt, normalization_method);
                kpruner.prune(&method.as_str(), &whitehash2, &mut writer, *partial_whitelist);
            }
            
            log::info!("Allelic and cross-allelic information written into `{}`", prunetable);

        }
        Some(("prune", sub_matches)) => {
            use rayon::prelude::*;
            use std::sync::{Arc, Mutex};
            let alleletable = sub_matches.get_one::<String>("ALLELETABLE").expect("required");
            let allele_strand_table = sub_matches.get_one::<String>("ALLELESTRANDTABLE").expect("required");
            let contacts = sub_matches.get_one::<String>("CONTACTS").expect("required");
            let prunetable = sub_matches.get_one::<String>("PRUNETABLE").expect("required");
            let count_re_opt: Option<String> = sub_matches
                .get_one::<String>("COUNTRE")
                .and_then(|s| if s == "none" || s == "-" { None } else { Some(s.to_string()) });
            let method = sub_matches.get_one::<String>("METHOD").expect("error");
            let normalization_method = sub_matches.get_one::<String>("NORMALIZATION_METHOD").expect("error");
            let whitelist = sub_matches.get_one::<String>("WHITELIST").expect("error");
            let first_cluster = sub_matches.get_one::<String>("FIRST_CLUSTER").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");

            assert !(method == "fast" || method == "precise" || method == "greedy", 
                     "method must be in ['fast', 'precise', 'greedy']");
            let method = method.as_str();

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let mut whitehash: HashSet<String> = HashSet::new();

            if whitelist != "none" {
                let reader = common_reader(&whitelist);
                let reader = BufReader::new(reader);
                for line in reader.lines() {
                    let line = line.unwrap();
                    whitehash.insert(line);
                }
            }
            let mut whitehash2: HashSet<&String> = HashSet::new();
            for x in whitehash.iter() {
                whitehash2.insert(x);
            }

            let first_cluster_hashmap = Arc::new(Mutex::new(HashMap::new()));
            if first_cluster != "none" {
                let reader = common_reader(&first_cluster);
                let reader = BufReader::new(reader);

                let first_cluster_hashmap = Arc::new(Mutex::new(HashMap::<String, HashSet<String>>::new()));
            }
            if first_cluster != "none" {
                let reader = common_reader(&first_cluster);
                let reader = BufReader::new(reader);
              
                let lines: Vec<String> = reader.lines().map(|line| line.unwrap()).collect();
                lines.into_par_iter().for_each(|line| {
                    // split tab
                    let mut line = line.split("\t");
                    let group = line.next().unwrap().to_string();
                    let _count = line.next().unwrap().to_string();
                    // split next with space
                    let contigs: HashSet<_> = line.next().unwrap().split(" ").collect();
                   
                    let mut contigs: HashSet<String> = contigs.into_par_iter().map(|x| x.to_string()).collect();
                    if !whitehash.is_empty() {
                        contigs = contigs.intersection(&whitehash).cloned().collect();
                    }
                    let mut first_cluster_hashmap = first_cluster_hashmap.lock().unwrap();
                    first_cluster_hashmap.insert(group, contigs);
                });
            }

            let first_cluster_hashmap = Arc::clone(&first_cluster_hashmap);
            let first_cluster_hashmap = first_cluster_hashmap.lock().unwrap();
            let first_cluster_hashmap = first_cluster_hashmap.clone();

            let mut writer = common_writer(&prunetable);
       
            if first_cluster != "none" {
            
                log::set_max_level(log::LevelFilter::Off);
                let mut pruners: HashMap<String, Pruner> = first_cluster_hashmap
                            .keys().map(|x| (x.to_string(), Pruner::new(
                                &alleletable,  &allele_strand_table, &contacts,
                                &count_re_opt,
                                normalization_method))).collect();
                
                for (k, v) in first_cluster_hashmap {
                    log::info!("Pruning cluster `{}`", k);
                    let pruner = pruners.get_mut(&k).unwrap();
                    let _whitehash: HashSet<String> = v.into_iter().collect();
                    let mut _whitehash2: HashSet<&String> = HashSet::new();
                    for x in _whitehash.iter() {
                        _whitehash2.insert(&x);
                    }
                    
                    match method {
                        "fast" => {
                            pruner.kprune(&_whitehash2, &method, &mut writer);
                        }
                        "precise" => {
                            pruner.kprune(&_whitehash2, &method, &mut writer);
                        }
                        "greedy" => {
                            pruner.prune(&_whitehash2, &mut writer);
                        }
                        _ => {
                            eprintln!("No such method.");
                        }
                    }
                  
                }
                log::set_max_level(log::LevelFilter::Info);
                
            } else {
                let mut pruner = Pruner::new(&alleletable, &allele_strand_table, &contacts, 
                                                    &count_re_opt,
                                                 normalization_method);
                match method {
                        "fast" => {
                            pruner.kprune(&whitehash2, &method, &mut writer);
                        }
                        "precise" => {
                            pruner.kprune(&whitehash2, &method, &mut writer);
                        }
                        "greedy" => {
                            pruner.prune(&whitehash2, &mut writer);
                        }
                        _ => {
                            eprintln!("No such method.");
                        }
                    }
                
            }
            
            log::info!("Allelic and cross-allelic information written into `{}`", prunetable);

        }
        Some(("splitbam", sub_matches)) => {
            let input_bam = sub_matches.get_one::<String>("BAM").expect("required");
            let output_prefix = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let record_num = sub_matches.get_one::<usize>("RECORD_NUM").expect("error");

            split_bam(&input_bam, &output_prefix, *record_num).unwrap();
        }
        Some(("mergeclm", sub_matches)) => {
            let inputs: Vec<_> = sub_matches.get_many::<String>("INPUTS").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mut clm_files = Vec::new();
            for input in inputs {
                let path = std::path::Path::new(input);
                if path.is_dir() {
                    for entry in std::fs::read_dir(path).unwrap() {
                        let entry = entry.unwrap();
                        let path_str = entry.path().to_str().unwrap().to_string();
                        if path_str.ends_with(".clm") || path_str.ends_with(".clm.gz") {
                            clm_files.push(path_str);
                        }
                    }
                } else if path.is_file() {
                    clm_files.push(input.to_string());
                } else {
                    log::warn!("Input `{}` is neither a file nor a directory, ignored.", input);
                }
            }

            if clm_files.is_empty() {
                log::error!("No clm files found in the given inputs.");
                std::process::exit(1);
            }

            merge_clm(clm_files, &output).unwrap();

        }
        Some(("splitclm", sub_matches)) => {
            let input_clm = sub_matches.get_one::<String>("CLM").expect("required");
            let cluster_file = sub_matches.get_one::<String>("CLUSTER").expect("required");
            let output_dir = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let clm = Clm::new(&input_clm);
            clm.split_clm(&cluster_file, &output_dir).unwrap();

        }
        Some(("splitcontacts", sub_matches)) => {
            let input_contacts = sub_matches.get_one::<String>("CONTACTS").expect("required");
            let cluster_file = sub_matches.get_one::<String>("CLUSTER").expect("required");
            let output_dir = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let _ = split_contacts_by_clusters(&input_contacts, &cluster_file, &output_dir);

        }
        Some(("splitfastq", sub_matches)) => {
            let input_fastq = sub_matches.get_one::<String>("FASTQ").expect("required");
            let output_prefix = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let record_num = sub_matches.get_one::<usize>("RECORD_NUM").expect("error");

            split_fastq(&input_fastq, &output_prefix, *record_num).unwrap();
        }
        Some(("extract-fasta", sub_matches)) => {
            let input_fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let clusters = sub_matches.get_one::<String>("CLUSTERS").expect("required");
            let trim_length = sub_matches.get_one::<usize>("TRIM_LENGTH").expect("error");

            let fasta = Fastx::new(&input_fasta);
            fasta.split_by_cluster(&clusters, *trim_length).unwrap();
        }
        Some(("slidefastq", sub_matches)) => {
            let input_fastq = sub_matches.get_one::<String>("FASTQ").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let window = sub_matches.get_one::<u64>("WINDOW").expect("error");
            let step = sub_matches.get_one::<u64>("STEP").expect("error");
            let min_length = sub_matches.get_one::<u64>("MIN_LENGTH").expect("error");
            let filetype = sub_matches.get_one::<String>("FILETYPE").expect("error").as_str();
            let coordinate_suffix = sub_matches.get_one::<bool>("COORDINATE").expect("error");
            let fa = Fastx::new(&input_fastq);
            let _ = fa.slide(&output, *window, *step, *min_length, &filetype, *coordinate_suffix);
        }
        Some(("slidefasta", sub_matches)) => {
            let input_fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let window = sub_matches.get_one::<u64>("WINDOW").expect("error");
            let step = sub_matches.get_one::<u64>("STEP").expect("error");
            
            let fa = Fastx::new(&input_fasta);
            let _ = fa.slidefasta(&output, *window, *step);
        }
        Some(("slide2raw", sub_matches)) => {
            let input_bam = sub_matches.get_one::<String>("BAM").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            slide2raw(&input_bam, &output, *threads);
            
        }
        Some(("simulator", sub_matches)) => {
            match sub_matches.subcommand() {
                Some(("split-ont", sub_sub_matches)) => {
                    let input_bam = sub_sub_matches.get_one::<String>("BAM").expect("required");
                    let output = sub_sub_matches.get_one::<String>("OUTPUT").expect("error");
                    let min_quality = sub_sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
                    simulation_from_split_read(&input_bam, &output, *min_quality);
                }
                Some(("porec", sub_sub_matches)) => {
                    let fasta = sub_sub_matches.get_one::<String>("FASTA").expect("required");
                    let vcf = sub_sub_matches.get_one::<String>("VCF").expect("required");
                    let bed = sub_sub_matches.get_one::<String>("BED").expect("required");
                    let output = sub_sub_matches.get_one::<String>("OUTPUT").expect("error");

                    simulate_porec(&fasta, &vcf, &bed, &output);
                }
                Some(("hic", sub_sub_matches)) => {
                    let fasta = sub_sub_matches.get_one::<String>("FASTA").expect("required");
                    let vcf = sub_sub_matches.get_one::<String>("VCF").expect("required");
                    let bam = sub_sub_matches.get_one::<String>("BAM").expect("required");
                    let min_mapq = sub_sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
                    let threads = sub_sub_matches.get_one::<usize>("THREADS").expect("error");
                    let output = sub_sub_matches.get_one::<String>("OUTPUT").expect("error");

                    simulate_hic(&fasta, &vcf, &bam, *min_mapq, *threads, &output);
                }
                _ => {
                    eprintln!("No such subcommand.");
                }
            }
            
        }
        Some(("kmer", sub_matches)) => {
            match sub_matches.subcommand() {
                Some(("count", sub_sub_matches)) => {
                    let input_fasta = sub_sub_matches.get_one::<String>("FASTA").expect("required");
                    let k = sub_sub_matches.get_one::<usize>("K").expect("error");
                    let output = sub_sub_matches.get_one::<String>("OUTPUT").expect("error");

                    let fasta = Fastx::new(&input_fasta);
                    let kmer_count = fasta.kmer_count(*k).unwrap();
                    let mut writer = common_writer(&output);
                    for (k, v) in kmer_count.iter() {
                        writer.write(format!("{}\t{}\n", k, v).as_bytes()).unwrap();
                    }

                }
                Some(("mask", sub_sub_matches)) => {
                    let input_fasta = sub_sub_matches.get_one::<String>("FASTA").expect("required");
                    let k = sub_sub_matches.get_one::<usize>("K").expect("error");
                    let ploidy = sub_sub_matches.get_one::<u64>("PLOIDY").expect("error");
                    let output = sub_sub_matches.get_one::<String>("OUTPUT").expect("error");
                    
                    let fasta = Fastx::new(&input_fasta);
                    fasta.mask_high_frequency_kmer(*k, *ploidy, output).unwrap();
                },
                Some(("position", sub_sub_matches)) => {
                    let input_fasta = sub_sub_matches.get_one::<String>("FASTA").expect("required");
                    let k = sub_sub_matches.get_one::<usize>("K").expect("required");
                    let kmer_list = sub_sub_matches.get_one::<String>("KMER_LIST").expect("required");
                    let output = sub_sub_matches.get_one::<String>("OUTPUT").expect("error");
                    
                    let fasta = Fastx::new(&input_fasta);
                    let _ = fasta.kmer_positions(*k, &kmer_list, &output);

                }
                _ => {
                    eprintln!("{:?}", sub_matches.subcommand());
                    eprintln!("No such subcommand.");
                }
            }
        }
        Some(("digest", sub_matches)) => {
            let input_fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let pattern = sub_matches.get_one::<String>("PATTERN").expect("error");
            let slope = sub_matches.get_one::<i64>("SLOPE").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let fasta = Fastx::new(&input_fasta);
            let positions = fasta.digest(&pattern, *slope).unwrap();
            let mut writer = common_writer(&output);
            for (k, v) in positions {
                for pos in v {
                    let pos_start = pos[0];
                    let pos_end = pos[1];
                    writer.write(format!("{}\t{}\t{}\n", k, pos_start, pos_end).as_bytes()).unwrap();
                }
                
            }

            log::info!("Digest results written to `{}`", output);
        }

        Some(("count_re", sub_matches)) => {
            let input_fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let pattern = sub_matches.get_one::<String>("PATTERN").expect("error");
            let min_re = sub_matches.get_one::<u64>("MIN_RE").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            // let fasta = Fastx::new(&input_fasta);
            // let contigsizes = fasta.get_chrom_size().unwrap();
            let fasta = Fastx::new(&input_fasta);
            let (counts, contigsizes) = fasta.count_re(&pattern).unwrap();
            let mut count_re = CountRE::new(&output);

            // counts + 1 
            let counts: IndexMap<String, u64> = counts.into_iter()
                                                    .map(|(k, v)| (k, v + 1))
                                                    .collect();
            // filter counts which are less than min re
            let counts: IndexMap<String, u64> = counts.into_iter()
                                                    .filter(|(_, v)| *v >= *min_re)
                                                    .collect();
            let contigsizes = contigsizes.into_iter()
                                        .filter(|(k, _)| counts.contains_key(k))
                                        .collect();
            count_re.from_hashmap(counts, contigsizes);
            count_re.write(&output);
        }
        Some(("cutsite", sub_matches)) => {
            let fastq = sub_matches.get_one::<String>("FASTQ").expect("required");
            let pattern = sub_matches.get_one::<String>("PATTERN").expect("error");
            let _output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            cut_site(fastq.to_string(), pattern.as_bytes(), _output.clone()).unwrap();
        }
        Some(("paf2depth", sub_matches)) => {
            let paf = sub_matches.get_one::<String>("PAF").expect("required");
            let chromsizes = sub_matches.get_one::<String>("CHROMSIZES").expect("required");
            let window_size = sub_matches.get_one::<usize>("WINSIZE").expect("error");
            let mut step_size = sub_matches.get_one::<usize>("STEPSIZE").expect("error");
            let min_mapq = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let secondary = sub_matches.get_one::<bool>("SECONDARY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let pt = PAFTable::new(&paf);
            if *step_size == 0 {
                step_size = window_size;
            } 
            pt.to_depth(&chromsizes, *window_size, *step_size, *min_mapq, *secondary, &output).unwrap();

        }
        Some(("paf2porec", sub_matches)) => {
            let paf = sub_matches.get_one::<String>("PAF").expect("required");
            let empty_string = String::new();
            let bed = sub_matches.get_one::<String>("BED").unwrap_or(&empty_string);
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let min_quality = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let min_identity = sub_matches.get_one::<f32>("MIN_IDENTITY").expect("error");
            let min_length = sub_matches.get_one::<u32>("MIN_LENGTH").expect("error");
            let max_edge = sub_matches.get_one::<u64>("MAX_EDGE").expect("error");
            // let min_order = sub_matches.get_one::<u32>("MIN_ORDER").expect("error");
            let max_order = sub_matches.get_one::<u32>("MAX_ORDER").expect("error");
            let secondary = sub_matches.get_one::<bool>("SECONDARY").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");

            let pt = PAFTable::new(&paf);

            pt.paf2table(bed, output, min_quality, min_identity, 
                            min_length, max_order, max_edge, *secondary, *threads).unwrap();

        }
        Some(("porec2pairs", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let chromsizes = sub_matches.get_one::<String>("CHROMSIZES").expect("required");

            if !Path::new(chromsizes).exists() {
                log::error!("Chromsizes file `{}` not found.", chromsizes);
                std::process::exit(1);
            }
            
            let output_raw = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let output = output_raw.strip_suffix("/").unwrap_or(output_raw).to_string();
            let chunksize = sub_matches.get_one::<usize>("CHUNKSIZE").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let min_order = sub_matches.get_one::<usize>("MIN_ORDER").expect("error");
            let max_order = sub_matches.get_one::<usize>("MAX_ORDER").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let mut prt = PoreCTable::new(&table);

            if output.ends_with(".pqs") {
                let output = if output == "-" {
                    let output = table.strip_suffix(".porec.gz").unwrap_or(&table).strip_suffix(".con.gz").unwrap_or(&table).strip_suffix(".concatemer.gz").unwrap_or(&table);
                    format!("{}.pqs", output)
                } else {
                    output.to_string()
                };
                
                if Path::new(&output).exists() {
                    log::info!("Output directory {} already exists. Removing it first.", output);
                    std::fs::remove_dir_all(&output).unwrap();
                }
                prt.to_pairs_pqs(&chromsizes, &output, *chunksize, *min_quality, *min_order, *max_order, *threads).unwrap();
            } else {
                prt.to_pairs(&chromsizes, &output, *min_quality, *min_order, *max_order).unwrap();
            }
            
        }
        Some(("porec-break", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let break_bed = sub_matches.get_one::<String>("BREAK_BED").expect("required");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mut prt = PoreCTable::new(&table);
            prt.break_contigs(&break_bed, &output, *threads);

        }

        Some(("porec-dup", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let collapsed_list = sub_matches.get_one::<String>("COLLAPSED").expect("required");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mut prt = PoreCTable::new(&table);
            prt.dup(&collapsed_list, 123, &output, *threads);

        }
        Some(("porec-merge", sub_matches)) => {
            let tables: Vec<_> = sub_matches.get_many::<String>("TABLES").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            merge_porec_tables(tables, output)
        }
        Some(("porec2reads", m)) => {
            // let table = m.get_one::<String>("TABLE").unwrap();
            // let fasta = m.get_one::<String>("FASTA").unwrap();
            // let out_prefix = m.get_one::<String>("OUTPUT").unwrap();
            // let length_arg = *m.get_one::<usize>("LENGTH").unwrap();
            
           
            // let read_len = if length_arg == 0 { None } else { Some(length_arg) };
            
            // let mut porec_table = PoreCTable::new(table);
            // porec_table.to_pe_pair_reads(fasta, out_prefix, read_len).unwrap();
            let table = m.get_one::<String>("TABLE").unwrap();
            let fasta = m.get_one::<String>("FASTA").unwrap();
            let out_prefix = m.get_one::<String>("OUTPUT").unwrap();
            let length_arg = *m.get_one::<usize>("LENGTH").unwrap();
            
           
            let read_len = if length_arg == 0 { None } else { Some(length_arg) };
            
            let read_len = if length_arg == 0 { None } else { Some(length_arg) };
                        
            let mut temp_files = Vec::new();

            let actual_table = if table.ends_with(".bam") {
                let tmp_paf = format!("{}.tmp.paf", out_prefix);
                let tmp_porec = format!("{}.tmp.porec.gz", out_prefix);
                let threads = *m.get_one::<usize>("THREADS").unwrap_or(&4);
                let min_quality = *m.get_one::<u8>("MIN_MAPQ").unwrap_or(&1);
                let min_identity = *m.get_one::<f32>("MIN_IDENTITY").unwrap_or(&0.0);
                let min_length = *m.get_one::<u32>("MIN_LENGTH").unwrap_or(&0);
                let max_edge = *m.get_one::<u64>("MAX_EDGE").unwrap_or(&100000);
                let max_order = *m.get_one::<u32>("MAX_ORDER").unwrap_or(&1000);
                let secondary = *m.get_one::<bool>("SECONDARY").unwrap_or(&false);
                std::fs::remove_file(&tmp_paf).ok();
                std::fs::remove_file(&tmp_porec).ok();

                temp_files.push(tmp_paf.clone());
                temp_files.push(tmp_porec.clone());

                // 1. BAM -> PAF
                bam2paf(table, &tmp_paf, threads, secondary);
                
                // 2. PAF -> PoreC
                let pt = PAFTable::new(&tmp_paf);
                let empty_string = String::new();
                pt.paf2table(&empty_string, &tmp_porec, &min_quality, &min_identity, 
                                &min_length, &max_order, &max_edge, secondary, threads).unwrap();

                tmp_porec
            } else if table.ends_with(".paf") || table.ends_with(".paf.gz") {
                let tmp_porec = format!("{}.tmp.porec.gz", out_prefix);
                let threads = *m.get_one::<usize>("THREADS").unwrap_or(&4);
                let min_quality = *m.get_one::<u8>("MIN_MAPQ").unwrap_or(&1);
                let min_identity = *m.get_one::<f32>("MIN_IDENTITY").unwrap_or(&0.0);
                let min_length = *m.get_one::<u32>("MIN_LENGTH").unwrap_or(&0);
                let max_edge = *m.get_one::<u64>("MAX_EDGE").unwrap_or(&100000);
                let max_order = *m.get_one::<u32>("MAX_ORDER").unwrap_or(&1000);
                let secondary = *m.get_one::<bool>("SECONDARY").unwrap_or(&false);
                std::fs::remove_file(&tmp_porec).ok();
                temp_files.push(tmp_porec.clone());

                let pt = PAFTable::new(table);
                let empty_string = String::new();
                pt.paf2table(&empty_string, &tmp_porec, &min_quality, &min_identity, 
                                &min_length, &max_order, &max_edge, secondary, threads).unwrap();

                tmp_porec
            } else {
                table.to_string()
            };

            let mut porec_table = PoreCTable::new(&actual_table);
            porec_table.to_pe_pair_reads(fasta, out_prefix, read_len).unwrap();

            for tmp in temp_files {
                std::fs::remove_file(tmp).ok();
            }
            
        }
        Some(("porec-intersect", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let bed = sub_matches.get_one::<String>("BED").expect("required");
            let invert = sub_matches.get_one::<bool>("INVERT").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();
            
            let mut prt = PoreCTable::new(&table);
            if *invert {
                log::info!("Invert the table.");
            }
            // prt.intersect(&bed, *invert, &output);
            prt.intersect_multi_threads(&bed, *invert, &output);
        }
              
        Some(("paf-downsample", sub_matches)) => {
            let paf = sub_matches.get_one::<String>("PAF").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mode = sub_matches.get_one::<String>("MODE").expect("error").as_str();
            let by = sub_matches.get_one::<String>("BY").expect("error").as_str();

            let pairs = *sub_matches.get_one::<u64>("PAIRS").expect("error");
            let reads = *sub_matches.get_one::<u64>("READS").expect("error");
            let frac = *sub_matches.get_one::<f64>("FRAC").expect("error");
            let bases = *sub_matches.get_one::<u64>("BASES").expect("error");
            let seed = *sub_matches.get_one::<u64>("SEED").expect("error");
            let min_quality = *sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let min_identity = *sub_matches.get_one::<f32>("MIN_IDENTITY").expect("error");
            let min_length = *sub_matches.get_one::<u32>("MIN_LENGTH").expect("error");
            let min_order = *sub_matches.get_one::<usize>("MIN_ORDER").expect("error");
            let max_order = *sub_matches.get_one::<usize>("MAX_ORDER").expect("error");
            let keep_comments = *sub_matches.get_one::<bool>("KEEP_COMMENTS").unwrap_or(&false);

            // translate args into Option targets (same as porec-downsample)
            let (target_reads, target_pairs, target_bases, frac_opt) = match mode {
                "pairs" => {
                    if pairs == 0 {
                        log::error!("--mode pairs requires --pairs > 0");
                        std::process::exit(1);
                    }
                    (None, Some(pairs), None, None)
                }
                "reads" => {
                    if reads == 0 {
                        log::error!("--mode reads requires --reads > 0");
                        std::process::exit(1);
                    }
                    (Some(reads), None, None, None)
                }
                "bases" => {
                    if bases == 0 {
                        log::error!("--mode bases requires --bases > 0");
                        std::process::exit(1);
                    }
                    (None, None, Some(bases), None)
                }
                "frac" => {
                    if !(0.0 < frac && frac <= 1.0) {
                        log::error!("--mode frac requires --frac in (0,1]");
                        std::process::exit(1);
                    }
                    (None, None, None, Some(frac))
                }
                _ => {
                    log::error!("Invalid --mode: {}", mode);
                    std::process::exit(1);
                }
            };
            let frac_by_pairs = match by {
                "pairs" => true,
                "reads" => false,
                _ => {
                    log::error!("Invalid --by: {}", by);
                    std::process::exit(1);
                }
            };

            let pt = PAFTable::new(&paf);
            pt.downsample_by_vpc(
                &output.to_string(),
                seed,
                min_quality,
                min_identity,
                min_length,
                min_order,
                max_order,
                target_reads,
                target_pairs,
                target_bases,
                frac_opt,
                frac_by_pairs,
                keep_comments,
            )
            .unwrap();
        }
        Some(("porec-downsample", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mode = sub_matches.get_one::<String>("MODE").expect("error").as_str();
            let by = sub_matches.get_one::<String>("BY").expect("error").as_str();

            let pairs = *sub_matches.get_one::<u64>("PAIRS").expect("error");
            let reads = *sub_matches.get_one::<u64>("READS").expect("error");
            let frac = *sub_matches.get_one::<f64>("FRAC").expect("error");

            let seed = *sub_matches.get_one::<u64>("SEED").expect("error");
            let min_mapq = *sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let min_order = *sub_matches.get_one::<usize>("MIN_ORDER").expect("error");
            let max_order = *sub_matches.get_one::<usize>("MAX_ORDER").expect("error");

            // translate args into Option targets
            let (target_reads, target_pairs, frac_opt) = match mode {
                "pairs" => {
                    if pairs == 0 {
                        log::error!("--mode pairs requires --pairs > 0");
                        std::process::exit(1);
                    }
                    (None, Some(pairs), None)
                }
                "reads" => {
                    if reads == 0 {
                        log::error!("--mode reads requires --reads > 0");
                        std::process::exit(1);
                    }
                    (Some(reads), None, None)
                }
                "frac" => {
                    if !(0.0 < frac && frac <= 1.0) {
                        log::error!("--mode frac requires --frac in (0,1]");
                        std::process::exit(1);
                    }
                    (None, None, Some(frac))
                }
                _ => {
                    log::error!("Invalid --mode: {}", mode);
                    std::process::exit(1);
                }
            };

            let frac_by_pairs = match by {
                "pairs" => true,
                "reads" => false,
                _ => {
                    log::error!("Invalid --by: {}", by);
                    std::process::exit(1);
                }
            };

            let mut prt = PoreCTable::new(&table);
            prt.downsample(
                &output.to_string(),
                seed,
                min_mapq,
                min_order,
                max_order,
                target_reads,
                target_pairs,
                frac_opt,
                frac_by_pairs,
            )
            .unwrap();
        }

        Some(("porec2depth", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let chromsizes = sub_matches.get_one::<String>("CHROMSIZES").expect("required");
            let window_size = sub_matches.get_one::<usize>("WINSIZE").expect("error");
            let mut step_size = *sub_matches.get_one::<usize>("STEPSIZE").expect("error");
            let min_mapq = *sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            if step_size == 0 {
                step_size = *window_size;
            }

            let mut prt = PoreCTable::new(&table);
            prt.to_depth(&chromsizes, *window_size, step_size, min_mapq, &output).unwrap();
        }
        Some(("paf2pairs", sub_matches)) => {
            let paf = sub_matches.get_one::<String>("PAF").expect("required");
            let chromsizes = sub_matches.get_one::<String>("CHROMSIZES").expect("required");
            let empty_string = String::new();
            let bed = sub_matches.get_one::<String>("BED").unwrap_or(&empty_string);
           
            let min_quality = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let min_identity = sub_matches.get_one::<f32>("MIN_IDENTITY").expect("error");
            let min_length = sub_matches.get_one::<u32>("MIN_LENGTH").expect("error");
            let max_edge = sub_matches.get_one::<u64>("MAX_EDGE").expect("error");
            let min_order = sub_matches.get_one::<usize>("MIN_ORDER").expect("error");
            let max_order = sub_matches.get_one::<u32>("MAX_ORDER").expect("error");
            let secondary = sub_matches.get_one::<bool>("SECONDARY").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let output_raw = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let output = output_raw.strip_suffix("/").unwrap_or(output_raw).to_string();

            let pt = PAFTable::new(&paf);
            let prefix = output.strip_suffix(".pairs").unwrap_or(&output);
            let prefix = prefix.strip_suffix(".pairs.gz").unwrap_or(prefix);
            let prefix: &str = prefix.strip_suffix(".pairs.pqs").unwrap_or(prefix);
            let table_output = format!("{}.porec.gz", prefix);
            pt.paf2table(&bed, &table_output, min_quality, min_identity, 
                            min_length, max_order, max_edge, *secondary, *threads).unwrap();
            let mut prt = PoreCTable::new(&table_output);
            if output.ends_with(".pqs") {
                let output = if output == "-" {
                    let output = table_output.strip_suffix(".porec.gz").unwrap_or(&table_output).strip_suffix(".con.gz").unwrap_or(&table_output).strip_suffix(".concatemer.gz").unwrap_or(&table_output);
                    format!("{}.pqs", output)
                } else {
                    output.to_string()
                };
                if Path::new(&output).exists() {
                    log::info!("Output directory {} already exists. Removing it first.", output);
                    std::fs::remove_dir_all(&output).unwrap();
                }
                prt.to_pairs_pqs(&chromsizes, &output, 1000000, *min_quality, *min_order, *max_order as usize, *threads).unwrap();
               
               
            } else {
                prt.to_pairs(&chromsizes, &output, *min_quality, *min_order, *max_order as usize).unwrap();
            }

        }
        Some(("pairs-downsample", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let percent = sub_matches.get_one::<f64>("PERCENT").expect("error");
            let number = sub_matches.get_one::<usize>("NUMBER").expect("error");
            let seed = sub_matches.get_one::<usize>("SEED").expect("error");  

            let min_mapq = *sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let threads = *sub_matches.get_one::<usize>("THREADS").expect("error");


            match (*percent, *number) {
                (0.0, 0) => {
                    log::error!("Either percent or number must be greater than 0.");
                    std::process::exit(1);
                }
                (0.0, _) => {
                    if Path::new(pairs).is_dir() {
                        let p = PQS::new(pairs);
                        match p.is_pqs() {
                            true => {
                                if Path::new(output).exists() {
                                    log::info!("Output directory {} already exists. Removing it first.", output);
                                    std::fs::remove_dir_all(output).unwrap();
                                }
                                // ✅ prefer using unified downsample for PQS (supports mapq + threads)
                                downsample_pqs(
                                    pairs,
                                    output,
                                    Some(*number as u64),
                                    None,
                                    *seed as u64,
                                    min_mapq,
                                    threads,
                                ).unwrap();
                            }
                            false => log::error!("The input directory is not a PQS directory."),
                        }
                    } else {
                        // Pairs text format path (no min_mapq support here yet)
                        Pairs::new(pairs).downsample(*number, *percent, *seed, output);
                    }
                }
                (_, _) => {
                    if Path::new(pairs).is_dir() {
                        let p = PQS::new(pairs);
                        match p.is_pqs() {
                            true => {
                                if Path::new(output).exists() {
                                    log::info!("Output directory {} already exists. Removing it first.", output);
                                    std::fs::remove_dir_all(output).unwrap();
                                }
                                downsample_pqs(
                                    pairs,
                                    output,
                                    None,
                                    Some(*percent),
                                    *seed as u64,
                                    min_mapq,
                                    threads,
                                ).unwrap();
                            }
                            false => log::error!("The input directory is not a PQS directory."),
                        }
                    } else {
                        Pairs::new(pairs).downsample(*number, *percent, *seed, output);
                    }
                }
            }

        }
        Some(("pairs-merge", sub_matches)) => {
            
            let files: Vec<_> = sub_matches.get_many::<String>("FILES").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            // get first file
            let first_file = files.first().unwrap();
     
            if Path::new(first_file).is_dir() {
                let p = PQS::new(&first_file);
                match p.is_pqs() {
                    true => {
                        if Path::new(&output).exists() {
                            log::info!("Output directory {} already exists. Removing it first.", output);
                            std::fs::remove_dir_all(&output).unwrap();
                        }
                        let _ = merge_pqs(files, &output);
                    },
                    false => {
                        log::error!("The first input directory is not a PQS directory.");
                    }
                }
            } else {
                merge_pairs(files, output)

            }
        }
        Some(("pairs-break", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let break_bed = sub_matches.get_one::<String>("BREAK_BED").expect("required");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            if Path::new(&pairs).is_dir() {
                let mut p = PQS::new(&pairs);
                match p.is_pqs() {
                    true => {
                        if Path::new(&output).exists() {
                            log::info!("Output directory {} already exists. Removing it first.", output);
                            std::fs::remove_dir_all(&output).unwrap();
                        }
                        let _ = p.break_contigs(&break_bed, &output);
                    },
                    false => {
                        log::error!("The input directory is not a PQS directory.");
                    }
                }
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.break_contigs(&break_bed, &output);
            }
        }
        Some(("pairs-dup", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let collapsed_list = sub_matches.get_one::<String>("COLLAPSED").expect("required");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            if Path::new(&pairs).is_dir() {
                let p = PQS::new(&pairs);
                match p.is_pqs() {
                    true => {
                        if Path::new(&output).exists() {
                            log::info!("Output directory {} already exists. Removing it first.", output);
                            std::fs::remove_dir_all(&output).unwrap();
                        }
                        let _ = p.dup(&collapsed_list, 123, &output);
                    },
                    false => {
                        log::error!("The input directory is not a PQS directory.");
                    }
                }
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.dup(&collapsed_list, 123, &output);
            }
        }
        Some(("pairs-filter", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let whitelist = sub_matches.get_one::<String>("WHITELIST").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let mut whitehash: HashSet<String> = HashSet::new();

            if whitelist != "none" {
                let reader = common_reader(&whitelist);
                let reader = BufReader::new(reader);
                for line in reader.lines() {
                    let line = line.unwrap();
                    whitehash.insert(line);
                }
            }

          
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            

            if Path::new(pairs).is_dir() {
               log::error!("The input directory is not a pairs or pairs.gz file.");
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.filter_by_mapq(*min_quality, &whitehash, &output);
            }
        }
        Some(("pairs-split", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let chunksize = sub_matches.get_one::<usize>("CHUNKSIZE").expect("error");
            let _output_prefix = sub_matches.get_one::<String>("OUTPUT").expect("error");

            if Path::new(&pairs).is_dir() {
                log::error!("The input directory is not a pairs or pairs.gz file.");
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.split(*chunksize);
            }
        }
        Some(("pairs2contacts", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let min_contacts = sub_matches.get_one::<u32>("MIN_CONTACTS").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let split_num = sub_matches.get_one::<u32>("SPLIT_NUM").expect("error");
            
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            if Path::new(&pairs).is_dir() {
                let p = PQS::new(&pairs);
                match p.is_pqs() {
                    true => {
                        if *split_num > 1 {
                            let _ = p.to_split_contacts(*min_contacts, *split_num, *min_quality, &output);
                        } else {
                            let _ = p.to_contacts(*min_contacts, *min_quality, &output);
                        }
                    },
                    false => {
                        log::error!("The input directory is not a PQS directory.");
                    }
                }
            } else {
                let mut _pairs = Pairs::new(&pairs);
                let contacts = if *split_num > 1 {
                    _pairs.to_split_contacts(*min_contacts, *split_num, *min_quality).unwrap()
                } else {
                    let mut _pairs = Pairs::new(&pairs);
                    _pairs.to_contacts(*min_contacts, *min_quality).unwrap()
                };
            
                contacts.write(&output);
            }
            log::info!("Contacts written to {}", output);
            
        }
        Some(("pairs2clm", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let binsize = sub_matches.get_one::<u32>("BINSIZE").expect("error");
            let min_contacts = sub_matches.get_one::<u32>("MIN_CONTACTS").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let no_output_split_contacts = sub_matches.get_one::<bool>("NO_OUTPUT_SPLIT_CONTACTS").expect("error");
            let output_depth = sub_matches.get_one::<bool>("OUTPUT_DEPTH").expect("error");
            // let low_memory = sub_matches.get_one::<bool>("LOW_MEMORY").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let disable_filter = sub_matches.get_one::<bool>("DISABLE_FILTER").expect("error");
            let max_depth_ratio = sub_matches.get_one::<f64>("MAX_DEPTH_RATIO").expect("error");
            let max_q0_ratio = sub_matches.get_one::<f64>("MAX_Q0_RATIO").expect("error");
            
            
            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let output_split_contacts = match no_output_split_contacts {
                true => false,
                false => true
            };
            
            // pairs is dir 
            if std::path::Path::new(&pairs).is_dir() {
                let p = PQS::new(&pairs);
                match p.is_pqs() {
                    true => { 
                        let _ = p.to_clm(*min_contacts, 
                                    *min_quality, &output, 
                                    output_split_contacts, 
                                    *output_depth, *binsize,
                                    *threads,
                                    *disable_filter,
                                    *max_depth_ratio,
                                    *max_q0_ratio,
                                    );
                    },
                    false => {
                        log::error!("The input directory is not a PQS directory.");
                    }
                }
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.to_clm(*min_contacts, *min_quality, 
                            &output, output_split_contacts, 
                            *output_depth, *binsize, *threads, true,
                            *disable_filter,
                            *max_depth_ratio
                        );
            }

        }
        Some(("pairs2depth", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let binsize = sub_matches.get_one::<u32>("BINSIZE").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            if Path::new(&pairs).is_dir() {
                let p = PQS::new(&pairs);
                match p.is_pqs() {
                    true => {
                        let _ = p.to_depth(*binsize, *min_quality, &output);
                    },
                    false => {
                        log::error!("The input directory is not a PQS directory.");
                    }
                }
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.to_depth(*binsize, *min_quality, &output);
            }
        }
        Some(("pairs2mnd", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            

            if Path::new(&pairs).is_dir() {
                let p = PQS::new(&pairs);
                match p.is_pqs() {
                    true => {
                        let _ = p.to_mnd(*min_quality, &output);
                    },
                    false => {
                        log::error!("The input directory is not a PQS directory.");
                    }
                }
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.to_mnd(*min_quality, &output).unwrap();
            }
        }
        Some(("pairs2pqs", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let chunksize = sub_matches.get_one::<usize>("CHUNKSIZE").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mut pairs = Pairs::new(&pairs);
            // if path of output existed, remove it first
            if Path::new(&output).exists() {
                log::info!("Output path {} existed, removing it first.", output);
                std::fs::remove_dir_all(&output).unwrap();
            }
            pairs.to_pqs(*chunksize, &output);
        }
        Some(("bam-chr2ctg", sub_m)) => {
            let input = sub_m.get_one::<String>("INPUT").unwrap();
            let bed = sub_m.get_one::<String>("BED").unwrap();
            let output = sub_m.get_one::<String>("OUTPUT").unwrap();
            let threads = *sub_m.get_one::<usize>("THREADS").unwrap();

            chr_to_ctg(input, bed, output, threads).unwrap();
        }
        Some(("porec-chr2ctg", sub_m)) => {
            let input = sub_m.get_one::<String>("INPUT").unwrap();
            let bed = sub_m.get_one::<String>("BED").unwrap();
            let output = sub_m.get_one::<String>("OUTPUT").unwrap();
            let threads = *sub_m.get_one::<usize>("THREADS").unwrap();

            let mut prt = PoreCTable::new(input);
            prt.chr_porec_to_contig_porec(bed, output, threads).unwrap();
        }
        Some(("pqs-chr2ctg", sub_m)) => {
            let input = sub_m.get_one::<String>("INPUT").unwrap();
            let bed = sub_m.get_one::<String>("BED").unwrap();
            let output = sub_m.get_one::<String>("OUTPUT").unwrap();
            let threads = *sub_m.get_one::<usize>("THREADS").unwrap();

            let _ = chr_pqs_to_contig_pqs(input, bed, output, threads);
        }
        Some(("pairs-prune", sub_matches)) => {
            let input = sub_matches.get_one::<String>("INPUT").expect("required");
            let prune_table = sub_matches.get_one::<String>("PRUNETABLE").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("required");
            let threads = *sub_matches.get_one::<usize>("THREADS").expect("error");
            if Path::new(&input).is_dir(){
                prune_pqs(input, prune_table, output, threads).unwrap();
            } else {
                prune_pairs(input, prune_table, output, threads).unwrap();
            }
            
        }
        Some(("bam-prune", sub_matches)) => {
            let input_bam = sub_matches.get_one::<String>("BAM").expect("required");
            let prune_table = sub_matches.get_one::<String>("PRUNETABLE").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("required");
            let threads = *sub_matches.get_one::<usize>("THREADS").expect("error");

            prune_bam(input_bam, prune_table, output, threads).unwrap();
        }
        Some(("pairs2bam", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");

            if Path::new(&pairs).is_dir() {
                // log::error!("The input is a directory, that is not a pairs or pairs.gz file.");
                let pqs = PQS::new(&pairs);
                let _ = pqs.to_bam(*min_quality, &output, *threads);
            } else {
                let mut pairs = Pairs::new(&pairs);
                pairs.to_bam(*min_quality, &output, *threads);
            }
        }

        Some(("pairs2porec", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let mapq = sub_matches.get_one::<u8>("MAPQ").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let empty_string = String::new();
            let enzyme_bed = sub_matches.get_one::<String>("ENZYME").unwrap_or(&empty_string);
            let min_edge_support = *sub_matches.get_one::<u32>("MIN_EDGE_SUPPORT").unwrap_or(&1);
            let min_clique_size = *sub_matches.get_one::<usize>("MIN_CLIQUE_SIZE").unwrap_or(&3);
            let max_comp_size = *sub_matches.get_one::<usize>("MAX_COMP_SIZE").unwrap_or(&1000);
            let threads = *sub_matches.get_one::<usize>("THREADS").unwrap_or(&4);

            if Path::new(&pairs).is_dir() {
                let pqs = PQS::new(&pairs);
                let _ = pqs.to_porec(*mapq, &output, enzyme_bed, min_edge_support, min_clique_size, max_comp_size, threads);
            } else {
                log::error!("The input is not a directory, that is not a pairs or pairs.gz file.");
            
            }
        }

        Some(("bam2paf", sub_matches)) => {
            let bam = sub_matches.get_one::<String>("BAM").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let is_secondary = sub_matches.get_one::<bool>("SECONDARY").expect("error");

            bam2paf(&bam, &output, *threads, *is_secondary);
        }

        
        Some(("bam2fastq", sub_matches)) => {
            let bams: Vec<_> = sub_matches.get_many::<String>("BAM").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");

            bam2fastq(&bams, &output, *threads);
        }

        Some(("bam2fasta", sub_matches)) => {
            let bam = sub_matches.get_one::<String>("BAM").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            bam2fasta(&bam, &output, *threads);
        }
        
        Some(("bamstat", sub_matches)) => {
            let bam: Vec<_> = sub_matches.get_many::<String>("BAM").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            bamstat(&bam,  &output, *threads);
        }
        Some(("hicbamstat", sub_matches)) => {
            let bam: Vec<_> = sub_matches.get_many::<String>("BAM").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            bamstat_hic(&bam,  &output, *threads);
        }
        Some(("porecbamstat", sub_matches)) => {
            let bam: Vec<_> = sub_matches.get_many::<String>("BAM").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            
            bamstat_porec(&bam,  &output, *threads);
        }
        Some(("bam2pairs", sub_matches)) => {
            let bam = sub_matches.get_one::<String>("BAM").expect("required");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
           
            if Path::new(&output).exists() {
                log::info!("Output path {} existed, removing it first.", output);
                std::fs::remove_dir_all(&output).unwrap();
            }

            if output.ends_with(".pqs") {
                bam2pqs(&bam, *min_quality, 1000000, &output, *threads).unwrap();
            } else {
                bam2pairs(&bam, *min_quality, &output, *threads);
            }
               
            

        }
        
        Some(("pairs-intersect", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let bed = sub_matches.get_one::<String>("BED").expect("required");
            let invert = sub_matches.get_one::<bool>("INVERT").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let edge_length = sub_matches.get_one::<u64>("EDGE_LENGTH").expect("error");
            let max_q0_ratio = sub_matches.get_one::<f64>("MAX_Q0_RATIO").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");


            if *invert {
                log::info!("Invert the table.");
            }

            if Path::new(&pairs).is_dir() {
                let p = PQS::new(&pairs);
                match p.is_pqs() {
                    true => {
                        
                        if output == "-" {
                            let output = pairs.strip_suffix(".pqs").unwrap_or(&pairs);
                            log::info!("Output path set to {}.", output);
                        }
                        if Path::new(&output).exists() {
                            log::info!("Output path {} existed, removing it first.", output);
                            std::fs::remove_dir_all(&output).unwrap();
                        }
                        let _ = p.intersect(&bed, *invert, *min_quality, *threads, *max_q0_ratio, &output);
                    },
                    false => {
                        log::error!("The input directory is not a PQS directory.");
                    }
                }
            } else {

                ThreadPoolBuilder::new()
                    .num_threads(*threads)
                    .build_global()
                    .unwrap();
                let mut pairs = Pairs::new(&pairs);
                pairs.intersect_multi_threads(&bed, *invert, *min_quality,
                                                 *edge_length, &output);
            }
        }

        Some(("chromsizes", sub_matches)) => {
            let fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let fasta = Fastx::new(&fasta);
            let mut writer = common_writer(&output);
            let chromsizes = fasta.get_chrom_size().unwrap();
            let mut chromsizes: Vec<_> = chromsizes.into_iter().collect();
            chromsizes.sort_by(|a, b| a.0.cmp(&b.0));
            for (k, v) in chromsizes {
                writer.write(format!("{}\t{}\n", k, v).as_bytes()).unwrap();
            }
        }

        Some(("prunepairs", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let prune = sub_matches.get_one::<String>("PRUNE").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mut pairs = Pairs::new(&pairs);
            let prunetable = PruneTable::new(&prune);
            let contigs:HashSet<ContigPair> = prunetable.contig_pairs().unwrap().into_iter().collect();

            pairs.remove_by_contig_pairs(contigs, &output).unwrap();

        }

        Some(("modbam2fq", sub_matches)) => {
            let input_bam = sub_matches.get_one::<String>("BAM").expect("required");
            let min_prob = sub_matches.get_one::<f32>("MIN_PROB").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            modbam2fastq(&input_bam, *min_prob, &output).unwrap();
        }
        Some(("phase-reads", sub_matches)) => {
            let input = sub_matches.get_one::<String>("INPUT").expect("required");
            let groups = sub_matches.get_one::<String>("GROUPS").expect("required");
            let format = sub_matches.get_one::<String>("FORMAT").expect("error").as_str();
            let min_quality = *sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let threads = *sub_matches.get_one::<usize>("THREADS").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mut actual_format = format.to_string();
            if actual_format == "auto" {
                if input.ends_with(".bam") || input.ends_with(".sam") {
                    actual_format = "bam".to_string();
                } else if input.ends_with(".paf") || input.ends_with(".paf.gz") {
                    actual_format = "paf".to_string();
                } else {
                    log::error!("Cannot determine format from file extension. Please specify --format.");
                    std::process::exit(1);
                }
            }

            match actual_format.as_str() {
                "bam" => {
                
                    cphasing::bam::phase_reads(input, groups, output, min_quality, threads);
                }
                "paf" => {
                    let pt = PAFTable::new(input);
                    pt.phase_reads(groups, output, min_quality).unwrap();
                }
                _ => {
                    log::error!("Unsupported format: {}", actual_format);
                    std::process::exit(1);
                }
            }
        }
        Some(("modfa", sub_matches)) => {
            let input_fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let bed = sub_matches.get_one::<String>("BED").expect("required");
            
            let min_frac = sub_matches.get_one::<f64>("MIN_FRAC").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            modify_fasta(&input_fasta, &bed, 
                            *min_frac, &output).unwrap();
        }

        Some(("optimize", sub_matches)) => {
            use rand::rngs::SmallRng;
            use rand::SeedableRng;
            use rand::prelude::SliceRandom;
            use hashbrown::HashMap;
        
            let split_contacts = sub_matches.get_one::<String>("SPLITCONTACTS").expect("required");
            let count_re = sub_matches.get_one::<String>("COUNTRE").expect("required");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let mutation_rate = sub_matches.get_one::<f64>("MUTATION").expect("error");

            let ngen = sub_matches.get_one::<usize>("NGEN").expect("error");
            let npop = sub_matches.get_one::<usize>("NPOP").expect("error");
            let resume = sub_matches.get_one::<bool>("RESUME").expect("error");
            let seed = sub_matches.get_one::<u64>("SEED").expect("error");
            let skip_ga = sub_matches.get_one::<bool>("SKIPGA").expect("error");
            let run_lkh = sub_matches.get_one::<bool>("RUNLKH").expect("error");

            let _output = Path::new(&count_re)
                            .file_stem()
                            .unwrap()
                            .to_str()
                            .unwrap();
            let output = Path::new(_output).with_extension("tour");

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let mut count_re = CountRE::new(&count_re);
            count_re.parse();

            let contigsizes = count_re.to_lengths();
            
            let initial_contigs = contigsizes.keys().cloned().collect::<Vec<String>>();
            
            let mut contigsizes_idx: IndexMap<usize, usize> = IndexMap::new();
            let mut contig2idx: HashMap<String, usize> = HashMap::new();
            let mut idx2contig: HashMap<usize, String> = HashMap::new();
            for (i, contig) in initial_contigs.iter().enumerate() {
                let size = contigsizes.get(contig).unwrap();
                contigsizes_idx.insert(i, (*size).try_into().unwrap());
                contig2idx.insert(contig.clone(), i);
                idx2contig.insert(i, contig.clone());
            }

            let white_list = contig2idx.keys().cloned().collect::<HashSet<String>>();
            let mut split_contacts = SplitContacts::read_from_file(&split_contacts, Some(&white_list)).unwrap();
          
            let contacts_map = split_contacts.to_contacts(&contig2idx);
            let num_contigs = contigsizes_idx.len();
            let mut contacts_vec: Vec<HashMap<usize, u32>> = vec![HashMap::new(); num_contigs];
            for (&(u, v), &weight) in contacts_map.iter() {
                
                contacts_vec[u].insert(v, weight as u32);
                contacts_vec[v].insert(u, weight as u32);
            }
            let signs = vec![true; contigsizes_idx.len()];
            let mut tour = Tour {
                contigs: contigsizes_idx.keys().cloned().collect(),
                signs: signs
            };

            if *resume && Path::new(&output).exists() {
                log::info!("Resume optimization from existing output directory: {}", output.display());
                // mv output to output.sav
               
                let input_tour = common_reader(&output.to_str().unwrap());
                // read the last line of the file
                let reader = BufReader::new(input_tour);
                let last_line = reader.lines().last().unwrap().unwrap();
                let _tour: Vec<String> = last_line
                                        .split_whitespace()
                                        .map(|x| x.parse::<String>().unwrap())
                                        .collect();
                
                // get signs and contigs from tour : tig1+, tig2-, ...
                let mut initial_tour: Vec<usize> = Vec::new();
                let mut signs = Vec::new();
                for contig in _tour.iter() {
                    let (contig, sign) = if contig.ends_with('+') {
                        (&contig[..contig.len()-1], true)
                    } else if contig.ends_with('-') {
                        (&contig[..contig.len()-1], false)
                    } else {
                        (contig.as_str(), true)
                    };
                    let contig_idx = contig2idx.get(contig).unwrap();
                    initial_tour.push(*contig_idx);
                    signs.push(sign);
                } 

                log::info!("Backup previous to {}.sav", output.display());
                std::fs::rename(&output, format!("{}.sav", output.display())).unwrap();

                
                tour = Tour {
                    contigs: initial_tour,
                    signs: vec![true; num_contigs],
                };
                
            } else {
                log::info!("Starting optimization with random tour.");
                let mut rng = SmallRng::seed_from_u64(*seed);
                tour.contigs.shuffle(&mut rng);
            }

            
            let mut tour = run_hybrid(
                &tour,
                contigsizes_idx, 
                contacts_vec, 
                *mutation_rate,
                *npop, 
                *ngen, 
                1,
                *run_lkh,
                !*skip_ga,
                *resume,
                *seed,
                threads.clone(),
                Some(&split_contacts),
                Some(&contig2idx),
            );
            // let matrix = ContactMatrix::new(&contigsizes_idx, contacts_vec);
            // let mut tour = run_lkh_optimizer_dual_node(
            //     &tour,
            //     contigsizes_idx.clone(), 
            //     &matrix,
            //     &split_contacts,
            //     &contig2idx,
            //     10,
            //     *seed,
            // );
            split_contacts.normalize_by_cis();
            let detailed_matrix = DetailedContactMatrix::new(split_contacts.to_detailed_contact_matrix(&contig2idx));

            robust_orient_contigs(&mut tour, &detailed_matrix, 10);

            let mut writer = common_writer(output.to_str().unwrap());

            for (i, &idx) in tour.contigs.iter().enumerate() {
                let name = &idx2contig[&idx];
                let sign = if tour.signs[i] { "+" } else { "-" };
                write!(writer, "{}{} ", name, sign).unwrap();
            }
            writeln!(writer).unwrap();

        }

        _ => {
            let input_arg = matches.subcommand().unwrap().0; 
            eprintln!("No such subcommand: {}", input_arg);
        },
    }
}

