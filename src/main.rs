#[macro_use] extern crate scan_fmt;
use indexmap::IndexMap;
use cphasing::aligner::read_bam;
use cphasing::bam::split_bam;
use cphasing::cli::cli;
use cphasing::core::{ 
    BaseTable,  common_reader, 
    common_writer, ContigPair,
    check_program};
use cphasing::count_re::CountRE;
// use cphasing::contacts::Contacts;
use cphasing::cutsite::cut_site;
use cphasing::fastx::{ Fastx, split_fastq };
use cphasing::methy::{ modbam2fastq, modify_fasta };
use cphasing::optimize::ContigScoreTable;
use cphasing::optimize::SimulatedAnnealing;
use cphasing::paf::PAFTable;
use cphasing::pairs::Pairs;
use cphasing::porec::{
        PoreCTable, merge_porec_tables };
use cphasing::kprune::{ PruneTable, KPruner };
use cphasing::prune::{ Pruner };
use cphasing::simulation::{ 
        simulation_from_split_read, simulate_porec };
use std::collections::{ HashMap, HashSet };
use std::io::Write; 
use chrono::Local;
use env_logger::Builder;
use log::LevelFilter;


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
            let min_mapq = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            check_program("minimap2");
            match import_format.as_str() {
                "bam" => {
                    read_bam(&input, *min_mapq, &output);
                }
                "paf" => {
                    read_paf(&input, *min_mapq, &output);
                }
                _ => {
                    eprintln!("No such format.");
                }
            }

        }
        Some(("kprune", sub_matches)) => {
            use rayon::ThreadPoolBuilder;
            use rayon::prelude::*;
            use std::io::BufReader;
            use std::io::BufRead;
            use std::sync::{Arc, Mutex};
            let alleletable = sub_matches.get_one::<String>("ALLELETABLE").expect("required");
            let contacts = sub_matches.get_one::<String>("CONTACTS").expect("required");
            // let count_re = sub_matches.get_one::<String>("COUNT_RE").expect("required");
            let prunetable = sub_matches.get_one::<String>("PRUNETABLE").expect("required");
            let method = sub_matches.get_one::<String>("METHOD").expect("error");
            let normalization_method = sub_matches.get_one::<String>("NORMALIZATION_METHOD").expect("error");
            let whitelist = sub_matches.get_one::<String>("WHITELIST").expect("error");
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

            let mut whitehash2: HashSet<&String> = HashSet::new();
            for x in whitehash.iter() {
                whitehash2.insert(x);
            }
          
            // let mut first_cluster_hashmap: HashMap<String, HashSet<String>> = HashMap::new();
            let first_cluster_hashmap = Arc::new(Mutex::new(HashMap::new()));
            if first_cluster != "none" {
                let reader = common_reader(&first_cluster);
                let reader = BufReader::new(reader);
                // for line in reader.lines() {
                //     let line = line.unwrap();
                //     // split tab
                //     let mut line = line.split("\t");
                //     let group = line.next().unwrap().to_string();
                //     let count = line.next().unwrap().to_string();
                //     // split next with space
                //     let contigs: HashSet<_> = line.next().unwrap().split(" ").collect();
                   
                //     let contigs: HashSet<String> = contigs.par_iter().map(|x| x.to_string()).collect();
                //     if whitehash.len() > 0 {
                //         let contigs: HashSet<String> = contigs.intersection(&whitehash).map(|x| x.to_string()).collect();
                //     }
                //     first_cluster_hashmap.insert(group, contigs);
                // }
                let lines: Vec<String> = reader.lines().map(|line| line.unwrap()).collect();
                lines.into_par_iter().for_each(|line| {
                    // split tab
                    let mut line = line.split("\t");
                    let group = line.next().unwrap().to_string();
                    let _count = line.next().unwrap().to_string();
                    // split next with space
                    let contigs: HashSet<_> = line.next().unwrap().split(" ").collect();
                   
                    let mut contigs: HashSet<String> = contigs.into_par_iter().map(|x| x.to_string()).collect();
                    if whitehash.len() > 0 {
                        contigs = contigs.intersection(&whitehash).map(|x| x.to_string()).collect();
                    }
                    let mut first_cluster_hashmap = first_cluster_hashmap.lock().unwrap();
                    first_cluster_hashmap.insert(group, contigs);
                });
            }

            // Arc::new(Mutex::new(HashMap::new())) to HashMap
            let first_cluster_hashmap = Arc::clone(&first_cluster_hashmap);
            let first_cluster_hashmap = first_cluster_hashmap.lock().unwrap();
            let first_cluster_hashmap = first_cluster_hashmap.clone();

            let mut writer = common_writer(&prunetable);
            // let mut wtr = csv::WriterBuilder::new()
            //     .delimiter(b'\t')
            //     .has_headers(false)
            //     .from_writer(writer);

            if first_cluster != "none" {
            
                log::set_max_level(log::LevelFilter::Off);
                let mut kpruners: HashMap<String, KPruner> = first_cluster_hashmap
                            .keys().map(|x| (x.to_string(), KPruner::new(
                                &alleletable, &contacts, &prunetable,
                                normalization_method))).collect();
                
                for (k, v) in first_cluster_hashmap {
                    log::info!("Pruning cluster `{}`", k);
                    let kpruner = kpruners.get_mut(&k).unwrap();
                    let _whitehash: HashSet<String> = v.into_iter().collect();
                    let mut _whitehash2: HashSet<&String> = HashSet::new();
                    for x in _whitehash.iter() {
                        _whitehash2.insert(&x);
                    }
                    kpruner.prune(&method.as_str(), &_whitehash2, &mut writer );
                }
                log::set_max_level(log::LevelFilter::Info);

                // write kpruners prunetable into a file
                let mut allelic_counts: u32 = 0;
                let mut cross_allelic_counts: u32 = 0;
                for (_, v) in kpruners {
                    allelic_counts += v.allelic_counts;
                    cross_allelic_counts += v.cross_allelic_counts;
                    // v.prunetable.write(&mut wtr);
                }
                log::info!("Allelic counts: {}", allelic_counts);
                log::info!("Cross-allelic counts: {}", cross_allelic_counts);
                
            } else {
                let mut kpruner = KPruner::new(&alleletable, &contacts, &prunetable, normalization_method);
                kpruner.prune(&method.as_str(), &whitehash2, &mut writer);
                // kpruner.prunetable.write(&mut wtr);
            }
            
            log::info!("Allelic and cross-allelic information written into `{}`", prunetable);

        }
        Some(("prune", sub_matches)) => {
            use rayon::ThreadPoolBuilder;
            use rayon::prelude::*;
            use std::io::BufReader;
            use std::io::BufRead;
            use std::sync::{Arc, Mutex};
            let alleletable = sub_matches.get_one::<String>("ALLELETABLE").expect("required");
            let allele_strand_table = sub_matches.get_one::<String>("ALLELESTRANDTABLE").expect("required");
            let contacts = sub_matches.get_one::<String>("CONTACTS").expect("required");
            let prunetable = sub_matches.get_one::<String>("PRUNETABLE").expect("required");
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
        Some(("splitfastq", sub_matches)) => {
            let input_fastq = sub_matches.get_one::<String>("FASTQ").expect("required");
            let output_prefix = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let record_num = sub_matches.get_one::<usize>("RECORD_NUM").expect("error");

            split_fastq(&input_fastq, &output_prefix, *record_num).unwrap();
        }
        Some(("slidefastq", sub_matches)) => {
            let input_fastq = sub_matches.get_one::<String>("FASTQ").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let window = sub_matches.get_one::<u64>("WINDOW").expect("error");
            let step = sub_matches.get_one::<u64>("STEP").expect("error");
            let min_length = sub_matches.get_one::<u64>("MIN_LENGTH").expect("error");

            let fa = Fastx::new(&input_fastq);
            let _ = fa.slide(&output, *window, *step, *min_length);
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
            
            let mut count_re = CountRE::new(&output);
            let fasta = Fastx::new(&input_fasta);
            let contigsizes = fasta.get_chrom_size().unwrap();
            let counts = fasta.count_re(&pattern).unwrap();

            // counts + 1 
            let counts: IndexMap<String, u64> = counts.into_iter()
                                                    .map(|(k, v)| (k, v + 1))
                                                    .collect();
            // filter counts which are less than min re
            let counts: IndexMap<String, u64> = counts.into_iter()
                                                    .filter(|(_, v)| *v >= *min_re)
                                                    .collect();
            count_re.from_hashmap(counts, contigsizes);
            count_re.write(&output);
        }
        Some(("cutsite", sub_matches)) => {
            let fastq = sub_matches.get_one::<String>("FASTQ").expect("required");
            let pattern = sub_matches.get_one::<String>("PATTERN").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            cut_site(fastq.to_string(), pattern.as_bytes(), "-".to_string()).unwrap();
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

            let pt = PAFTable::new(&paf);

            pt.paf2table(bed, output, min_quality, min_identity, 
                            min_length, max_order, max_edge).unwrap();

        }
        Some(("porec2pairs", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let chromsizes = sub_matches.get_one::<String>("CHROMSIZES").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let min_order = sub_matches.get_one::<usize>("MIN_ORDER").expect("error");
            let max_order = sub_matches.get_one::<usize>("MAX_ORDER").expect("error");
            let prt = PoreCTable::new(&table);

            prt.to_pairs(&chromsizes, &output, *min_quality, *min_order, *max_order).unwrap();
        }
        Some(("porec-merge", sub_matches)) => {
            let tables: Vec<_> = sub_matches.get_many::<String>("TABLES").expect("required").collect();
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            merge_porec_tables(tables, output)
        }
        Some(("porec-intersect", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let bed = sub_matches.get_one::<String>("BED").expect("required");
            let invert = sub_matches.get_one::<bool>("INVERT").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let mut prt = PoreCTable::new(&table);
            if *invert {
                log::info!("Invert the table.");
            }
            prt.intersect(&bed, *invert, &output);
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
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let pt = PAFTable::new(&paf);
            
            let prefix = output.strip_suffix(".pairs").unwrap_or(&output);
            let prefix = prefix.strip_suffix(".pairs.gz").unwrap_or(prefix);
            let table_output = format!("{}.porec.gz", prefix);
            pt.paf2table(&bed, &table_output, min_quality, min_identity, 
                            min_length, max_order, max_edge).unwrap();
            let prt = PoreCTable::new(&table_output);
            prt.to_pairs(&chromsizes, &output, *min_quality, *min_order, *max_order as usize).unwrap();

        }
        Some(("pairs-filter", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let mut pairs = Pairs::new(&pairs);

            pairs.filter_by_mapq(*min_quality, &output);
        }
        Some(("pairs2contacts", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let min_contacts = sub_matches.get_one::<u32>("MIN_CONTACTS").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let split_num = sub_matches.get_one::<u32>("SPLIT_NUM").expect("error");
            let mut pairs = Pairs::new(&pairs);
            let contacts = if *split_num > 1 {
                pairs.to_split_contacts(*min_contacts, *split_num, *min_quality).unwrap()
            } else {
                pairs.to_contacts(*min_contacts, *min_quality).unwrap()
            };
           
            contacts.write(&output);
            log::info!("Contacts written to {}", output);
            
        }
        Some(("pairs2clm", sub_matches)) => {
            use rayon::ThreadPoolBuilder;
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            let min_contacts = sub_matches.get_one::<u32>("MIN_CONTACTS").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let no_output_split_contacts = sub_matches.get_one::<bool>("NO_OUTPUT_SPLIT_CONTACTS").expect("error");
            let low_memory = sub_matches.get_one::<bool>("LOW_MEMORY").expect("error");
            let threads = sub_matches.get_one::<usize>("THREADS").expect("error");
            let mut pairs = Pairs::new(&pairs);

            ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            let output_split_contacts = match no_output_split_contacts {
                true => false,
                false => true
            };

            pairs.to_clm(*min_contacts, *min_quality, &output, output_split_contacts, *low_memory);
            // let contacts = Contacts::from_clm(&output);
            // contacts.write(&contacts.file);
        }
        Some(("pairs2mnd", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let mut pairs = Pairs::new(&pairs);

            pairs.to_mnd(*min_quality, &output).unwrap();
        }
        
        Some(("pairs2bam", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let mut pairs = Pairs::new(&pairs);

            pairs.to_bam(*min_quality, &output);
        }
        
        Some(("pairs-intersect", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let bed = sub_matches.get_one::<String>("BED").expect("required");
            let invert = sub_matches.get_one::<bool>("INVERT").expect("error");
            let min_quality = sub_matches.get_one::<u8>("MIN_QUALITY").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let mut pairs = Pairs::new(&pairs);

            if *invert {
                log::info!("Invert the table.");
            }
            pairs.intersect(&bed, *invert, *min_quality, &output);
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

        Some(("modfa", sub_matches)) => {
            let input_fasta = sub_matches.get_one::<String>("FASTA").expect("required");
            let bed = sub_matches.get_one::<String>("BED").expect("required");
            let min_score = sub_matches.get_one::<u32>("MIN_SCORE").expect("error");
            let min_frac = sub_matches.get_one::<f32>("MIN_FRAC").expect("error");
            let bed_fmt = sub_matches.get_one::<String>("BED_FORMAT").expect("error");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            modify_fasta(&input_fasta, &bed, *min_score, 
                            *min_frac, &bed_fmt, &output).unwrap();
        }

        Some(("optimize", sub_matches)) => {
            let score = sub_matches.get_one::<String>("SCORE").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let cst = ContigScoreTable::new(&score);
            let co = cst.read();
            let mut sa = SimulatedAnnealing::new(co, 2000.0, 0.9999, 0.01, 100000000);
            let best = sa.run();

            best.save(&output);

        }
        _ => {
            eprintln!("{:?}", matches.subcommand());
            eprintln!("No such subcommand.");
        },
    }
}
