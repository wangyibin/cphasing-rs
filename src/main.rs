use cphasing::cli::cli;
use cphasing::core::{ BaseTable, common_writer, ContigPair };
use cphasing::cutsite::cut_site;
use cphasing::fastx::Fastx;
use cphasing::paf::PAFTable;
use cphasing::pairs::Pairs;
use cphasing::porec::PoreCTable;
use cphasing::prune::PruneTable;
use std::collections::HashSet;
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
        Some(("cutsite", sub_matches)) => {
            let fastq = sub_matches.get_one::<String>("FASTQ").expect("required");
            let pattern = sub_matches.get_one::<String>("PATTERN").expect("error");

            cut_site(fastq.to_string(), pattern.as_bytes(), "-".to_string());
        }
        Some(("paf2table", sub_matches)) => {
            let paf = sub_matches.get_one::<String>("PAF").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let min_quality = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let min_identity = sub_matches.get_one::<f32>("MIN_IDENTITY").expect("error");
            let min_length = sub_matches.get_one::<u32>("MIN_LENGTH").expect("error");
            let max_order = sub_matches.get_one::<u32>("MAX_ORDER").expect("error");

            let pt = PAFTable::new(&paf);

            pt.paf2table(output, min_quality, min_identity, min_length, max_order);

        }
        Some(("porec2pairs", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let chromsizes = sub_matches.get_one::<String>("CHROMSIZES").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let prt = PoreCTable::new(&table);

            prt.to_pairs(&chromsizes, &output);
        }
        Some(("pairs2mnd", sub_matches)) => {
            let pairs = sub_matches.get_one::<String>("PAIRS").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let mut pairs = Pairs::new(&pairs);

            pairs.to_mnd(&output);
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

            pairs.remove_by_contig_pairs(contigs, &output);

        }
        _ => {
            eprintln!("No such subcommand.");
        },
    }
}
