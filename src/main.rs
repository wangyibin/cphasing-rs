use cphasing::cli::cli;
use cphasing::core::BaseTable;
use cphasing::cutsite::cut_site;
use cphasing::paf::PAFTable;
use cphasing::porec::PoreCTable;
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

            cut_site(fastq.to_string(), pattern.to_string(), "-".to_string());
        }
        Some(("paf2table", sub_matches)) => {
            let paf = sub_matches.get_one::<String>("PAF").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");
            
            let min_quality = sub_matches.get_one::<u8>("MIN_MAPQ").expect("error");
            let min_identity = sub_matches.get_one::<f32>("MIN_IDENTITY").expect("error");
            let min_length = sub_matches.get_one::<u32>("MIN_LENGTH").expect("error");

            let pt = PAFTable::new(&paf);

            pt.paf2table(output, min_quality, min_identity, min_length);

        }
        Some(("porec2pairs", sub_matches)) => {
            let table = sub_matches.get_one::<String>("TABLE").expect("required");
            let chromsizes = sub_matches.get_one::<String>("CHROMSIZES").expect("required");
            let output = sub_matches.get_one::<String>("OUTPUT").expect("error");

            let prt = PoreCTable::new(&table);

            prt.to_pairs(&chromsizes, &output);
        }
        _ => unreachable!(),
    }
}
