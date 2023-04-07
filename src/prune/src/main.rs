#! [allow (dead_code)]
use clap::{ arg, Arg, Command, value_parser };
use std::error::Error;
use std::io;
use std::process;
use std::path::PathBuf;

use prune::core::BaseTable;
use prune::core::{AlleleTable, PairTable, pruner};


fn main() {
    let matches = Command::new("prune")
        .version("0.1")
        .author("Yibin Wang")
        .about("prune allelic signal")
        // .arg(Arg::new("verbose")
        //     .short('v')
        //     .help("verbosity level"),)
        .arg(arg!(--allele <PATH> "Path to allele table.")
                .short('a')
                .required(true)
                .value_parser(value_parser!(String)))
        .arg(arg!(--pairtable <PATH> "Path to pair table.")
                .short('p')
                .required(true)
                .value_parser(value_parser!(String)))
        .arg(arg!(--count <PATH< "Path to count REs.")
                .short('c')
                .required(true)
                .value_parser(value_parser!(String)))
        .get_matches();
    
    let allele = matches.get_one::<String>("allele").unwrap();
    let count_re = matches.get_one::<String>("count").unwrap();
    let pairtable = matches.get_one::<String>("pairtable").unwrap();
    
    let at = AlleleTable::new(&allele);
    let pt = PairTable::new(&pairtable);
    match pruner(at, pt) {
        Ok(()) => eprintln!("Success."),
        Err(e) => eprintln!("Failaled to write data: {}", {e}),
    }
}







// #[cfg(test)]
// mod tests {
//     use super::*;
// }
