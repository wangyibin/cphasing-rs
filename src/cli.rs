use std::ffi::OsString;
use std::path::PathBuf;

use clap::{arg, Arg, Command, Subcommand, value_parser};

const VERSION: &str = "0.0.6";

pub fn cli() -> Command {
    Command::new("cphasing")
        .about("Phasing and scaffolding based on Pore-C or Hi-C data")
        .subcommand_required(true)
        .version(VERSION)
        .arg_required_else_help(true)
        .allow_external_subcommands(true)
        .subcommand(
            Command::new("aligner")
                .about("align methylation reads")
                .arg(arg!(<FASTA> "fasta"))
                .arg(arg!(<BAM> "align bam from `dorado`"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("10"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                .arg_required_else_help(true),
        ).subcommand(
            Command::new("cutsite")
                .about("cut restriction site on pore-c reads")
                .arg(arg!(<FASTQ> "pore-c reads"))
                .arg(
                    Arg::new("PATTERN")
                        .long("pattern")
                        .value_parser(value_parser!(String))
                        .default_value("GATCGATC"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("paf2table")
                .about("convert paf to pore-c table")
                .arg(arg!(<PAF> "paf"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                .arg(
                    Arg::new("MIN_MAPQ")
                        .long("min-mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1"))
                .arg(
                    Arg::new("MIN_IDENTITY")
                        .long("min-identity")
                        .short('p')
                        .value_parser(value_parser!(f32))
                        .default_value("0.75"))
                .arg(
                    Arg::new("MIN_LENGTH")
                        .long("min-length")
                        .short('l')
                        .value_parser(value_parser!(u32))
                        .default_value("30"))
                .arg(
                    Arg::new("MAX_ORDER")
                        .long("max-order")
                        .short('m')
                        .value_parser(value_parser!(u32))
                        .default_value("50")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("porec2pairs")
                .about("convert pore-c table to pairs")
                .arg(arg!(<TABLE> "pore-c table"))
                .arg(arg!(<CHROMSIZES> "chromsizes"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
        )
        .subcommand(
            Command::new("pairs2mnd")
                .about("convert pairs to mnd file")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
        )
        .subcommand(
            Command::new("chromsizes")
                .about("generate chromsizes file")
                .arg(arg!(<FASTA> "fasta"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))

        )
        .subcommand(
            Command::new("prunepairs")
                .about("prune pairs")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(arg!(<PRUNE> "Prune contigs"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                
        )
        .subcommand(
            Command::new("modbam2fq")
                .about("convert modified bam to fastq with modified base")
                .arg(arg!(<BAM> "modified bam"))
                .arg(
                    Arg::new("MIN_PROB")
                        .long("min-prob")
                        .short('p')
                        .value_parser(value_parser!(f32))
                        .default_value("0.5"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
        )
        .subcommand(
            Command::new("modfa")
                .about("modify fasta by bedMethy file")
                .arg(arg!(<FASTA> "fasta"))
                .arg(arg!(<BED> "bed"))
                .arg(
                    Arg::new("MIN_SCORE")
                        .long("min-score")
                        .short('s')
                        .value_parser(value_parser!(u32))
                        .default_value("1"))
                .arg(
                    Arg::new("MIN_FRAC")
                        .long("min-frac")
                        .short('f')
                        .value_parser(value_parser!(f32))
                        .default_value("0.5"))
                .arg(
                    Arg::new("BED_FORMAT")
                        .long("bed-fmt")
                        .short('f')
                        .default_value("bedMethy")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
        )
        .subcommand(
            Command::new("optimize")
            .about("optimize contigs")
            .arg(arg!(<SCORE> "score"))
            .arg(
                Arg::new("OUTPUT")
                    .long("output")
                    .short('o')
                    .value_parser(value_parser!(String))
                    .default_value("-"))
      )
        
}
