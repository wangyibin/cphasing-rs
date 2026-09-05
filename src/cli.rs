use std::ffi::OsString;
use std::io::IsTerminal;
use std::path::PathBuf;

use clap::{
    Arg, ArgAction, ColorChoice, Command, Subcommand, arg,
    builder::{
        ArgPredicate, Styles,
        styling::{AnsiColor, Effects},
    },
    value_parser,
};
use cool2mcool::mcool::BASE_RESOLUTION;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const COOL2MCOOL_DEFAULT_RESOLUTIONS: &str =
    "2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000";
const COOL2MCOOL_DEFAULT_KR_MIN_RESOLUTION: &str = "5000";

#[allow(dead_code)]
fn non_negative_f64(value: &str) -> Result<f64, String> {
    let parsed = value.parse::<f64>().map_err(|error| error.to_string())?;
    if !parsed.is_finite() || parsed < 0.0 {
        return Err("value must be finite and non-negative".to_string());
    }
    Ok(parsed)
}

#[allow(dead_code)]
fn unit_interval_f64(value: &str) -> Result<f64, String> {
    let parsed = non_negative_f64(value)?;
    if parsed > 1.0 {
        return Err("value must be between 0 and 1".to_string());
    }
    Ok(parsed)
}

#[allow(dead_code)]
fn positive_unit_interval_f64(value: &str) -> Result<f64, String> {
    let parsed = unit_interval_f64(value)?;
    if parsed == 0.0 {
        return Err("value must be greater than 0".to_string());
    }
    Ok(parsed)
}

#[allow(dead_code)]
fn positive_usize(value: &str) -> Result<usize, String> {
    let parsed = value.parse::<usize>().map_err(|error| error.to_string())?;
    if parsed == 0 {
        return Err("value must be greater than 0".to_string());
    }
    Ok(parsed)
}

fn orientation_window(value: &str) -> Result<usize, String> {
    let parsed = positive_usize(value)?;
    if parsed > 16 {
        return Err("orientation window must be between 1 and 16".to_string());
    }
    Ok(parsed)
}

fn disabled_or_block_span(value: &str) -> Result<usize, String> {
    let parsed = value.parse::<usize>().map_err(|error| error.to_string())?;
    if parsed == 1 {
        return Err("block span must be 0 (disabled) or at least 2".to_string());
    }
    Ok(parsed)
}

fn cool2mcool_resolution(value: &str) -> Result<u64, String> {
    let resolution = value.parse::<u64>().map_err(|error| {
        format!("expected a positive resolution in base pairs, got {value:?}: {error}")
    })?;
    if resolution == 0 || resolution % BASE_RESOLUTION != 0 {
        return Err(format!(
            "resolution must be a positive multiple of {BASE_RESOLUTION} bp"
        ));
    }
    Ok(resolution)
}

fn cool2mcool_compression_level(value: &str) -> Result<u8, String> {
    let level = value.parse::<u8>().map_err(|error| {
        format!("expected a gzip compression level from 1 to 9, got {value:?}: {error}")
    })?;
    if (1..=9).contains(&level) {
        Ok(level)
    } else {
        Err("gzip compression level must be between 1 and 9".to_string())
    }
}

const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Yellow.on_default());

const ROOT_HELP_TEMPLATE: &str = "\
{before-help}{about-with-newline}
{usage-heading} {usage}{after-help}{options}";

struct HelpGroup {
    heading: &'static str,
    commands: &'static [&'static str],
}

const HELP_GROUPS: &[HelpGroup] = &[
    HelpGroup {
        heading: "PAF",
        commands: &[
            "methalign",
            "paf2concat",
            "paf2depth",
            "paf2pairs",
            "paf-downsample",
        ],
    },
    HelpGroup {
        heading: "CONCATEMER",
        commands: &[
            "concat2pqs",
            "concat2pairs",
            "concat2depth",
            "concat2reads",
            "concat-split",
            "concat-merge",
            "concat-intersect",
            "concat-downsample",
            "concat-break",
            "concat-dup",
            "concat-chr2ctg",
        ],
    },
    HelpGroup {
        heading: "PAIRS",
        commands: &[
            "pairs-downsample",
            "pairs-merge",
            "pairs-filter",
            "pairs-break",
            "pairs-dup",
            "pairs-split",
            "pairs2contacts",
            "pairs2clm",
            "pairs2depth",
            "pairs2mnd",
            "pairs2bam",
            "pairs2porec",
            "pairs-chr2ctg",
            "pairs-prune",
            "pairs-intersect",
            "splitcontacts",
            "cool2mcool",
        ],
    },
    HelpGroup {
        heading: "BAM",
        commands: &[
            "bam2paf",
            "bam2pairs",
            "bam-chr2ctg",
            "bam-prune",
            "bam2fastq",
            "bam2fasta",
            "modbam2fq",
            "bamstat",
            "splitbam",
        ],
    },
    HelpGroup {
        heading: "ASSEMBLY & GRAPH",
        commands: &[
            "kprune",
            "clm",
            "mergeclm",
            "splitclm",
            "gfa",
            "optimize",
            "phase-reads",
        ],
    },
    HelpGroup {
        heading: "SEQUENCE, RESTRICTION ENZYME & SIMULATION",
        commands: &[
            "splitfastq",
            "extract-fasta",
            "slidefastq",
            "slidefasta",
            "slide2raw",
            "digest",
            "count_re",
            "cutsite",
            "chromsizes",
            "modfa",
            "simulator",
        ],
    },
];

fn use_colored_help() -> bool {
    std::env::var_os("NO_COLOR").is_none() && std::io::stdout().is_terminal()
}

fn append_styled(output: &mut String, text: &str, style: &str, use_color: bool) {
    if use_color {
        output.push_str(style);
    }
    output.push_str(text);
    if use_color {
        output.push_str("\x1b[0m");
    }
}

fn append_command_group(
    output: &mut String,
    heading: &str,
    commands: &mut Vec<&Command>,
    use_color: bool,
) {
    if commands.is_empty() {
        return;
    }

    commands.sort_by_key(|command| (command.get_display_order(), command.get_name()));
    append_styled(
        output,
        &format!("{heading}:\n"),
        "\x1b[1;32m",
        use_color,
    );

    let width = commands
        .iter()
        .map(|command| command.get_name().len())
        .max()
        .unwrap_or_default();
    for command in commands {
        output.push_str("  ");
        append_styled(output, command.get_name(), "\x1b[1;36m", use_color);
        output.push_str(&" ".repeat(width - command.get_name().len() + 2));
        output.push_str(
            &command
                .get_about()
                .map(|about| about.to_string())
                .unwrap_or_default(),
        );
        output.push('\n');
    }
    output.push('\n');
}

fn grouped_subcommands_help(command: &Command) -> String {
    let use_color = use_colored_help();
    let visible_commands = command
        .get_subcommands()
        .filter(|subcommand| !subcommand.is_hide_set())
        .collect::<Vec<_>>();
    let mut output = String::new();

    for group in HELP_GROUPS {
        let mut commands = visible_commands
            .iter()
            .copied()
            .filter(|command| group.commands.contains(&command.get_name()))
            .collect::<Vec<_>>();
        append_command_group(&mut output, group.heading, &mut commands, use_color);
    }

    let mut ungrouped = visible_commands
        .into_iter()
        .filter(|command| {
            !HELP_GROUPS
                .iter()
                .any(|group| group.commands.contains(&command.get_name()))
        })
        .collect::<Vec<_>>();
    append_command_group(&mut output, "OTHER COMMANDS", &mut ungrouped, use_color);
    append_styled(&mut output, "OPTIONS:\n", "\x1b[1;32m", use_color);

    output
}

pub fn cli() -> Command {
    let command = Command::new("cphasing")
        .color(ColorChoice::Auto)
        .about("Phasing and scaffolding based on Pore-C or Hi-C data")
        .subcommand_required(true)
        .version(VERSION)
        .styles(STYLES)
        .arg_required_else_help(true)
        .allow_external_subcommands(true)
        .subcommand(
            Command::new("aligner")
                .hide(true)
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
                    Arg::new("MIN_PROB")
                        .long("min-prob")
                        .short('p')
                        .value_parser(value_parser!(f32))
                        .default_value("0.75"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        ).subcommand(
            Command::new("realign")
                .hide(true)
                .about("rescue secondary alignment by high-quality alignments")
                .arg(arg!(<PAF> "paf from minimap2 with secondary"))
                .arg(
                    Arg::new("MIN_MAPQ")
                        .long("min-mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1"))
                .arg(
                    Arg::new("FORMAT")
                        .long("format")
                        .short('f')
                        .value_parser(value_parser!(String))
                        .default_value("paf")
                        .help("input format")
                        )
                .arg(
                    Arg::new("CONTACTS")
                        .long("contacts")
                        .short('c')
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        .help("contacts file for pore-c or hi-c data")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("methalign")
                .display_order(10)
                     .about("refine alignments by methylation or without methylation")
                        .arg(arg!(<BAM> "align bam from `dorado` "))
                        .arg(
                            Arg::new("FASTA")
                                .long("fasta")
                                .short('f')
                                .value_parser(value_parser!(String))
                                .default_value("none")
                                .help("fasta file for reference, if not set, only use sequence matching"))
                        .arg(
                            Arg::new("BEDGRAPH")
                                .long("bed")
                                .short('b')
                                .value_parser(value_parser!(String))
                                .default_value("none")
                                .help("bed file for regions to refine, if not set, use whole genome"))
                        .arg(
                            Arg::new("MATCH_SCORE")
                                .long("match-score")
                                .short('A')
                                .value_parser(value_parser!(i32))
                                .default_value("0")
                                .help("match score of methylation sites")
                        )
                        .arg(
                            Arg::new("REF_PENALTY")
                                .long("ref-penalty")
                                .value_parser(value_parser!(i32))
                                .default_value("2")
                                .help("penalty for mismatch methylation on reference")
                            )
                        .arg(
                            Arg::new("READ_PENALTY")
                                .long("read-penalty")
                                .value_parser(value_parser!(i32))
                                .default_value("2")
                                .help("penalty for mismatch methylation on read")
                            )
                        .arg(
                            Arg::new("REF_PROB_CUTOFF")
                                .long("ref-prob-cutoff")
                                .short('r')
                                .value_parser(value_parser!(f64))
                                .default_value("50.0")
                                .help("cutoff of probability for reference methylation, 0-100"))
                        .arg(
                            Arg::new("PROB_CUTOFF")
                                .long("prob-cutoff")
                                .short('c')
                                .value_parser(value_parser!(u8))
                                .default_value("128")
                                .help("cutoff of probability for read methylation, 0-255"))
                        .arg(
                            Arg::new("DESIGNATE_MAPQ")
                                .long("designate-mapq")
                                .short('q')
                                .value_parser(value_parser!(u8))
                                .default_value("2")
                                .help("designate mapq for refined alignments"))
                        .arg(
                            Arg::new("IS_SET_Y")
                                .long("is-set-y")
                                .short('y')
                                .action(ArgAction::SetTrue)
                                .default_value("false")
                                .help("the aligner was set -Y for supplementary alignments")
                        )
                        .arg(
                            Arg::new("CPG")
                                .long("cpg")
                                .action(ArgAction::SetTrue)
                                .default_value("false")
                                .help("use cpg methylation only, default use all 5mC sites")
                        )
                        .arg(
                            Arg::new("OUTPUT_SECONDARY")
                                .long("output-secondary")
                                .short('s')
                                .action(ArgAction::SetTrue)
                                .default_value("false")
                                .help("output secondary alignments")
                        )
                        .arg(
                            Arg::new("THREADS")
                                .long("threads")
                                .short('t')
                                .value_parser(value_parser!(usize))
                                .default_value("8"))
                        .arg(
                            Arg::new("OUTPUT")
                                .long("output")
                                .short('o')
                                .value_parser(value_parser!(String))
                                .default_value("-")
                                .help("output file, default is stdout"))
                        .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("cool2mcool")
                .display_order(57)
                .about("generate a normalized multi-resolution Cooler file from a fixed 1 kb COOL input")
                .arg(
                    Arg::new("INPUT")
                        .value_name("INPUT.cool")
                        .value_parser(value_parser!(PathBuf))
                        .required(true)
                        .help("input fixed-resolution, symmetric-upper COOL file at 1 kb"),
                )
                .arg(
                    Arg::new("OUTPUT")
                        .value_name("OUTPUT.mcool")
                        .value_parser(value_parser!(PathBuf))
                        .required(true)
                        .help("destination MCOOL file"),
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .value_name("N")
                        .value_parser(positive_usize)
                        .help("Rust aggregation worker count (default: available cores, capped at 8)"),
                )
                .arg(
                    Arg::new("LEVEL_PARALLELISM")
                        .long("level-parallelism")
                        .value_name("N")
                        .value_parser(positive_usize)
                        .default_value("2")
                        .help("concurrent resolution levels (at most --threads)"),
                )
                .arg(
                    Arg::new("AGGREGATION_MODE")
                        .long("aggregation-mode")
                        .value_name("MODE")
                        .value_parser(["pyramid", "direct"])
                        .default_value("pyramid")
                        .help("reuse parent levels (pyramid), or aggregate every level from 1 kb (direct)"),
                )
                .arg(
                    Arg::new("COMPRESSION_LEVEL")
                        .long("compression-level")
                        .value_name("LEVEL")
                        .value_parser(cool2mcool_compression_level)
                        .default_value("1")
                        .help("gzip compression level for generated datasets (1-9)"),
                )
                .arg(
                    Arg::new("RESOLUTIONS")
                        .long("resolutions")
                        .value_name("BP[,BP...]")
                        .value_delimiter(',')
                        .value_parser(cool2mcool_resolution)
                        .default_value(COOL2MCOOL_DEFAULT_RESOLUTIONS)
                        .help("comma-separated output resolutions in base pairs"),
                )
                .arg(
                    Arg::new("KR_MIN_RESOLUTION")
                        .long("kr-min-resolution")
                        .value_name("BP")
                        .value_parser(cool2mcool_resolution)
                        .default_value(COOL2MCOOL_DEFAULT_KR_MIN_RESOLUTION)
                        .help("compute and store KR at resolutions greater than or equal to this value"),
                )
                .arg(
                    Arg::new("FORCE")
                        .long("force")
                        .action(ArgAction::SetTrue)
                        .help("atomically replace an existing output after successful generation"),
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("alleles")
                .alias("allele")
                .hide(true)
                .about("Identify allelic contig pairs with parallel partig-style self comparison")
                .arg(arg!(<FASTA> "fasta"))
                .arg(
                    Arg::new("K")
                        .long("kmer-size")
                        .short('k')
                        .value_parser(value_parser!(usize))
                        .default_value("19"))
                .arg(
                    Arg::new("W")
                        .long("window-size")
                        .short('w')
                        .value_parser(value_parser!(usize))
                        .default_value("19"))
                .arg(
                    Arg::new("M")
                        .long("minimum-similarity")
                        .short('m')
                        .value_parser(value_parser!(f64))
                        .default_value("0.85"))
                .arg(
                    Arg::new("MAX_OCCURRENCE")
                        .long("max-occurrence")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("100")
                        .help("ignore minimizers occurring more than this many times"))
                .arg(
                    Arg::new("MIN_CHAIN")
                        .long("min-chain")
                        .short('n')
                        .value_parser(value_parser!(usize))
                        .default_value("5")
                        .help("minimum number of collinear minimizers"))
                .arg(
                    Arg::new("DIFF_THRESHOLD")
                        .long("diff-threshold")
                        .short('d')
                        .value_parser(value_parser!(f64))
                        .default_value("0.1")
                        .help("retain matches within this fraction of the best chain"))
                .arg(
                    Arg::new("TRIM_LENGTH")
                        .long("trim-length")
                        .value_parser(value_parser!(usize))
                        .default_value("0")
                        .help("trim this many bases from both ends when contig length is greater than three times this value"))
                .arg(
                    Arg::new("SPLIT_REGIONS")
                        .long("split-regions")
                        .value_parser(value_parser!(String))
                        .help("four-column TSV: split name, source contig, 0-based start, end"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("kprune")
                .display_order(70)
                .about("Identify the allelic and cross-allelic contig pairs by allele table")
                .arg(arg!(<ALLELETABLE> "allele table"))
                .arg(arg!(<CONTACTS> "contacts"))
                .arg(arg!(<PRUNETABLE> "output path of prune table"))
                .arg(
                    Arg::new("COUNTRE")
                        .long("count-re")
                        .short('c')
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        .help("restriction enzyme count file")
                )
                .arg(
                    Arg::new("METHOD")
                        .long("method")
                        .short('m')
                        .value_parser(value_parser!(String))
                        .default_value("precise")
                        .help("method of prune: [fast, precise, greedy]")
                        )
                .arg(
                    Arg::new("NORMALIZATION_METHOD")
                        .long("normalization-method")
                        .short('n')
                        .value_parser(value_parser!(String))
                        .default_value("cis")
                        .help("normalization method: cis, cis_unique, none")
                )
                .arg(
                    Arg::new("FIRST_CLUSTER")
                        .long("first-cluster")
                        .short('f')
                        .help("first cluster from hyperpartition")
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        )
                .arg(
                    Arg::new("WHITELIST")
                        .long("whitelist")
                        .short('w')
                        .help("whitelist file, only keep the contig pairs in whitelist")
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        )
                .arg(
                    Arg::new("PARTIAL_WHITELIST")
                        .long("partial-whitelist")
                        .short('p')
                        .help("partial whitelist file, only keep the contig pairs with one contig in partial whitelist")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("prune")
                .hide(true)
                .about("Identity the allelic and cross-allelic contig pairs by raw allele table")
                .arg(arg!(<ALLELETABLE> "raw allele table"))
                .arg(arg!(<ALLELESTRANDTABLE> "allele table"))
                .arg(arg!(<CONTACTS> "contacts"))
                .arg(arg!(<PRUNETABLE> "output path of prune table"))
                .arg(
                    Arg::new("COUNTRE")
                        .long("count-re")
                        .short('c')
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        .help("restriction enzyme count file")
                )
                .arg(
                    Arg::new("METHOD")
                        .long("method")
                        .short('m')
                        .value_parser(value_parser!(String))
                        .default_value("fast")
                        .help("method of prune: [fast, precise, greedy]")
                        )
                .arg(
                    Arg::new("NORMALIZATION_METHOD")
                        .long("normalization-method")
                        .short('n')
                        .value_parser(value_parser!(String))
                        .default_value("cis")
                        .help("normalization method, cis, cis_unique, none")
                ).arg(
                    Arg::new("FIRST_CLUSTER")
                        .long("first-cluster")
                        .short('f')
                        .help("first cluster from hyperpartition")
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        )
                .arg(
                    Arg::new("WHITELIST")
                        .long("whitelist")
                        .short('w')
                        .help("whitelist file, only keep the contigs in whitelist")
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("splitbam")
                .display_order(69)
                .about("split bam by record number")
                .alias("split-bam")
                .arg(arg!(<BAM> "align bam from `dorado`"))
                .arg(
                    Arg::new("RECORD_NUM")
                        .long("record-num")
                        .short('n')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("output.split"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("clm")
                .display_order(71)
                .about("Read, write, and convert CLM/CLMB files.")
                .subcommand(
                    Command::new("convert")
                        .about("Convert between text CLM/CLM.GZ and binary CLMB.")
                        .arg(Arg::new("INPUT").required(true).help("input CLM, CLM.GZ, or CLMB"))
                        .arg(
                            Arg::new("OUTPUT")
                                .long("output")
                                .short('o')
                                .required(true)
                                .help("output .clmb, .clm, or .clm.gz"),
                        )
                        .arg(
                            Arg::new("BLOCK_MIB")
                                .long("block-mib")
                                .value_parser(value_parser!(usize))
                                .default_value("8")
                                .help("target uncompressed CLMB block size in MiB"),
                        )
                        .arg_required_else_help(true),
                )
                .subcommand_required(true)
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("mergeclm")
                .display_order(72)
                .alias("clm-merge")
                .alias("merge-clm")
                .about("merge clm files")
                .arg(
                    Arg::new("INPUTS")
                        .action(ArgAction::Set)
                        .num_args(1..)
                        .required(true)
                        .help("the dir contains clm files or multiple clm files")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout")
                )
                .arg_required_else_help(true),

        )
        .subcommand(
            Command::new("splitclm")
                .display_order(73)
                .alias("split-clm")
                .alias("clm-split")
                .about("Split CLM or CLMB by the cluster file.")
                .arg(arg!(<CLM> "input CLM, CLM.GZ, or CLMB"))
                .arg(arg!(<CLUSTER> "cluster"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .help("output dir, default splitclm")
                        .value_parser(value_parser!(String))
                        .default_value("./"))
                        .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("splitcontacts")
                .display_order(56)
                .alias("split-contacts")
                .about("Split contacts by the cluster file.")
                .arg(arg!(<CONTACTS> "contacts file"))
                .arg(arg!(<CLUSTER> "cluster"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .help("output dir, default splitclm")
                        .value_parser(value_parser!(String))
                        .default_value("./"))
                .arg_required_else_help(true)
                
        )
        .subcommand(
            Command::new("gfa")
                .display_order(74)
                .about("Process GFA-derived inputs for scaffolding.")
                .subcommand(
                    Command::new("aggregate-contacts")
                        .about("Aggregate split contacts into normalized GFA-end evidence.")
                        .arg(Arg::new("SEGMENTS").long("segments").required(true).help("allowed segment lengths TSV"))
                        .arg(Arg::new("CONTACTS").long("contacts").required(true).help("input split contacts"))
                        .arg(Arg::new("GFA_LINKS").long("gfa-links").required(true).help("GFA physical links to summarize"))
                        .arg(Arg::new("OUTPUT").long("output").required(true).help("output normalized end-contact TSV"))
                        .arg(
                            Arg::new("THREADS")
                                .long("threads")
                                .short('t')
                                .value_parser(value_parser!(usize))
                                .default_value("8")
                                .help("parallel contact parsing and aggregation workers")
                        )
                        .arg_required_else_help(true)
                )
                .subcommand(
                    Command::new("contract-inputs")
                        .about("Contract split contacts and CLM records onto selected GFA blocks.")
                        .arg(Arg::new("MAPPING").long("mapping").required(true).help("precomputed GFA block contraction table"))
                        .arg(Arg::new("CONTACTS").long("contacts").required(true).help("input split contacts"))
                        .arg(Arg::new("CONTACTS_OUTPUT").long("contacts-output").required(true).help("output contracted split contacts"))
                        .arg(Arg::new("CLM").long("clm").help("optional input CLM"))
                        .arg(Arg::new("CLM_OUTPUT").long("clm-output").help("output contracted CLM"))
                        .arg(Arg::new("CLUSTERS").long("clusters").help("cluster table used to route contracted CLM records directly by group"))
                        .arg(Arg::new("CLM_OUTPUT_DIR").long("clm-output-dir").help("directory for directly split <group>.clm files"))
                        .arg(Arg::new("TMP_DIR").long("tmp-dir").help("temporary directory"))
                        .arg(
                            Arg::new("SORT_BUFFER")
                                .long("sort-buffer")
                                .value_parser(value_parser!(String))
                                .default_value("64M")
                                .help("maximum memory used by external sort")
                        )
                        .arg(
                            Arg::new("THREADS")
                                .long("threads")
                                .short('t')
                                .value_parser(value_parser!(usize))
                                .default_value("8")
                                .help("parallel sort and CLM block workers")
                        )
                        .arg_required_else_help(true)
                )
                .subcommand_required(true)
                .arg_required_else_help(true)
        )
        .subcommand(
            Command::new("splitfastq")
                .display_order(80)
                .about("split fastq by record number")
                .alias("splitfq")
                .arg(arg!(<FASTQ> "fastq"))
                .arg(
                    Arg::new("RECORD_NUM")
                        .long("record-num")
                        .short('n')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("output.split"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("extract-fasta")
                .display_order(81)
                .about("extract fasta by a cluster file")
                .alias("split-fasta-by-cluster")
                .alias("split-fasta-by-clusters")
                .arg(arg!(<FASTA> "fasta"))
                .arg(arg!(<CLUSTERS> "clusters file"))
                .arg(
                    Arg::new("TRIM_LENGTH")
                        .long("trim-length")
                        .short('l')
                        .value_parser(value_parser!(usize))
                        .default_value("25000"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("slidefastq")
                .display_order(82)
                .about("slide fastq")
                .alias("slidefq")
                .arg(arg!(<FASTQ> "fastq"))
                .arg(
                    Arg::new("WINDOW")
                        .long("window")
                        .short('w')
                        .value_parser(value_parser!(u64))
                        .default_value("5000"))
                .arg(
                    Arg::new("STEP")
                        .long("step")
                        .short('s')
                        .value_parser(value_parser!(u64))
                        .default_value("0"))
                .arg(
                    Arg::new("MIN_LENGTH")
                        .long("min-lenth")
                        .short('l')
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                )
                .arg(
                    Arg::new("FILETYPE")
                        .long("filetype")
                        .short('f')
                        .value_parser(value_parser!(String))
                        .default_value("auto")
                        .help("filtype of sequences, auto dont support for stream input")
                )
                .arg(
                    Arg::new("COORDINATE")
                        .long("coordinate")
                        .short('c')
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("add coordinate suffix of record id")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("slidefasta")
                .display_order(83)
                .about("slide fasta")
                .alias("slidefa")
                .arg(arg!(<FASTA> "fasta"))
                .arg(
                    Arg::new("WINDOW")
                        .long("window")
                        .short('w')
                        .value_parser(value_parser!(u64))
                        .default_value("10000"))
                .arg(
                    Arg::new("STEP")
                        .long("step")
                        .short('s')
                        .value_parser(value_parser!(u64))
                        .default_value("0"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("slide2raw")
                .display_order(84)
                .about("Convert slided read id to raw in bam file")
                .arg(arg!(<BAM> "slided read mapping bam"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("simulator")
                .display_order(91)
                .about("simulating test data")
                .arg_required_else_help(true)
                .allow_external_subcommands(true)
                .subcommand(
                    Command::new("split-ont")
                        .about("simulate short from long read")
                        .arg(arg!(<BAM> "raw bam with MM/ML tags"))
                        .arg(
                            Arg::new("MIN_QUALITY")
                                .long("min-quality")
                                .short('q')
                                .value_parser(value_parser!(u8))
                                .default_value("40"))
                        .arg(
                            Arg::new("OUTPUT")
                                .long("output")
                                .short('o')
                                .value_parser(value_parser!(String))
                                .default_value("-")
                        .help("output file, default is stdout"))
                        .arg_required_else_help(true),
                ).subcommand(
                    Command::new("concat")
                        .about("simulate pore-c reads")
                        .alias("porec")
                        .arg(arg!(<FASTA> "fasta"))
                        .arg(arg!(<VCF> "vcf"))
                        .arg(arg!(<BED> "bed"))
                        .arg(
                            Arg::new("OUTPUT")
                                .long("output")
                                .short('o')
                                .value_parser(value_parser!(String))
                                .default_value("-")
                        .help("output file, default is stdout"))
                        .arg_required_else_help(true),
                        )   
                    .subcommand(
                        Command::new("hic")
                            .about("simulate hic reads")
                            .arg(arg!(<FASTA> "fasta"))
                            .arg(arg!(<VCF> "vcf"))
                            .arg(arg!(<BAM> "bam"))
                            .arg(
                                Arg::new("MIN_MAPQ")
                                    .long("--min-mapq")
                                    .short('q')
                                    .value_parser(value_parser!(u8))
                                    .default_value("1")
                            )
                            .arg(
                                Arg::new("THREADS")
                                    .long("threads")
                                    .short('t')
                                    .value_parser(value_parser!(usize))
                                    .default_value("4")
                            )
                            .arg(
                                Arg::new("OUTPUT")
                                    .long("output")
                                    .short('o')
                                    .value_parser(value_parser!(String))
                                    .default_value("-")
                            .help("output file, default is stdout"))
                            .arg_required_else_help(true),
                            )
        )
        .subcommand(
            Command::new("kmer")
                .display_order(85)
                .hide(true)
                .about("some kmer operations")
                .arg_required_else_help(true)
                .allow_external_subcommands(true)
                .subcommand(
                    Command::new("count")
                        .about("count kmer")
                        .arg(arg!(<FASTA> "fasta"))
                        .arg(
                            Arg::new("K")
                                .long("k")
                                .short('k')
                                .value_parser(value_parser!(usize))
                                .default_value("19"))
                        .arg(
                            Arg::new("OUTPUT")
                                .long("output")
                                .short('o')
                                .value_parser(value_parser!(String))
                                .help("output file, default is stdout")
                                .default_value("-")
                        ).arg_required_else_help(true),
                )
                .subcommand(                
                    Command::new("mask")
                        .about("mask high frequency kmer")
                        .arg(arg!(<FASTA> "fasta"))
                        .arg(
                            Arg::new("K")
                                .long("k")
                                .short('k')
                                .value_parser(value_parser!(usize))
                                .default_value("19"))
                        .arg(
                            Arg::new("PLOIDY")
                                .long("ploidy")
                                .short('p')
                                .value_parser(value_parser!(u64))
                                .default_value("12")
                        )
                        .arg(
                            Arg::new("OUTPUT")
                                .long("output")
                                .short('o')
                                .value_parser(value_parser!(String))
                                .help("output file, default is stdout")
                                .default_value("-")
                        )
                        
                        .arg_required_else_help(true),
                )
                .subcommand(
                Command::new("position")
                    .about("export kmer positions to bed")
                    .arg_required_else_help(true)
                    .arg(arg!(<FASTA> "fasta"))
                    .arg(arg!(<KMER_LIST> "kmer list"))
                    .arg(
                        Arg::new("K")
                            .long("k")
                            .short('k')
                            .value_parser(value_parser!(usize))
                            .default_value("19")
                    )
                    .arg(
                        Arg::new("OUTPUT")
                            .long("output")
                            .short('o')
                            .value_parser(value_parser!(String))
                            .help("output file, default is stdout")
                            .default_value("-")
                    ).arg_required_else_help(true),

                    )
        )
        .subcommand(
            Command::new("digest")
                .display_order(86)
                .about("digest genome by restriction enzyme, output restriction sites")
                .arg(arg!(<FASTA> "fasta"))
                .arg(
                    Arg::new("PATTERN")
                        .long("pattern")
                        .short('p')
                        .value_parser(value_parser!(String))
                        .default_value("GATC")
                        .help("restriction enzyme pattern, multiple pattern use comma to seperate"))
                .arg(
                        Arg::new("SLOPE")
                            .long("slope")
                            .short('s')
                            .value_parser(value_parser!(i64))
                            .default_value("30"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("count_re")
                .display_order(87)
                .alias("countre")
                .about("count restriction enzyme sites, only support single enzyme")
                .arg(arg!(<FASTA> "fasta"))
                .arg(
                    Arg::new("PATTERN")
                        .long("pattern")
                        .short('p')
                        .value_parser(value_parser!(String))
                        .default_value("AAGCTT")
                        .help("restriction enzyme pattern, multiple pattern use comma to seperate"))
                .arg(
                    Arg::new("MIN_RE")
                        .long("min-re")
                        .short('r')
                        .value_parser(value_parser!(u64))
                        .default_value("1"))
                    
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("cutsite")
                .display_order(88)
                .about("cut restriction site on pore-c reads")
                .arg(arg!(<FASTQ> "pore-c reads"))
                .arg(
                    Arg::new("PATTERN")
                        .long("pattern")
                        .short('p')
                        .value_parser(value_parser!(String))
                        .default_value("GATC"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("paf2depth")
                .display_order(12)
                .about("Calculate depth from paf file")
                .arg(arg!(<PAF> "paf"))
                .arg(arg!(<CHROMSIZES> "chromsizes"))
                .arg(
                    Arg::new("WINSIZE")
                        .long("winsize")
                        .short('w')
                        .value_parser(value_parser!(usize))
                        .default_value("10000")
                )
                .arg(
                    Arg::new("STEPSIZE")
                        .long("stepsize")
                        .short('s')
                        .value_parser(value_parser!(usize))
                        .default_value("0")
                )
                .arg(
                    Arg::new("MIN_MAPQ")
                        .long("min-mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0")
                )
                .arg(
                    Arg::new("SECONDARY")
                        .long("secondary")
                        .action(ArgAction::SetTrue)
                        .help("include secondary alignments")
                        .value_parser(value_parser!(bool))
                        .default_value("false")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),

        )
        .subcommand(
            Command::new("concat2depth")
                .display_order(23)
                .about("Calculate depth from pore-c table file")
                .alias("porec2depth")
                .alias("porec-depth")
                .arg(arg!(<TABLE> "pore-c table"))
                .arg(arg!(<CHROMSIZES> "chromsizes"))
                .arg(
                    Arg::new("WINSIZE")
                        .long("winsize")
                        .short('w')
                        .value_parser(value_parser!(usize))
                        .default_value("10000")
                )
                .arg(
                    Arg::new("STEPSIZE")
                        .long("stepsize")
                        .short('s')
                        .value_parser(value_parser!(usize))
                        .default_value("0")
                )
                .arg(
                    Arg::new("MIN_MAPQ")
                        .long("min-mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("bam-chr2ctg")
                .display_order(62)
                .about("Convert a chromosome-level BAM file to contig-level BAM using a BED mapping")
                .arg(
                    Arg::new("INPUT")
                        .long("input")
                        .short('i')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Input chromosome-level BAM file")
                )
                .arg(
                    Arg::new("BED")
                        .long("bed")
                        .short('b')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("4-column BED: chrom start end contig_name (0-based, half-open)")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Output contig-level BAM file")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("Number of threads")
                )
                .arg_required_else_help(true)
        )
        .subcommand(
            Command::new("paf2concat")
                .display_order(11)
                .about("convert PAF to text Pore-C or alignment-level concat PQS")
                .alias("paf2concatemer")
                .alias("paf2porec")
                .alias("paf2pcon")
                .alias("paf2table")
                .arg(arg!(<PAF> "paf"))
                .arg(
                    Arg::new("BED")
                        .long("bed")
                        .short('b')
                        .value_parser(value_parser!(String))
                        .default_value(""))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output path; .concat.pqs/.porec.pqs selects PQS; .concat/.porec selects text; adding .gz compresses text"))
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
                        .default_value("0.8"))
                .arg(
                    Arg::new("MIN_LENGTH")
                        .long("min-length")
                        .short('l')
                        .value_parser(value_parser!(u32))
                        .default_value("150"))
                .arg(
                    Arg::new("MAX_EDGE")
                        .long("min-edge")
                        .short('e')
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("remove the alignments located in the edge of contigs")
                        
                )
                // .arg(
                //     Arg::new("MIN_ORDER")
                //         .long("min-order")
                //         .short('m')
                //         .value_parser(value_parser!(usize))
                //         .default_value("2")
                // )
                .arg(
                    Arg::new("MAX_ORDER")
                        .long("max-order")
                        .short('M')
                        .value_parser(value_parser!(u32))
                        .default_value("50")
                )
                .arg(
                    Arg::new("SECONDARY")
                        .long("secondary")
                        .action(ArgAction::SetTrue)
                        .help("include secondary alignments")
                        .value_parser(value_parser!(bool))
                        .default_value("false")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat2pqs")
                .display_order(21)
                .about("convert an alignment-level Pore-C table to concat.pqs")
                .alias("porec2pqs")
                .alias("concatemer2pqs")
                .alias("con2pqs")
                .arg(arg!(<TABLE> "Pore-C alignment table, optionally gzip-compressed"))
                .arg(arg!(<CHROMSIZES> "two-column contig sizes file"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output concat.pqs directory; derived from TABLE by default")
                )
                .arg(
                    Arg::new("CHUNKSIZE")
                        .long("chunksize")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000")
                        .help("target alignment rows per shard; read_idx groups are never split")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("parallel parser and Parquet writer threads")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat-split")
                .display_order(25)
                .alias("porec-split")
                .alias("split-porec")
                .about("split a Pore-C table into complete-read concat.pqs shards")
                .arg(arg!(<TABLE> "Pore-C table or concat.pqs directory"))
                .arg(
                    Arg::new("CHROMSIZES")
                        .long("chromsizes")
                        .value_parser(value_parser!(String))
                        .help("required for text input; concat.pqs uses its _contigsizes")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("output .concat.pqs/.porec.pqs directory")
                )
                .arg(
                    Arg::new("CHUNKSIZE")
                        .long("chunksize")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000")
                        .help("target rows per shard; never splits read_idx")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("parallel parser and Parquet writer threads")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat2pairs")
                .display_order(22)
                .about("convert concatemer (con) table to pairs")
                .alias("porec2pairs")
                .alias("concatemer2pairs")
                .alias("con2pairs")
                .alias("pore2pairs")
                .alias("porec2pair")
                .arg(arg!(<TABLE> "pore-c table"))
                .arg(arg!(<CHROMSIZES> "chromsizes"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg(
                    Arg::new("MIN_MAPQ")
                        .long("min-mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1")
                        .help("minimum mapping quality of the alignment")
                    )
                .arg(
                    Arg::new("MIN_ORDER")
                        .long("min-order")
                        .value_parser(value_parser!(usize))
                        .default_value("2")
                        .help("min order of the concatemers")
                )
                .arg(
                    Arg::new("MAX_ORDER")
                        .long("max-order")
                        .value_parser(value_parser!(usize))
                        .default_value("50")
                        .help("max order of the concatemers")
                )
                .arg(
                    Arg::new("CHUNKSIZE")
                        .long("chunksize")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000")
                        .help("chunksize of the pqs output, unused for pairs or pairs.gz output")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat-break")
                .display_order(29)
                .about("Break contigs at break points.")
                .alias("porec-break")
                .alias("concatemer-break")
                .alias("con-break")
                .arg(arg!(<TABLE> "pore-c table"))
                .arg(arg!(<BREAK_BED> "break points in bed format"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("number of threads")
                )
                .arg(
                    Arg::new("CHUNKSIZE")
                        .long("chunksize")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000")
                        .help("target alignment rows per shard for concat.pqs output")
                )
                .arg(
                    Arg::new("CHROMSIZES")
                        .long("chromsizes")
                        .value_parser(value_parser!(String))
                        .help("contig sizes required when text input is written as concat.pqs")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),

        )
        .subcommand(
            Command::new("concat-dup")
                .display_order(30)
                .about("Break contigs at break points.")
                .alias("porec-dup")
                .alias("concatemer-dup")
                .alias("con-dup")
                .arg(arg!(<TABLE> "pore-c table"))
                .arg(arg!(<COLLAPSED> "collapsed contigs list, two columns with raw contigs and dup contigs."))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                )
                .arg(
                    Arg::new("CHUNKSIZE")
                        .long("chunksize")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000")
                        .help("target alignment rows per shard for concat.pqs output")
                )
                .arg(
                    Arg::new("CHROMSIZES")
                        .long("chromsizes")
                        .value_parser(value_parser!(String))
                        .help("contig sizes required when text input is written as concat.pqs")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .help("materialize a duplicated output; without this option only input PQS/cn.info is updated"))
                .arg_required_else_help(true),

        )
        .subcommand(
            Command::new("concat-merge")
                .display_order(26)
                .about("merge Pore-C text tables or natively merge concat PQS directories")
                .alias("porec-merge")
                .alias("concatemer-merge")
                .alias("con-merge")
                .arg(
                    Arg::new("TABLES")
                        .action(ArgAction::Set)
                        .num_args(0..)
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("parallel Parquet workers for native PQS output")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat2reads")
                .display_order(24)
                .alias("porec2reads")
                .alias("concatemer2reads")
                .alias("porec2read")
                .about("Simulate PE reads from pore-c table and reference genome")
                .arg(arg!(<TABLE> "pore-c table, or bam/paf file (will convert automatically)"))
                .arg(arg!(<FASTA> "reference fasta"))
                .arg(
                    Arg::new("LENGTH")
                        .long("length")
                        .short('l')
                        .value_parser(value_parser!(usize))
                        .default_value("0")
                        .help("Simulated read length. Use 0 for full length of the fragment")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("simulated_reads")
                        .help("Output file prefix for _R1.fa and _R2.fa")
                )
                .arg(
                    Arg::new("MIN_MAPQ")
                        .long("min-mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1")
                        .help("Minimum mapping quality (Only used for BAM/PAF conversion)")
                )
                .arg(
                    Arg::new("MIN_IDENTITY")
                        .long("min-identity")
                        .short('p')
                        .value_parser(value_parser!(f32))
                        .default_value("0.0")
                        .help("Minimum identity (Only used for BAM/PAF conversion)")
                )
                .arg(
                    Arg::new("MIN_LENGTH")
                        .long("min-length")
                        .value_parser(value_parser!(u32))
                        .default_value("0")
                        .help("Minimum alignment length (Only used for BAM/PAF conversion)")
                )
                .arg(
                    Arg::new("MAX_EDGE")
                        .long("max-edge")
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("remove alignments located in the edge of contigs (Only used for BAM/PAF conversion)")
                )
                .arg(
                    Arg::new("MAX_ORDER")
                        .long("max-order")
                        .value_parser(value_parser!(u32))
                        .default_value("50")
                        .help("Max order of the concatemers (Only used for BAM/PAF conversion)")
                )
                .arg(
                    Arg::new("SECONDARY")
                        .long("secondary")
                        .action(ArgAction::SetTrue)
                        .help("include secondary alignments (Only used for BAM/PAF conversion)")
                        .value_parser(value_parser!(bool))
                        .default_value("false")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4")
                        .help("Number of threads (Only used for BAM/PAF conversion)")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat-intersect")
                .display_order(27)
                .about("According a bed file to intersection a concatemer (con) table.")
                .alias("porec-intersect")
                .alias("concatemer-intersect")
                .alias("con-intersect")
                .arg(arg!(<TABLE> "pore-c table"))
                .arg(arg!(<BED> "3-columns bed file"))
                .arg(
                    Arg::new("INVERT")
                        .long("invert")
                        .short('v')
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("invert the selection")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout")
                    
                )
                .arg_required_else_help(true),
        )

        .subcommand(
            Command::new("paf-downsample")
                .display_order(14)
                .alias("down-paf")
                .alias("downsample-paf")
                .about("Downsample PAF by virtual pairs, reads, or fraction (random by query/read; VPC=order)")
                .arg(arg!(<PAF> "paf"))
                .arg(
                    Arg::new("MODE")
                        .long("mode")
                        .value_parser(["pairs", "reads", "frac", "bases"])
                        .default_value("pairs")
                        .help("Downsample mode: pairs (target virtual interactions), reads (target reads), frac (fraction), bases (target total bases)")
                )
                .arg(
                    Arg::new("PAIRS")
                        .long("pairs")
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("Target virtual pair interactions to keep (only when --mode pairs)")
                )
                .arg(
                    Arg::new("READS")
                        .long("reads")
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("Target reads/queries to keep (only when --mode reads)")
                )
                .arg(
                    Arg::new("FRAC")
                        .long("frac")
                        .value_parser(value_parser!(f64))
                        .default_value("0.0")
                        .help("Fraction to keep in (0,1] (only when --mode frac)")
                )
                .arg(
                    Arg::new("BASES")
                        .long("bases")
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("Target total bases to keep (only when --mode bases)")
                )
                .arg(
                    Arg::new("BY")
                        .long("by")
                        .value_parser(["reads", "pairs"])
                        .default_value("pairs")
                        .help("For --mode frac: compute fraction by reads or by virtual pairs")
                )
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0")
                        .help("Minimum MAPQ to count an alignment piece into order (VPC)")
                )
                .arg(
                    Arg::new("MIN_IDENTITY")
                        .long("min-identity")
                        .short('p')
                        .value_parser(value_parser!(f32))
                        .default_value("0")
                        .help("Minimum identity to count an alignment piece into order (VPC)")
                )
                .arg(
                    Arg::new("MIN_LENGTH")
                        .long("min-length")
                        .short('l')
                        .value_parser(value_parser!(u32))
                        .default_value("0")
                        .help("Minimum aligned length to count an alignment piece into order (VPC)")
                )
                .arg(
                    Arg::new("MIN_ORDER")
                        .long("min-order")
                        .value_parser(value_parser!(usize))
                        .default_value("1")
                        .help("Min order/VPC of reads to be eligible (inclusive)")
                )
                .arg(
                    Arg::new("MAX_ORDER")
                        .long("max-order")
                        .value_parser(value_parser!(usize))
                        .default_value("50")
                        .help("Max order/VPC of reads to be eligible (exclusive)")
                )
                .arg(
                    Arg::new("SEED")
                        .long("seed")
                        .short('s')
                        .value_parser(value_parser!(u64))
                        .default_value("42")
                        .help("Random seed")
                )
                .arg(
                    Arg::new("KEEP_COMMENTS")
                        .long("keep-comments")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Keep comment/header lines starting with #")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("Output PAF, default stdout")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat-downsample")
                .display_order(28)
                .alias("porec-downsample")
                .alias("down-porec")
                .alias("downsample-porec")
                .about("Downsample pore-c table by virtual pairs, reads, or fraction (random by read)")
                .arg(arg!(<TABLE> "pore-c table"))
                .arg(
                    Arg::new("MODE")
                        .long("mode")
                        .value_parser(["pairs", "reads", "frac"])
                        .default_value("pairs")
                        .help("Downsample mode: pairs (target virtual interactions), reads (target reads), frac (fraction)")
                )
                .arg(
                    Arg::new("PAIRS")
                        .long("pairs")
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("Target virtual pair interactions to keep (only when --mode pairs)")
                )
                .arg(
                    Arg::new("READS")
                        .long("reads")
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("Target reads to keep (only when --mode reads)")
                )
                .arg(
                    Arg::new("FRAC")
                        .long("frac")
                        .value_parser(value_parser!(f64))
                        .default_value("0.0")
                        .help("Fraction to keep in (0,1] (only when --mode frac)")
                )
                .arg(
                    Arg::new("BY")
                        .long("by")
                        .value_parser(["reads", "pairs"])
                        .default_value("pairs")
                        .help("For --mode frac: compute fraction by reads or by virtual pairs")
                )
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1")
                        .help("Minimum MAPQ to count an alignment piece into order")
                )
                .arg(
                    Arg::new("MIN_ORDER")
                        .long("min-order")
                        .value_parser(value_parser!(usize))
                        .default_value("2")
                        .help("Min order of the concatemers (inclusive)")
                )
                .arg(
                    Arg::new("MAX_ORDER")
                        .long("max-order")
                        .value_parser(value_parser!(usize))
                        .default_value("50")
                        .help("Max order of the concatemers (exclusive)")
                )
                .arg(
                    Arg::new("SEED")
                        .long("seed")
                        .short('s')
                        .value_parser(value_parser!(u64))
                        .default_value("42")
                        .help("Random seed")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("Parallel workers for native concat PQS input/output")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("Output pore-c table, default stdout")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("paf2pairs")
                .display_order(13)
                .about("convert paf to pairs")
                .arg(arg!(<PAF> "paf"))
                .arg(arg!(<CHROMSIZES> "chromsizes"))
                .arg(
                    Arg::new("BED")
                        .long("bed")
                        .short('b')
                        .value_parser(value_parser!(String))
                        .default_value(""))
                .arg(
                    Arg::new("MIN_MAPQ")
                        .long("min-mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1")
                        .help("minimum mapping quality of the alignment")
                    )
                .arg(
                    Arg::new("MIN_IDENTITY")
                        .long("min-identity")
                        .short('p')
                        .value_parser(value_parser!(f32))
                        .default_value("0.8")
                        .help("minimum identity of the alignment")
                    )
                .arg(
                    Arg::new("MIN_LENGTH")
                        .long("min-length")
                        .short('l')
                        .value_parser(value_parser!(u32))
                        .default_value("150")
                        .help("minimum length of the alignment")
                    )
                .arg(
                    Arg::new("MAX_EDGE")
                        .long("min-edge")
                        .short('e')
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .help("remove the alignments located in the edge of contigs")
                    )
                .arg(
                    Arg::new("MIN_ORDER")
                        .long("min-order")
                        .short('m')
                        .value_parser(value_parser!(usize))
                        .default_value("2")
                        .help("max order of the concatemers")
                    )
                .arg(
                    Arg::new("MAX_ORDER")
                        .long("max-order")
                        .short('M')
                        .value_parser(value_parser!(u32))
                        .default_value("50")
                        .help("max order of the concatemers")
                    )
                .arg(
                    Arg::new("SECONDARY")
                        .long("secondary")
                        .action(ArgAction::SetTrue)
                        .help("include secondary alignments")
                        .value_parser(value_parser!(bool))
                        .default_value("false")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs-downsample")
                .display_order(40)
                .alias("down-pairs")
                .alias("downsample-pairs")
                .alias("downpairs")
                .about("downsample pairs")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0"))
                .arg(
                    Arg::new("NUMBER")
                        .long("number")
                        .short('n')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000"))
            
                .arg(
                    Arg::new("PERCENT")
                        .long("percent")
                        .short('p')
                        .value_parser(value_parser!(f64))
                        .default_value("0.0"))
                .arg(
                    Arg::new("SEED")
                        .long("seed")
                        .short('s')
                        .value_parser(value_parser!(usize))
                        .default_value("42"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs-merge")
                .display_order(41)
                .alias("pairsmerge")
                .alias("merge-pairs")
                .about("merge multiple pairs to one")
                .arg(
                    Arg::new("FILES")
                        .action(ArgAction::Set)
                        .num_args(0..)
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs-filter")
                .display_order(42)
                .alias("filterpairs")
                .alias("filter-pairs")
                .about("filter pairs by mapq")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1"))
                .arg(
                    Arg::new("WHITELIST")
                        .long("whitelist")
                        .short('w')
                        .help("whitelist file, only keep the contigs in whitelist")
                        .value_parser(value_parser!(String))
                        .default_value("none")
                        )
                .arg(Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs-break")
                .display_order(43)
                .about("Break contigs at chimeric points follwed a bed")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(arg!(<BREAK_BED> "break contigs with a bed format"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs-dup")
                .display_order(44)
                .about("dup collapsed contigs by collapsed rescue")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(arg!(<COLLAPSED> "collapsed contigs with raw contig and duplicated contigs"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .help("materialize a duplicated output; without this option only input PQS/cn.info is updated"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs-split")
                .display_order(45)
                .about("split pairs by chunksize")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("CHUNKSIZE")
                        .long("chunksize")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("10000000"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("output.split"))
                .arg_required_else_help(true),

        )
        .subcommand(
            Command::new("pairs2contacts")
                .display_order(46)
                .alias("pairs2contact")
                .about("calculate the contacts between contigs")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(Arg::new("MIN_CONTACTS")
                        .long("min-contacts")
                        .short('c')
                        .value_parser(value_parser!(u32))
                        .default_value("1"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1"))
                .arg(
                    Arg::new("SPLIT_NUM")
                        .long("split-num")
                        .short('n')
                        .value_parser(value_parser!(u32))
                        .default_value("1"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs2clm")
                .display_order(47)
                .about("convert pairs to clm file, which is used for `allhic optimize`")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(Arg::new("MIN_CONTACTS")
                        .long("min-contacts")
                        .short('c')
                        .value_parser(value_parser!(u32))
                        .default_value("1"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0"))
                .arg(
                    Arg::new("NO_OUTPUT_SPLIT_CONTACTS")
                        .long("no-output-split")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Dont output split contacts."))
                .arg(
                    Arg::new("OUTPUT_DEPTH")
                        .short('d')
                        .long("output-depth")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("The path of ooutput depth."))
                .arg(
                    Arg::new("BINSIZE")
                        .long("binsize")
                        .short('b')
                        .value_parser(value_parser!(u32))
                        .default_value("10000")
                )
                .arg(
                    Arg::new("DISABLE_FILTER")
                        .long("disable-filter")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Disable filter bins with extremely high contact depth.")
                )
                .arg(
                    Arg::new("MAX_DEPTH_RATIO")
                        .long("max-depth-ratio")
                        .value_parser(value_parser!(f64))
                        .default_value("5.0")
                        .help("The ratio threshold for high depth filtering (e.g., 3.0 means > 3 * mean_depth).")
                )
                .arg(
                    Arg::new("MAX_Q0_RATIO")
                        .long("max-q0-ratio")
                        .value_parser(value_parser!(f64))
                        .default_value("0.0")
                        .help("The ratio threshold for contacts (MAPQ=0) (e.g., 3.0 means Q0_data > 3 * Q1_data ), only effect for pairs.pqs and mapq=0.")
                )
                .arg(
                    Arg::new("USE_CN")
                        .long("use-cn")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Use optional PAIRS/cn.info copy numbers when generating PQS CLM outputs."),
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4"))
                .arg(Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("output.clmb")
                        .help("output CLMB by default; use .clm or .clm.gz for legacy text"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs2depth")
                .display_order(48)
                .about("convert pairs to depth file")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("BINSIZE")
                        .long("binsize")
                        .short('b')
                        .value_parser(value_parser!(u32))
                        .default_value("10000")
                )
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs2mnd")
                .display_order(49)
                .about("convert pairs to mnd file")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("IGNORE_CN")
                        .long("ignore-cn")
                        .action(ArgAction::SetTrue)
                        .help("ignore optional pairs.pqs/cn.info and keep original contig names"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        ).subcommand(
            Command::new("pairs2bam")
                .display_order(50)
                .about("convert pairs to pseudo bam file")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs2pqs")
                .display_order(51)
                .hide(true)
                .about("convert pairs to pairs.pqs file")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(
                    Arg::new("CHUNKSIZE")
                        .long("chunksize")
                        .short('c')
                        .value_parser(value_parser!(usize))
                        .default_value("1000000")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("output.pqs")
                        .help("output dir, default is output.pqs"))
                .arg_required_else_help(true),
                    
        )
        .subcommand(
            Command::new("pairs2porec")
                .display_order(52)
                .about("convert pairs to pore-c table")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(arg!(<ENZYME> "enzyme bed file"))
                .arg(
                    Arg::new("MAPQ")
                        .long("mapq")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0")
                        .help("minimum mapping quality of the alignment")
                )
                .arg(
                    Arg::new("MIN_EDGE_SUPPORT")
                        .long("min-edge-support")
                        .value_parser(value_parser!(u32))
                        .default_value("1")
                        .help("Minimum edge support for graph construction")
                )
                .arg(
                    Arg::new("MIN_CLIQUE_SIZE")
                        .long("min-clique-size")
                        .value_parser(value_parser!(usize))
                        .default_value("3")
                        .help("Minimum clique size for outputting components")
                )
                .arg(
                    Arg::new("MAX_COMP_SIZE")
                        .long("max-comp-size")
                        .value_parser(value_parser!(usize))
                        .default_value("1000")
                        .help("Maximum component size. Larger components will be skipped")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4")
                        .help("Number of threads")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is output.porec"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concat-chr2ctg")
                .display_order(31)
                .about("Convert a chromosome-level Pore-C table to contig-level Pore-C table using a BED mapping")
                .alias("porec-chr2ctg")
                .arg(
                    Arg::new("INPUT")
                        .long("input")
                        .short('i')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Input chromosome-level Pore-C table")
                )
                .arg(
                    Arg::new("BED")
                        .long("bed")
                        .short('b')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("4-column BED: chrom start end contig_name (0-based, half-open)")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Output contig-level Pore-C table")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("Number of threads")
                )
                .arg_required_else_help(true)
        )
        .subcommand(
            Command::new("pairs-chr2ctg")
                .display_order(53)
            .about("Convert a chromosome-level pairs PQS directory to contig-level pairs PQS; text .pairs/.pairs.gz inputs are not supported")
            .arg(
                Arg::new("INPUT")
                    .long("input")
                    .short('i')
                    .value_parser(value_parser!(String))
                    .required(true)
                    .help("Input chromosome-level pairs PQS directory; text .pairs and .pairs.gz files are not supported")
            )
            .arg(
                Arg::new("BED")
                    .long("bed")
                    .short('b')
                    .value_parser(value_parser!(String))
                    .required(true)
                    .help("4-column BED: chrom start end contig_name (0-based, half-open). Used to map chr positions to contigs")
            )
            .arg(
                Arg::new("OUTPUT")
                    .long("output")
                    .short('o')
                    .value_parser(value_parser!(String))
                    .required(true)
                    .help("Output contig-level pairs PQS directory")
            )
            .arg(
                Arg::new("THREADS")
                    .long("threads")
                    .short('t')
                    .value_parser(value_parser!(usize))
                    .default_value("8")
                    .help("Number of threads")
            )
            .arg(
                Arg::new("DROP_UNMAPPED")
                    .long("drop-unmapped")
                    .action(ArgAction::SetTrue)
                    .default_value("false")
                    .help("Drop contacts whose positions do not fall into any BED contig interval (default: true behavior). This flag is reserved for future; currently unmapped are dropped")
            )
            .arg_required_else_help(true) 

        )
        .subcommand(
            Command::new("pairs-prune")
                .display_order(54)
                .alias("pqs-prune")
                .about("Prune contacts in pairs.pqs format based on a prune table")
                .arg(
                    Arg::new("INPUT")
                        .long("input")
                        .short('i')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Input pairs.pqs directory (contains q0/q1 parquet)")
                )
                .arg(
                    Arg::new("PRUNETABLE")
                        .long("prune-table")
                        .short('p')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Prune table file (contains contig pairs to prune)")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Output pruned PQS directory")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("Number of threads to use")
                )
                .arg_required_else_help(true)
        )
        .subcommand(
            Command::new("bam-prune")
                .display_order(63)
                .about("Prune contacts in BAM format based on a prune table")
                .arg(
                    Arg::new("BAM")
                        .long("bam")
                        .short('b')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Input BAM file")
                )
                .arg(
                    Arg::new("PRUNETABLE")
                        .long("prune-table")
                        .short('p')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Prune table file (contains contig pairs to prune)")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .required(true)
                        .help("Output pruned BAM file")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8")
                        .help("Number of threads to use to process BAM file")
                )
                .arg_required_else_help(true)
        )
        .subcommand(
            Command::new("bam2paf")
                .display_order(60)
            .about("convert dorado align bam to paf")
            .arg(arg!(<BAM> "bam file should be sorted by read name, dont sort it by coordinate"))
            .arg(
                Arg::new("SECONDARY")
                    .long("secondary")
                    .short('s')
                    .value_parser(value_parser!(bool))
                    .action(ArgAction::SetTrue)
                    .default_value("false")
                    .help("Is output secondary, default is false")
            )
            .arg(
                Arg::new("THREADS")
                    .long("threads")
                    .short('t')
                    .value_parser(value_parser!(usize))
                    .default_value("8"))
            .arg(
                Arg::new("OUTPUT")
                    .long("output")
                    .short('o')
                    .value_parser(value_parser!(String))
                    .default_value("-")
                    .help("output file, default is stdout"))
            .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("bam2fastq")
                .display_order(64)
                .about("convert bam to fastq file")
                .arg(
                    Arg::new("BAM")
                        .action(ArgAction::Set)
                        .num_args(1..)
                        .required(true)
                        .help("one or more bam files sorted by query name")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )

        .subcommand(
            Command::new("bam2fasta")
                .display_order(65)
                .about("convert bam to fasta file")
                .arg(arg!(<BAM> "bam"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )

        .subcommand(
            Command::new("bamstat")
                .display_order(67)
                .about("stat the sequences in bam file")
                .arg(
                    Arg::new("BAM")
                        .action(ArgAction::Set)
                        .num_args(0..)
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )

        .subcommand(
            Command::new("hicbamstat")
                .about("stat the hic bam")
                .hide(true)
                .arg(
                    Arg::new("BAM")
                        .action(ArgAction::Set)
                        .num_args(0..)
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("concatbamstat")
                .display_order(32)
                .about("stat the porec bam")
                .alias("porecbamstat")
                .hide(true)
                .arg(
                    Arg::new("BAM")
                        .action(ArgAction::Set)
                        .num_args(0..)
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("bam2pairs")
                .display_order(61)
                .about("convert read align bam to pairs, bam should be sorted by read name, dont sort it by coordinate")
                .arg(arg!(<BAM> "bam"))
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0"))
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("pairs-intersect")
                .display_order(55)
                .about("According a bed file to intersection a pairs file.")
                .arg(arg!(<PAIRS> "pairs"))
                .arg(arg!(<BED> "3-columns bed file"))
                .arg(
                    Arg::new("INVERT")
                        .long("invert")
                        .short('v')
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("invert the selection")
                )
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("1"))
                .arg(
                    Arg::new("EDGE_LENGTH")
                        .long("edge-length")
                        .short('e')
                        .value_parser(value_parser!(u64))
                        .default_value("0")
                        .hide(true)
                        .help("remove the alignments located in the edge of contigs")
                )
                .arg(
                    Arg::new("MAX_Q0_RATIO")
                        .long("max-q0-ratio")
                        .value_parser(value_parser!(f64))
                        .default_value("0.0")
                        .help("The ratio threshold for contacts (MAPQ=0) (e.g., 3.0 means Q0_data > 3 * Q1_data ), only effect for pairs.pqs and mapq=0.")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("8"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout")
                    
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("chromsizes")
                .display_order(89)
                .about("generate chromsizes file")
                .alias("contigsizes")
                .arg(arg!(<FASTA> "fasta"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),

        )
        .subcommand(
            Command::new("prunepairs")
                .about("prune pairs")
                .hide(true)
                .arg(arg!(<PAIRS> "pairs"))
                .arg(arg!(<PRUNE> "Prune contigs"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
                
        )
        .subcommand(
            Command::new("phase-reads")
                .display_order(76)
                .alias("phasereads")
                .about("Phase reads based on contig grouping information")
                .arg(arg!(<INPUT> "alignment file (paf or name-sorted bam)"))
                .arg(arg!(<GROUPS> "contigs group file, two columns: contig_name and group_name"))
                .arg(
                    Arg::new("FORMAT")
                        .long("format")
                        .short('f')
                        .value_parser(["auto", "paf", "bam"])
                        .default_value("auto")
                        .help("Format of the input alignment file")
                )
                .arg(
                    Arg::new("MIN_QUALITY")
                        .long("min-quality")
                        .short('q')
                        .value_parser(value_parser!(u8))
                        .default_value("0")
                        .help("Minimum mapping quality to retain an alignment")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("4")
                        .help("Number of threads (mainly for BAM parsing)")
                )
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("Output file, default is stdout")
                )
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("modbam2fq")
                .display_order(66)
                .about("convert modified bam to fastq with modified base")
                .arg(arg!(<BAM> "modified bam"))
                .arg(
                    Arg::new("MIN_PROB")
                        .long("min-prob")
                        .short('p')
                        .value_parser(value_parser!(f32))
                        .default_value("0.75"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("modfa")
                .display_order(90)
                .about("modify fasta by bedMethy file")
                .arg(arg!(<FASTA> "fasta"))
                .arg(arg!(<BED> "bed"))
                .arg(
                    Arg::new("MIN_FRAC")
                        .long("min-frac")
                        .short('f')
                        .value_parser(value_parser!(f64))
                        .default_value("0.5"))
                .arg(
                    Arg::new("OUTPUT")
                        .long("output")
                        .short('o')
                        .value_parser(value_parser!(String))
                        .default_value("-")
                        .help("output file, default is stdout"))
                .arg_required_else_help(true),
        )
        .subcommand(
            Command::new("optimize")
                .display_order(75)
                .about("optimize contigs order and orientation.")
                .arg(arg!(<COUNTRE> "count RE file of single cluster"))
                .arg(arg!(<CLMB> "CLMB contact-distance file for ALLHiC ordering and orientation"))
                .arg(
                    Arg::new("MUTATION")
                        .long("mutapb")
                        .alias("mutabp")
                        .short('m')
                        .value_parser(value_parser!(f64))
                        .default_value("0.2")
                        .help("mutation probability in GA")
                )
                .arg(
                    Arg::new("NGEN")
                        .long("ngen")
                        .short('n')
                        .value_parser(value_parser!(usize))
                        .default_value("5000")
                        .help("Maximum number of generations in GA")
                )
                .arg(
                    Arg::new("NPOP")
                        .long("npop")
                        .short('p')
                        .value_parser(value_parser!(usize))
                        .default_value("100")
                        .help("Population size in GA")
                )
                .arg(
                    Arg::new("RESUME")
                        .long("resume")
                        .short('r')
                        .action(ArgAction::SetTrue)
                        .help("resume from previous run")
                        .value_parser(value_parser!(bool))
                        .default_value("false")
                )
                .arg(
                    Arg::new("INITIALIZER")
                        .long("initializer")
                        .value_parser(["random", "seriation", "end-tsp", "end-greedy", "end-hierarchical", "end-beam"])
                        .default_value_if("SPLIT_CONTACTS", ArgPredicate::IsPresent, "end-hierarchical")
                        .default_value("random")
                        .help("de-novo ordering initializer; defaults to end-hierarchical when --split-contacts is provided")
                )
                .arg(
                    Arg::new("SPLIT_CONTACTS")
                        .long("split-contacts")
                        .value_parser(value_parser!(String))
                        .required_if_eq_any([("INITIALIZER", "end-tsp"), ("INITIALIZER", "end-greedy"), ("INITIALIZER", "end-hierarchical"), ("INITIALIZER", "end-beam")])
                        .help("half-contig contacts used by endpoint initializers")
                )
                .arg(
                    Arg::new("SEED")
                        .long("seed")
                        .short('s')
                        .value_parser(value_parser!(u64))
                        .default_value("42")
                        .help("random seed of GA")

                )
                .arg(
                    Arg::new("SKIPGA")
                        .long("skipGA")
                        .alias("skipga")
                        .action(ArgAction::SetTrue)
                        .help("skip genetic algorithm and")
                        .value_parser(value_parser!(bool))
                        .default_value("false")
                )
                .arg(
                    Arg::new("LOGDIST")
                        .long("logDist")
                        .alias("log-dist")
                        .action(ArgAction::SetTrue)
                        .help("use ALLHiC's links * log(distance) ordering objective")
                        .value_parser(value_parser!(bool))
                        .default_value("false")
                )
                .arg(
                    Arg::new("LENGTH_TIERED")
                        .long("length-tiered")
                        .alias("length-tiered-objective")
                        .action(ArgAction::SetTrue)
                        .conflicts_with("LOGDIST")
                        .help("use the experimental 70/20/10 anchor/interval/fragment ordering objective")
                )
                .arg(
                    Arg::new("ENDPOINT_MULTISCALE")
                        .long("endpoint-multiscale")
                        .action(ArgAction::SetTrue)
                        .requires("SPLIT_CONTACTS")
                        .conflicts_with_all(["LOGDIST", "LENGTH_TIERED"])
                        .help("use the experimental signed endpoint multi-scale anchor objective")
                )
                .arg(
                    Arg::new("NO_BACKBONE")
                        .long("no-backbone")
                        .action(ArgAction::SetTrue)
                        .help("disable high-confidence path-block initialization")
                )
                .arg(
                    Arg::new("ORIENTATION_METHOD")
                        .long("orientation-method")
                        .value_parser(["banded-legacy", "robust", "banded", "banded-contact", "intervening", "legacy"])
                        .default_value("banded-legacy")
                        .help("orientation solver: historical banded DP and signed-block refinement (banded-legacy, default), conservative banded DP (banded), experimental endpoint evidence (robust), experimental contact-only margin (banded-contact), intervening-gap refinement, or legacy ALLHiC behavior")
                )
                .arg(
                    Arg::new("ORIENTATION_WINDOW")
                        .long("orientation-window")
                        .value_name("CONTIGS")
                        .value_parser(orientation_window)
                        .default_value("3")
                        .help("local rank window for orientation evidence and signed block refinement (1-16)")
                )
                .arg(
                    Arg::new("ORIENTATION_PAIR_WEIGHT")
                        .long("orientation-pair-weight")
                        .value_parser(["links", "sqrt-links", "equal-pair"])
                        .default_value("sqrt-links")
                        .help("normalization applied to each contig-pair orientation score")
                )
                .arg(
                    Arg::new("ORIENTATION_MIN_LINKS")
                        .long("orientation-min-links")
                        .value_name("N")
                        .value_parser(positive_usize)
                        .default_value("3")
                        .help("minimum retained links required for a contig pair to orient")
                )
                .arg(
                    Arg::new("ORIENTATION_TRUST_INPUT")
                        .long("orientation-trust-input")
                        .action(ArgAction::SetTrue)
                        .requires("RESUME")
                        .help("apply an input-sign prior to an explicitly trusted resumed tour in robust mode")
                )
                .arg(
                    Arg::new("ORIENTATION_AUDIT")
                        .long("orientation-audit")
                        .value_name("TSV")
                        .help("write per-contig endpoint evidence diagnostics in robust mode (not calibrated probabilities)")
                )
                .arg(
                    Arg::new("ORIENTATION_PRIOR")
                        .long("orientation-prior")
                        .value_name("STRENGTH")
                        .value_parser(non_negative_f64)
                        .default_value("0.05")
                        .help("dimensionless input-sign penalty in banded modes; in robust mode applied only with --orientation-trust-input")
                )
                .arg(
                    Arg::new("ORIENTATION_MIN_CONFIDENCE")
                        .long("orientation-min-confidence")
                        .value_name("FRACTION")
                        .hide_default_value(true)
                        .value_parser(unit_interval_f64)
                        .default_value_if("ORIENTATION_METHOD", "banded-legacy", "0")
                        .default_value_if("ORIENTATION_METHOD", "robust", "0.2")
                        .default_value("0.95")
                        .help("minimum normalized max-marginal effect (banded: 0.95) or endpoint effect (robust: 0.2); banded-legacy requires 0 (disabled); not a calibrated probability")
                )
                .arg(
                    Arg::new("ORIENTATION_MAX_FLIP_BP_FRACTION")
                        .long("orientation-max-flip-bp-fraction")
                        .value_name("FRACTION")
                        .hide_default_value(true)
                        .value_parser(positive_unit_interval_f64)
                        .default_value_if("ORIENTATION_METHOD", "banded-legacy", "1")
                        .default_value("0.05")
                        .help("largest consecutive changed-sign run as a fraction of scaffold bp; banded-legacy requires 1 (unrestricted); other modes default to 0.05")
                )
                .arg(
                    Arg::new("ORIENTATION_BLOCK_SPAN")
                        .long("orientation-block-span")
                        .value_name("CONTIGS")
                        .hide_default_value(true)
                        .value_parser(disabled_or_block_span)
                        .default_value_if("ORIENTATION_METHOD", "banded-legacy", "32")
                        .default_value("0")
                        .help("largest reverse-complement block; banded-legacy defaults to 32, other modes to 0 (disabled); robust tests blocks before single signs, banded modes afterwards")
                )
                .arg(
                    Arg::new("ORIENTATION_BLOCK_MAX_BP_FRACTION")
                        .long("orientation-block-max-bp-fraction")
                        .value_name("FRACTION")
                        .hide_default_value(true)
                        .value_parser(positive_unit_interval_f64)
                        .default_value_if("ORIENTATION_METHOD", "banded-legacy", "1")
                        .default_value("0.05")
                        .help("largest reverse-complement block as a fraction of scaffold bp; banded-legacy requires 1 (unrestricted); conservative modes default to 0.05")
                )
                .arg(
                    Arg::new("ORIENTATION_BLOCK_PASSES")
                        .long("orientation-block-passes")
                        .value_name("N")
                        .value_parser(positive_usize)
                        .default_value("4")
                        .help("maximum signed block-refinement sweeps")
                )
                .arg(
                    Arg::new("ORIENTATION_BLOCK_MIN_GAIN")
                        .long("orientation-block-min-gain")
                        .value_name("FRACTION")
                        .hide_default_value(true)
                        .value_parser(non_negative_f64)
                        .default_value_if("ORIENTATION_METHOD", "banded-legacy", "0.0001")
                        .default_value("0.05")
                        .help("minimum relative local-band gain for banded-legacy (0.0001), or independent improvement at both boundaries for conservative banded modes (0.05)")
                )
                .arg(
                    Arg::new("THREADS")
                        .long("threads")
                        .short('t')
                        .value_parser(value_parser!(usize))
                        .default_value("10")
                        .help("number of threads"))
                .arg_required_else_help(true),
        );

    let grouped_help = grouped_subcommands_help(&command);
    command
        .help_template(ROOT_HELP_TEMPLATE)
        .after_help(grouped_help)
}
