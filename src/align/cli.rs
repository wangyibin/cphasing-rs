
use clap::{arg, Arg, ArgAction,
    builder::{
        styling::{AnsiColor, Effects},
        Styles,
    },
    Command,
    value_parser};

const VERSION: &str = env!("CARGO_PKG_VERSION");

const STYLES: Styles = Styles::styled()
.header(AnsiColor::Green.on_default().effects(Effects::BOLD))
.usage(AnsiColor::Green.on_default().effects(Effects::BOLD))
.literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
.placeholder(AnsiColor::Yellow.on_default());

fn parse_si_size_u64(s: &str) -> Result<u64, String> {
let s = s.trim();
if s.is_empty() {
return Err("Empty input".to_string());
}

let last = s.chars().last().unwrap();
let (num_str, multiplier) = match last {
'G' | 'g' => (&s[..s.len() - 1], 1_000_000_000.0),
'M' | 'm' => (&s[..s.len() - 1], 1_000_000.0),
'K' | 'k' => (&s[..s.len() - 1], 1_000.0),
_ => (s, 1.0),
};

let num = num_str.parse::<f64>()
.map_err(|_| format!("Invalid number format: {}", s))?;

Ok((num * multiplier) as u64)
}

fn parse_si_size_i64(s: &str) -> Result<i64, String> {
let s = s.trim();
if s.is_empty() {
return Err("Empty input".to_string());
}

let last = s.chars().last().unwrap();
let (num_str, multiplier) = match last {
'G' | 'g' => (&s[..s.len() - 1], 1_000_000_000.0),
'M' | 'm' => (&s[..s.len() - 1], 1_000_000.0),
'K' | 'k' => (&s[..s.len() - 1], 1_000.0),
_ => (s, 1.0),
};

let num = num_str.parse::<f64>()
.map_err(|_| format!("Invalid number format: {}", s))?;

Ok((num * multiplier) as i64)
}

fn parse_si_size_i32(s: &str) -> Result<i32, String> {
let s = s.trim();
if s.is_empty() {
return Err("Empty input".to_string());
}

let last = s.chars().last().unwrap();
let (num_str, multiplier) = match last {
'G' | 'g' => (&s[..s.len() - 1], 1_000_000_000.0),
'M' | 'm' => (&s[..s.len() - 1], 1_000_000.0),
'K' | 'k' => (&s[..s.len() - 1], 1_000.0),
_ => (s, 1.0),
};

let num = num_str.parse::<f64>()
.map_err(|_| format!("Invalid number format: {}", s))?;

Ok((num * multiplier) as i32)
}

fn parse_float_pair(s: &str) -> Result<(f32, Option<f32>), String> {
let s = s.trim();
if s.is_empty() {
return Err("Empty input".to_string());
}

if let Some((first, second)) = s.split_once(',') {
let v1 = first.trim().parse::<f32>()
    .map_err(|_| format!("Invalid float format: {}", first))?;
let v2 = second.trim().parse::<f32>()
    .map_err(|_| format!("Invalid float format: {}", second))?;
Ok((v1, Some(v2)))
} else {
let v1 = s.parse::<f32>()
    .map_err(|_| format!("Invalid float format: {}", s))?;
Ok((v1, None))
}
}

fn parse_int_pair(s: &str) -> Result<(i32, Option<i32>), String> {
let s = s.trim();
if s.is_empty() {
return Err("Empty input".to_string());
}

if let Some((first, second)) = s.split_once(',') {
let v1 = first.trim().parse::<i32>()
    .map_err(|_| format!("Invalid integer format: {}", first))?;
let v2 = second.trim().parse::<i32>()
    .map_err(|_| format!("Invalid integer format: {}", second))?;
Ok((v1, Some(v2)))
} else {
let v1 = s.parse::<i32>()
    .map_err(|_| format!("Invalid integer format: {}", s))?;
Ok((v1, None))
}
}


pub fn cli() -> Command {
Command::new("align")
.version(VERSION)
.styles(STYLES)
.about("Align Pore-C/CiFi reads with in-memory methylation refinement and context rescue")
.arg(
    arg!(<reference> "reference FASTA" )
)
.arg(
    Arg::new("query")
        .help("query sequences with BAM or fastq(.gz)/fasta(.gz), multiple files allowed. Use '-' for stdin. Stdin and pipe are only supported for BAM input and single file.")
        .required(true)
        .num_args(1..)
)
.next_help_heading("Indexing")
.arg(
    Arg::new("hpc")
        .short('H')
        .hide(true)
        .help("Use homopolymer-compressed k-mer (always true for map-pb/ont)")
        .action(ArgAction::SetTrue),
)
.arg(
    Arg::new("kmer")
        .short('k')
        .help("k-mer size (no larger than 28) [15]")
        .value_parser(value_parser!(i16))
        .value_name("INT")
)
.arg(
    Arg::new("window")
        .short('w')
        .help("minimizer window size [5]")
        .value_parser(value_parser!(i16))
        .value_name("INT")
)
.arg(
    Arg::new("batch_size")
    .short('I')
    .help("split index for every ~NUM input bases [16G]")
    .value_parser(parse_si_size_u64)
    .value_name("NUM")
    .default_value("16G"),

)
.next_help_heading("Mapping")
.arg(
    Arg::new("mid_occ_frac")
        .short('f')
        .help("If fraction, ignore top FLOAT fraction of most frequent minimizers [0.0002]. If integer, ignore minimizers occuring more than INT1 times. INT2 is only effective in the --sr or -xsr mode, which sets the threshold for a second round of seeding. [0.00002]")
        .value_parser(parse_float_pair)
        .value_name("FLOAT|INT1[,INT2]")

)
.arg(
    Arg::new("bounds_of_occurrence")
        .short('U')
        .help("Lower and upper bounds of k-mer occurrences [10,1000000]. The final k-mer occurrence threshold is max{INT1, min{INT2, -f}}. This option prevents excessively small or large -f estimated from the input reference.")
        .value_parser(parse_int_pair)
        .value_name("INT1,[INT2]")
)
.arg(
    Arg::new("max_gap")
        .short('g')
        .help("stop chain enlongation if there are no minimizers in INT-bp [5000]")
        .value_parser(parse_si_size_i32)
        .value_name("INT")
)
.arg(
    Arg::new("max_gap_ref")
        .short('G')
        .help("max intron length (effective with -xsplice; changing -r) [200k]")
        .value_parser(parse_si_size_i32)
        .value_name("INT")
)
.arg(
    Arg::new("max_frag_len")
        .short('F')
        .help("max fragment length (effective with -xsr or in the fragment mode) [800]")
        .value_parser(parse_si_size_i32)
        .value_name("INT")
)
.arg(
    Arg::new("mask_level")
        .short('M')
        .hide(true)
        .help("Mark as secondary a chain that overlaps with a better chain by FLOAT or more of the shorter chain [0.5]")
        .value_parser(value_parser!(f32))
        .value_name("FLOAT")
)
.arg(
    Arg::new("bw")
        .short('r')
        .help("chaining/alignment bandwidth and long-join bandwidth [500,20000]")
        .value_parser(parse_int_pair)
        .value_name("INT,[INT]"),
)
.arg(
    Arg::new("min_cnt")
        .short('n')
        .help("minimal number of minimizers on a chain [3]")
        .value_parser(value_parser!(i32))
        .value_name("INT")
)
.arg(
    Arg::new("min_chain_score")
        .short('m')
        .help("minimal chaining score (matching bases minus log gap penalty) [40]")
        .value_parser(value_parser!(i32))
        .value_name("INT")
)
.arg(
    Arg::new("pri_ratio")
        .short('p')
        .help("Min secondary-to-primary score ratio [0.8]")
        .value_parser(value_parser!(f32))
        .value_name("FLOAT"),
)
.arg(
    Arg::new("best_n")
        .short('N')
        .help("Retain at most N secondary alignments [5]")
        .value_parser(value_parser!(i32))
        .value_name("INT"),
)

.next_help_heading("Alignments")
.arg(
    Arg::new("matching_score")
        .short('A')
        .help("matching score [2]")
        .value_parser(value_parser!(i32))
        .value_name("INT"),
)
.arg(
    Arg::new("mismatch_penalty")
        .short('B')
        .help("mismatch penalty (larger value for lower divergence) [4]")
        .value_parser(value_parser!(i32))
        .value_name("INT"),
)
.arg(
    Arg::new("gap_open")
        .short('O')
        .help("gap open penalty [4,24]")
        .value_parser(parse_int_pair)
        .value_name("INT,[INT]"),
)
.arg(
    Arg::new("gap_extension")
        .short('E')
        .help("gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]")
        .value_parser(parse_int_pair)
        .value_name("INT,[INT]"),
)
.arg(
    Arg::new("z_drop")
        .short('z')
        .help("Z-drop score and inversion Z-drop score [400,200]")
        .value_parser(parse_int_pair)
        .value_name("INT,[INT]"),
)
.arg(
    Arg::new("min_dp_max")
        .short('s')
        .help("minimal peak DP alignment score [80]")
        .value_parser(value_parser!(i32))
        .value_name("INT")
)

.next_help_heading("Input/Output")
.arg(
    Arg::new("output_cigar")
        .short('c')
        .long("cigar")
        .help("output CIGAR in PAF")
        .action(ArgAction::SetTrue)
)
.arg(
    Arg::new("eqx")
        .long("eqx")
        .help("write CIGAR with =/X operators")
        .action(ArgAction::SetTrue)
)
.arg(
    Arg::new("cs")
        .long("cs")
        .num_args(0..=1)
        .require_equals(true)
        .default_missing_value("short")
        .value_parser(["short", "long", "none"])
        .help("output the cs tag: --cs, --cs=short, --cs=long, or --cs=none")
        .value_name("STR")
)
.arg(
    Arg::new("output_bam")
        .short('a')
        .help("output in the BAM format (PAF by default)")
        .action(ArgAction::SetTrue)
)
.arg(
    Arg::new("output")
        .short('o')
        .help("output file path, PAF or BAM formats [stdout]")
        .value_parser(value_parser!(String))
        .default_value("-")
        .value_name("FILE"),
)
.arg(
    Arg::new("soft_clip")
        .short('Y')
        .help("use soft clipping for supplementary alignments")
        .action(ArgAction::SetTrue),

)
.arg(
    Arg::new("seed")
        .long("seed")
        .help("Integer seed for randomizing equally best hits. Minimap2 hashes INT and read name when choosing between equally best hits. [11]" )
        .value_parser(value_parser!(i32))
        .value_name("INT")
        .default_value("11")
)
.arg(
    Arg::new("max_qlen")
        .long("max-qlen")
        .help("skip reads longer than INT [0 (disabled)]")
        .value_parser(value_parser!(i32))
        .value_name("INT")
)
.arg(
    Arg::new("secondary")
        .long("secondary")
        .help("Whether to output secondary alignments" )
        .value_parser(["yes", "no"])
        .default_value("no")
)
.arg(
    Arg::new("threads")
        .short('t')
        .help("number of threads")
        .value_parser(value_parser!(usize))
        .value_name("INT")
        .default_value("8"),
)
.arg(
    Arg::new("mini_batch_size")
        .short('K')
        .help("minibatch size for mapping")
        .value_parser(parse_si_size_i64)
        .value_name("STR")
        .default_value("500M"),

)
.next_help_heading("Gap Rescue for Pore-C")
.arg(
    Arg::new("porec_gap_rescue")
        .long("gap-rescue")
        .alias("porec-gap-rescue")
        .help("Enable a second sensitive k/w pass to rescue uncovered query gaps in Pore-C reads")
        .action(ArgAction::SetTrue),
)
.arg(
    Arg::new("rescue_k")
        .long("rescue-k")
        .help("k-mer size for Pore-C gap rescue sensitive pass [11]")
        .value_parser(value_parser!(i16))
        .value_name("INT")
        .default_value("11"),
)
.arg(
    Arg::new("rescue_w")
        .long("rescue-w")
        .help("minimizer window size for Pore-C gap rescue sensitive pass [5]")
        .value_parser(value_parser!(i16))
        .value_name("INT")
        .default_value("5"),
)
.arg(
    Arg::new("min_gap_len")
        .long("min-gap-len")
        .help("minimum uncovered query gap length to trigger sensitive rescue [100]")
        .value_parser(value_parser!(usize))
        .value_name("INT")
        .default_value("100"),
)
.next_help_heading("Alignment Rescue for Multi-mapping Reads")
.arg(
    Arg::new("no_realign")
        .long("no-realign")
        .conflicts_with("rescue_mode")
        .help("Disable local alignment MAPQ rescue for multi-mapping reads")
        .action(ArgAction::SetTrue),
)
.arg(
    Arg::new("rescue_mode")
        .long("realign")
        .alias("rescue-mode")
        .value_name("STR")
        .value_parser([
            "precise",
            "sensitive",
            // "homeolog",
            // "homoeolog",
            // "homeolog2",
            // "homoeolog2",
            // "zero-homeolog",
            // "zero-homoeolog",
            // "zero",
        ])
        .default_value("precise")
        .help("Rescue mode:\n - precise: use robust rescue algorithm (rescue_robust)\n - sensitive: use sensitive rescue algorithm (rescue)")
        // "\n - homeolog: resolve near-tie A/B/C/D candidates using read-level anchors\n - homeolog2: resolve near-tie A/B/C/D candidates using only other-interval support\n - zero-homeolog: conservatively rescue all-zero-MAPQ homeolog reads")
)
.arg(
    Arg::new("mapq_rescue")
        .long("mapq-rescue")
        .help("MAPQ threshold for realignment rescue")
        .value_parser(value_parser!(u8))
        .value_name("INT")
        .default_value("1"),
)
.arg(
    Arg::new("porec_mapq_calibrate")
        .long("porec-mapq-calibrate")
        .hide(true)
        .help("Recalibrate MAPQ for Pore-C segments using query-interval competition and restriction-site support")
        .action(ArgAction::SetTrue),
)
.arg(
    Arg::new("candidate_probability")
        .hide(true)
        .long("candidate-probability")
        .help("Annotate Pore-C/CiFi candidate alignments with posterior-like probability tags cp/ac/am")
        .action(ArgAction::SetTrue),
)
.arg(
    Arg::new("graph_assignment")
        .long("graph-assignment")
        .hide(true)
        .help("Annotate Pore-C/CiFi multi-fragment reads with graph-based assignment tags ga/gp/gs/gc")
        .action(ArgAction::SetTrue),
)
.arg(
    Arg::new("re_site")
        .long("re-site")
        .help("restriction enzyme motif for Pore-C reads, comma-separated for multiple motifs (e.g. AAGCTT)")
        .value_parser(value_parser!(String))
        .value_name("STR")
)
.next_help_heading("Methylation refinement")
.arg(Arg::new("meth_bed").long("meth-bed").value_name("BEDGRAPH")
    .help("Reference methylation bedGraph; requires unaligned BAM with MM/ML tags"))
.arg(Arg::new("meth_match").long("meth-match-score").value_parser(value_parser!(i32)).default_value("0"))
.arg(Arg::new("meth_ref_penalty").long("meth-ref-penalty").value_parser(value_parser!(i32)).default_value("2"))
.arg(Arg::new("meth_read_penalty").long("meth-read-penalty").value_parser(value_parser!(i32)).default_value("2"))
.arg(Arg::new("meth_ref_cutoff").long("meth-ref-prob-cutoff").value_parser(value_parser!(f64)).default_value("50"))
.arg(Arg::new("meth_cutoff").long("meth-prob-cutoff").value_parser(value_parser!(u8)).default_value("128"))
.arg(Arg::new("meth_mapq").long("meth-designate-mapq").value_parser(value_parser!(u8)).default_value("2"))
.arg(Arg::new("meth_cpg").long("meth-cpg").action(ArgAction::SetTrue))
.next_help_heading("Presets")
.arg(
    Arg::new("preset")
        .short('x')
        .value_name("STR")
        .help(
r##"-
- porec/porec:map-ont/poerc:lr:hq - Pore-C reads, porec equal to porec:lr:hq
- cifi  - CiFi reads
- hic   - Hi-C/OmniC short reads
- lr:hq - accurate long reads (error rate <1%) against a reference genome
- splice/splice:hq - spliced alignment for long reads/accurate long reads
- asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
- sr - short reads against a reference
- map-pb/map-hifi/map-ont - CLR/HiFi/Nanopore vs reference mapping
- ava-pb/ava-ont - PacBio CLR/Nanopore read overlap
"##
    )
        .value_parser([
            "porec", "porec:lr:hq", "porec:map-ont", "cifi", "hic", "lr:hq", "lr:hqae", "map-hifi", "map-pb", "map-ont", "asm5", "asm10", "asm20", "sr", "splice", "splice:hq"
        ])
        .default_value("porec"),
)
.arg_required_else_help(true)

}
