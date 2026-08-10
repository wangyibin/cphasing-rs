use cphasing::alleles::{AllelesFasta, AllelesOptions, Anchor, Anchor32};
use cphasing::core::BaseTable;
use rayon::ThreadPoolBuilder;
use std::collections::HashSet;
use std::fs;

const SEQUENCE: &str = "TTTCGACAGTTCTCCCTGGCACCTCTGAAAGCTTTCCTGGTTTGATTGTTGAAAGTCTTAGGGCTCAACTTGGTCAGCCCTTCTTCATGGAAATTGTTATGACCATGTGTTGGTCCATCTGGATGATGCGCAATGATGTCATTTTCAAAGGTTTAC";

#[test]
fn identical_contigs_match_partig_counts_and_similarity() {
    assert_eq!(std::mem::size_of::<Anchor>(), 24);
    assert_eq!(std::mem::size_of::<Anchor32>(), 16);
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("identical.fa");
    let output = directory.path().join("identical.allele.table");
    fs::write(&fasta, format!(">s1\n{SEQUENCE}\n>s2\n{SEQUENCE}\n")).unwrap();

    let mut alleles = AllelesFasta::new(&fasta.to_string_lossy().into_owned());
    alleles
        .run_with_options(
            AllelesOptions {
                min_similarity: 0.0,
                diff_threshold: 0.0,
                ..AllelesOptions::default()
            },
            &output.to_string_lossy(),
        )
        .unwrap();

    let text = fs::read_to_string(output).unwrap();
    let lines = text.lines().collect::<Vec<_>>();
    assert_eq!(lines[0], "#s1 156 12 0");
    assert_eq!(lines[1], "#s2 156 12 0");
    assert_eq!(lines.len(), 4);

    for (index, line) in lines[2..].iter().enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[0], (index + 1).to_string());
        assert_eq!(fields[4..7], ["12", "12", "12"]);
        assert_eq!(fields[7], "1.0");
        assert_eq!(fields[8], "1");
    }
}

#[test]
fn reverse_complement_contigs_emit_identical_symmetric_reverse_matches() {
    let reverse_complement = SEQUENCE
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!(),
        })
        .collect::<String>();
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("reverse.fa");
    let output = directory.path().join("reverse.allele.table");
    fs::write(
        &fasta,
        format!(">z_reverse\n{reverse_complement}\n>a_forward\n{SEQUENCE}\n"),
    )
    .unwrap();

    let mut alleles = AllelesFasta::new(&fasta.to_string_lossy().into_owned());
    alleles
        .run_with_options(
            AllelesOptions {
                min_similarity: 0.0,
                diff_threshold: 0.0,
                ..AllelesOptions::default()
            },
            &output.to_string_lossy(),
        )
        .unwrap();

    let text = fs::read_to_string(output).unwrap();
    let matches = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(matches.len(), 2);
    assert_eq!(matches[0][2..4], ["z_reverse", "a_forward"]);
    assert_eq!(matches[1][2..4], ["a_forward", "z_reverse"]);
    assert_eq!(matches[0][4..8], matches[1][4..8]);
    assert_eq!(matches[0][8], "-1");
    assert_eq!(matches[1][8], "-1");
}

#[test]
fn internal_trimming_matches_a_pretrimmed_fasta() {
    const TRIM: usize = 10;
    let directory = tempfile::tempdir().unwrap();
    let untrimmed_fasta = directory.path().join("untrimmed.fa");
    let trimmed_fasta = directory.path().join("trimmed.fa");
    let internal_output = directory.path().join("internal.allele.table");
    let pretrimmed_output = directory.path().join("pretrimmed.allele.table");
    let padded = format!("{}{SEQUENCE}{}", "A".repeat(TRIM), "C".repeat(TRIM));
    fs::write(&untrimmed_fasta, format!(">s1\n{padded}\n>s2\n{padded}\n")).unwrap();
    fs::write(
        &trimmed_fasta,
        format!(">s1\n{SEQUENCE}\n>s2\n{SEQUENCE}\n"),
    )
    .unwrap();

    let mut internal = AllelesFasta::new(&untrimmed_fasta.to_string_lossy().into_owned());
    internal
        .run_with_options(
            AllelesOptions {
                trim_length: TRIM,
                min_similarity: 0.0,
                diff_threshold: 0.0,
                ..AllelesOptions::default()
            },
            &internal_output.to_string_lossy(),
        )
        .unwrap();

    let mut pretrimmed = AllelesFasta::new(&trimmed_fasta.to_string_lossy().into_owned());
    pretrimmed
        .run_with_options(
            AllelesOptions {
                min_similarity: 0.0,
                diff_threshold: 0.0,
                ..AllelesOptions::default()
            },
            &pretrimmed_output.to_string_lossy(),
        )
        .unwrap();

    assert_eq!(
        fs::read(internal_output).unwrap(),
        fs::read(pretrimmed_output).unwrap()
    );
}

#[test]
fn split_regions_match_a_materialized_split_fasta() {
    let directory = tempfile::tempdir().unwrap();
    let source_fasta = directory.path().join("source.fa");
    let materialized_fasta = directory.path().join("materialized.split.fa");
    let regions = directory.path().join("split.regions.tsv");
    let native_output = directory.path().join("native.allele.table");
    let materialized_output = directory.path().join("materialized.allele.table");
    fs::write(
        &source_fasta,
        format!(">ctg1 description\n{SEQUENCE}{SEQUENCE}\n>ctg2\n{SEQUENCE}\n"),
    )
    .unwrap();
    // Region rows need not be grouped by source. Output order follows the
    // source FASTA, while preserving row order within each source contig.
    fs::write(
        &regions,
        format!(
            "ctg2|0_156\tctg2\t0\t156\nctg1|156_312\tctg1\t156\t312\nctg1|0_156\tctg1\t0\t156\n"
        ),
    )
    .unwrap();
    fs::write(
        &materialized_fasta,
        format!(">ctg1|156_312\n{SEQUENCE}\n>ctg1|0_156\n{SEQUENCE}\n>ctg2|0_156\n{SEQUENCE}\n"),
    )
    .unwrap();

    let options = AllelesOptions {
        min_similarity: 0.0,
        diff_threshold: 0.0,
        ..AllelesOptions::default()
    };
    let mut native = AllelesFasta::new(&source_fasta.to_string_lossy().into_owned());
    native.set_split_regions(Some(regions.to_string_lossy().into_owned()));
    native
        .run_with_options(options, &native_output.to_string_lossy())
        .unwrap();

    let mut materialized = AllelesFasta::new(&materialized_fasta.to_string_lossy().into_owned());
    materialized
        .run_with_options(options, &materialized_output.to_string_lossy())
        .unwrap();

    assert_eq!(
        fs::read(native_output).unwrap(),
        fs::read(materialized_output).unwrap()
    );
}

#[test]
fn split_regions_report_missing_source_contigs() {
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("source.fa");
    let regions = directory.path().join("split.regions.tsv");
    let output = directory.path().join("native.allele.table");
    fs::write(&fasta, format!(">ctg1\n{SEQUENCE}\n")).unwrap();
    fs::write(&regions, "missing|0_10\tmissing\t0\t10\n").unwrap();

    let mut alleles = AllelesFasta::new(&fasta.to_string_lossy().into_owned());
    alleles.set_split_regions(Some(regions.to_string_lossy().into_owned()));
    let error = alleles
        .run_with_options(AllelesOptions::default(), &output.to_string_lossy())
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("split-regions source contigs were not found in FASTA: missing")
    );
}

#[test]
fn radix_partition_preserves_every_directional_match() {
    const CONTIGS: usize = 110;
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("radix.fa");
    let output = directory.path().join("radix.allele.table");
    let input = (0..CONTIGS)
        .map(|index| format!(">s{index:03}\n{SEQUENCE}\n"))
        .collect::<String>();
    fs::write(&fasta, input).unwrap();

    let mut alleles = AllelesFasta::new(&fasta.to_string_lossy().into_owned());
    alleles
        .run_with_options(
            AllelesOptions {
                min_similarity: 0.0,
                max_occurrence: 128,
                diff_threshold: 0.0,
                ..AllelesOptions::default()
            },
            &output.to_string_lossy(),
        )
        .unwrap();

    let text = fs::read_to_string(output).unwrap();
    let mut pairs = HashSet::new();
    let mut headers = 0;
    for line in text.lines() {
        if line.starts_with('#') {
            headers += 1;
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[4..7], ["12", "12", "12"]);
        assert_eq!(fields[7].parse::<f64>().unwrap(), 1.0);
        assert_eq!(fields[8], "1");
        assert!(pairs.insert((fields[2].to_owned(), fields[3].to_owned())));
    }

    assert_eq!(headers, CONTIGS);
    assert_eq!(pairs.len(), CONTIGS * (CONTIGS - 1));
}

#[test]
fn pipelined_sketch_preserves_input_order_across_thread_counts() {
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("unbalanced.fa");
    let serial_output = directory.path().join("serial.allele.table");
    let parallel_output = directory.path().join("parallel.allele.table");
    let slow = SEQUENCE.repeat(128);
    let middle = SEQUENCE.repeat(3);
    fs::write(
        &fasta,
        format!(
            ">slow_first description\n{slow}\n>fast_second\n{SEQUENCE}\n>middle_third\n{middle}\n"
        ),
    )
    .unwrap();

    let options = AllelesOptions {
        min_similarity: 0.0,
        diff_threshold: 0.0,
        ..AllelesOptions::default()
    };
    for (threads, output) in [(1, &serial_output), (4, &parallel_output)] {
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let mut alleles = AllelesFasta::new(&fasta.to_string_lossy().into_owned());
        pool.install(|| {
            alleles
                .run_with_options(options, &output.to_string_lossy())
                .unwrap()
        });
    }

    let serial = fs::read(&serial_output).unwrap();
    let parallel = fs::read(&parallel_output).unwrap();
    assert_eq!(serial, parallel);
    let headers = std::str::from_utf8(&parallel)
        .unwrap()
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| line.split_ascii_whitespace().next().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(headers, ["#slow_first", "#fast_second", "#middle_third"]);
}

#[test]
fn trim_threshold_and_saturating_length_are_unchanged() {
    const TRIM: usize = 4;
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("trim-boundary.fa");
    fs::write(
        &fasta,
        ">exact description\nACGTACGTACGT\n>over\nACGTACGTACGTA\n",
    )
    .unwrap();

    let mut alleles = AllelesFasta::new(&fasta.to_string_lossy().into_owned());
    let sequences = alleles.seqs(TRIM).unwrap();
    assert_eq!(sequences[0], b"ACGTACGTACGT");
    assert_eq!(sequences[1], b"ACGTA");
    assert_eq!(alleles.contigs, ["exact", "over"]);
    assert_eq!(alleles.contig_lengths["exact"], 12);
    assert_eq!(alleles.contig_lengths["over"], 5);

    let untrimmed = alleles.seqs(usize::MAX).unwrap();
    assert_eq!(untrimmed[0], b"ACGTACGTACGT");
    assert_eq!(untrimmed[1], b"ACGTACGTACGTA");
}

#[test]
fn pipelined_parse_error_does_not_commit_partial_metadata() {
    let directory = tempfile::tempdir().unwrap();
    let fasta = directory.path().join("invalid.fa");
    let output = directory.path().join("invalid.allele.table");
    let mut input = format!(">valid\n{SEQUENCE}\n>").into_bytes();
    input.push(0xff);
    input.extend_from_slice(b"\nACGT\n");
    fs::write(&fasta, input).unwrap();

    let mut alleles = AllelesFasta::new(&fasta.to_string_lossy().into_owned());
    let error = alleles
        .run_with_options(AllelesOptions::default(), &output.to_string_lossy())
        .unwrap_err();
    assert!(error.to_string().contains("record 2"), "{error:#}");
    assert!(alleles.contigs.is_empty());
    assert!(alleles.contig_lengths.is_empty());
    assert!(!output.exists());
}
