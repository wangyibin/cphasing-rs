use std::fs::{self, File};
use std::io::Read;
use std::process::Command;

use cphasing::clm::ClmbReader;
use flate2::read::GzDecoder;
use polars::prelude::*;

fn run_pairs2clm(input: &std::path::Path, output: &std::path::Path, use_cn: bool, threads: usize) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"));
    command
        .arg("pairs2clm")
        .arg(input)
        .arg("--min-contacts")
        .arg("1")
        .arg("--no-output-split")
        .arg("--disable-filter")
        .arg("--threads")
        .arg(threads.to_string())
        .arg("--output")
        .arg(output);
    if use_cn {
        command.arg("--use-cn");
    }
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn pairs2clm_writes_clmb_directly() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.clmb");
    let input = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("test/test.pairs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2clm")
        .arg(input)
        .arg("--min-contacts")
        .arg("1")
        .arg("--no-output-split")
        .arg("--disable-filter")
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let mut reader = ClmbReader::open(output).unwrap();
    let mut records = 0usize;
    let mut distances = 0usize;
    while let Some(block) = reader.next_block().unwrap() {
        for record in block {
            records += 1;
            distances += record.distances.len();
        }
    }
    assert!(records > 0);
    assert_eq!(records % 4, 0);
    assert!(distances >= records);
}

#[test]
fn pairs2clm_preserves_coordinates_beyond_u32() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("large.pairs");
    let output = directory.path().join("large.clmb");
    fs::write(
        &input,
        "## pairs format 1.0\n\
#shape: upper triangle\n\
#chromsize: longA 5000000000\n\
#chromsize: longB 6000000000\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n\
r1\tlongA\t4500000000\tlongB\t5500000000\t+\t-\n",
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2clm")
        .arg(&input)
        .arg("--min-contacts")
        .arg("1")
        .arg("--no-output-split")
        .arg("--disable-filter")
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut reader = ClmbReader::open(output).unwrap();
    let mut maximum = 0u64;
    while let Some(block) = reader.next_block().unwrap() {
        for record in block {
            maximum = maximum.max(record.distances.into_iter().max().unwrap());
        }
    }
    assert!(maximum > u32::MAX as u64);
}

#[test]
fn pairs2clm_reads_u64_pqs_coordinates_without_narrowing() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("large.pairs.pqs");
    let output = directory.path().join("large.clmb");
    fs::create_dir_all(input.join("q0")).unwrap();
    fs::create_dir_all(input.join("q1")).unwrap();
    fs::write(
        input.join("_contigsizes"),
        "longA\t5000000000\nlongB\t6000000000\n",
    )
    .unwrap();
    fs::write(
        input.join("_metadata_counts"),
        "q0_records\t1\nq1_records\t1\n",
    )
    .unwrap();

    let categorical = |name: &str, values: &[&str]| {
        Series::new(name.into(), values)
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
            .unwrap()
    };
    let mut frame = DataFrame::new(vec![
        categorical("chrom1", &["longA"]).into(),
        categorical("chrom2", &["longB"]).into(),
        Series::new("pos1".into(), [4_500_000_000u64]).into(),
        Series::new("pos2".into(), [5_500_000_000u64]).into(),
        Series::new("mapq".into(), [60u8]).into(),
    ])
    .unwrap();
    ParquetWriter::new(File::create(input.join("q0/0.parquet")).unwrap())
        .finish(&mut frame)
        .unwrap();
    ParquetWriter::new(File::create(input.join("q1/0.parquet")).unwrap())
        .finish(&mut frame)
        .unwrap();

    run_pairs2clm(&input, &output, false, 1);
    let mut reader = ClmbReader::open(output).unwrap();
    let mut maximum = 0u64;
    while let Some(block) = reader.next_block().unwrap() {
        for record in block {
            maximum = maximum.max(record.distances.into_iter().max().unwrap());
        }
    }
    assert!(maximum > u32::MAX as u64);
}

#[test]
fn pairs2clm_writes_clmb_from_pqs_shards() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.pairs");
    let pqs = directory.path().join("input.pairs.pqs");
    let output = directory.path().join("output.clmb");
    fs::write(
        &input,
        "## pairs format 1.0\n\
#shape: upper triangle\n\
#chromsize: A 1000\n\
#chromsize: B 1200\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t10\tB\t20\t+\t-\t60\n\
r2\tA\t30\tB\t40\t-\t+\t60\n\
r3\tA\t50\tB\t60\t+\t+\t0\n",
    )
    .unwrap();
    let convert = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(input)
        .arg("--chunksize")
        .arg("100")
        .arg("--output")
        .arg(&pqs)
        .output()
        .unwrap();
    assert!(
        convert.status.success(),
        "{}",
        String::from_utf8_lossy(&convert.stderr)
    );

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2clm")
        .arg(&pqs)
        .arg("--min-contacts")
        .arg("1")
        .arg("--no-output-split")
        .arg("--disable-filter")
        .arg("--threads")
        .arg("4")
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let mut reader = ClmbReader::open(output).unwrap();
    let mut records = 0usize;
    let mut distances = 0usize;
    while let Some(block) = reader.next_block().unwrap() {
        for record in block {
            records += 1;
            distances += record.distances.len();
        }
    }
    assert!(records > 0);
    assert_eq!(records % 4, 0);
    assert!(distances >= records);
}

#[test]
fn pairs2clm_applies_optional_cn_info_only_at_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.pairs");
    let pqs = directory.path().join("input.pairs.pqs");
    fs::write(
        &input,
        "## pairs format 1.0\n\
#shape: upper triangle\n\
#chromsize: A 1000\n\
#chromsize: B 1200\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t10\tB\t20\t+\t-\t60\n\
r2\tA\t30\tB\t50\t+\t-\t60\n\
r3\tA\t70\tB\t110\t+\t-\t60\n\
r4\tA\t130\tB\t170\t+\t-\t60\n\
r5\tA\t190\tB\t230\t+\t-\t60\n\
r6\tA\t290\tB\t310\t+\t-\t60\n\
r7\tA\t370\tB\t410\t+\t-\t60\n",
    )
    .unwrap();
    let convert = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(&input)
        .arg("--chunksize")
        .arg("2")
        .arg("--output")
        .arg(&pqs)
        .output()
        .unwrap();
    assert!(
        convert.status.success(),
        "{}",
        String::from_utf8_lossy(&convert.stderr)
    );
    fs::write(pqs.join("cn.info"), "A\t2\nB\t2\n").unwrap();

    let raw_output = directory.path().join("raw.clm");
    run_pairs2clm(&pqs, &raw_output, false, 1);
    let raw = fs::read_to_string(&raw_output).unwrap();
    assert!(!raw.contains("A_d2"));
    assert_eq!(raw.lines().count(), 4);

    let cn_output_1 = directory.path().join("cn.t1.clm");
    let cn_output_4 = directory.path().join("cn.t4.clm");
    run_pairs2clm(&pqs, &cn_output_1, true, 1);
    run_pairs2clm(&pqs, &cn_output_4, true, 4);
    let cn_text = fs::read_to_string(&cn_output_1).unwrap();
    assert_eq!(cn_text, fs::read_to_string(&cn_output_4).unwrap());
    assert!(cn_text.contains("A_d2"));
    assert!(cn_text.contains("B_d2"));
    assert_eq!(cn_text.lines().count(), 16);

    let mut total_counts = 0usize;
    let mut total_distances = 0usize;
    for line in cn_text.lines() {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), 3);
        let count = fields[1].parse::<usize>().unwrap();
        let distance_count = fields[2].split_whitespace().count();
        assert_eq!(count, distance_count);
        assert!(matches!(count, 1 | 2));
        total_counts += count;
        total_distances += distance_count;
    }
    assert_eq!(total_counts, 7 * 4);
    assert_eq!(total_distances, 7 * 4);

    let clmb_output = directory.path().join("cn.clmb");
    run_pairs2clm(&pqs, &clmb_output, true, 2);
    let mut reader = ClmbReader::open(&clmb_output).unwrap();
    assert!(reader.header.contigs.iter().any(|contig| contig == "A_d2"));
    assert!(reader.header.contigs.iter().any(|contig| contig == "B_d2"));
    let mut records = 0usize;
    let mut distances = 0usize;
    while let Some(block) = reader.next_block().unwrap() {
        for record in block {
            records += 1;
            assert!(matches!(record.distances.len(), 1 | 2));
            distances += record.distances.len();
        }
    }
    assert_eq!(records, 16);
    assert_eq!(distances, 7 * 4);

    let auxiliary_output = directory.path().join("aux.clm");
    let auxiliary = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2clm")
        .arg(&pqs)
        .arg("--use-cn")
        .arg("--min-contacts")
        .arg("1")
        .arg("--disable-filter")
        .arg("--output-depth")
        .arg("--output")
        .arg(&auxiliary_output)
        .output()
        .unwrap();
    assert!(
        auxiliary.status.success(),
        "{}",
        String::from_utf8_lossy(&auxiliary.stderr)
    );
    let mut split_text = String::new();
    GzDecoder::new(fs::File::open(directory.path().join("aux.split.contacts.gz")).unwrap())
        .read_to_string(&mut split_text)
        .unwrap();
    assert!(split_text.contains("A_d2_0"));
    assert!(split_text.contains("B_d2_0"));
    assert_eq!(
        split_text
            .lines()
            .map(|line| line.split('\t').nth(2).unwrap().parse::<usize>().unwrap())
            .sum::<usize>(),
        7
    );
    let depth_text = fs::read_to_string(directory.path().join("aux.depth")).unwrap();
    assert!(depth_text.lines().any(|line| line.starts_with("A_d2\t")));
    assert!(depth_text.lines().any(|line| line.starts_with("B_d2\t")));
}

#[test]
fn pairs2clm_use_cn_treats_missing_cn_info_as_cn_one() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.pairs");
    let pqs = directory.path().join("input.pairs.pqs");
    let output = directory.path().join("output.clm");
    fs::write(
        &input,
        "## pairs format 1.0\n\
#shape: upper triangle\n\
#chromsize: A 1000\n\
#chromsize: B 1200\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t10\tB\t20\t+\t-\t60\n",
    )
    .unwrap();
    let convert = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(&input)
        .arg("--output")
        .arg(&pqs)
        .output()
        .unwrap();
    assert!(convert.status.success());

    run_pairs2clm(&pqs, &output, true, 1);
    let text = fs::read_to_string(output).unwrap();
    assert_eq!(text.lines().count(), 4);
    assert!(!text.contains("_d2"));
}
