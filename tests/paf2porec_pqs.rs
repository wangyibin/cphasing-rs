use std::fs::{self, File};
use std::io::Write;
use std::process::{Command, Stdio};
use std::thread;
use std::time::Duration;

use flate2::Compression;
use flate2::write::GzEncoder;
use polars::prelude::*;

fn run_paf2porec(paf: &std::path::Path, output: &std::path::Path) {
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("paf2porec")
        .arg(paf)
        .arg("--output")
        .arg(output)
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn porec_chr2ctg_preserves_zero_based_target_start() {
    let directory = tempfile::Builder::new()
        .prefix("porec-chr2ctg-zero-start-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let paf = directory.path().join("input.paf");
    fs::write(
        &paf,
        "read1\t1000\t0\t300\t+\tA\t1000\t0\t300\t290\t300\t60\ttp:A:P\n\
read1\t1000\t310\t610\t+\tB\t1000\t100\t400\t290\t300\t60\ttp:A:P\n",
    )
    .unwrap();
    let mapping = directory.path().join("mapping.bed");
    fs::write(&mapping, "A\t0\t500\tA1\nA\t500\t1000\tA2\nB\t0\t500\tB1\n").unwrap();

    let text_input = directory.path().join("input.concat");
    let text_output = directory.path().join("converted.concat");
    run_paf2porec(&paf, &text_input);
    let text_result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-chr2ctg")
        .arg("--input")
        .arg(&text_input)
        .arg("--bed")
        .arg(&mapping)
        .arg("--output")
        .arg(&text_output)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        text_result.status.success(),
        "{}",
        String::from_utf8_lossy(&text_result.stderr)
    );
    let text = fs::read_to_string(&text_output).unwrap();
    let records = text
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(records.len(), 2);
    assert_eq!(&records[0][5..8], &["A1", "0", "300"]);

    let boundary_input = directory.path().join("boundary.concat");
    let boundary_output = directory.path().join("boundary-converted.concat");
    fs::write(
        &boundary_input,
        "0\t1000\t0\t1\t+\tA\t500\t501\t60\t1\tPASS\n",
    )
    .unwrap();
    let boundary_result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-chr2ctg")
        .arg("--input")
        .arg(&boundary_input)
        .arg("--bed")
        .arg(&mapping)
        .arg("--output")
        .arg(&boundary_output)
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        boundary_result.status.success(),
        "{}",
        String::from_utf8_lossy(&boundary_result.stderr)
    );
    let boundary = fs::read_to_string(&boundary_output).unwrap();
    let boundary_fields = boundary.trim_end().split('\t').collect::<Vec<_>>();
    assert_eq!(&boundary_fields[5..8], &["A2", "0", "1"]);

    let pqs_input = directory.path().join("input.concat.pqs");
    let pqs_output = directory.path().join("converted.concat.pqs");
    run_paf2porec(&paf, &pqs_input);
    let pqs_result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-chr2ctg")
        .arg("--input")
        .arg(&pqs_input)
        .arg("--bed")
        .arg(&mapping)
        .arg("--output")
        .arg(&pqs_output)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        pqs_result.status.success(),
        "{}",
        String::from_utf8_lossy(&pqs_result.stderr)
    );
    let frame = ParquetReader::new(File::open(pqs_output.join("q0/0.parquet")).unwrap())
        .finish()
        .unwrap();
    assert_eq!(frame.height(), 2);
    assert_eq!(
        frame.column("start").unwrap().u64().unwrap().get(0),
        Some(0)
    );
    let chrom = frame
        .column("chrom")
        .unwrap()
        .cast(&DataType::String)
        .unwrap();
    assert_eq!(chrom.str().unwrap().get(0), Some("A1"));
}

#[test]
fn paf2porec_suffixes_and_downstream_reader_support_text_and_concat_pqs() {
    let directory = tempfile::Builder::new()
        .prefix("paf2porec-pqs-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let paf = directory.path().join("input.paf");
    fs::write(
        &paf,
        "read1\t1000\t0\t300\t+\tA\t1000\t10\t310\t290\t300\t60\ttp:A:P\n\
read1\t1000\t310\t610\t-\tB\t2000\t20\t320\t285\t300\t50\ttp:A:P\n\
read2\t900\t0\t250\t+\tA\t1000\t100\t350\t240\t250\t40\ttp:A:P\n\
read2\t900\t260\t560\t+\tC\t3000\t200\t500\t290\t300\t30\ttp:A:P\n",
    )
    .unwrap();

    let plain_text = directory.path().join("output.concat");
    let text = directory.path().join("output.concat.gz");
    let concat_pqs = directory.path().join("output.concat.pqs");
    let porec_pqs = directory.path().join("alias.porec.pqs");
    run_paf2porec(&paf, &plain_text);
    run_paf2porec(&paf, &text);
    run_paf2porec(&paf, &concat_pqs);
    run_paf2porec(&paf, &porec_pqs);
    assert_eq!(fs::read_to_string(&plain_text).unwrap().lines().count(), 4);

    for output in [&concat_pqs, &porec_pqs] {
        assert!(output.join("_metadata").is_file());
        assert!(output.join("_contigsizes").is_file());
        assert!(!std::path::Path::new(&format!("{}.tmp.porec.gz", output.display())).exists());
        assert!(
            fs::read_to_string(output.join("_metadata"))
                .unwrap()
                .contains("'format': 'concat'")
        );
        let frame = ParquetReader::new(File::open(output.join("q0/0.parquet")).unwrap())
            .finish()
            .unwrap();
        assert_eq!(frame.height(), 4);
        assert_eq!(frame.column("read_idx").unwrap().dtype(), &DataType::UInt64);
        assert_eq!(frame.column("start").unwrap().dtype(), &DataType::UInt64);
    }
    assert_eq!(
        fs::read_to_string(concat_pqs.join("_contigsizes")).unwrap(),
        "A\t1000\nB\t2000\nC\t3000\n"
    );
    let cloned_pqs = directory.path().join("cloned.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pqs")
        .arg(&concat_pqs)
        .arg(concat_pqs.join("_contigsizes"))
        .arg("--output")
        .arg(&cloned_pqs)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        fs::read_to_string(concat_pqs.join("_metadata_counts")).unwrap(),
        fs::read_to_string(cloned_pqs.join("_metadata_counts")).unwrap()
    );
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        assert_eq!(
            fs::metadata(concat_pqs.join("q0/0.parquet")).unwrap().ino(),
            fs::metadata(cloned_pqs.join("q0/0.parquet")).unwrap().ino()
        );
    }

    let text_depth = directory.path().join("text.depth");
    let pqs_depth = directory.path().join("pqs.depth");
    for (input, output) in [(&text, &text_depth), (&concat_pqs, &pqs_depth)] {
        let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .arg("porec2depth")
            .arg(input)
            .arg(concat_pqs.join("_contigsizes"))
            .arg("--winsize")
            .arg("100")
            .arg("--stepsize")
            .arg("50")
            .arg("--threads")
            .arg("2")
            .arg("--output")
            .arg(output)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    assert_eq!(
        fs::read_to_string(&text_depth).unwrap(),
        fs::read_to_string(&pqs_depth).unwrap()
    );
    assert_eq!(
        fs::read_to_string(&pqs_depth)
            .unwrap()
            .lines()
            .take(8)
            .collect::<Vec<_>>(),
        vec![
            "A\t0\t100\t0.900",
            "A\t50\t150\t1.500",
            "A\t100\t200\t2.000",
            "A\t150\t250\t2.000",
            "A\t200\t300\t2.000",
            "A\t250\t350\t1.600",
            "A\t300\t400\t0.600",
            "A\t350\t450\t0.000",
        ]
    );
    let filtered_depth = directory.path().join("pqs.q45.depth");
    let filtered_result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2depth")
        .arg(&concat_pqs)
        .arg(concat_pqs.join("_contigsizes"))
        .arg("--winsize")
        .arg("100")
        .arg("--stepsize")
        .arg("50")
        .arg("--min-mapq")
        .arg("45")
        .arg("--threads")
        .arg("2")
        .arg("--output")
        .arg(&filtered_depth)
        .output()
        .unwrap();
    assert!(
        filtered_result.status.success(),
        "{}",
        String::from_utf8_lossy(&filtered_result.stderr)
    );
    assert_eq!(
        fs::read_to_string(&filtered_depth)
            .unwrap()
            .lines()
            .take(8)
            .collect::<Vec<_>>(),
        vec![
            "A\t0\t100\t0.900",
            "A\t50\t150\t1.000",
            "A\t100\t200\t1.000",
            "A\t150\t250\t1.000",
            "A\t200\t300\t1.000",
            "A\t250\t350\t0.600",
            "A\t300\t400\t0.100",
            "A\t350\t450\t0.000",
        ]
    );
    let porec_metadata_path = porec_pqs.join("_metadata");
    let porec_metadata = fs::read_to_string(&porec_metadata_path)
        .unwrap()
        .replace("'format': 'concat'", "'format': 'porec'");
    fs::write(&porec_metadata_path, porec_metadata).unwrap();

    let text_pairs = directory.path().join("text.pairs");
    let pqs_pairs = directory.path().join("pqs.pairs");
    for (input, output) in [(&text, &text_pairs), (&concat_pqs, &pqs_pairs)] {
        let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .arg("porec2pairs")
            .arg(input)
            .arg(concat_pqs.join("_contigsizes"))
            .arg("--output")
            .arg(output)
            .arg("--threads")
            .arg("1")
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    assert_eq!(
        fs::read_to_string(&text_pairs).unwrap(),
        fs::read_to_string(&pqs_pairs).unwrap()
    );

    let pairs_pqs = directory.path().join("from-concat.pairs.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pairs")
        .arg(&porec_pqs)
        .arg(concat_pqs.join("_contigsizes"))
        .arg("--output")
        .arg(&pairs_pqs)
        .arg("--threads")
        .arg("1")
        .arg("--chunksize")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        fs::read_to_string(pairs_pqs.join("_metadata"))
            .unwrap()
            .contains("'format': 'pairs'")
    );
    assert!(pairs_pqs.join("q0/0.parquet").is_file());
    assert!(
        fs::read_to_string(pairs_pqs.join("_metadata_counts"))
            .unwrap()
            .contains("q0_records\t2\n")
    );

    let break_bed = directory.path().join("break.bed");
    fs::write(concat_pqs.join("cn.info"), "A\t3\nB\t2\n").unwrap();
    fs::write(&break_bed, "A\t0\t500\tA_1\n").unwrap();
    let broken_pqs = directory.path().join("broken.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-break")
        .arg(&concat_pqs)
        .arg(&break_bed)
        .arg("--output")
        .arg(&broken_pqs)
        .arg("--chunksize")
        .arg("1")
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        fs::read_to_string(broken_pqs.join("_contigsizes"))
            .unwrap()
            .contains("A_1\t501\n")
    );
    assert_eq!(
        fs::read_to_string(concat_pqs.join("_metadata_counts")).unwrap(),
        fs::read_to_string(broken_pqs.join("_metadata_counts")).unwrap()
    );
    assert_eq!(
        fs::read_to_string(broken_pqs.join("cn.info")).unwrap(),
        "A_1\t3\nB\t2\n"
    );

    let collapsed = directory.path().join("collapsed.tsv");
    fs::write(&collapsed, "C\tC_d2\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-dup")
        .arg(&broken_pqs)
        .arg(&collapsed)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(result.stdout.is_empty());
    assert_eq!(
        fs::read_to_string(broken_pqs.join("cn.info")).unwrap(),
        "A_1\t3\nB\t2\nC\t2\n"
    );
    let duplicated_pqs = directory.path().join("duplicated.porec.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-dup")
        .arg(&broken_pqs)
        .arg(&collapsed)
        .arg("--output")
        .arg(&duplicated_pqs)
        .arg("--chunksize")
        .arg("1")
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        fs::read_to_string(duplicated_pqs.join("_contigsizes"))
            .unwrap()
            .contains("C_d2\t3000\n")
    );
    assert_eq!(
        fs::read_to_string(broken_pqs.join("_metadata")).unwrap(),
        fs::read_to_string(duplicated_pqs.join("_metadata")).unwrap()
    );
    assert_eq!(
        fs::read_to_string(duplicated_pqs.join("cn.info")).unwrap(),
        "A_1\t3\nB\t2\nC\t1\nC_d2\t1\n"
    );
    let mut seen_reads = std::collections::HashSet::new();
    let mut record_count = 0;
    for entry in fs::read_dir(duplicated_pqs.join("q0")).unwrap() {
        let frame = ParquetReader::new(File::open(entry.unwrap().path()).unwrap())
            .finish()
            .unwrap();
        record_count += frame.height();
        let shard_reads = frame
            .column("read_idx")
            .unwrap()
            .u64()
            .unwrap()
            .into_no_null_iter()
            .collect::<std::collections::HashSet<_>>();
        assert!(
            shard_reads
                .iter()
                .all(|read_idx| seen_reads.insert(*read_idx))
        );
    }
    assert_eq!(record_count, 4);
    assert_eq!(seen_reads.len(), 2);
    let duplicated_q0 =
        ParquetReader::new(File::open(duplicated_pqs.join("q0/0.parquet")).unwrap())
            .finish()
            .unwrap();
    let duplicated_q1 =
        ParquetReader::new(File::open(duplicated_pqs.join("q1/0.parquet")).unwrap())
            .finish()
            .unwrap();
    assert_eq!(
        duplicated_q0.column("chrom").unwrap(),
        duplicated_q1.column("chrom").unwrap()
    );

    let native_merge = directory.path().join("native-merge.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-merge")
        .arg(&concat_pqs)
        .arg(&porec_pqs)
        .arg("--output")
        .arg(&native_merge)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let counts = fs::read_to_string(native_merge.join("_metadata_counts")).unwrap();
    assert!(counts.contains("q0_records\t8"));
    assert!(counts.contains("q0_concats\t4"));
    assert!(
        fs::read_to_string(native_merge.join("_metadata"))
            .unwrap()
            .contains("'read_idx_scope': 'shard'")
    );
    assert_eq!(
        fs::read_to_string(native_merge.join("cn.info")).unwrap(),
        "A\t3\nB\t2\n"
    );
    let mut native_q0 = fs::read_dir(native_merge.join("q0"))
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .collect::<Vec<_>>();
    native_q0.sort();
    assert_eq!(native_q0.len(), 2);
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        assert_eq!(
            fs::metadata(concat_pqs.join("q0/0.parquet")).unwrap().ino(),
            fs::metadata(&native_q0[0]).unwrap().ino()
        );
        assert_eq!(
            fs::metadata(porec_pqs.join("q0/0.parquet")).unwrap().ino(),
            fs::metadata(&native_q0[1]).unwrap().ino()
        );
    }

    let linked_pairs = directory.path().join("linked-merge.pairs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pairs")
        .arg(&native_merge)
        .arg(native_merge.join("_contigsizes"))
        .arg("--output")
        .arg(&linked_pairs)
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let linked_materialized = directory.path().join("linked-remapped.porec");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-merge")
        .arg(&native_merge)
        .arg("--output")
        .arg(&linked_materialized)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let linked_materialized = fs::read_to_string(linked_materialized).unwrap();
    assert_eq!(linked_materialized.lines().count(), 8);
    assert_eq!(
        linked_materialized
            .lines()
            .map(|line| line.split('\t').next().unwrap())
            .collect::<std::collections::HashSet<_>>()
            .len(),
        4
    );

    let merged = directory.path().join("merged.porec");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-merge")
        .arg(&text)
        .arg(&porec_pqs)
        .arg("--output")
        .arg(&merged)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(fs::read_to_string(merged).unwrap().lines().count(), 8);
}

#[test]
fn native_concat_pqs_transform_commands_preserve_complete_reads() {
    let directory = tempfile::Builder::new()
        .prefix("porec-pqs-transforms-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let paf = directory.path().join("input.paf");
    fs::write(
        &paf,
        "read1\t1000\t0\t300\t+\tA\t1000\t10\t310\t290\t300\t60\ttp:A:P\n\
read1\t1000\t310\t610\t-\tB\t2000\t20\t320\t285\t300\t50\ttp:A:P\n\
read2\t900\t0\t250\t+\tA\t1000\t100\t350\t240\t250\t40\ttp:A:P\n\
read2\t900\t260\t560\t+\tC\t3000\t200\t500\t290\t300\t30\ttp:A:P\n",
    )
    .unwrap();
    let input = directory.path().join("input.concat.pqs");
    run_paf2porec(&paf, &input);
    fs::write(input.join("cn.info"), "A\t2\n").unwrap();

    let split = directory.path().join("split.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-split")
        .arg(&input)
        .arg("--output")
        .arg(&split)
        .arg("--chunksize")
        .arg("1")
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        fs::read_to_string(input.join("_metadata_counts")).unwrap(),
        fs::read_to_string(split.join("_metadata_counts")).unwrap()
    );
    assert_eq!(fs::read_to_string(split.join("cn.info")).unwrap(), "A\t2\n");
    let mut split_reads = std::collections::HashSet::new();
    let split_shards = fs::read_dir(split.join("q0")).unwrap().collect::<Vec<_>>();
    assert_eq!(split_shards.len(), 2);
    for shard in split_shards {
        let frame = ParquetReader::new(File::open(shard.unwrap().path()).unwrap())
            .finish()
            .unwrap();
        let reads = frame
            .column("read_idx")
            .unwrap()
            .u64()
            .unwrap()
            .into_no_null_iter()
            .collect::<std::collections::HashSet<_>>();
        assert_eq!(reads.len(), 1);
        assert!(split_reads.insert(*reads.iter().next().unwrap()));
    }

    let regions = directory.path().join("regions.bed");
    fs::write(&regions, "A\t0\t400\n").unwrap();
    let intersected = directory.path().join("intersected.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-intersect")
        .arg(&input)
        .arg(&regions)
        .arg("--output")
        .arg(&intersected)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let counts = fs::read_to_string(intersected.join("_metadata_counts")).unwrap();
    assert!(counts.contains("q0_records\t2\n"));
    assert!(counts.contains("q0_concats\t2\n"));
    assert_eq!(
        fs::read_to_string(intersected.join("cn.info")).unwrap(),
        "A\t2\n"
    );

    let mapping = directory.path().join("mapping.bed");
    fs::write(
        &mapping,
        "A\t0\t500\tA1\nA\t500\t1000\tA2\nB\t0\t1000\tB1\nC\t0\t1000\tC1\n",
    )
    .unwrap();
    let converted = directory.path().join("converted.porec.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-chr2ctg")
        .arg("--input")
        .arg(&input)
        .arg("--bed")
        .arg(&mapping)
        .arg("--output")
        .arg(&converted)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        fs::read_to_string(converted.join("_contigsizes")).unwrap(),
        "A1\t500\nA2\t500\nB1\t1000\nC1\t1000\n"
    );
    let converted_frame = ParquetReader::new(File::open(converted.join("q0/0.parquet")).unwrap())
        .finish()
        .unwrap();
    assert_eq!(converted_frame.height(), 4);
    assert_eq!(
        fs::read_to_string(converted.join("cn.info")).unwrap(),
        "A1\t2\nA2\t2\n"
    );

    let downsampled = directory.path().join("downsampled.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-downsample")
        .arg(&input)
        .arg("--mode")
        .arg("reads")
        .arg("--reads")
        .arg("1")
        .arg("--output")
        .arg(&downsampled)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let counts = fs::read_to_string(downsampled.join("_metadata_counts")).unwrap();
    assert!(counts.contains("q0_records\t2\n"));
    assert!(counts.contains("q0_concats\t1\n"));
    assert_eq!(
        fs::read_to_string(downsampled.join("cn.info")).unwrap(),
        "A\t2\n"
    );

    let linked = directory.path().join("linked.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-merge")
        .arg(&input)
        .arg(&input)
        .arg("--output")
        .arg(&linked)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        fs::read_to_string(linked.join("cn.info")).unwrap(),
        "A\t2\n"
    );
    let linked_split = directory.path().join("linked-split.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-split")
        .arg(&linked)
        .arg("--output")
        .arg(&linked_split)
        .arg("--chunksize")
        .arg("100")
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(fs::read_dir(linked_split.join("q0")).unwrap().count(), 2);
    assert!(
        fs::read_to_string(linked_split.join("_metadata"))
            .unwrap()
            .contains("'read_idx_scope': 'shard'")
    );
    let linked_downsampled = directory.path().join("linked-downsampled.concat.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-downsample")
        .arg(&linked)
        .arg("--mode")
        .arg("reads")
        .arg("--reads")
        .arg("1")
        .arg("--output")
        .arg(&linked_downsampled)
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let counts = fs::read_to_string(linked_downsampled.join("_metadata_counts")).unwrap();
    assert!(counts.contains("q0_records\t2\n"));
    assert!(counts.contains("q0_concats\t1\n"));

    let fasta = directory.path().join("reference.fa");
    fs::write(
        &fasta,
        format!(
            ">A\n{}\n>B\n{}\n>C\n{}\n",
            "A".repeat(1000),
            "C".repeat(2000),
            "G".repeat(3000)
        ),
    )
    .unwrap();
    let reads_prefix = directory.path().join("reads");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2reads")
        .arg(&input)
        .arg(&fasta)
        .arg("--output")
        .arg(&reads_prefix)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        directory
            .path()
            .join("reads_R1.fa.gz")
            .metadata()
            .unwrap()
            .len()
            > 0
    );
    assert!(
        directory
            .path()
            .join("reads_R2.fa.gz")
            .metadata()
            .unwrap()
            .len()
            > 0
    );
}

#[test]
fn paf2porec_short_stdin_reads_preserve_complete_read_groups() {
    let directory = tempfile::Builder::new()
        .prefix("paf2porec-short-stdin-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let records = [
        "read1\t1000\t0\t300\t+\tA\t1000\t10\t310\t290\t300\t60\ttp:A:P\n",
        "read1\t1000\t310\t610\t-\tB\t2000\t20\t320\t285\t300\t50\ttp:A:P\n",
        "read2\t900\t0\t250\t+\tA\t1000\t100\t350\t240\t250\t40\ttp:A:P\n",
        "read2\t900\t260\t560\t+\tC\t3000\t200\t500\t290\t300\t30\ttp:A:P\n",
    ];

    let regular_paf = directory.path().join("regular.paf");
    fs::write(&regular_paf, records.concat()).unwrap();
    let regular_output = directory.path().join("regular.concat");
    run_paf2porec(&regular_paf, &regular_output);

    let streamed_output = directory.path().join("streamed.concat");
    let mut child = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("paf2porec")
        .arg("-")
        .arg("--output")
        .arg(&streamed_output)
        .arg("--threads")
        .arg("4")
        .stdin(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let mut stdin = child.stdin.take().unwrap();
    for record in records {
        stdin.write_all(record.as_bytes()).unwrap();
        stdin.flush().unwrap();
        thread::sleep(Duration::from_millis(20));
    }
    drop(stdin);
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    assert_eq!(
        fs::read_to_string(&streamed_output).unwrap(),
        fs::read_to_string(&regular_output).unwrap()
    );
    assert_eq!(
        fs::read_to_string(directory.path().join("streamed.read.summary")).unwrap(),
        fs::read_to_string(directory.path().join("regular.read.summary")).unwrap()
    );
}

#[test]
fn paf2porec_keeps_a_read_group_crossing_the_target_batch_boundary() {
    let directory = tempfile::Builder::new()
        .prefix("paf2porec-batch-boundary-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let target_boundary = 2 * 1024 * 1024;
    let mut paf_text = String::with_capacity(target_boundary + 1024);
    let mut filler_index = 0usize;
    while paf_text.len() + 160 < target_boundary {
        paf_text.push_str(&format!(
            "filler{filler_index}\t300\t0\t300\t+\tA\t1000\t10\t310\t290\t300\t60\ttp:A:P\n"
        ));
        filler_index += 1;
    }

    let crossing_qname = "x".repeat(target_boundary - paf_text.len() + 16);
    paf_text.push_str(&format!(
        "{crossing_qname}\t1000\t0\t300\t+\tA\t1000\t10\t310\t290\t300\t60\ttp:A:P\n"
    ));
    paf_text.push_str(&format!(
        "{crossing_qname}\t1000\t310\t610\t-\tB\t2000\t20\t320\t285\t300\t50\ttp:A:P\n"
    ));
    paf_text.push_str("sentinel\t300\t0\t300\t+\tC\t3000\t30\t330\t290\t300\t60\ttp:A:P\n");
    assert!(paf_text.len() > target_boundary);
    let paf_text = paf_text.repeat(3);

    let paf = directory.path().join("input.paf");
    fs::write(&paf, &paf_text).unwrap();
    let output = directory.path().join("output.concat");
    run_paf2porec(&paf, &output);

    let regular_text = fs::read_to_string(&output).unwrap();
    let rows = regular_text.lines().collect::<Vec<_>>();
    assert_eq!(rows.len(), 6);
    assert_eq!(
        rows[0].split('\t').next().unwrap(),
        rows[1].split('\t').next().unwrap()
    );

    let streamed_output = directory.path().join("streamed.concat");
    let mut child = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("paf2porec")
        .arg("-")
        .arg("--output")
        .arg(&streamed_output)
        .arg("--threads")
        .arg("4")
        .stdin(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(paf_text.as_bytes())
        .unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let streamed_text = fs::read_to_string(streamed_output).unwrap();
    let mut regular_fields = regular_text
        .lines()
        .map(|line| line.split_once('\t').unwrap().1)
        .collect::<Vec<_>>();
    let mut streamed_fields = streamed_text
        .lines()
        .map(|line| line.split_once('\t').unwrap().1)
        .collect::<Vec<_>>();
    regular_fields.sort_unstable();
    streamed_fields.sort_unstable();
    assert_eq!(streamed_fields, regular_fields);
    assert_eq!(
        fs::read_to_string(directory.path().join("streamed.read.summary")).unwrap(),
        fs::read_to_string(directory.path().join("output.read.summary")).unwrap()
    );

    #[cfg(unix)]
    {
        let fifo = directory.path().join("input-fifo.paf");
        let mkfifo = Command::new("mkfifo").arg(&fifo).output().unwrap();
        assert!(
            mkfifo.status.success(),
            "{}",
            String::from_utf8_lossy(&mkfifo.stderr)
        );
        let fifo_output = directory.path().join("fifo.concat");
        let child = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .arg("paf2porec")
            .arg(&fifo)
            .arg("--output")
            .arg(&fifo_output)
            .arg("--threads")
            .arg("4")
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let mut fifo_writer = fs::OpenOptions::new().write(true).open(&fifo).unwrap();
        fifo_writer.write_all(paf_text.as_bytes()).unwrap();
        drop(fifo_writer);
        let result = child.wait_with_output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let fifo_text = fs::read_to_string(fifo_output).unwrap();
        let mut fifo_fields = fifo_text
            .lines()
            .map(|line| line.split_once('\t').unwrap().1)
            .collect::<Vec<_>>();
        fifo_fields.sort_unstable();
        assert_eq!(fifo_fields, regular_fields);
        assert_eq!(
            fs::read_to_string(directory.path().join("fifo.read.summary")).unwrap(),
            fs::read_to_string(directory.path().join("output.read.summary")).unwrap()
        );
    }
}

#[test]
fn paf2porec_rapidgzip_backend_matches_flate2_for_ordinary_gzip() {
    let rapidgzip_available = Command::new("rapidgzip").arg("--version").output().is_ok();

    let directory = tempfile::Builder::new()
        .prefix("paf2porec-rapidgzip-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let records = "read1\t1000\t0\t300\t+\tA\t1000\t10\t310\t290\t300\t60\ttp:A:P\n\
read1\t1000\t310\t610\t-\tB\t2000\t20\t320\t285\t300\t50\ttp:A:P\n\
read2\t900\t0\t250\t+\tA\t1000\t100\t350\t240\t250\t40\ttp:A:P\n\
read2\t900\t260\t560\t+\tC\t3000\t200\t500\t290\t300\t30\ttp:A:P\n";
    let gzip_path = directory.path().join("input.paf.gz");
    let mut encoder = GzEncoder::new(File::create(&gzip_path).unwrap(), Compression::default());
    encoder.write_all(records.as_bytes()).unwrap();
    encoder.finish().unwrap();

    let run_backend = |backend: &str, output: &std::path::Path| {
        Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .arg("paf2porec")
            .arg(&gzip_path)
            .arg("--output")
            .arg(output)
            .arg("--threads")
            .arg("2")
            .env("CPHASING_GZIP_BACKEND", backend)
            .env("CPHASING_IO_THREADS", "2")
            .output()
            .unwrap()
    };

    let flate2_output = directory.path().join("flate2.concat");
    let flate2_result = run_backend("flate2", &flate2_output);
    assert!(
        flate2_result.status.success(),
        "{}",
        String::from_utf8_lossy(&flate2_result.stderr)
    );

    let fallback_output = directory.path().join("fallback.concat");
    let fallback_result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("paf2porec")
        .arg(&gzip_path)
        .arg("--output")
        .arg(&fallback_output)
        .arg("--threads")
        .arg("2")
        .env("CPHASING_GZIP_BACKEND", "rapidgzip")
        .env(
            "CPHASING_RAPIDGZIP",
            directory.path().join("missing-rapidgzip"),
        )
        .env("CPHASING_IO_THREADS", "2")
        .env("PATH", "")
        .output()
        .unwrap();
    assert!(
        fallback_result.status.success(),
        "{}",
        String::from_utf8_lossy(&fallback_result.stderr)
    );
    assert!(String::from_utf8_lossy(&fallback_result.stderr).contains("was not found"));
    assert_eq!(
        fs::read_to_string(&fallback_output).unwrap(),
        fs::read_to_string(&flate2_output).unwrap()
    );
    assert_eq!(
        fs::read_to_string(directory.path().join("fallback.read.summary")).unwrap(),
        fs::read_to_string(directory.path().join("flate2.read.summary")).unwrap()
    );

    if !rapidgzip_available {
        return;
    }

    let rapidgzip_output = directory.path().join("rapidgzip.concat");
    let rapidgzip_result = run_backend("rapidgzip", &rapidgzip_output);
    assert!(
        rapidgzip_result.status.success(),
        "{}",
        String::from_utf8_lossy(&rapidgzip_result.stderr)
    );
    assert!(String::from_utf8_lossy(&rapidgzip_result.stderr).contains("using rapidgzip"));
    assert_eq!(
        fs::read_to_string(&rapidgzip_output).unwrap(),
        fs::read_to_string(&flate2_output).unwrap()
    );
    assert_eq!(
        fs::read_to_string(directory.path().join("rapidgzip.read.summary")).unwrap(),
        fs::read_to_string(directory.path().join("flate2.read.summary")).unwrap()
    );
}

#[test]
fn porec2pairs_compact_parser_skips_unused_columns_and_preserves_pair_fields() {
    let directory = tempfile::Builder::new()
        .prefix("porec2pairs-compact-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let porec = directory.path().join("input.porec");
    fs::write(
        &porec,
        "1\tbad-length\tbad-start\tbad-end\t+\tB\t10\t30\t60\tbad-identity\tignored\n\
1\tbad-length\tbad-start\tbad-end\t-\tA\t40\t60\t50\tbad-identity\tignored\n\
2\tbad-length\tbad-start\tbad-end\t+\tA\t0\t10\t0\tbad-identity\tignored\n",
    )
    .unwrap();
    let chromsizes = directory.path().join("contigsizes");
    fs::write(&chromsizes, "A\t100\nB\t100\n").unwrap();
    let pairs = directory.path().join("output.pairs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pairs")
        .arg(&porec)
        .arg(&chromsizes)
        .arg("--output")
        .arg(&pairs)
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let data_lines = fs::read_to_string(&pairs)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_string)
        .collect::<Vec<_>>();
    assert_eq!(data_lines, ["1\tA\t50\tB\t20\t-\t+\t50"]);
    assert_eq!(
        fs::read_to_string(directory.path().join("output.concatemer.summary")).unwrap(),
        "2\t1"
    );
}

#[test]
fn porec2pairs_pqs_compact_parser_keeps_complete_chunks_and_unique_pair_ids() {
    let directory = tempfile::Builder::new()
        .prefix("porec2pairs-pqs-compact-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let porec = directory.path().join("input.porec");
    fs::write(
        &porec,
        "1\tbad-length\tbad-start\tbad-end\t+\tD\t30\t50\t60\tbad-identity\tignored\n\
1\tbad-length\tbad-start\tbad-end\t-\tA\t0\t20\t60\tbad-identity\tignored\n\
1\tbad-length\tbad-start\tbad-end\t+\tC\t20\t40\t60\tbad-identity\tignored\n\
1\tbad-length\tbad-start\tbad-end\t-\tB\t10\t30\t60\tbad-identity\tignored\n\
2\tbad-length\tbad-start\tbad-end\t+\tC\t50\t70\t0\tbad-identity\tignored\n\
2\tbad-length\tbad-start\tbad-end\t-\tA\t40\t60\t50\tbad-identity\tignored\n\
2\tbad-length\tbad-start\tbad-end\t+\tB\t45\t65\t40\tbad-identity\tignored\n",
    )
    .unwrap();
    let chromsizes = directory.path().join("contigsizes");
    fs::write(&chromsizes, "A\t100\nB\t100\nC\t100\nD\t100\n").unwrap();
    let pairs_pqs = directory.path().join("output.pairs.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pairs")
        .arg(&porec)
        .arg(&chromsizes)
        .arg("--output")
        .arg(&pairs_pqs)
        .arg("--chunksize")
        .arg("2")
        .arg("--threads")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    assert_eq!(
        fs::read_to_string(pairs_pqs.join("_metadata_counts")).unwrap(),
        "q0_records\t9\nq1_records\t7\n"
    );
    let summary =
        fs::read_to_string(directory.path().join("output.pairs.concatemer.summary")).unwrap();
    assert!(summary.lines().any(|line| line == "3\t1"));
    assert!(summary.lines().any(|line| line == "4\t1"));

    let mut q0_files = fs::read_dir(pairs_pqs.join("q0"))
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().and_then(|value| value.to_str()) == Some("parquet"))
        .collect::<Vec<_>>();
    q0_files.sort();
    assert_eq!(q0_files.len(), 2);

    let mut pair_ids = Vec::new();
    for path in q0_files {
        let frame = ParquetReader::new(File::open(path).unwrap())
            .finish()
            .unwrap();
        pair_ids.extend(
            frame
                .column("read_idx")
                .unwrap()
                .as_materialized_series()
                .str()
                .unwrap()
                .into_no_null_iter()
                .map(|value| value.parse::<u64>().unwrap()),
        );
    }
    pair_ids.sort_unstable();
    assert_eq!(pair_ids, (1..=9).collect::<Vec<_>>());
}
