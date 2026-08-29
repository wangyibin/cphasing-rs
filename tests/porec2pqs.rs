use std::collections::HashMap;
use std::fs::{self, File};
use std::process::Command;
use std::sync::{Arc, Mutex};

use cphasing::porec::{PoreCRecord, write_concat_pqs_record_batches};
use polars::prelude::*;

fn parquet_files(path: &std::path::Path) -> Vec<std::path::PathBuf> {
    let mut files = fs::read_dir(path)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "parquet")
        })
        .collect::<Vec<_>>();
    files.sort();
    files
}

#[test]
fn porec2pqs_preserves_full_records_and_never_splits_a_read_id() {
    let directory = tempfile::Builder::new()
        .prefix("porec2pqs-test-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let input = directory.path().join("input.porec");
    let chromsizes = directory.path().join("chromsizes.txt");
    let output = directory.path().join("input.concat.pqs");
    fs::write(&chromsizes, "A\t1000\nB\t2000\n").unwrap();
    fs::write(
        &input,
        "# read_idx read_length read_start read_end strand chrom start end mapq identity filter\n\
10\t1000\t0\t100\t+\tA\t10\t110\t60\t0.99\tpass\n\
10\t1000\t110\t210\t-\tB\t20\t120\t0\t0.95\tlow_mapq\n\
11\t900\t0\t90\t+\tA\t30\t120\t1\t0.91\tpass\n\
11\t900\t100\t190\t+\tB\t40\t130\t2\t0.92\tpass\n\
11\t900\t200\t290\t-\tA\t50\t140\t3\t0.93\tpass\n\
12\t800\t0\t80\t+\tB\t60\t140\t0\t0.90\tlow_mapq\n",
    )
    .unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pqs")
        .arg(&input)
        .arg(&chromsizes)
        .arg("--output")
        .arg(&output)
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

    let q0_files = parquet_files(&output.join("q0"));
    assert_eq!(q0_files.len(), 3);
    let expected_columns = [
        "read_idx",
        "read_length",
        "read_start",
        "read_end",
        "strand",
        "chrom",
        "start",
        "end",
        "mapping_quality",
        "identity",
        "filter_reason",
    ];
    let mut read_to_shard = HashMap::<u64, usize>::new();
    let mut q0_rows = 0usize;
    for (shard_id, path) in q0_files.iter().enumerate() {
        let frame = ParquetReader::new(File::open(path).unwrap())
            .finish()
            .unwrap();
        assert_eq!(frame.get_column_names_str(), expected_columns);
        q0_rows += frame.height();
        for read_idx in frame
            .column("read_idx")
            .unwrap()
            .as_materialized_series()
            .u64()
            .unwrap()
            .into_no_null_iter()
        {
            assert_eq!(
                *read_to_shard.entry(read_idx).or_insert(shard_id),
                shard_id,
                "read_idx {read_idx} was split across Parquet shards"
            );
        }
    }
    assert_eq!(q0_rows, 6);
    assert_eq!(read_to_shard.len(), 3);

    let q1_rows = parquet_files(&output.join("q1"))
        .into_iter()
        .map(|path| {
            let frame = ParquetReader::new(File::open(path).unwrap())
                .finish()
                .unwrap();
            assert!(
                frame
                    .column("mapping_quality")
                    .unwrap()
                    .as_materialized_series()
                    .u8()
                    .unwrap()
                    .into_no_null_iter()
                    .all(|mapq| mapq >= 1)
            );
            frame.height()
        })
        .sum::<usize>();
    assert_eq!(q1_rows, 4);

    let metadata = fs::read_to_string(output.join("_metadata")).unwrap();
    assert!(metadata.contains("'format': 'concat'"));
    assert!(metadata.contains("'complete_groups': True"));
    let counts = fs::read_to_string(output.join("_metadata_counts")).unwrap();
    assert!(counts.contains("q0_records\t6"));
    assert!(counts.contains("q1_records\t4"));
    assert!(counts.contains("q0_concats\t3"));
}

#[test]
fn porec2pqs_accepts_out_of_order_read_ids_without_splitting_groups() {
    let directory = tempfile::Builder::new()
        .prefix("porec2pqs-unsorted-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let input = directory.path().join("unsorted.porec");
    let chromsizes = directory.path().join("chromsizes.txt");
    let output = directory.path().join("unsorted.concat.pqs");
    fs::write(&chromsizes, "A\t1000\n").unwrap();
    fs::write(
        &input,
        "20\t100\t0\t10\t+\tA\t1\t11\t1\t0.9\tpass\n\
20\t100\t10\t20\t+\tA\t11\t21\t1\t0.9\tpass\n\
3\t100\t0\t10\t+\tA\t21\t31\t1\t0.9\tpass\n\
3\t100\t10\t20\t+\tA\t31\t41\t1\t0.9\tpass\n",
    )
    .unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pqs")
        .arg(&input)
        .arg(&chromsizes)
        .arg("--output")
        .arg(&output)
        .arg("--chunksize")
        .arg("2")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );

    let q0_files = parquet_files(&output.join("q0"));
    assert_eq!(q0_files.len(), 2);
    let read_ids = q0_files
        .iter()
        .map(|path| {
            let frame = ParquetReader::new(File::open(path).unwrap())
                .finish()
                .unwrap();
            let ids = frame
                .column("read_idx")
                .unwrap()
                .as_materialized_series()
                .u64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>();
            assert!(ids.iter().all(|read_idx| *read_idx == ids[0]));
            ids[0]
        })
        .collect::<Vec<_>>();
    assert_eq!(read_ids, vec![20, 3]);
}

#[test]
fn porec2pqs_parallel_parser_reports_bad_rows_and_removes_partial_output() {
    let directory = tempfile::Builder::new()
        .prefix("porec2pqs-malformed-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let input = directory.path().join("malformed.porec");
    let chromsizes = directory.path().join("chromsizes.txt");
    let output = directory.path().join("malformed.concat.pqs");
    fs::write(&chromsizes, "A\t1000\n").unwrap();
    fs::write(
        &input,
        "1\t100\t0\t10\t+\tA\t1\t11\t1\t0.9\tpass\n\
2\tmissing-fields\n",
    )
    .unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec2pqs")
        .arg(&input)
        .arg(&chromsizes)
        .arg("--output")
        .arg(&output)
        .arg("--threads")
        .arg("4")
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("expected 11 tab-separated fields"));
    assert!(!output.exists());
}

#[test]
fn direct_record_batches_keep_complete_reads_across_batch_boundaries() {
    let directory = tempfile::Builder::new()
        .prefix("direct-concat-pqs-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let output = directory.path().join("direct.concat.pqs");
    let record = |read_idx, start| PoreCRecord {
        read_idx,
        query_length: 100,
        query_start: start,
        query_end: start + 10,
        query_strand: '+',
        target: "A".to_string(),
        target_start: start as u64,
        target_end: (start + 10) as u64,
        mapq: 10,
        identity: 0.95,
        filter_reason: "pass".to_string(),
    };
    let batches = vec![
        vec![record(20, 0), record(20, 10)],
        vec![record(3, 0), record(3, 10), record(3, 20), record(4, 0)],
    ];
    let target_sizes = Arc::new(Mutex::new(HashMap::from([("A".to_string(), 1000)])));
    write_concat_pqs_record_batches(batches, target_sizes, output.to_str().unwrap(), 2, 2).unwrap();

    let files = parquet_files(&output.join("q0"));
    assert_eq!(files.len(), 3);
    let ids = files
        .iter()
        .map(|path| {
            let frame = ParquetReader::new(File::open(path).unwrap())
                .finish()
                .unwrap();
            let values = frame
                .column("read_idx")
                .unwrap()
                .as_materialized_series()
                .u64()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>();
            assert!(values.iter().all(|value| *value == values[0]));
            values[0]
        })
        .collect::<Vec<_>>();
    assert_eq!(ids, vec![20, 3, 4]);
}
