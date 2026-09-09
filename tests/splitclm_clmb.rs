use std::fs;

use cphasing::clm::{Clm, ClmbReader, ClmbRecord, ClmbWriter, convert_clm, encode_endpoint};
use cphasing::core::BaseTable;

fn read_records(path: impl AsRef<std::path::Path>) -> (Vec<String>, Vec<ClmbRecord>) {
    let mut reader = ClmbReader::open(path).unwrap();
    let contigs = reader.header.contigs.clone();
    let mut records = Vec::new();
    while let Some(block) = reader.next_block().unwrap() {
        records.extend(block);
    }
    (contigs, records)
}

#[test]
fn splitclm_splits_clmb_without_losing_binary_records() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.clmb");
    let clusters = directory.path().join("clusters.txt");
    let output = directory.path().join("split");
    let contigs = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    let group1_record = ClmbRecord {
        endpoint1: encode_endpoint(0, 1).unwrap(),
        endpoint2: encode_endpoint(1, 0).unwrap(),
        distances: vec![2, u32::MAX as u64 + 7],
    };
    let group2_record = ClmbRecord {
        endpoint1: encode_endpoint(1, 1).unwrap(),
        endpoint2: encode_endpoint(2, 1).unwrap(),
        distances: vec![9, 11],
    };
    let cross_group_record = ClmbRecord {
        endpoint1: encode_endpoint(0, 0).unwrap(),
        endpoint2: encode_endpoint(2, 0).unwrap(),
        distances: vec![13],
    };
    let mut writer = ClmbWriter::create_synchronous(&input, &contigs, 64, None, None).unwrap();
    for record in [&group1_record, &group2_record, &cross_group_record] {
        writer
            .write_record(record.endpoint1, record.endpoint2, &record.distances)
            .unwrap();
    }
    writer.finish().unwrap();
    fs::write(&clusters, "group1 A B\ngroup2 B C\nempty D\n").unwrap();

    Clm::new(&input.to_string_lossy().into_owned())
        .split_clm(
            &clusters.to_string_lossy().into_owned(),
            &output.to_string_lossy().into_owned(),
        )
        .unwrap();

    let (group1_contigs, group1_records) = read_records(output.join("group1.clmb"));
    let (group2_contigs, group2_records) = read_records(output.join("group2.clmb"));
    let (empty_contigs, empty_records) = read_records(output.join("empty.clmb"));
    assert_eq!(group1_contigs, contigs);
    assert_eq!(group2_contigs, contigs);
    assert_eq!(empty_contigs, contigs);
    assert_eq!(group1_records, vec![group1_record]);
    assert_eq!(group2_records, vec![group2_record]);
    assert!(empty_records.is_empty());
}

#[test]
fn splitclm_preserves_or_overrides_each_input_format() {
    let directory = tempfile::tempdir().unwrap();
    let text = directory.path().join("source.clm");
    let clusters = directory.path().join("clusters.txt");
    fs::write(&text, "A- B+\t2\t2 4294967302\nB- C-\t2\t9 11\nA+ C+\t1\t13\n").unwrap();
    fs::write(&clusters, "group1 A B\ngroup2 B C\nempty D\n").unwrap();
    for input_format in ["clm", "clm.gz", "clmb"] {
        let input = directory.path().join(format!("input.{input_format}"));
        convert_clm(text.to_str().unwrap(), input.to_str().unwrap(), 64).unwrap();
        for requested in ["auto", "clm", "clm.gz", "clmb"] {
            let output_format = if requested == "auto" { input_format } else { requested };
            let output = directory.path().join(format!("split-{input_format}-{requested}"));
            let status = std::process::Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
                .args(["splitclm", input.to_str().unwrap(), clusters.to_str().unwrap(),
                    "-o", output.to_str().unwrap()])
                .args(if requested == "auto" { vec![] } else { vec!["--output-format", requested] })
                .status().unwrap();
            assert!(status.success(), "{input_format} -> {requested}");
            for (group, expected) in [
                ("group1", "A- B+\t2\t2 4294967302\n"),
                ("group2", "B- C-\t2\t9 11\n"),
                ("empty", ""),
            ] {
                let result = output.join(format!("{group}.{output_format}"));
                if output_format == "clm.gz" {
                    assert!(fs::read(&result).unwrap().starts_with(&[0x1f, 0x8b]));
                }
                let decoded = directory.path().join("decoded.clm");
                convert_clm(result.to_str().unwrap(), decoded.to_str().unwrap(), 64).unwrap();
                assert_eq!(fs::read_to_string(decoded).unwrap(), expected,
                    "{input_format} -> {requested}, {group}");
            }
        }
    }
}

#[test]
fn splitclm_rejects_invalid_output_format() {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .args(["splitclm", "missing.clm", "missing.clusters", "--output-format", "bam"])
        .output().unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("invalid value"));
}
