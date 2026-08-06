use std::fs;

use cphasing::clm::{Clm, ClmbReader, ClmbRecord, ClmbWriter, encode_endpoint};
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
