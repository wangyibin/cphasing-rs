use std::fs;

use cphasing::clm::{
    ClmbReader, ClmbRecord, ClmbWriter, clm_output_prefix, convert_clm, encode_endpoint,
    is_clmb_file,
};

#[test]
fn clmb_round_trip_uses_adaptive_distance_width() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("adaptive.clmb");
    let contigs = vec!["ctgA".to_string(), "ctgB".to_string()];
    let expected = vec![
        ClmbRecord {
            endpoint1: encode_endpoint(0, 0).unwrap(),
            endpoint2: encode_endpoint(1, 1).unwrap(),
            distances: vec![2, 10, 100],
        },
        ClmbRecord {
            endpoint1: encode_endpoint(0, 1).unwrap(),
            endpoint2: encode_endpoint(1, 0).unwrap(),
            distances: vec![7, u32::MAX as u64 + 11],
        },
    ];
    let mut writer = ClmbWriter::create(&path, &contigs, 20, Some(2), Some(5)).unwrap();
    for record in &expected {
        writer
            .write_record(record.endpoint1, record.endpoint2, &record.distances)
            .unwrap();
    }
    writer.finish().unwrap();
    assert!(is_clmb_file(&path).unwrap());

    let mut reader = ClmbReader::open(&path).unwrap();
    assert_eq!(reader.header.contigs, contigs);
    let mut observed = Vec::new();
    while let Some(block) = reader.next_block().unwrap() {
        observed.extend(block);
    }
    assert_eq!(observed, expected);
}

#[test]
fn clmb_rejects_truncated_block() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("truncated.clmb");
    let contigs = vec!["ctgA".to_string(), "ctgB".to_string()];
    let mut writer = ClmbWriter::create(&path, &contigs, 1024, None, None).unwrap();
    writer.write_record(0, 2, &[2, 3, 4]).unwrap();
    writer.finish().unwrap();
    let length = fs::metadata(&path).unwrap().len();
    fs::OpenOptions::new()
        .write(true)
        .open(&path)
        .unwrap()
        .set_len(length - 1)
        .unwrap();
    let mut reader = ClmbReader::open(&path).unwrap();
    assert!(reader.next_block().is_err());
}

#[test]
fn text_clm_round_trips_through_clmb_without_narrowing_distances() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("input.clm");
    let binary = directory.path().join("output.clmb");
    let output = directory.path().join("roundtrip.clm");
    let expected = format!(
        "ctgA+ ctgB-\t3\t2 10 {}\nctgA- ctgB+\t1\t7\n",
        u32::MAX as u64 + 99
    );
    fs::write(&input, &expected).unwrap();
    convert_clm(input.to_str().unwrap(), binary.to_str().unwrap(), 32).unwrap();
    convert_clm(binary.to_str().unwrap(), output.to_str().unwrap(), 32).unwrap();
    assert_eq!(fs::read_to_string(output).unwrap(), expected);
}

#[test]
fn clmb_auxiliary_outputs_drop_the_binary_suffix() {
    assert_eq!(clm_output_prefix("input.clmb"), "input");
    assert_eq!(clm_output_prefix("input.clm.gz"), "input");
    assert_eq!(clm_output_prefix("input.clm"), "input");
    assert_eq!(clm_output_prefix("input.gz"), "input");
}
