use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::Path;
use std::process::Command;

use polars::prelude::*;

fn categorical(name: &str, values: &[&str]) -> Series {
    Series::new(name.into(), values)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
        .unwrap()
}

fn run(command: &mut Command) {
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

fn pairs_rows(path: &Path) -> BTreeMap<String, (String, u32, String, u32, u8)> {
    let frame = ParquetReader::new(File::open(path).unwrap())
        .finish()
        .unwrap();
    let read_idx = frame.column("read_idx").unwrap().str().unwrap();
    let chrom1 = frame
        .column("chrom1")
        .unwrap()
        .cast(&DataType::String)
        .unwrap();
    let chrom1 = chrom1.str().unwrap();
    let chrom2 = frame
        .column("chrom2")
        .unwrap()
        .cast(&DataType::String)
        .unwrap();
    let chrom2 = chrom2.str().unwrap();
    let pos1 = frame.column("pos1").unwrap().u32().unwrap();
    let pos2 = frame.column("pos2").unwrap().u32().unwrap();
    let mapq = frame.column("mapq").unwrap().u8().unwrap();

    (0..frame.height())
        .map(|row| {
            (
                read_idx.get(row).unwrap().to_string(),
                (
                    chrom1.get(row).unwrap().to_string(),
                    pos1.get(row).unwrap(),
                    chrom2.get(row).unwrap().to_string(),
                    pos2.get(row).unwrap(),
                    mapq.get(row).unwrap(),
                ),
            )
        })
        .collect()
}

#[test]
fn pairs_break_remaps_typed_coordinate_columns() {
    let directory = tempfile::Builder::new()
        .prefix("pairs-break-typed-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let pairs = directory.path().join("input.pairs");
    let input = directory.path().join("input.pairs.pqs");
    let output = directory.path().join("output.pairs.pqs");
    let bed = directory.path().join("break.bed");
    fs::write(
        &pairs,
        "## pairs format 1.0\n\
#chromsize: A 1000\n\
#chromsize: B 1000\n\
#chromsize: C 1000\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t100\tB\t200\t+\t-\t60\n\
r2\tA\t700\tB\t800\t-\t+\t0\n\
r3\tB\t300\tA\t600\t+\t-\t60\n\
r4\tB\t400\tC\t500\t-\t+\t60\n\
r5\tA\t499\tB\t250\t+\t-\t60\n",
    )
    .unwrap();
    fs::write(&bed, "A\t0\t499\tA_1\nA\t500\t999\tA_2\n").unwrap();

    run(Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(&pairs)
        .arg("--output")
        .arg(&input));

    // Keep q1 intentionally different from a fresh mapq filter over q0. This
    // verifies that pairs-break transforms the existing q1 shards directly.
    let q1_path = input.join("q1/0.parquet");
    let q1 = ParquetReader::new(File::open(&q1_path).unwrap())
        .finish()
        .unwrap();
    let mut q1 = q1
        .lazy()
        .filter(col("read_idx").eq(lit("r4")))
        .collect()
        .unwrap();
    ParquetWriter::new(File::create(&q1_path).unwrap())
        .finish(&mut q1)
        .unwrap();

    run(Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs-break")
        .arg(&input)
        .arg(&bed)
        .arg("--output")
        .arg(&output)
        .arg("--threads")
        .arg("2"));

    assert_eq!(
        pairs_rows(&output.join("q0/0.parquet")),
        BTreeMap::from([
            ("r1".into(), ("A_1".into(), 101, "B".into(), 200, 60)),
            ("r2".into(), ("A_2".into(), 201, "B".into(), 800, 0)),
            ("r3".into(), ("A_2".into(), 101, "B".into(), 300, 60)),
            ("r4".into(), ("B".into(), 400, "C".into(), 500, 60)),
            ("r5".into(), ("A".into(), 499, "B".into(), 250, 60)),
        ])
    );
    assert_eq!(
        pairs_rows(&output.join("q1/0.parquet")),
        BTreeMap::from([("r4".into(), ("B".into(), 400, "C".into(), 500, 60))])
    );
}

#[test]
fn porec_break_preserves_u32_positions() {
    let directory = tempfile::Builder::new()
        .prefix("porec-break-u32-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let input = directory.path().join("input.concat.pqs");
    let output = directory.path().join("output.concat.pqs");
    fs::create_dir_all(input.join("q0")).unwrap();
    fs::create_dir_all(input.join("q1")).unwrap();
    fs::write(input.join("_metadata"), "{'format': 'concat'}\n").unwrap();
    fs::write(input.join("_contigsizes"), "A\t1000\nB\t1000\n").unwrap();
    fs::write(
        input.join("_metadata_counts"),
        "q0_records\t4\nq1_records\t3\nq0_concats\t4\nq1_concats\t3\n",
    )
    .unwrap();
    fs::write(input.join("_readme"), "test concat PQS\n").unwrap();
    let mut frame = DataFrame::new(vec![
        Series::new("read_idx".into(), [1u64, 2, 3, 4]).into(),
        Series::new("read_length".into(), [1000u32; 4]).into(),
        Series::new("read_start".into(), [0u32; 4]).into(),
        Series::new("read_end".into(), [100u32; 4]).into(),
        categorical("strand", &["+", "-", "+", "+"]).into(),
        categorical("chrom", &["A", "A", "B", "A"]).into(),
        Series::new("start".into(), [100u32, 600, 300, 450]).into(),
        Series::new("end".into(), [200u32, 700, 400, 550]).into(),
        Series::new("mapping_quality".into(), [60u8, 0, 60, 60]).into(),
        Series::new("identity".into(), [0.9f32, 0.8, 0.7, 0.95]).into(),
        categorical("filter_reason", &["pass", "low_mapq", "pass", "pass"]).into(),
    ])
    .unwrap();
    ParquetWriter::new(File::create(input.join("q0/0.parquet")).unwrap())
        .finish(&mut frame)
        .unwrap();

    let bed = directory.path().join("break.bed");
    fs::write(&bed, "A\t0\t499\tA_1\nA\t500\t999\tA_2\n").unwrap();
    run(Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("porec-break")
        .arg(&input)
        .arg(&bed)
        .arg("--output")
        .arg(&output)
        .arg("--threads")
        .arg("2"));

    let q0 = ParquetReader::new(File::open(output.join("q0/0.parquet")).unwrap())
        .finish()
        .unwrap();
    assert_eq!(q0.column("start").unwrap().dtype(), &DataType::UInt32);
    assert_eq!(
        q0.column("start")
            .unwrap()
            .u32()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        [101, 101, 300, 450]
    );
    assert_eq!(
        q0.column("end")
            .unwrap()
            .u32()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        [201, 201, 400, 550]
    );
    let chrom = q0.column("chrom").unwrap().cast(&DataType::String).unwrap();
    assert_eq!(
        chrom.str().unwrap().into_no_null_iter().collect::<Vec<_>>(),
        ["A_1", "A_2", "B", "A"]
    );
    let q1 = ParquetReader::new(File::open(output.join("q1/0.parquet")).unwrap())
        .finish()
        .unwrap();
    assert_eq!(q1.height(), 3);
}
