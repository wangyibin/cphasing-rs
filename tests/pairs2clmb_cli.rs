use std::fs;
use std::process::Command;

use cphasing::clm::ClmbReader;

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
