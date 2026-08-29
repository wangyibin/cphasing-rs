use std::fs;
use std::process::Command;

use polars::prelude::*;

#[test]
fn pairs_chr2ctg_remaps_cn_info_to_bed_contigs() {
    let directory = tempfile::tempdir().unwrap();
    let pairs = directory.path().join("input.pairs");
    let input = directory.path().join("input.pairs.pqs");
    fs::write(
        &pairs,
        "## pairs format 1.0\n\
#shape: upper triangle\n\
#chromsize: A 1000\n\
#chromsize: B 1000\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t100\tB\t200\t+\t-\t60\n\
r2\tA\t700\tB\t800\t+\t-\t60\n",
    )
    .unwrap();
    let converted = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(&pairs)
        .arg("--output")
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        converted.status.success(),
        "{}",
        String::from_utf8_lossy(&converted.stderr)
    );
    fs::write(input.join("cn.info"), "A\t2\nB\t3\n").unwrap();

    let bed = directory.path().join("mapping.bed");
    fs::write(&bed, "A\t0\t500\tA1\nA\t500\t1000\tA2\nB\t0\t1000\tB1\n").unwrap();
    let output = directory.path().join("output.pairs.pqs");
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs-chr2ctg")
        .arg("--input")
        .arg(&input)
        .arg("--bed")
        .arg(&bed)
        .arg("--output")
        .arg(&output)
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
        fs::read_to_string(output.join("_contigsizes")).unwrap(),
        "A1\t500\nA2\t500\nB1\t1000\n"
    );
    assert_eq!(
        fs::read_to_string(output.join("cn.info")).unwrap(),
        "A1\t2\nA2\t2\nB1\t3\n"
    );

    let frame = ParquetReader::new(fs::File::open(output.join("q0/0.parquet")).unwrap())
        .finish()
        .unwrap();
    let chrom1 = frame
        .column("chrom1")
        .unwrap()
        .cast(&DataType::String)
        .unwrap();
    assert_eq!(
        chrom1
            .str()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        vec!["A1", "A2"]
    );
}

#[test]
fn pairs_chr2ctg_rejects_text_pairs_and_documents_the_limitation() {
    let directory = tempfile::tempdir().unwrap();
    let bed = directory.path().join("mapping.bed");
    fs::write(&bed, "A\t0\t1000\tA1\n").unwrap();

    for suffix in ["pairs", "pairs.gz"] {
        let input = directory.path().join(format!("input.{suffix}"));
        fs::write(&input, "not a pairs PQS directory\n").unwrap();
        let output = directory.path().join(format!("{suffix}.output.pairs.pqs"));
        let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .arg("pairs-chr2ctg")
            .arg("--input")
            .arg(&input)
            .arg("--bed")
            .arg(&bed)
            .arg("--output")
            .arg(&output)
            .output()
            .unwrap();
        assert!(!result.status.success());
        let stderr = String::from_utf8_lossy(&result.stderr);
        assert!(
            stderr.contains("only supports a pairs PQS directory"),
            "{stderr}"
        );
        assert!(stderr.contains(".pairs.gz"), "{stderr}");
        assert!(!output.exists());
    }

    let help = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs-chr2ctg")
        .arg("--help")
        .output()
        .unwrap();
    assert!(help.status.success());
    let stdout = String::from_utf8_lossy(&help.stdout);
    assert!(stdout.contains(".pairs/.pairs.gz"), "{stdout}");
}
