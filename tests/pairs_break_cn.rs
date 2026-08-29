use std::fs;
use std::process::Command;

#[test]
fn pairs_break_moves_source_cn_to_each_fragment() {
    let directory = tempfile::tempdir().unwrap();
    let pairs = directory.path().join("input.pairs");
    let pqs = directory.path().join("input.pairs.pqs");
    fs::write(
        &pairs,
        "## pairs format 1.0\n\
#chromsize: A 1000\n\
#chromsize: B 1000\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t100\tB\t200\t+\t-\t60\n\
r2\tA\t700\tB\t800\t+\t-\t60\n",
    )
    .unwrap();
    let convert = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(&pairs)
        .arg("--output")
        .arg(&pqs)
        .output()
        .unwrap();
    assert!(convert.status.success());
    fs::write(pqs.join("cn.info"), "A\t3\nB\t2\n").unwrap();
    let bed = directory.path().join("break.bed");
    fs::write(&bed, "A\t0\t499\tA_1\nA\t500\t999\tA_2\n").unwrap();
    let output = directory.path().join("broken.pairs.pqs");

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs-break")
        .arg(&pqs)
        .arg(&bed)
        .arg("--output")
        .arg(&output)
        .arg("--threads")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        fs::read_to_string(output.join("cn.info")).unwrap(),
        "A_1\t3\nA_2\t3\nB\t2\n"
    );
}
