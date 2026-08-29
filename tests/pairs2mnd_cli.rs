use std::collections::HashSet;
use std::fs;
use std::path::Path;
use std::process::Command;

fn run_pairs2mnd(input: &Path, output: &Path, ignore_cn: bool) {
    let mut command = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"));
    command
        .arg("pairs2mnd")
        .arg(input)
        .arg("--output")
        .arg(output);
    if ignore_cn {
        command.arg("--ignore-cn");
    }
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn pairs2mnd_uses_cn_info_by_default_at_output() {
    let directory = tempfile::tempdir().unwrap();
    let pairs = directory.path().join("input.pairs");
    let pqs = directory.path().join("input.pairs.pqs");
    let mut text = String::from(
        "## pairs format 1.0\n\
#shape: upper triangle\n\
#chromsize: A 1000\n\
#chromsize: B 1000\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n",
    );
    for index in 0..40 {
        text.push_str(&format!("cross{index}\tA\t10\tB\t20\t+\t-\t60\n"));
    }
    for index in 0..10 {
        text.push_str(&format!("cis{index}\tA\t30\tA\t40\t+\t-\t60\n"));
    }
    fs::write(&pairs, text).unwrap();
    let convert = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(&pairs)
        .arg("--chunksize")
        .arg("7")
        .arg("--output")
        .arg(&pqs)
        .output()
        .unwrap();
    assert!(convert.status.success());
    fs::write(pqs.join("cn.info"), "A\t2\nB\t2\n").unwrap();

    let output1 = directory.path().join("cn1.mnd");
    let output2 = directory.path().join("cn2.mnd");
    run_pairs2mnd(&pqs, &output1, false);
    run_pairs2mnd(&pqs, &output2, false);
    let cn_text = fs::read_to_string(&output1).unwrap();
    assert_eq!(cn_text, fs::read_to_string(&output2).unwrap());
    assert_eq!(cn_text.lines().count(), 50);
    assert!(cn_text.contains("A_d2") || cn_text.contains("B_d2"));

    let mut cross_combinations = HashSet::new();
    for line in cn_text.lines() {
        let fields = line.split_whitespace().collect::<Vec<_>>();
        assert_eq!(fields.len(), 16);
        if fields[1].starts_with('A') && fields[5].starts_with('A') {
            assert_eq!(fields[1], fields[5]);
        } else {
            cross_combinations.insert((fields[1].to_string(), fields[5].to_string()));
        }
    }
    assert_eq!(cross_combinations.len(), 4);

    let raw_output = directory.path().join("raw.mnd");
    run_pairs2mnd(&pqs, &raw_output, true);
    let raw_text = fs::read_to_string(raw_output).unwrap();
    assert_eq!(raw_text.lines().count(), 50);
    assert!(!raw_text.contains("_d2"));
}

#[test]
fn pairs2mnd_treats_missing_cn_info_as_cn_one() {
    let directory = tempfile::tempdir().unwrap();
    let pairs = directory.path().join("input.pairs");
    let pqs = directory.path().join("input.pairs.pqs");
    fs::write(
        &pairs,
        "## pairs format 1.0\n\
#chromsize: A 100\n\
#chromsize: B 100\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t10\tB\t20\t+\t-\t60\n",
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

    let output = directory.path().join("output.mnd");
    run_pairs2mnd(&pqs, &output, false);
    let text = fs::read_to_string(output).unwrap();
    assert!(text.contains(" A "));
    assert!(text.contains(" B "));
    assert!(!text.contains("_d2"));
}
