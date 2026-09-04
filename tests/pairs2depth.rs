use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use std::process::Command;

fn read_depth(path: &Path) -> BTreeMap<(String, u32, u32), u32> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            (
                (
                    fields[0].to_string(),
                    fields[1].parse().unwrap(),
                    fields[2].parse().unwrap(),
                ),
                fields[3].parse().unwrap(),
            )
        })
        .collect()
}

fn run_pairs2depth(input: &Path, output: &Path, min_quality: u8, threads: usize) {
    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2depth")
        .arg(input)
        .arg("--binsize")
        .arg("250")
        .arg("--min-quality")
        .arg(min_quality.to_string())
        .arg("--threads")
        .arg(threads.to_string())
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

fn duplicate_parquet_shards(pqs: &Path, quality: &str) {
    let directory = pqs.join(quality);
    let shards = fs::read_dir(&directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            matches!(
                path.extension().and_then(|value| value.to_str()),
                Some("pq" | "parquet")
            )
        })
        .collect::<Vec<_>>();
    assert!(!shards.is_empty());

    for (index, shard) in shards.iter().enumerate() {
        let extension = shard.extension().unwrap().to_str().unwrap();
        fs::copy(
            shard,
            directory.join(format!("duplicate-{index}.{extension}")),
        )
        .unwrap();
    }
}

#[test]
fn pqs_pairs2depth_counts_endpoints_without_storing_positions() {
    let directory = tempfile::tempdir().unwrap();
    let pairs = directory.path().join("input.pairs");
    let pqs = directory.path().join("input.pairs.pqs");
    fs::write(
        &pairs,
        "## pairs format 1.0\n\
#shape: upper triangle\n\
#chromsize: A 1000\n\
#chromsize: B 750\n\
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 mapq\n\
r1\tA\t100\tB\t200\t+\t-\t0\n\
r2\tA\t249\tA\t250\t+\t-\t1\n\
r3\tA\t999\tB\t750\t+\t-\t10\n\
r4\tA\t500\tB\t501\t+\t-\t5\n",
    )
    .unwrap();
    let converted = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("pairs2pqs")
        .arg(&pairs)
        .arg("--output")
        .arg(&pqs)
        .output()
        .unwrap();
    assert!(
        converted.status.success(),
        "{}",
        String::from_utf8_lossy(&converted.stderr)
    );

    let q0_output = directory.path().join("q0.depth");
    run_pairs2depth(&pqs, &q0_output, 0, 2);
    let q0_depth = read_depth(&q0_output);
    assert_eq!(
        q0_depth,
        BTreeMap::from([
            (("A".to_string(), 0, 250), 2),
            (("A".to_string(), 250, 500), 1),
            (("A".to_string(), 500, 750), 1),
            (("A".to_string(), 750, 1000), 1),
            (("B".to_string(), 0, 250), 1),
            (("B".to_string(), 250, 500), 0),
            (("B".to_string(), 500, 750), 1),
        ])
    );

    let q5_output = directory.path().join("q5.depth");
    run_pairs2depth(&pqs, &q5_output, 5, 2);
    assert_eq!(
        read_depth(&q5_output),
        BTreeMap::from([
            (("A".to_string(), 0, 250), 0),
            (("A".to_string(), 250, 500), 0),
            (("A".to_string(), 500, 750), 1),
            (("A".to_string(), 750, 1000), 1),
            (("B".to_string(), 0, 250), 0),
            (("B".to_string(), 250, 500), 0),
            (("B".to_string(), 500, 750), 1),
        ])
    );

    duplicate_parquet_shards(&pqs, "q0");
    let single_thread_output = directory.path().join("q0.single-thread.depth");
    let parallel_output = directory.path().join("q0.parallel.depth");
    run_pairs2depth(&pqs, &single_thread_output, 0, 1);
    run_pairs2depth(&pqs, &parallel_output, 0, 2);

    let doubled_depth = q0_depth
        .into_iter()
        .map(|(bin, count)| (bin, count * 2))
        .collect::<BTreeMap<_, _>>();
    assert_eq!(read_depth(&single_thread_output), doubled_depth);
    assert_eq!(read_depth(&parallel_output), doubled_depth);
}
