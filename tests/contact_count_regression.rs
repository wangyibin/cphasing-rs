use cphasing::{core::BaseTable, porec::PoreCTable, pqs::PQS};
use polars::prelude::*;
use std::fs::{self, File};

fn temp() -> tempfile::TempDir {
    tempfile::Builder::new()
        .prefix("contact-regression-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap()
}

#[test]
fn porec_intersection_keeps_tail_and_honors_invert() {
    let dir = temp();
    let input = dir.path().join("input.porec");
    let bed = dir.path().join("regions.bed");
    fs::write(&bed, "A\t0\t100\n").unwrap();
    for threads in [1, 3] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for count in [0, 1, 9_999, 10_000, 10_001] {
            let rows: Vec<_> = (0..count)
                .map(|i| {
                    format!(
                        "{i}\t0\t0\t0\t0\t{}\t10\t20",
                        if i % 2 == 0 { "A" } else { "B" }
                    )
                })
                .collect();
            fs::write(
                &input,
                if rows.is_empty() {
                    String::new()
                } else {
                    rows.join("\n") + "\n"
                },
            )
            .unwrap();
            for invert in [false, true] {
                let output = dir.path().join("out.porec").to_string_lossy().into_owned();
                pool.install(|| {
                    PoreCTable::new(&input.to_string_lossy().into_owned()).intersect_multi_threads(
                        &bed.to_string_lossy().into_owned(),
                        invert,
                        &output,
                    )
                });
                let expected: Vec<_> = rows
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| (i % 2 == 0) ^ invert)
                    .map(|(_, s)| s.as_str())
                    .collect();
                let actual = fs::read_to_string(&output).unwrap();
                assert_eq!(
                    actual.lines().collect::<Vec<_>>(),
                    expected,
                    "rows={count}, threads={threads}, invert={invert}"
                );
            }
        }
    }
}

fn write_shard(path: &std::path::Path, rows: &[(&str, &str, u32, u32, u8)]) {
    let categorical = |name: &str, values: Vec<&str>| {
        Series::new(name.into(), values)
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
            .unwrap()
            .into()
    };
    let mut df = DataFrame::new(vec![
        categorical("chrom1", rows.iter().map(|r| r.0).collect()),
        categorical("chrom2", rows.iter().map(|r| r.1).collect()),
        Series::new("pos1".into(), rows.iter().map(|r| r.2).collect::<Vec<_>>()).into(),
        Series::new("pos2".into(), rows.iter().map(|r| r.3).collect::<Vec<_>>()).into(),
        Series::new("mapq".into(), rows.iter().map(|r| r.4).collect::<Vec<_>>()).into(),
    ])
    .unwrap();
    ParquetWriter::new(File::create(path).unwrap())
        .finish(&mut df)
        .unwrap();
}

#[test]
fn pqs_counts_match_reference_across_shards_filters_and_threads() {
    let dir = temp();
    fs::create_dir(dir.path().join("q0")).unwrap();
    fs::create_dir(dir.path().join("q1")).unwrap();
    fs::write(dir.path().join("_contigsizes"), "A\t100\nB\t3\n").unwrap();
    let q0 = vec![
        ("A", "B", 0, 0, 0),
        ("A", "B", 49, 1, 1),
        ("A", "B", 50, 2, 60),
        ("B", "A", 3, 100, 60),
        ("A", "A", 100, 0, 20),
        ("A", "B", 50, 2, 60),
    ];
    // q1 is deliberately not equivalent to filtering q0.
    let q1 = vec![q0[2], q0[1], q0[5], q0[3]];
    for (name, rows) in [("q0", &q0), ("q1", &q1)] {
        for (index, chunk) in rows.chunks(2).enumerate() {
            write_shard(&dir.path().join(format!("{name}/{index}.parquet")), chunk);
        }
        write_shard(&dir.path().join(format!("{name}/empty.parquet")), &[]);
    }
    let pqs = PQS::new(&dir.path().to_string_lossy().into_owned());
    for threads in [1, 3] {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        for quality in [0, 1, 20, 61] {
            for split in [None, Some(1), Some(2), Some(5)] {
                for threshold in [1, 2, 10] {
                    let mut expected = std::collections::BTreeMap::new();
                    for &(a, b, p1, p2, q) in if quality == 0 { &q0 } else { &q1 } {
                        if quality > 1 && q < quality {
                            continue;
                        }
                        let label = |name: &str, pos: u32| match split {
                            None => name.to_owned(),
                            Some(n) => {
                                let size = if name == "A" { 100u32 } else { 3 };
                                format!("{name}_{}", (pos / (size / n).max(1)).min(n - 1))
                            }
                        };
                        *expected.entry((label(a, p1), label(b, p2))).or_insert(0u64) += 1;
                    }
                    let mut expected: Vec<_> = expected
                        .into_iter()
                        .filter(|(_, c)| *c >= threshold as u64)
                        .map(|((a, b), c)| format!("{a}\t{b}\t{c}"))
                        .collect();
                    expected.sort();
                    let output = dir.path().join("counts.tsv").to_string_lossy().into_owned();
                    pool.install(|| match split {
                        Some(n) => pqs.to_split_contacts(threshold, n, quality, &output),
                        None => pqs.to_contacts(threshold, quality, &output),
                    })
                    .unwrap();
                    let actual = fs::read_to_string(output).unwrap();
                    let mut actual: Vec<_> = actual.lines().map(str::to_owned).collect();
                    actual.sort();
                    assert_eq!(
                        actual, expected,
                        "threads={threads} quality={quality} split={split:?}"
                    );
                }
            }
        }
    }
    for threads in ["1", "3"] {
        let output = dir.path().join("cli-counts.tsv");
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .arg("pairs2contacts")
            .arg(dir.path())
            .args(["--threads", threads, "--min-quality", "1"])
            .arg("--output")
            .arg(&output)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(fs::read_to_string(output).unwrap(), "A\tB\t3\nB\tA\t1\n");
    }
    assert!(
        pqs.to_split_contacts(
            1,
            0,
            0,
            &dir.path().join("invalid").to_string_lossy().into_owned()
        )
        .is_err()
    );
    fs::write(dir.path().join("q0/bad.parquet"), "invalid parquet").unwrap();
    assert!(
        pqs.to_contacts(
            1,
            0,
            &dir.path().join("bad-output").to_string_lossy().into_owned()
        )
        .is_err()
    );
}
