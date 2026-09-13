use cphasing::{core::BaseTable, pqs::PQS};
use polars::prelude::*;
use polars::enable_string_cache;
use std::fs::{self, File};
use std::path::Path;

fn shard(path: &Path, n: usize, offset: u32) {
    let mut frame = DataFrame::new(vec![
        Series::new(
            "read_idx".into(),
            (0..n).map(|i| format!("read-{i}")).collect::<Vec<_>>(),
        )
        .into(),
        Series::new("chrom1".into(), vec!["A"; n])
            .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
            .unwrap()
            .into(),
        Series::new(
            "chrom2".into(),
            (0..n)
                .map(|i| if i % 2 == 0 { "A" } else { "B" })
                .collect::<Vec<_>>(),
        )
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
        .unwrap()
        .into(),
        Series::new(
            "pos1".into(),
            (0..n).map(|i| offset + i as u32).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "pos2".into(),
            (0..n).map(|i| offset + i as u32 + 1).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "mapq".into(),
            (0..n).map(|i| [0u8, 1, 20, 60][i % 4]).collect::<Vec<_>>(),
        )
        .into(),
    ])
    .unwrap();
    ParquetWriter::new(File::create(path).unwrap())
        .finish(&mut frame)
        .unwrap();
}

#[test]
fn mnd_stream_preserves_order_filters_cn_and_empty_shards() {
    enable_string_cache();
    let dir = tempfile::tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    for q in ["q0", "q1"] {
        fs::create_dir(dir.path().join(q)).unwrap();
        for i in 0..11 {
            shard(
                &dir.path().join(q).join(format!("{i}.parquet")),
                if i == 4 { 0 } else { 8 },
                if q == "q0" { i * 10 } else { 1000 + i * 10 },
            );
        }
    }
    let pqs = PQS::new(&dir.path().to_string_lossy().into_owned());
    for quality in [0, 1, 20, 61] {
        let q = if quality == 0 { "q0" } else { "q1" };
        let mut expected = String::new();
        for entry in walkdir::WalkDir::new(dir.path().join(q))
            .into_iter()
            .filter_map(Result::ok)
            .filter(|e| e.file_type().is_file())
        {
            let i: u32 = entry
                .path()
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .parse()
                .unwrap();
            if i == 4 {
                continue;
            }
            let offset = if q == "q0" { i * 10 } else { 1000 + i * 10 };
            for row in 0..8 {
                let mapq = [0u8, 1, 20, 60][row % 4];
                if quality > 1 && mapq < quality {
                    continue;
                }
                let chrom2 = if row % 2 == 0 { "A" } else { "B" };
                expected.push_str(&format!(
                    "0 A {} 0 0 {chrom2} {} 1 {mapq} - - {mapq} - - - -\n",
                    offset + row as u32,
                    offset + row as u32 + 1
                ));
            }
        }
        for threads in [1, 4] {
            let output = dir.path().join("out.mnd").to_string_lossy().into_owned();
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| pqs.to_mnd(quality, &output, false))
                .unwrap();
            assert_eq!(fs::read_to_string(&output).unwrap(), expected);
        }
    }
    fs::write(dir.path().join("cn.info"), "A\t2\nB\t3\n").unwrap();
    let mut previous = None;
    for threads in [1, 4] {
        let output = dir.path().join("cn.mnd").to_string_lossy().into_owned();
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| pqs.to_mnd(1, &output, true))
            .unwrap();
        let actual = fs::read_to_string(output).unwrap();
        for line in actual.lines() {
            let f: Vec<_> = line.split_whitespace().collect();
            assert!(f[1].starts_with('A'));
            if f[5].starts_with('A') {
                assert_eq!(f[1], f[5]);
            }
        }
        if let Some(previous) = previous {
            assert_eq!(actual, previous);
        }
        previous = Some(actual);
    }
    let empty = dir.path().join("empty");
    fs::create_dir_all(empty.join("q1")).unwrap();
    let output = dir.path().join("empty.mnd").to_string_lossy().into_owned();
    PQS::new(&empty.to_string_lossy().into_owned())
        .to_mnd(1, &output, false)
        .unwrap();
    assert_eq!(fs::metadata(output).unwrap().len(), 0);
}

#[test]
fn mnd_stream_returns_read_and_write_errors_without_hanging() {
    let dir = tempfile::tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    fs::create_dir(dir.path().join("q1")).unwrap();
    fs::write(dir.path().join("q1/bad.parquet"), "invalid parquet").unwrap();
    let pqs = PQS::new(&dir.path().to_string_lossy().into_owned());
    let error = pqs
        .to_mnd(
            1,
            &dir.path().join("out").to_string_lossy().into_owned(),
            false,
        )
        .unwrap_err();
    assert!(format!("{error:#}").contains("bad.parquet"));
    fs::remove_file(dir.path().join("q1/bad.parquet")).unwrap();
    for i in 0..12 {
        shard(&dir.path().join(format!("q1/{i}.parquet")), 1000, i * 1000);
    }
    #[cfg(target_os = "linux")]
    assert!(pqs.to_mnd(1, &"/dev/full".to_string(), false).is_err());
}

#[test]
#[ignore = "explicit bounded synthetic performance fixture"]
fn create_mnd_benchmark_fixture() {
    enable_string_cache();
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("pqs-mnd-stream-20260912/fixture.pqs");
    fs::create_dir_all(dir.join("q1")).unwrap();
    fs::create_dir_all(dir.join("q0")).unwrap();
    fs::write(dir.join("_contigsizes"), "A\t1000000\nB\t1000000\n").unwrap();
    fs::write(dir.join("cn.info"), "A\t2\nB\t3\n").unwrap();
    for i in 0..32 {
        shard(&dir.join(format!("q1/{i}.parquet")), 25000, i * 25000);
    }
}
