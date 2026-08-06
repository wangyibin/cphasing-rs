use std::fs::{self, File};
use std::path::Path;

use cphasing::core::BaseTable;
use cphasing::pqs::PQS;
use polars::prelude::*;

fn categorical(name: &str, values: &[&str]) -> Series {
    Series::new(name.into(), values)
        .cast(&DataType::Categorical(None, CategoricalOrdering::Physical))
        .unwrap()
}

fn write_shard(path: &Path, reads: &[&str], chroms: &[&str], mapq: &[u8]) {
    let mut df = DataFrame::new(vec![
        Series::new("read_idx".into(), reads).into(),
        categorical("chrom1", chroms).into(),
        Series::new("mapq".into(), mapq).into(),
        categorical("chrom2", chroms).into(),
    ])
    .unwrap();
    ParquetWriter::new(File::create(path).unwrap())
        .finish(&mut df)
        .unwrap();
}

fn read_chroms(path: &Path) -> Vec<String> {
    let df = LazyFrame::scan_parquet(path, ScanArgsParquet::default())
        .unwrap()
        .collect()
        .unwrap();
    df.column("chrom1")
        .unwrap()
        .categorical()
        .unwrap()
        .iter_str()
        .map(|value| value.unwrap().to_string())
        .collect()
}

#[test]
fn pairs_dup_keeps_quality_subsets_consistent_and_links_unaffected_shards() {
    let directory = tempfile::Builder::new()
        .prefix("pairs-dup-test-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR"))
        .unwrap();
    let input = directory.path().join("input.pqs");
    let output = directory.path().join("output.pqs");
    fs::create_dir_all(input.join("q0")).unwrap();
    fs::create_dir_all(input.join("q1")).unwrap();
    fs::write(input.join("_contigsizes"), "collapsed\t100\nother\t100\n").unwrap();
    fs::write(
        input.join("_metadata_counts"),
        "q0_records\t6\nq1_records\t5\n",
    )
    .unwrap();
    fs::write(input.join("_metadata"), "test metadata\n").unwrap();
    fs::write(input.join("_readme"), "test readme\n").unwrap();

    write_shard(
        &input.join("q0/0.parquet"),
        &["r0", "r1", "r2", "r3"],
        &["collapsed", "other", "collapsed", "collapsed"],
        &[0, 1, 1, 1],
    );
    write_shard(
        &input.join("q1/0.parquet"),
        &["r1", "r2", "r3"],
        &["other", "collapsed", "collapsed"],
        &[1, 1, 1],
    );
    write_shard(
        &input.join("q0/1.parquet"),
        &["u0", "u1"],
        &["other", "other"],
        &[0, 1],
    );
    write_shard(&input.join("q1/1.parquet"), &["u1"], &["other"], &[1]);
    let collapsed = directory.path().join("collapsed.tsv");
    fs::write(&collapsed, "collapsed\tcollapsed_d2\n").unwrap();

    PQS::new(&input.to_string_lossy().to_string())
        .dup(
            &collapsed.to_string_lossy().to_string(),
            123,
            &output.to_string_lossy().to_string(),
        )
        .unwrap();

    let q0 = read_chroms(&output.join("q0/0.parquet"));
    let q1 = read_chroms(&output.join("q1/0.parquet"));
    assert_eq!(q0[1..], q1);
    assert_eq!(q0[1], "other");
    assert_eq!(
        fs::read(input.join("q0/1.parquet")).unwrap(),
        fs::read(output.join("q0/1.parquet")).unwrap(),
    );

    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        assert_eq!(
            fs::metadata(input.join("q0/1.parquet")).unwrap().ino(),
            fs::metadata(output.join("q0/1.parquet")).unwrap().ino(),
        );
    }
}
