use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::time::Instant;
use anyhow::Result as anyResult;
use cphasing::core::common_reader;
use cphasing::methalign::{parse_bedgraph, MethConfig, MethRefiner};

// Frozen original loader: used only for exact compatibility checks.
fn legacy_parse_bedgraph(
    bedgraph: &String,
    cov_cutoff: f64,
) -> anyResult<HashMap<String, HashSet<i64>>> {
    let mut fh = common_reader(bedgraph);
    let mut map: HashMap<String, HashSet<i64>> = HashMap::new();
    let mut line = String::new();

    loop {
        line.clear();
        let n = fh.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.is_empty() {
            continue;
        }

        if line.as_bytes().first() == Some(&b'#') {
            continue;
        }

        let mut it = line.split_ascii_whitespace();
        let ctg = match it.next() {
            Some(x) => x,
            None => continue,
        };
        let start_str = match it.next() {
            Some(x) => x,
            None => continue,
        };
        let _end = it.next();
        let cov_str = match it.next() {
            Some(x) => x,
            None => continue,
        };

        let cov: f64 = match cov_str.parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if cov < cov_cutoff {
            continue;
        }
        let s: i64 = match start_str.parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        map.entry(ctg.to_owned()).or_default().insert(s);
    }

    log::info!(
        "Load {} methylation sites from bedgraph {} after filter by ref_prob_cutoff {}",
        map.values().map(|s| s.len()).sum::<usize>(),
        bedgraph,
        cov_cutoff
    );
    Ok(map)
}


#[test]
fn repeated_contigs_duplicates_and_filtering_preserve_sets() {
    let dir = tempfile::Builder::new().prefix(".meth-load-test-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    let path = dir.path().join("sites.bg");
    std::fs::write(&path, "# comment\n\nchrA\t10\t11\t100\nchrA 10 99 70\r\nchrA 20 21 49.9\nchrB 5 6 50\nchrB bad 7 100\nchrA 11 12 100\nchrC 1 2 bad\nchrB 6 8 90 extra\nchrA -1 0 100\nchrD 7 8 NaN\nchrD 8 9 inf\nchrA 12 13\n").unwrap();
    let name = path.to_str().unwrap().to_owned();
    let actual = parse_bedgraph(&name, 50.0).unwrap();
    let expected = HashMap::from([
        ("chrA".to_owned(), HashSet::from([-1, 10, 11])),
        ("chrB".to_owned(), HashSet::from([5, 6])),
        ("chrD".to_owned(), HashSet::from([7, 8])),
    ]);
    assert_eq!(actual, expected);
    assert_eq!(actual, legacy_parse_bedgraph(&name, 50.0).unwrap());
    std::fs::write(&path, "# no retained sites\nchrA 1 2 0\n").unwrap();
    assert!(parse_bedgraph(&name, 50.0).unwrap().is_empty());
}

#[test]
#[ignore = "explicit full bedGraph benchmark; set CPHASING_METH_BENCH_BED"]
fn benchmark_real_bedgraph() {
    let _ = env_logger::Builder::new().filter_level(log::LevelFilter::Info).try_init();
    let bed = std::env::var("CPHASING_METH_BENCH_BED").expect("named benchmark input required");
    if let Ok(reference) = std::env::var("CPHASING_METH_BENCH_REF") {
        let mut reader = needletail::parse_fastx_file(reference).unwrap();
        let mut sequences = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            sequences.push((String::from_utf8_lossy(record.id()).into_owned(), record.seq().to_vec()));
        }
        let start = Instant::now();
        let refiner = MethRefiner::new(MethConfig { bed, match_score: 0, ref_penalty: 2,
            read_penalty: 2, ref_prob_cutoff: 50.0, prob_cutoff: 128, designate_mapq: 2,
            cpg: true }, &sequences).unwrap();
        println!("METH_BENCH reference_case_seconds={:.9}", start.elapsed().as_secs_f64());
        std::hint::black_box(&refiner);
    } else {
        let start = Instant::now();
        let sites = parse_bedgraph(&bed, 50.0).unwrap();
        println!("METH_BENCH load_seconds={:.9} contigs={} sites={}",
            start.elapsed().as_secs_f64(), sites.len(), sites.values().map(HashSet::len).sum::<usize>());
        if std::env::var("CPHASING_METH_COMPARE_LEGACY").as_deref() == Ok("1") {
            let legacy = legacy_parse_bedgraph(&bed, 50.0).unwrap();
            assert!(sites == legacy, "all contig names and complete position sets must agree");
            println!("METH_BENCH exact_legacy_comparison=passed");
        }
        std::hint::black_box(&sites);
    }
}
