use std::fs;
use std::process::Command;

use cphasing::clm::{ClmbReader, ClmbWriter, encode_endpoint};

fn clmb_as_text(path: impl AsRef<std::path::Path>) -> String {
    let mut reader = ClmbReader::open(path).unwrap();
    let contigs = reader.header.contigs.clone();
    let mut output = String::new();
    while let Some(block) = reader.next_block().unwrap() {
        for record in block {
            let orientation1 = if record.orientation1() == 0 { '+' } else { '-' };
            let orientation2 = if record.orientation2() == 0 { '+' } else { '-' };
            let distances = record
                .distances
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(" ");
            output.push_str(&format!(
                "{}{} {}{}\t{}\t{}\n",
                contigs[record.contig1() as usize],
                orientation1,
                contigs[record.contig2() as usize],
                orientation2,
                record.distances.len(),
                distances
            ));
        }
    }
    output
}

#[test]
fn gfa_contacts_cli_aggregates_end_evidence() {
    let directory = tempfile::tempdir().unwrap();
    let segments = directory.path().join("segments.tsv");
    let contacts = directory.path().join("contacts.tsv");
    let links = directory.path().join("links.tsv");
    let output = directory.path().join("output.tsv");
    fs::write(&segments, "A\t100\nB\t200\nC\t100\n").unwrap();
    fs::write(&contacts, "A_1\tB_0\t2\nB_0\tA_1\t3\nA_1\tC_0\t10\n").unwrap();
    fs::write(&links, "A\tR\tB\tL\n").unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("gfa")
        .arg("aggregate-contacts")
        .arg("--segments")
        .arg(&segments)
        .arg("--contacts")
        .arg(&contacts)
        .arg("--gfa-links")
        .arg(&links)
        .arg("--output")
        .arg(&output)
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        fs::read_to_string(output).unwrap(),
        "#Segment1\tEnd1\tSegment2\tEnd2\tCount\tScore\tCompetitor1\tCompetitor2\nA\tR\tB\tL\t5\t0.001\t0.004\t0\n"
    );
}

#[test]
fn gfa_contract_inputs_cli_updates_block_distances() {
    let directory = tempfile::tempdir().unwrap();
    let mapping = directory.path().join("mapping.tsv");
    let contacts = directory.path().join("input.contacts");
    let clm = directory.path().join("input.clm");
    let output_contacts = directory.path().join("output.contacts");
    let output_clm = directory.path().join("output.clmb");
    fs::write(
        &mapping,
        "#Segment\tUnit\tHalf0\tHalf1\tPlusOrientation\tPlusLeftOffset\tPlusRightOffset\tMinusOrientation\tMinusLeftOffset\tMinusRightOffset\n\
A\tblock\tblock_0\tblock_0\t+\t135\t100\t-\t100\t135\n\
B\tblock\tblock_1\tblock_1\t+\t0\t95\t-\t95\t0\n\
C\tC\tC_0\tC_1\t+\t0\t7\t-\t0\t11\n",
    )
    .unwrap();
    fs::write(&contacts, "A_0\tC_1\t2.5\nA_0\tC_1\t3.5\nA_1\tB_0\t20\n").unwrap();
    fs::write(
        &clm,
        "A+ C+\t2\t10 20\nA+ C-\t1\t10\nB- C+\t1\t30\nB- C-\t1\t30\nA+ B+\t1\t5\n",
    )
    .unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("gfa")
        .arg("contract-inputs")
        .arg("--mapping")
        .arg(&mapping)
        .arg("--contacts")
        .arg(&contacts)
        .arg("--contacts-output")
        .arg(&output_contacts)
        .arg("--clm")
        .arg(&clm)
        .arg("--clm-output")
        .arg(&output_clm)
        .arg("--tmp-dir")
        .arg(directory.path())
        .arg("--sort-buffer")
        .arg("1M")
        .arg("--threads")
        .arg("4")
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        String::from_utf8_lossy(&result.stderr)
            .contains("GFA preprocessing timing: output sorting:")
    );
    assert_eq!(
        fs::read_to_string(output_contacts).unwrap(),
        "C_1\tblock_0\t6\n"
    );
    assert_eq!(
        clmb_as_text(output_clm),
        "block+ C+\t2\t152 162\nblock+ C-\t1\t156\nblock- C+\t1\t132\nblock- C-\t1\t136\n"
    );
}

#[test]
fn gfa_contract_inputs_cli_routes_clm_directly_to_cluster_files() {
    let directory = tempfile::tempdir().unwrap();
    let mapping = directory.path().join("mapping.tsv");
    let contacts = directory.path().join("input.contacts");
    let clm = directory.path().join("input.clm");
    let clusters = directory.path().join("clusters.txt");
    let output_contacts = directory.path().join("output.contacts");
    let output_dir = directory.path().join("split-clm");
    fs::write(
        &mapping,
        "A\tblock\tblock_0\tblock_0\t+\t135\t100\t-\t100\t135\n\
B\tblock\tblock_1\tblock_1\t+\t0\t95\t-\t95\t0\n\
C\tC\tC_0\tC_1\t+\t0\t7\t-\t0\t11\n\
D\tD\tD_0\tD_1\t+\t0\t0\t-\t0\t0\n",
    )
    .unwrap();
    fs::write(&contacts, "A_0\tC_1\t2\n").unwrap();
    fs::write(
        &clm,
        "A+ C+\t2\t10 20\nB- C+\t1\t30\nA+ D+\t1\t99\nA+ B+\t1\t5\n",
    )
    .unwrap();
    fs::write(
        &clusters,
        "group1\t2\tblock C\ngroup2\t2\tblock D\ngroup3\t0\n",
    )
    .unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("gfa")
        .arg("contract-inputs")
        .arg("--mapping")
        .arg(&mapping)
        .arg("--contacts")
        .arg(&contacts)
        .arg("--contacts-output")
        .arg(&output_contacts)
        .arg("--clm")
        .arg(&clm)
        .arg("--clusters")
        .arg(&clusters)
        .arg("--clm-output-dir")
        .arg(&output_dir)
        .arg("--tmp-dir")
        .arg(directory.path())
        .arg("--sort-buffer")
        .arg("1M")
        .arg("--threads")
        .arg("4")
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let group1_reader = ClmbReader::open(output_dir.join("group1.clmb")).unwrap();
    assert_eq!(group1_reader.header.target_block_size, 256 * 1024);
    assert_eq!(
        clmb_as_text(output_dir.join("group1.clmb")),
        "block+ C+\t2\t152 162\nblock- C+\t1\t132\n"
    );
    assert_eq!(
        clmb_as_text(output_dir.join("group2.clmb")),
        "block+ D+\t1\t234\n"
    );
    assert_eq!(clmb_as_text(output_dir.join("group3.clmb")), "");
}

#[test]
fn gfa_contract_inputs_cli_remaps_clmb_blocks_in_parallel() {
    let directory = tempfile::tempdir().unwrap();
    let mapping = directory.path().join("mapping.tsv");
    let contacts = directory.path().join("input.contacts");
    let clmb = directory.path().join("input.clmb");
    let clusters = directory.path().join("clusters.txt");
    let output_contacts = directory.path().join("output.contacts");
    let output_dir = directory.path().join("split-clm");
    fs::write(
        &mapping,
        "A\tblock\tblock_0\tblock_0\t+\t135\t100\t-\t100\t135\n\
B\tblock\tblock_1\tblock_1\t+\t0\t95\t-\t95\t0\n\
C\tC\tC_0\tC_1\t+\t0\t7\t-\t0\t11\n\
D\tD\tD_0\tD_1\t+\t0\t0\t-\t0\t0\n",
    )
    .unwrap();
    fs::write(&contacts, "A_0\tC_1\t2\n").unwrap();
    fs::write(
        &clusters,
        "group1\t2\tblock C\ngroup2\t2\tblock D\ngroup3\t0\n",
    )
    .unwrap();
    let contigs = vec!["A".into(), "B".into(), "C".into(), "D".into()];
    let mut writer = ClmbWriter::create_synchronous(&clmb, &contigs, 20, None, None).unwrap();
    writer
        .write_record(
            encode_endpoint(0, 0).unwrap(),
            encode_endpoint(2, 0).unwrap(),
            &[10, 20],
        )
        .unwrap();
    writer
        .write_record(
            encode_endpoint(1, 1).unwrap(),
            encode_endpoint(2, 0).unwrap(),
            &[30],
        )
        .unwrap();
    writer
        .write_record(
            encode_endpoint(0, 0).unwrap(),
            encode_endpoint(3, 0).unwrap(),
            &[99],
        )
        .unwrap();
    writer
        .write_record(
            encode_endpoint(0, 0).unwrap(),
            encode_endpoint(1, 0).unwrap(),
            &[5],
        )
        .unwrap();
    writer.finish().unwrap();

    let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .arg("gfa")
        .arg("contract-inputs")
        .arg("--mapping")
        .arg(&mapping)
        .arg("--contacts")
        .arg(&contacts)
        .arg("--contacts-output")
        .arg(&output_contacts)
        .arg("--clm")
        .arg(&clmb)
        .arg("--clusters")
        .arg(&clusters)
        .arg("--clm-output-dir")
        .arg(&output_dir)
        .arg("--threads")
        .arg("4")
        .output()
        .unwrap();

    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("GFA preprocessing timing: contract CLM:")
    );
    assert_eq!(
        clmb_as_text(output_dir.join("group1.clmb")),
        "block+ C+\t2\t152 162\nblock- C+\t1\t132\n"
    );
    assert_eq!(
        clmb_as_text(output_dir.join("group2.clmb")),
        "block+ D+\t1\t234\n"
    );
    assert_eq!(clmb_as_text(output_dir.join("group3.clmb")), "");
}
