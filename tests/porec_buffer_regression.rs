use cphasing::{core::BaseTable, porec::PoreCTable};
use std::fs;

#[test]
fn moved_buffers_preserve_break_and_chr_mapping_across_batches() {
    let dir = tempfile::tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    let bed = dir.path().join("map.bed");
    fs::write(&bed, "A\t0\t100\tcut\n").unwrap();
    let input = dir.path().join("in.porec");
    for n in [0, 1, 4999, 5000, 5001, 15001] {
        let rows: Vec<_> = (0..n)
            .map(|i| {
                format!(
                    "r{i}\t0\t0\t0\t0\t{}\t10\t20\t60\tx\ty",
                    if i % 3 == 0 { "B" } else { "A" }
                )
            })
            .collect();
        fs::write(
            &input,
            if n == 0 {
                String::new()
            } else {
                rows.join("\n") + "\n"
            },
        )
        .unwrap();
        for threads in [1, 4] {
            for convert in [false, true] {
                let output = dir.path().join("out.porec").to_string_lossy().into_owned();
                let mut table = PoreCTable::new(&input.to_string_lossy().into_owned());
                if convert {
                    table
                        .chr_porec_to_contig_porec(
                            &bed.to_string_lossy().into_owned(),
                            &output,
                            threads,
                        )
                        .unwrap();
                } else {
                    table
                        .break_contigs(
                            &bed.to_string_lossy().into_owned(),
                            &output,
                            threads,
                            5000,
                            None,
                        )
                        .unwrap();
                }
                let mut expected = String::new();
                for (i, row) in rows.iter().enumerate() {
                    if i % 3 == 0 {
                        if !convert {
                            expected.push_str(row);
                            expected.push('\n');
                        }
                    } else {
                        expected.push_str(&row.replace(
                            "\tA\t10\t20\t",
                            if convert {
                                "\tcut\t10\t20\t"
                            } else {
                                "\tcut\t11\t21\t"
                            },
                        ));
                        expected.push('\n');
                    }
                }
                assert_eq!(
                    fs::read_to_string(output).unwrap(),
                    expected,
                    "n={n}, threads={threads}, convert={convert}"
                );
            }
        }
    }
}

#[test]
fn moved_coitree_buffers_preserve_empty_chunks_invert_and_tail() {
    let dir = tempfile::tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    let bed = dir.path().join("regions.bed");
    fs::write(&bed, "A\t0\t100\n").unwrap();
    let input = dir.path().join("in.porec");
    let rows: Vec<_> = (0..10001)
        .map(|i| {
            format!(
                "r{i}\t0\t0\t0\t0\t{}\t10\t20\t60\tx",
                if i < 10000 { "B" } else { "A" }
            )
        })
        .collect();
    fs::write(&input, rows.join("\n") + "\n").unwrap();
    for invert in [false, true] {
        let output = dir.path().join("out.porec").to_string_lossy().into_owned();
        PoreCTable::new(&input.to_string_lossy().into_owned()).intersect_multi_threads_coitree(
            &bed.to_string_lossy().into_owned(),
            invert,
            &output,
        );
        let expected = rows
            .iter()
            .enumerate()
            .filter(|(i, _)| (*i == 10000) ^ invert)
            .map(|(_, r)| format!("{r}\n\n"))
            .collect::<String>();
        assert_eq!(fs::read_to_string(output).unwrap(), expected);
    }
}
