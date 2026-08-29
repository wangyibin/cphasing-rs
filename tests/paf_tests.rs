use std::fs;

use cphasing::core::BaseTable;
use cphasing::paf::PAFTable;

fn run_depth(
    paf_path: &std::path::Path,
    chromsizes: &std::path::Path,
    output: &std::path::Path,
    min_mapq: u8,
    secondary: bool,
) {
    let paf = PAFTable::new(&paf_path.to_string_lossy().into_owned());
    paf.to_depth(
        &chromsizes.to_string_lossy().into_owned(),
        100,
        50,
        min_mapq,
        secondary,
        &output.to_string_lossy().into_owned(),
    )
    .unwrap();
}

#[test]
fn paf2depth_integrates_overlapping_windows_and_filters_alignments() {
    let directory = tempfile::tempdir().unwrap();
    let paf_path = directory.path().join("input.paf");
    let chromsizes = directory.path().join("input.chromsizes");
    fs::write(&chromsizes, "A\t1000\nB\t500\n").unwrap();
    fs::write(
        &paf_path,
        "r1\t1000\t0\t300\t+\tA\t1000\t10\t310\t290\t300\t60\ttp:A:P\n\
r2\t900\t0\t250\t+\tA\t1000\t100\t350\t240\t250\t40\ttp:A:P\n\
r3\t800\t0\t200\t+\tA\t1000\t300\t500\t190\t200\t50\ttp:A:S\n\
r4\t500\t0\t100\t+\tB\t500\t0\t100\t95\t100\t50\ttp:A:P\n",
    )
    .unwrap();

    let primary_output = directory.path().join("primary.depth");
    run_depth(&paf_path, &chromsizes, &primary_output, 0, false);
    assert_eq!(
        fs::read_to_string(&primary_output)
            .unwrap()
            .lines()
            .take(8)
            .collect::<Vec<_>>(),
        vec![
            "A\t0\t100\t0.900",
            "A\t50\t150\t1.500",
            "A\t100\t200\t2.000",
            "A\t150\t250\t2.000",
            "A\t200\t300\t2.000",
            "A\t250\t350\t1.600",
            "A\t300\t400\t0.600",
            "A\t350\t450\t0.000",
        ]
    );

    let filtered_output = directory.path().join("filtered.depth");
    run_depth(&paf_path, &chromsizes, &filtered_output, 45, false);
    assert_eq!(
        fs::read_to_string(&filtered_output)
            .unwrap()
            .lines()
            .take(8)
            .collect::<Vec<_>>(),
        vec![
            "A\t0\t100\t0.900",
            "A\t50\t150\t1.000",
            "A\t100\t200\t1.000",
            "A\t150\t250\t1.000",
            "A\t200\t300\t1.000",
            "A\t250\t350\t0.600",
            "A\t300\t400\t0.100",
            "A\t350\t450\t0.000",
        ]
    );

    let secondary_output = directory.path().join("secondary.depth");
    run_depth(&paf_path, &chromsizes, &secondary_output, 45, true);
    assert_eq!(
        fs::read_to_string(&secondary_output)
            .unwrap()
            .lines()
            .skip(5)
            .take(3)
            .collect::<Vec<_>>(),
        vec![
            "A\t250\t350\t1.100",
            "A\t300\t400\t1.100",
            "A\t350\t450\t1.000",
        ]
    );
}
