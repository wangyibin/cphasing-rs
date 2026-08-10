use cphasing::core::common_reader;
use cphasing::splitcontacts::split_contacts_by_clusters;
use std::fs;
use std::io::BufRead;
use std::path::Path;

fn read_lines(path: &Path) -> Vec<String> {
    common_reader(path.to_str().unwrap())
        .lines()
        .collect::<Result<Vec<_>, _>>()
        .unwrap()
}

#[test]
fn splits_contacts_into_single_threaded_gzip_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let clusters = directory.path().join("clusters.txt");
    let contacts = directory.path().join("contacts.tsv");
    fs::write(&clusters, "group1 A B\ngroup2 C D\n").unwrap();
    fs::write(&contacts, "A_0\tB_1\t5\nC_1\tD_0\t7\nA_1\tC_0\t11\n").unwrap();

    split_contacts_by_clusters(
        &contacts.to_string_lossy().into_owned(),
        &clusters.to_string_lossy().into_owned(),
        &directory.path().to_string_lossy().into_owned(),
    )
    .unwrap();

    assert_eq!(
        read_lines(&directory.path().join("group1.split.contacts.gz")),
        vec!["A_0\tB_1\t5"]
    );
    assert_eq!(
        read_lines(&directory.path().join("group2.split.contacts.gz")),
        vec!["C_1\tD_0\t7"]
    );
}
