use cphasing::sketch::{MinimizerData, MinimizerInfo, sketch};

const SEQUENCE: &[u8] = b"TTTCGACAGTTCTCCCTGGCACCTCTGAAAGCTTTCCTGGTTTGATTGTTGAAAGTCTTAGGGCTCAACTTGGTCAGCCCTTCTTCATGGAAATTGTTATGACCATGTGTTGGTCCATCTGGATGATGCGCAATGATGTCATTTTCAAAGGTTTAC";

#[test]
fn minimizer_layout_stays_compact() {
    assert_eq!(std::mem::size_of::<MinimizerInfo>(), 16);
    assert_eq!(std::mem::size_of::<MinimizerData>(), 24);
}

#[test]
fn partig_minimizer_count_is_reproduced() {
    let minimizers = sketch(SEQUENCE, 7, 19, 19);

    // partig_pthreads_large_contig emits 12 minimizers for this sequence.
    assert_eq!(minimizers.len(), 12);
    assert!(minimizers.iter().all(|item| item.info.rid == 7));
    assert!(
        minimizers
            .windows(2)
            .all(|pair| pair[0].info.pos < pair[1].info.pos)
    );
}

#[test]
fn sketch_handles_short_lowercase_and_ambiguous_sequences() {
    assert!(sketch(b"ACGT", 0, 19, 19).is_empty());
    assert_eq!(
        sketch(SEQUENCE, 0, 19, 19),
        sketch(&SEQUENCE.to_ascii_lowercase(), 0, 19, 19)
    );

    let mut with_gap = SEQUENCE.to_vec();
    with_gap[70] = b'N';
    assert!(
        sketch(&with_gap, 0, 19, 19)
            .iter()
            .all(|item| item.info.pos < 70 || item.info.pos >= 89)
    );
}
