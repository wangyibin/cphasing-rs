use rust_htslib::bam::{self, Read, record::{Aux, Cigar, CigarString}, header::HeaderRecord};
use std::process::Command;

fn run_case(mm: &str, ml: &[u8], scores: &[Option<i16>], reverse_candidate: bool) -> Vec<bam::Record> {
    let dir = tempfile::Builder::new().prefix(".meth-rescue-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    let fa = dir.path().join("ref.fa");
    let bed = dir.path().join("meth.bg");
    std::fs::write(&fa, format!(">a\n{}\n>b\n{}\n", "CG".repeat(50), "CG".repeat(50))).unwrap();
    std::fs::write(&bed, "a\t0\t1\t100\n").unwrap();
    let mut header = bam::Header::new();
    for name in ["a", "b"] { header.push_record(HeaderRecord::new(b"SQ").push_tag(b"SN", name).push_tag(b"LN", 100)); }
    let input = dir.path().join("input.bam");
    let output = dir.path().join("output.bam");
    {
        let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
        for (i, score) in scores.iter().enumerate() {
            let mut r = bam::Record::new();
            let seq: &[u8] = if i == 0 { b"CGCG" } else { b"" };
            let qual: &[u8] = if i == 0 { &[10,20,30,40] } else { &[] };
            r.set(b"read", Some(&CigarString(vec![Cigar::Match(4)])), seq, qual);
            r.set_tid(if i == 0 {0} else {1}); r.set_pos(0); r.set_mapq(0);
            r.set_flags(if i == 0 {0} else {256 | if reverse_candidate {16} else {0}});
            if let Some(score) = score { r.push_aux(b"AS", Aux::I16(*score)).unwrap(); }
            r.push_aux(b"tp", Aux::Char(if i == 0 {b'P'} else {b'S'})).unwrap();
            r.push_aux(b"SA", Aux::String("a,99,+,4M,60,0;")).unwrap();
            if i == 0 {
                r.push_aux(b"MM", Aux::String(mm)).unwrap();
                r.push_aux(b"ML", Aux::ArrayU8(ml.into())).unwrap();
                r.push_aux(b"MN", Aux::I32(4)).unwrap();
            }
            writer.write(&r).unwrap();
        }
    }
    let mut command = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"));
    command.args(["methalign", input.to_str().unwrap(), "-t", "1", "-s", "-o", output.to_str().unwrap()]);
    if !mm.is_empty() { command.args(["-f", fa.to_str().unwrap(), "-b", bed.to_str().unwrap(), "--cpg"]); }
    let result = command.output().unwrap();
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    bam::Reader::from_path(output).unwrap().records().map(Result::unwrap).collect()
}

#[test]
fn unknown_skipped_bases_do_not_rescue_but_implicit_low_calls_do() {
    for mm in ["C+m?;", "C+m?,1;"] {
        let r = run_case(mm, &[255], &[Some(100),Some(100)], false);
        assert_eq!(r[0].tid(), 0); assert_eq!(r[0].mapq(), 0);
        assert_eq!(r[0].aux(b"m0").unwrap(), r[1].aux(b"m0").unwrap());
    }
    for mm in ["C+m;", "C+m.;", "C+m?,0;"] {
        let r = run_case(mm, &[0], &[Some(100),Some(100)], false);
        assert_eq!(r[0].tid(), 1); assert_eq!(r[0].mapq(), 2);
    }
}

#[test]
fn promoted_primary_recovers_sequence_qualities_modifications_and_tags() {
    let r = run_case("C+m?,0;", &[0], &[Some(100),Some(100)], true);
    assert_eq!(r[0].tid(), 1); assert_eq!(r[0].flags(), 16);
    assert_eq!(r[0].seq().as_bytes(), b"CGCG"); assert_eq!(r[0].qual(), &[40,30,20,10]);
    assert_eq!(r[0].aux(b"MM").unwrap(), Aux::String("C+m?,0;"));
    assert!(r[0].aux(b"ML").is_ok()); assert!(r[0].aux(b"MN").is_ok());
    assert_eq!(r[0].aux(b"tp").unwrap(), Aux::Char(b'P'));
    assert_eq!(r[1].aux(b"tp").unwrap(), Aux::Char(b'S'));
    assert!(r.iter().all(|r| r.aux(b"SA").is_err()));
}

#[test]
fn missing_scores_cannot_shift_indices_and_negative_scores_are_supported() {
    for mm in ["", "C+m?,0;"] {
        let r = run_case(mm, &[0], &[None,Some(100),Some(90)], false);
        assert_eq!(r[0].tid(), 0); assert_eq!(r[0].mapq(), 0); assert!(r[0].aux(b"RF").is_err());
    }
    for (primary, secondary) in [(-10, -20), (-300, -400)] {
        let r = run_case("", &[], &[Some(primary),Some(secondary)], false);
        assert_eq!(r[0].mapq(), 2);
    }
}

#[test]
fn multicode_ml_stride_preserves_methylation_probability() {
    // m=0, h=255 at the first C, then m=0,h=255 at the second C.
    let single = run_case("C+m?,0,0;", &[0,0], &[Some(100),Some(100)], false);
    for (mm, ml) in [("C+mh?,0,0;", & [0,255,0,255][..]),
        ("C+hm?,0,0;", &[255,0,255,0][..]),
        ("C+h?,0;C+m?,0,0;", &[255,0,0][..])] {
        let multi = run_case(mm, ml, &[Some(100),Some(100)], false);
        for (a,b) in multi.iter().zip(single.iter()) {
            assert_eq!((a.tid(),a.mapq(),a.flags()), (b.tid(),b.mapq(),b.flags()));
            assert_eq!(a.aux(b"m0").unwrap(), b.aux(b"m0").unwrap());
        }
    }
}

#[test]
fn an_explicit_call_without_ml_is_not_unmethylated_evidence() {
    let r = run_case("C+m?,0;", &[], &[Some(100),Some(100)], false);
    assert_eq!(r[0].tid(), 0); assert_eq!(r[0].mapq(), 0);
}
