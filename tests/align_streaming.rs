use rust_htslib::bam::{self, Read, record::Aux};
use std::path::Path;
use std::process::Command;

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
        .args(args).env("CPHASING_IO_THREADS", "1").output().unwrap()
}

fn success(args: &[&str]) {
    let output = run(args);
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
}

fn records(path: &Path) -> Vec<bam::Record> {
    bam::Reader::from_path(path).unwrap().records().map(Result::unwrap).collect()
}

fn sequence() -> Vec<u8> {
    let mut state = 421u64;
    (0..2400).map(|_| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        b"ACGT"[(state & 3) as usize]
    }).collect()
}

#[test]
fn align_help_and_option_validation() {
    success(&["align", "--help"]);
    let output = run(&["align", "missing.fa", "reads.fq", "--realign", "unsupported"]);
    assert!(!output.status.success());
    let output = run(&["align", "missing.fa", "reads.fq", "--meth-bed", "ref.bed"]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("requires unaligned BAM/SAM"));
}

#[test]
fn cached_headers_preserve_contigs_across_reads_threads_and_batches() {
    let dir = tempfile::Builder::new().prefix(".align-header-test-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    let reference = dir.path().join("ref.fa");
    let input = dir.path().join("reads.fq");
    let first = sequence();
    let second: Vec<u8> = first.iter().copied().rev().collect();
    let targets = [first, second];
    std::fs::write(&reference, format!(">contig_a\n{}\n>contig_b\n{}\n",
        std::str::from_utf8(&targets[0]).unwrap(),
        std::str::from_utf8(&targets[1]).unwrap())).unwrap();
    let mut fastq = String::new();
    for i in 0..128 {
        let mut bases = targets[i % 2][400..1600].to_vec();
        if i % 4 >= 2 { bases = bio::alphabets::dna::revcomp(&bases); }
        fastq.push_str(&format!("@read{i}\n{}\n+\n{}\n",
            std::str::from_utf8(&bases).unwrap(), "I".repeat(bases.len())));
    }
    std::fs::write(&input, fastq).unwrap();
    let mut expected = None;
    for (threads, batch) in [("1", "50M"), ("4", "50M"), ("4", "8K")] {
        let output = dir.path().join(format!("{threads}-{batch}.bam"));
        let result = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"))
            .args(["align", reference.to_str().unwrap(), input.to_str().unwrap(),
                "-x", "porec:map-ont", "-t", threads, "-K", batch, "-p", "0.8",
                "--no-realign", "-a", "-o", output.to_str().unwrap()])
            .env("CPHASING_IO_THREADS", "1").env("CPHASING_ALIGN_TIMING", "1")
            .output().unwrap();
        let log = String::from_utf8_lossy(&result.stderr);
        assert!(result.status.success(), "{log}");
        let calls = |stage: &str| -> usize {
            log.lines().find(|line| line.contains(&format!("stage={stage} ")))
                .unwrap().split("calls=").nth(1).unwrap().parse().unwrap()
        };
        assert_eq!(calls("rammap_align_and_format"), 128);
        assert!(calls("sam_parse") >= 128);
        assert!(calls("header_create") > 0 && calls("header_create") <= 128);
        if threads == "1" {
            assert!(calls("header_create") < 128, "the single-thread state must reuse headers");
        }
        let output_records = records(&output);
        assert_eq!(output_records.len(), 128);
        for (i, record) in output_records.iter().enumerate() {
            assert_eq!(record.qname(), format!("read{i}").as_bytes());
            assert_eq!(record.tid(), (i % 2) as i32);
            assert_eq!(record.pos(), 400);
            assert_eq!(record.is_reverse(), i % 4 >= 2);
            assert!(!record.is_unmapped() && !record.is_secondary() && !record.is_supplementary());
            assert_eq!(record.seq().as_bytes(), targets[i % 2][400..1600]);
            assert_eq!(record.qual(), &[40; 1200]);
        }
        let identities: Vec<_> = output_records.iter().map(|r|
            (r.tid(), r.pos(), r.flags(), r.mapq(), r.cigar().to_string())).collect();
        if let Some(ref previous) = expected { assert_eq!(&identities, previous); }
        else { expected = Some(identities); }
    }
}

#[test]
fn reference_ownership_and_timing_preserve_full_read_records() {
    let dir = tempfile::Builder::new().prefix(".align-ownership-test-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    let reference = dir.path().join("ref.fa");
    let input = dir.path().join("reads.bam");
    let bed = dir.path().join("meth.bg");
    let seq = sequence();
    // Owned reference normalization must preserve lower-case FASTA behavior.
    std::fs::write(&reference, format!(">ref\n{}\n",
        String::from_utf8(seq.clone()).unwrap().to_ascii_lowercase())).unwrap();
    std::fs::write(&bed, "ref\t0\t1\t100\n").unwrap();
    let forward = seq[400..1600].to_vec();
    let reverse = bio::alphabets::dna::revcomp(&forward);
    let sequences = [forward, reverse, vec![b'N'; 1200]];
    let mut originals = Vec::new();
    {
        let mut writer = bam::Writer::from_path(&input, &bam::Header::new(), bam::Format::Bam).unwrap();
        for (i, bases) in sequences.iter().enumerate() {
            let mut rec = bam::Record::new();
            let quality: Vec<_> = (0..bases.len()).map(|p| 10 + (p % 30) as u8).collect();
            rec.set(format!("read{i}").as_bytes(), None, bases, &quality);
            rec.set_flags(4);
            rec.set_tid(-1);
            rec.set_pos(-1);
            let nc = bases.iter().filter(|&&b| b == b'C').count();
            let mm = format!("C+m{};", ",0".repeat(nc));
            let ml: Vec<u8> = (0..nc).map(|p| 180 + (p % 70) as u8).collect();
            rec.push_aux(b"MM", Aux::String(&mm)).unwrap();
            rec.push_aux(b"ML", Aux::ArrayU8(ml.as_slice().into())).unwrap();
            rec.push_aux(b"MN", Aux::I32(bases.len() as i32)).unwrap();
            rec.push_aux(b"RX", Aux::String("keep-this-tag")).unwrap();
            writer.write(&rec).unwrap();
            originals.push(rec);
        }
    }
    let mut expected = None;
    for mode in ["sequence", "gap", "methylation", "methylation-gap"] {
        let output = dir.path().join(format!("{mode}.bam"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_cphasing-rs"));
        command.args(["align", reference.to_str().unwrap(), input.to_str().unwrap(),
            "-x", "porec:map-ont", "-t", "1", "-p", "0.8", "--no-realign", "-a", "-o", output.to_str().unwrap()])
            .env("CPHASING_IO_THREADS", "1");
        if mode.contains("gap") { command.arg("--gap-rescue"); }
        if mode.contains("methylation") {
            command.args(["--meth-bed", bed.to_str().unwrap(), "--meth-cpg"]);
        }
        if mode != "sequence" { command.env("CPHASING_ALIGN_TIMING", "1"); }
        let result = command.output().unwrap();
        assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
        if mode != "sequence" {
            assert!(String::from_utf8_lossy(&result.stderr).contains("ALIGN_TIMING stage=restore_primary"));
        }
        let output_records = records(&output);
        assert_eq!(output_records.len(), 3);
        for (record, original) in output_records.iter().zip(&originals) {
            assert_eq!(record.qname(), original.qname());
            let mut bases = original.seq().as_bytes();
            let mut quality = original.qual().to_vec();
            if record.is_reverse() {
                bases = bio::alphabets::dna::revcomp(&bases);
                quality.reverse();
            }
            assert_eq!(record.seq().as_bytes(), bases);
            assert_eq!(record.qual(), quality);
            for tag in [b"MM", b"ML", b"MN", b"RX"] {
                assert_eq!(record.aux(tag).unwrap(), original.aux(tag).unwrap());
            }
        }
        assert!(!output_records[0].is_reverse());
        assert!(output_records[1].is_reverse());
        assert!(output_records[2].is_unmapped());
        let identities: Vec<_> = output_records.iter().map(|r|
            (r.tid(), r.pos(), r.flags(), r.mapq(), r.cigar().to_string())).collect();
        if let Some(ref previous) = expected { assert_eq!(&identities, previous); }
        else { expected = Some(identities); }
    }
}

#[test]
fn raw_modbam_to_refined_bam_and_paf_without_intermediate_alignment() {
    // Keep the small, generated fixture inside the repository.
    let dir = tempfile::Builder::new().prefix(".align-test-")
        .tempdir_in(env!("CARGO_MANIFEST_DIR")).unwrap();
    let reference = dir.path().join("ref.fa");
    let reads = dir.path().join("reads.bam");
    let baseline = dir.path().join("baseline.bam");
    let refined = dir.path().join("refined.bam");
    let paf = dir.path().join("refined.paf");
    let bed = dir.path().join("ref.bed");
    let seq = sequence();
    let text = std::str::from_utf8(&seq).unwrap();
    std::fs::write(&reference, format!(">copy1\n{text}\n>copy2\n{text}\n")).unwrap();
    let read_seq = &seq[400..1600];
    let n_c = read_seq.iter().filter(|&&b| b == b'C').count();
    let mm = format!("C+m{};", ",0".repeat(n_c));
    let ml = vec![255u8; n_c];
    let mut read = bam::Record::new();
    read.set(b"molecule", None, read_seq, &vec![30; read_seq.len()]);
    read.set_flags(4);
    read.set_tid(-1);
    read.set_pos(-1);
    read.push_aux(b"MM", Aux::String(&mm)).unwrap();
    read.push_aux(b"ML", Aux::ArrayU8(ml.as_slice().into())).unwrap();
    {
        let mut writer = bam::Writer::from_path(&reads, &bam::Header::new(), bam::Format::Bam).unwrap();
        writer.write(&read).unwrap();
    }
    let rf = reference.to_str().unwrap();
    let rd = reads.to_str().unwrap();
    success(&["align", rf, rd, "-x", "porec", "-t", "1", "--no-realign",
        "--secondary", "yes", "-a", "-o", baseline.to_str().unwrap()]);
    let before = records(&baseline);
    let primary = before.iter().find(|r| !r.is_secondary() && !r.is_supplementary()).unwrap();
    assert!(before.iter().any(|r| r.is_secondary()));
    let desired_tid = 1 - primary.tid();
    let desired_name = if desired_tid == 0 { "copy1" } else { "copy2" };
    let bed_text: String = seq.iter().enumerate().filter(|(_, b)| **b == b'C')
        .map(|(p, _)| format!("{desired_name}\t{p}\t{}\t100\n", p + 1)).collect();
    std::fs::write(&bed, bed_text).unwrap();
    success(&["align", rf, rd, "-x", "porec", "-t", "1", "--realign", "precise",
        "--meth-bed", bed.to_str().unwrap(), "-a", "-o", refined.to_str().unwrap()]);
    let after = records(&refined);
    assert_eq!(after.len(), 1, "secondary filtering must occur after methylation scoring");
    assert_eq!(after[0].tid(), desired_tid);
    assert_eq!(after[0].mapq(), 2);
    assert!(after[0].aux(b"s0").is_ok());
    assert!(after[0].aux(b"m0").is_ok());
    assert_eq!(after[0].seq_len(), read_seq.len());
    assert!(after[0].aux(b"MM").is_ok());
    assert_eq!(after[0].basemods_iter().unwrap().count(), n_c);
    success(&["align", rf, rd, "-x", "porec", "-t", "1", "--realign", "sensitive",
        "--meth-bed", bed.to_str().unwrap(), "-c", "-o", paf.to_str().unwrap()]);
    let text = std::fs::read_to_string(&paf).unwrap();
    let cols: Vec<_> = text.lines().next().unwrap().split('\t').collect();
    assert_eq!(cols[5], desired_name);
    assert!(text.contains("\ts0:i:"));

    // Methylation alone must retain alternatives internally even with --no-realign.
    success(&["align", rf, rd, "-x", "porec", "-t", "2", "--no-realign",
        "--meth-bed", bed.to_str().unwrap(), "-a", "-o", refined.to_str().unwrap()]);
    assert_eq!(records(&refined)[0].tid(), desired_tid);

    // The compatibility command uses the same core on an existing alignment.
    let offline = dir.path().join("offline.bam");
    success(&["methalign", baseline.to_str().unwrap(), "-f", rf, "-b", bed.to_str().unwrap(),
        "-t", "1", "-o", offline.to_str().unwrap()]);
    assert_eq!(records(&offline)[0].tid(), desired_tid);

    // Sequence-only FASTQ and BAM stdout are also supported.
    let fastq = dir.path().join("reads.fq");
    std::fs::write(&fastq, format!("@molecule\n{}\n+\n{}\n",
        std::str::from_utf8(read_seq).unwrap(), "I".repeat(read_seq.len()))).unwrap();
    let output = run(&["align", rf, fastq.to_str().unwrap(), "-x", "cifi", "-t", "1",
        "--no-realign", "-a", "-o", "-"]);
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    assert!(output.stdout.starts_with(&[0x1f, 0x8b]));

    let bad = dir.path().join("missing-tags.bam");
    {
        read.remove_aux(b"MM").unwrap();
        let mut writer = bam::Writer::from_path(&bad, &bam::Header::new(), bam::Format::Bam).unwrap();
        writer.write(&read).unwrap();
    }
    let output = run(&["align", rf, bad.to_str().unwrap(), "-t", "1", "--meth-bed",
        bed.to_str().unwrap(), "-o", dir.path().join("bad.paf").to_str().unwrap()]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("missing MM/ML"));
}

mod methylation_cases {
    use super::*;
    use cphasing::methalign::{MethConfig, MethRefiner, get_as};
    use cphasing::align::engine::AlignmentUnit;

    fn original() -> bam::Record {
        let mut rec = bam::Record::new();
        rec.set(b"read", None, b"ACGTCG", &[30; 6]);
        rec.set_flags(4);
        rec.push_aux(b"MM", Aux::String("C+m,0,0;")).unwrap();
        rec.push_aux(b"ML", Aux::ArrayU8((&[255u8, 255][..]).into())).unwrap();
        rec
    }

    fn candidate(tid: i32, flags: u16) -> bam::Record {
        let mut rec = bam::Record::new();
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(6)]);
        rec.set(b"read", Some(&cigar), b"ACGTCG", &[30; 6]);
        rec.set_tid(tid);
        rec.set_pos(0);
        rec.set_flags(flags);
        rec.set_mapq(0);
        rec.push_aux(b"AS", Aux::I32(100)).unwrap();
        rec
    }

    fn refiner(seq: &[u8], bed: &str) -> MethRefiner {
        use std::io::Write;
        let mut file = tempfile::Builder::new().prefix(".meth-test-")
            .tempfile_in(env!("CARGO_MANIFEST_DIR")).unwrap();
        file.write_all(bed.as_bytes()).unwrap();
        MethRefiner::new(MethConfig {
            bed: file.path().to_str().unwrap().to_string(), match_score: 0,
            ref_penalty: 2, read_penalty: 2, ref_prob_cutoff: 50.0,
            prob_cutoff: 128, designate_mapq: 2, cpg: true,
        }, &[("a".into(), seq.to_vec()), ("b".into(), seq.to_vec())]).unwrap()
    }

    #[test]
    fn reverse_supplementary_promotion_preserves_segment_role() {
        let refiner = refiner(b"CGACGT", "b\t0\t1\t100\nb\t3\t4\t100\n");
        let result = refiner.refine(vec![candidate(0, 0x810), candidate(1, 0x110)],
            &original(), &["a".into(), "b".into()]);
        assert_eq!(result[0].tid(), 1);
        assert!(result[0].is_reverse());
        assert!(result[0].is_supplementary());
        assert!(!result[0].is_secondary());
        assert!(result[1].is_secondary());
        assert!(!result[1].is_supplementary());
        assert_eq!(get_as(&result[0]), Some(100));
        assert_eq!(get_as(&result[1]), Some(96));
    }

    #[test]
    fn equal_methylation_evidence_does_not_promote_a_tie() {
        let refiner = refiner(b"ACGTCG", "a\t1\t2\t100\na\t4\t5\t100\nb\t1\t2\t100\nb\t4\t5\t100\n");
        let result = refiner.refine(vec![candidate(0, 0), candidate(1, 0x100)],
            &original(), &["a".into(), "b".into()]);
        assert_eq!(result[0].tid(), 0);
        assert_eq!(result[0].mapq(), 0);
    }

    #[test]
    fn precise_and_sensitive_rescue_use_other_fragment_anchor() {
        for sensitive in [false, true] {
            let mut anchor = candidate(0, 0);
            anchor.set_mapq(60);
            let ambiguous = candidate(1, 0x800);
            let alternative = candidate(0, 0x100);
            let mut unit = AlignmentUnit::from_records(vec![anchor, ambiguous, alternative]);
            if sensitive { unit.rescue(1, &[b"a", b"b"]); }
            else { unit.rescue_robust(1, &[b"a", b"b"]); }
            let result = unit.into_records(true);
            assert_eq!(result.len(), 3);
            assert_eq!(result[2].tid(), 1, "displaced primary must remain an alternative");
            assert!(result[2].is_secondary());
            assert_eq!(result[1].tid(), 0);
            assert!(result[1].mapq() > 0);
            assert!(result[1].is_supplementary());
        }
    }
}
