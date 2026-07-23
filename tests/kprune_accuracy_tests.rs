use cphasing::cli;
use cphasing::contacts::{ContactEvidence, ContactRecord, Contacts2};
use cphasing::kprune::{
    GroupDecision, KPruner, PruneThresholds, all_informative_groups_cross_allelic,
    classify_selected_alternative,
};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{self, Write};
use std::sync::{Arc, Mutex};

#[test]
fn aggregates_duplicate_contacts_before_normalization() {
    let contacts = Contacts2 {
        file: String::new(),
        records: vec![
            ContactRecord {
                chrom1: "a".into(),
                chrom2: "b".into(),
                count: 2.0,
            },
            ContactRecord {
                chrom1: "b".into(),
                chrom2: "a".into(),
                count: 3.0,
            },
        ],
    };

    let evidence = contacts.to_evidence_data(&HashMap::new(), &"none".to_string(), &None, None);
    let value = evidence.values().next().unwrap();

    assert_eq!(evidence.len(), 1);
    assert_eq!(value.raw_count, 5.0);
    assert_eq!(value.score, 5.0);
}

#[test]
fn missing_alternative_is_insufficient_evidence() {
    assert_eq!(
        classify_selected_alternative(
            Some(ContactEvidence {
                raw_count: 5.0,
                score: 10.0,
            }),
            None,
            false,
            PruneThresholds::default(),
        ),
        GroupDecision::InsufficientEvidence,
    );
}

#[test]
fn tied_or_weak_alternative_is_insufficient_evidence() {
    for selected_score in [10.0, 10.5] {
        assert_eq!(
            classify_selected_alternative(
                Some(ContactEvidence {
                    raw_count: 5.0,
                    score: 10.0,
                }),
                Some(ContactEvidence {
                    raw_count: 5.0,
                    score: selected_score,
                }),
                false,
                PruneThresholds::default(),
            ),
            GroupDecision::InsufficientEvidence,
        );
    }
}

#[test]
fn supported_stronger_alternative_votes_cross_allelic() {
    assert_eq!(
        classify_selected_alternative(
            Some(ContactEvidence {
                raw_count: 2.0,
                score: 10.0,
            }),
            Some(ContactEvidence {
                raw_count: 3.0,
                score: 12.0,
            }),
            false,
            PruneThresholds::default(),
        ),
        GroupDecision::CrossAllelic,
    );
}

#[test]
fn conflicting_group_votes_retain_pair() {
    assert!(!all_informative_groups_cross_allelic([
        GroupDecision::CrossAllelic,
        GroupDecision::Compatible,
    ]));
    assert!(!all_informative_groups_cross_allelic([
        GroupDecision::InsufficientEvidence,
    ]));
}

#[derive(Clone, Default)]
struct SharedBuffer(Arc<Mutex<Vec<u8>>>);

impl Write for SharedBuffer {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(bytes);
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn run_kprune_filtered(whitelist_names: &[&str], partial_whitelist: bool) -> String {
    let directory = tempfile::tempdir().unwrap();
    let allele_path = directory.path().join("accuracy.allele.table");
    let contacts_path = directory.path().join("accuracy.contacts");
    fs::write(
        &allele_path,
        concat!(
            "#a 1000 100 100\n#a2 1000 100 100\n#b 1000 100 100\n#b2 1000 100 100\n",
            "#c 1000 100 100\n#c2 1000 100 100\n#d 1000 100 100\n#d2 1000 100 100\n",
            "#e 1000 100 100\n#e2 1000 100 100\n#f 1000 100 100\n#f2 1000 100 100\n",
            "0\t1\ta\ta2\t100\t100\t90\t0.9\t1\n",
            "2\t3\tb\tb2\t100\t100\t90\t0.9\t1\n",
            "4\t5\tc\tc2\t100\t100\t90\t0.9\t1\n",
            "6\t7\td\td2\t100\t100\t90\t0.9\t1\n",
            "8\t9\te\te2\t100\t100\t90\t0.9\t1\n",
            "10\t11\tf\tf2\t100\t100\t90\t0.9\t1\n",
        ),
    )
    .unwrap();
    fs::write(
        &contacts_path,
        concat!(
            "a\ta2\t1\nb\tb2\t1\nc\tc2\t1\nd\td2\t1\ne\te2\t1\nf\tf2\t1\n",
            "a\tb\t4\nb\ta\t6\na\tb2\t12\na2\tb\t12\na2\tb2\t10\n",
            "c\td\t10\nc\td2\t10\nc2\td\t10\nc2\td2\t10\n",
            "e\tf\t10\ne\tf2\t10.5\ne2\tf\t10.5\ne2\tf2\t10\n",
        ),
    )
    .unwrap();

    let allele_path = allele_path.to_string_lossy().to_string();
    let contacts_path = contacts_path.to_string_lossy().to_string();
    let output_path = "-".to_string();
    let normalization = "none".to_string();
    let sink = SharedBuffer::default();
    let mut writer: Box<dyn Write + Send> = Box::new(sink.clone());
    let mut pruner = KPruner::new(
        &allele_path,
        &contacts_path,
        &output_path,
        &None,
        &normalization,
    );
    pruner.thresholds = PruneThresholds::default();
    let whitelist: HashSet<String> = whitelist_names
        .iter()
        .map(|name| name.to_string())
        .collect();
    let whitelist_refs: HashSet<&String> = whitelist.iter().collect();
    pruner.prune("precise", &whitelist_refs, &mut writer, partial_whitelist);
    drop(writer);

    let bytes = sink.0.lock().unwrap().clone();
    String::from_utf8(bytes).unwrap()
}

fn run_kprune() -> String {
    run_kprune_filtered(&[], false)
}

#[test]
fn kprune_uses_margin_and_emits_stable_output() {
    let first = run_kprune();
    let second = run_kprune();

    assert_eq!(first, second);
    assert!(first.lines().all(|line| line.split('\t').count() == 7));
    assert!(first.lines().any(|line| line == "a\tb\t0\t0\t0\t0\t1"));
    assert!(!first.lines().any(|line| line == "c\td\t0\t0\t0\t0\t1"));
    assert!(!first.lines().any(|line| line == "e\tf\t0\t0\t0\t0\t1"));
    let mut sorted = first.lines().collect::<Vec<_>>();
    sorted.sort();
    assert_eq!(first.lines().collect::<Vec<_>>(), sorted);
}

#[test]
fn partial_whitelist_applies_to_direct_and_inferred_pairs() {
    let output = run_kprune_filtered(&["a"], true);

    assert!(output.lines().all(|line| {
        let mut fields = line.split('\t');
        fields.next() == Some("a") || fields.next() == Some("a")
    }));
}

#[test]
fn cli_parses_finite_non_negative_prune_thresholds() {
    let matches = cli::cli()
        .try_get_matches_from(["cphasing", "kprune", "a", "c", "o"])
        .unwrap();
    let (_, subcommand) = matches.subcommand().unwrap();
    assert_eq!(subcommand.get_one::<f64>("MIN_CONTACTS"), Some(&1.0));
    assert_eq!(subcommand.get_one::<f64>("MIN_MARGIN"), Some(&0.10));

    for argument in ["--min-contacts=-1", "--min-margin=NaN", "--min-margin=inf"] {
        assert!(
            cli::cli()
                .try_get_matches_from(["cphasing", "kprune", "a", "c", "o", argument])
                .is_err()
        );
    }
}
