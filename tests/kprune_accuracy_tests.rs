use cphasing::cli;
use cphasing::contacts::{ContactEvidence, ContactRecord, Contacts2};
use cphasing::kprune::{
    GroupDecision, GroupVote, KPruner, PruneThresholds, StabilityConfig,
    all_informative_groups_cross_allelic, classify_candidate_by_bidirectional_gap,
    classify_candidate_by_global_gap, classify_candidate_by_stability,
    classify_selected_alternative, group_decisions_indicate_cross_allelic,
    solve_maximum_bipartite_matching,
};
use ordered_float::OrderedFloat;
use pathfinding::matrix::Matrix;
use pathfinding::prelude::kuhn_munkres;
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
fn two_row_matching_matches_kuhn_munkres_for_scores_and_ties() {
    for columns in 2..=5 {
        let cases = 3usize.pow((2 * columns) as u32);
        for mut encoded in 0..cases {
            let matrix = Matrix::from_fn(2, columns, |_| {
                let value = encoded % 3;
                encoded /= 3;
                OrderedFloat(value as f64)
            });

            let actual = solve_maximum_bipartite_matching(&matrix);
            let expected = kuhn_munkres(&matrix);
            assert_eq!(actual, expected, "matrix={matrix:?}");
        }
    }
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
fn global_gap_detects_candidate_displaced_by_stronger_matching() {
    let matrix = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(0.0)],
        vec![OrderedFloat(10.0), OrderedFloat(0.0)],
    ])
    .unwrap();

    assert_eq!(
        classify_candidate_by_global_gap(&matrix, 0, 0, 5.0, PruneThresholds::default(),),
        GroupDecision::CrossAllelic,
    );
}

#[test]
fn global_gap_retains_candidate_in_an_optimal_matching() {
    let matrix = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(0.0)],
        vec![OrderedFloat(0.0), OrderedFloat(4.0)],
    ])
    .unwrap();

    assert_eq!(
        classify_candidate_by_global_gap(&matrix, 0, 0, 5.0, PruneThresholds::default(),),
        GroupDecision::Compatible,
    );
}

#[test]
fn bidirectional_gap_reuses_global_gap_for_unselected_candidate() {
    let normalized = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(0.0)],
        vec![OrderedFloat(10.0), OrderedFloat(0.0)],
    ])
    .unwrap();
    let raw = normalized.clone();

    assert_eq!(
        classify_candidate_by_bidirectional_gap(
            &normalized,
            &raw,
            0,
            0,
            5.0,
            PruneThresholds::default(),
            0.50,
        ),
        GroupDecision::CrossAllelic,
    );
}

#[test]
fn bidirectional_gap_recovers_weak_selected_candidate_rejected_by_raw_contacts() {
    let normalized = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(4.9)],
        vec![OrderedFloat(4.9), OrderedFloat(5.0)],
    ])
    .unwrap();
    let raw = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(10.0)],
        vec![OrderedFloat(10.0), OrderedFloat(5.0)],
    ])
    .unwrap();

    assert_eq!(
        classify_candidate_by_bidirectional_gap(
            &normalized,
            &raw,
            0,
            0,
            5.0,
            PruneThresholds::default(),
            0.05,
        ),
        GroupDecision::CrossAllelic,
    );
}

#[test]
fn bidirectional_gap_keeps_selected_candidate_supported_by_raw_contacts() {
    let normalized = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(4.9)],
        vec![OrderedFloat(4.9), OrderedFloat(5.0)],
    ])
    .unwrap();
    let raw = Matrix::from_rows(vec![
        vec![OrderedFloat(10.0), OrderedFloat(5.0)],
        vec![OrderedFloat(5.0), OrderedFloat(10.0)],
    ])
    .unwrap();

    assert_eq!(
        classify_candidate_by_bidirectional_gap(
            &normalized,
            &raw,
            0,
            0,
            10.0,
            PruneThresholds::default(),
            0.02,
        ),
        GroupDecision::Compatible,
    );
}

#[test]
fn bidirectional_gap_is_not_diluted_by_unrelated_large_assignments() {
    let normalized = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(4.9), OrderedFloat(0.0)],
        vec![OrderedFloat(4.9), OrderedFloat(5.0), OrderedFloat(0.0)],
        vec![OrderedFloat(0.0), OrderedFloat(0.0), OrderedFloat(100.0)],
    ])
    .unwrap();
    let raw = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(10.0), OrderedFloat(0.0)],
        vec![OrderedFloat(10.0), OrderedFloat(5.0), OrderedFloat(0.0)],
        vec![OrderedFloat(0.0), OrderedFloat(0.0), OrderedFloat(100.0)],
    ])
    .unwrap();

    assert_eq!(
        classify_candidate_by_bidirectional_gap(
            &normalized,
            &raw,
            0,
            0,
            5.0,
            PruneThresholds::default(),
            0.02,
        ),
        GroupDecision::Compatible,
    );
}

#[test]
fn stability_vote_detects_consistent_candidate_displacement() {
    let matrix = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(0.0)],
        vec![OrderedFloat(10.0), OrderedFloat(0.0)],
    ])
    .unwrap();
    let config = StabilityConfig {
        replicates: 8,
        jitter: 0.1,
        vote_threshold: 0.5,
    };

    assert_eq!(
        classify_candidate_by_stability(&matrix, 0, 0, 5.0, PruneThresholds::default(), config,),
        GroupDecision::CrossAllelic,
    );
}

#[test]
fn stability_vote_is_deterministic_for_stable_candidate() {
    let matrix = Matrix::from_rows(vec![
        vec![OrderedFloat(5.0), OrderedFloat(0.0)],
        vec![OrderedFloat(0.0), OrderedFloat(4.0)],
    ])
    .unwrap();
    let config = StabilityConfig {
        replicates: 8,
        jitter: 0.2,
        vote_threshold: 0.5,
    };

    let first =
        classify_candidate_by_stability(&matrix, 0, 0, 5.0, PruneThresholds::default(), config);
    let second =
        classify_candidate_by_stability(&matrix, 0, 0, 5.0, PruneThresholds::default(), config);

    assert_eq!(first, GroupDecision::Compatible);
    assert_eq!(first, second);
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

#[test]
fn any_cross_vote_prunes_despite_compatible_groups() {
    let decisions = [GroupDecision::CrossAllelic, GroupDecision::Compatible];

    assert!(!group_decisions_indicate_cross_allelic(
        decisions,
        GroupVote::Conservative,
    ));
    assert!(group_decisions_indicate_cross_allelic(
        decisions,
        GroupVote::AnyCross,
    ));
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
    assert_eq!(subcommand.get_one::<f64>("MIN_CONTACTS"), Some(&0.0));
    assert_eq!(subcommand.get_one::<f64>("MIN_MARGIN"), Some(&0.0));

    for argument in ["--min-contacts=-1", "--min-margin=NaN", "--min-margin=inf"] {
        assert!(
            cli::cli()
                .try_get_matches_from(["cphasing", "kprune", "a", "c", "o", argument])
                .is_err()
        );
    }
}

#[test]
fn cli_parses_group_vote_strategy() {
    let matches = cli::cli()
        .try_get_matches_from(["cphasing", "kprune", "a", "c", "o"])
        .unwrap();
    let (_, subcommand) = matches.subcommand().unwrap();
    assert_eq!(
        subcommand.get_one::<String>("GROUP_VOTE"),
        Some(&"any-cross".to_string()),
    );

    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "kprune",
            "a",
            "c",
            "o",
            "--group-vote",
            "any-cross",
        ])
        .unwrap();
    let (_, subcommand) = matches.subcommand().unwrap();
    assert_eq!(
        subcommand.get_one::<String>("GROUP_VOTE"),
        Some(&"any-cross".to_string()),
    );
}

#[test]
fn cli_parses_decision_mode() {
    let matches = cli::cli()
        .try_get_matches_from(["cphasing", "kprune", "a", "c", "o"])
        .unwrap();
    let (_, subcommand) = matches.subcommand().unwrap();
    assert_eq!(
        subcommand.get_one::<String>("DECISION_MODE"),
        Some(&"bidirectional-gap".to_string()),
    );
    assert_eq!(subcommand.get_one::<f64>("MAX_SELECTED_GAP"), Some(&0.50),);

    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "kprune",
            "a",
            "c",
            "o",
            "--decision-mode",
            "global-gap",
        ])
        .unwrap();
    let (_, subcommand) = matches.subcommand().unwrap();
    assert_eq!(
        subcommand.get_one::<String>("DECISION_MODE"),
        Some(&"global-gap".to_string()),
    );

    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "kprune",
            "a",
            "c",
            "o",
            "--decision-mode",
            "bidirectional-gap",
            "--max-selected-gap",
            "0.02",
        ])
        .unwrap();
    let (_, subcommand) = matches.subcommand().unwrap();
    assert_eq!(
        subcommand.get_one::<String>("DECISION_MODE"),
        Some(&"bidirectional-gap".to_string()),
    );
    assert_eq!(subcommand.get_one::<f64>("MAX_SELECTED_GAP"), Some(&0.02),);
}

#[test]
fn cli_parses_stability_voting_parameters() {
    let matches = cli::cli()
        .try_get_matches_from([
            "cphasing",
            "kprune",
            "a",
            "c",
            "o",
            "--decision-mode",
            "stability",
        ])
        .unwrap();
    let (_, subcommand) = matches.subcommand().unwrap();

    assert_eq!(
        subcommand.get_one::<String>("DECISION_MODE"),
        Some(&"stability".to_string()),
    );
    assert_eq!(
        subcommand.get_one::<usize>("STABILITY_REPLICATES"),
        Some(&16),
    );
    assert_eq!(subcommand.get_one::<f64>("STABILITY_JITTER"), Some(&0.10),);
    assert_eq!(
        subcommand.get_one::<f64>("STABILITY_VOTE_THRESHOLD"),
        Some(&0.50),
    );
}
