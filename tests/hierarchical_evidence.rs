use cphasing::core::ContigPair;
use cphasing::optimize::{BackboneConfig, BackbonePlan};
use cphasing::order::{Tour, run_hierarchical_end_initializer_with_evidence};
use cphasing::splitcontacts::SplitContacts;
use hashbrown::HashMap;
use indexmap::IndexMap;
use std::collections::HashSet;

#[test]
fn hierarchical_evidence_preserves_exact_raw_support_at_threshold() {
    // The half-length product is 1.5 * 50.5 = 75.75. Reconstructing the
    // support from (5 / 75.75) * 75.75 rounds below the exact threshold.
    let lengths = IndexMap::from([(0, 3), (1, 101)]);
    let names = HashMap::from([("a".to_string(), 0), ("b".to_string(), 1)]);
    let contacts = SplitContacts {
        file: "raw-support-boundary.split.contacts".to_string(),
        contigs: HashSet::from(["a".to_string(), "b".to_string()]),
        data: HashMap::from([(
            ContigPair::new("a".to_string(), "b".to_string()),
            vec![0.0, 0.0, 5.0, 0.0],
        )]),
    };
    let initial = Tour {
        contigs: vec![0, 1],
        signs: vec![true, true],
    };

    let result =
        run_hierarchical_end_initializer_with_evidence(&initial, &lengths, &contacts, &names);

    assert_eq!(result.joins.len(), 1);
    let join = result.joins[0];
    assert_eq!(join.raw_support.to_bits(), 5.0_f64.to_bits());
    assert_ne!(
        (join.normalized_score * 75.75).to_bits(),
        5.0_f64.to_bits(),
        "the fixture must expose the lossy normalized-score reconstruction"
    );

    let plan = BackbonePlan::from_hierarchical_joins(
        &result.tour.contigs,
        &result.joins,
        &BackboneConfig {
            min_reduction_fraction: 0.0,
            ..BackboneConfig::default()
        },
    )
    .unwrap();
    assert_eq!(plan.accepted_edges(), 1);
}

#[test]
fn batched_hierarchical_merges_emit_one_exact_final_adjacency_each() {
    let n = 513usize;
    let lengths = (0..n).map(|id| (id, 100)).collect::<IndexMap<_, _>>();
    let names = (0..n)
        .map(|id| (format!("c{id}"), id))
        .collect::<HashMap<_, _>>();
    let data = (0..n - 1)
        .map(|id| {
            (
                ContigPair::new(format!("c{id}"), format!("c{}", id + 1)),
                vec![0.0, 0.0, 10.0, 0.0],
            )
        })
        .collect::<HashMap<_, _>>();
    let contacts = SplitContacts {
        file: "batched-chain.split.contacts".to_string(),
        contigs: (0..n).map(|id| format!("c{id}")).collect(),
        data,
    };
    let initial = Tour {
        contigs: (0..n).collect(),
        signs: vec![true; n],
    };

    let result =
        run_hierarchical_end_initializer_with_evidence(&initial, &lengths, &contacts, &names);

    assert_eq!(result.joins.len(), n - 1);
    assert!(
        result
            .joins
            .iter()
            .all(|join| join.raw_support.to_bits() == 10.0_f64.to_bits())
    );
    let final_adjacencies = result
        .tour
        .contigs
        .windows(2)
        .map(|pair| {
            if pair[0] < pair[1] {
                (pair[0], pair[1])
            } else {
                (pair[1], pair[0])
            }
        })
        .collect::<HashSet<_>>();
    assert!(
        result
            .joins
            .iter()
            .all(|join| final_adjacencies.contains(&(join.left_contig, join.right_contig)))
    );
}
