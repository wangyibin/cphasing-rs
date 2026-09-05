use cphasing::clm::{ClmbWriter, encode_endpoint};
use cphasing::core::ContigPair;
use cphasing::optimize::{
    AllhicProblem, BackboneConfig, BackbonePlan, BandedOrientationConfig, OptimizeConfig,
    OptimizeProblem, OrderObjective, OrientationGapModel, OrientationPairWeight,
    OrientationProblem, OrientedDistanceRecord, SignedBlockRefineConfig, UnitGene, mutate_allhic,
    optimize_order, optimize_order_with_hierarchical_joins,
};
use cphasing::order::{
    HierarchicalJoinEvidence, Tour, run_beam_end_initializer, run_hierarchical_end_initializer,
    run_hierarchical_end_initializer_with_evidence,
};
use cphasing::splitcontacts::SplitContacts;
use hashbrown::HashMap;
use indexmap::IndexMap;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rayon::ThreadPoolBuilder;
use std::collections::HashSet;

fn problem() -> OptimizeProblem {
    let lengths = IndexMap::from([(0, 100), (1, 200), (2, 300), (3, 400)]);
    let mut contacts = vec![HashMap::new(); 4];
    for (u, v, links) in [(0, 2, 5), (2, 3, 7), (0, 1, 11), (1, 3, 13)] {
        contacts[u].insert(v, links);
        contacts[v].insert(u, links);
    }
    OptimizeProblem::from_sparse_contacts(&lengths, &contacts).unwrap()
}

fn path_block_problem() -> OptimizeProblem {
    let lengths = IndexMap::from([(0, 100), (1, 100), (2, 100), (3, 100), (4, 100), (5, 100)]);
    let mut contacts = vec![HashMap::new(); lengths.len()];
    for (u, v, links) in [(0, 1, 100), (1, 2, 90), (3, 4, 80)] {
        contacts[u].insert(v, links);
        contacts[v].insert(u, links);
    }
    OptimizeProblem::from_sparse_contacts(&lengths, &contacts).unwrap()
}

fn compact_config(backbone: BackboneConfig) -> OptimizeConfig {
    OptimizeConfig {
        population_size: 12,
        stale_generations: 8,
        max_generations: 40,
        mutation_probability: 0.5,
        seed: 42,
        phases: 2,
        report_interval: 0,
        backbone,
    }
}

fn initial_tour(contig_count: usize) -> Tour<usize> {
    Tour {
        contigs: (0..contig_count).rev().collect(),
        signs: (0..contig_count).map(|index| index % 2 == 0).collect(),
    }
}

fn hierarchical_join(
    left_contig: usize,
    right_contig: usize,
    raw_support: f64,
    confidence: f64,
    reciprocal_confident: bool,
) -> HierarchicalJoinEvidence {
    HierarchicalJoinEvidence {
        left_contig,
        right_contig,
        normalized_score: 1.0,
        raw_support,
        confidence,
        reciprocal_confident,
    }
}

#[test]
fn seriation_seed_improves_a_scrambled_path_and_is_deterministic() {
    let problem = fragmented_path_problem(24);
    let fallback = (0..24)
        .step_by(2)
        .chain((1..24).step_by(2))
        .collect::<Vec<_>>();

    let first = problem.seriation_seed(&fallback).unwrap();
    let second = problem.seriation_seed(&fallback).unwrap();

    assert_eq!(first.order, second.order);
    assert_eq!(first.final_score.to_bits(), second.final_score.to_bits());
    assert!(first.final_score < first.fallback_score);
    assert_eq!(first.order.len(), 24);
    let mut sorted = first.order.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, (0..24).collect::<Vec<_>>());
}

fn fragmented_path_problem(contig_count: usize) -> OptimizeProblem {
    let lengths = (0..contig_count)
        .map(|id| (id, 100))
        .collect::<IndexMap<_, _>>();
    let mut contacts = vec![HashMap::new(); contig_count];
    for u in 0..contig_count {
        for v in (u + 1)..contig_count {
            let links = if v == u + 1 { 10_000 } else { 1 };
            contacts[u].insert(v, links);
            contacts[v].insert(u, links);
        }
    }
    OptimizeProblem::from_sparse_contacts(&lengths, &contacts).unwrap()
}

#[test]
fn hierarchical_end_initializer_recovers_a_supported_path_without_ga() {
    let lengths = IndexMap::from([(0, 100), (1, 100), (2, 100)]);
    let names = HashMap::from([
        ("a".to_string(), 0),
        ("b".to_string(), 1),
        ("c".to_string(), 2),
    ]);
    let contacts = SplitContacts {
        file: "synthetic.split.contacts".to_string(),
        contigs: HashSet::from(["a".to_string(), "b".to_string(), "c".to_string()]),
        data: HashMap::from([
            (
                ContigPair::new("a".to_string(), "b".to_string()),
                vec![0.0, 0.0, 100.0, 0.0],
            ),
            (
                ContigPair::new("b".to_string(), "c".to_string()),
                vec![0.0, 0.0, 100.0, 0.0],
            ),
            (
                ContigPair::new("a".to_string(), "c".to_string()),
                vec![1.0; 4],
            ),
        ]),
    };
    let initial = Tour {
        contigs: vec![2, 0, 1],
        signs: vec![true; 3],
    };

    let legacy = run_hierarchical_end_initializer(&initial, &lengths, &contacts, &names);
    let result =
        run_hierarchical_end_initializer_with_evidence(&initial, &lengths, &contacts, &names);
    let adjacencies = result
        .tour
        .contigs
        .windows(2)
        .map(|pair| (pair[0].min(pair[1]), pair[0].max(pair[1])))
        .collect::<HashSet<_>>();

    assert_eq!(adjacencies, HashSet::from([(0, 1), (1, 2)]));
    assert_eq!(result.tour.contigs, legacy.contigs);
    assert_eq!(result.tour.signs, legacy.signs);
    assert_eq!(result.tour.contigs.len(), result.tour.signs.len());
    assert_eq!(result.joins.len(), result.tour.contigs.len() - 1);
    for join in &result.joins {
        assert!(join.left_contig < join.right_contig);
        assert!(adjacencies.contains(&(join.left_contig, join.right_contig)));
        assert!(join.normalized_score > 0.0);
        assert!(join.raw_support > 0.0);
        assert!(join.confidence >= 1.0);
        assert!(join.reciprocal_confident);
    }
}

#[test]
fn beam_end_initializer_escapes_a_greedy_endpoint_trap() {
    let names = HashMap::from([
        ("a".to_string(), 0),
        ("b".to_string(), 1),
        ("c".to_string(), 2),
        ("d".to_string(), 3),
    ]);
    let contacts = SplitContacts {
        file: "synthetic.split.contacts".to_string(),
        contigs: HashSet::from([
            "a".to_string(),
            "b".to_string(),
            "c".to_string(),
            "d".to_string(),
        ]),
        data: HashMap::from([
            (
                ContigPair::new("a".to_string(), "b".to_string()),
                vec![0.0, 0.0, 10.0, 0.0],
            ),
            (
                ContigPair::new("a".to_string(), "c".to_string()),
                vec![0.0, 0.0, 6.0, 0.0],
            ),
            (
                ContigPair::new("b".to_string(), "d".to_string()),
                vec![0.0, 6.0, 0.0, 0.0],
            ),
            (
                ContigPair::new("c".to_string(), "d".to_string()),
                vec![0.0, 0.0, 1.0, 0.0],
            ),
        ]),
    };
    let initial = Tour {
        contigs: vec![0, 1, 2, 3],
        signs: vec![true; 4],
    };

    let result = run_beam_end_initializer(&initial, &contacts, &names, 16);
    let adjacencies = result
        .contigs
        .windows(2)
        .map(|pair| (pair[0].min(pair[1]), pair[0].max(pair[1])))
        .collect::<HashSet<_>>();

    assert_eq!(adjacencies, HashSet::from([(0, 2), (1, 3), (2, 3)]));
    assert_eq!(result.contigs.len(), result.signs.len());
}

fn assert_same_optimization(
    left: &cphasing::optimize::OptimizeResult,
    right: &cphasing::optimize::OptimizeResult,
) {
    assert_eq!(left.tour.contigs, right.tour.contigs);
    assert_eq!(left.tour.signs, right.tour.signs);
    assert_eq!(left.initial_score.to_bits(), right.initial_score.to_bits());
    assert_eq!(left.final_score.to_bits(), right.final_score.to_bits());
    assert_eq!(left.generations, right.generations);
    assert_eq!(left.reports.len(), right.reports.len());
    for (left_report, right_report) in left.reports.iter().zip(&right.reports) {
        assert_eq!(left_report.phase, right_report.phase);
        assert_eq!(left_report.generation, right_report.generation);
        assert_eq!(
            left_report.best_score.to_bits(),
            right_report.best_score.to_bits()
        );
    }

    match (&left.backbone, &right.backbone) {
        (Some(left_backbone), Some(right_backbone)) => {
            assert_eq!(left_backbone.used, right_backbone.used);
            assert_eq!(left_backbone.contig_count, right_backbone.contig_count);
            assert_eq!(left_backbone.unit_count, right_backbone.unit_count);
            assert_eq!(left_backbone.block_count, right_backbone.block_count);
            assert_eq!(left_backbone.accepted_edges, right_backbone.accepted_edges);
            assert_eq!(
                left_backbone.seed_score.map(f64::to_bits),
                right_backbone.seed_score.map(f64::to_bits)
            );
            assert_eq!(
                left_backbone.coarse_score.map(f64::to_bits),
                right_backbone.coarse_score.map(f64::to_bits)
            );
        }
        (None, None) => {}
        _ => panic!("backbone reports differ"),
    }
}

#[test]
fn score_matches_allhic_dense_reference() {
    let problem = problem().with_objective(OrderObjective::LogDistance);
    let expected =
        5.0 * 200.0_f64.ln() + 7.0 * 450.0_f64.ln() + 11.0 * 550.0_f64.ln() + 13.0 * 300.0_f64.ln();
    assert!((problem.evaluate(&[2, 0, 3, 1]) - expected).abs() < 1e-10);
}

#[test]
fn default_score_matches_allhic_reciprocal_reference() {
    let problem = problem();
    let expected = -5.0 / 200.0 - 7.0 / 450.0 - 11.0 / 550.0 - 13.0 / 300.0;
    assert!((problem.evaluate(&[2, 0, 3, 1]) - expected).abs() < 1e-10);
}

#[test]
fn length_tiered_objective_prioritizes_the_large_contig_order() {
    let lengths = IndexMap::from([
        (0, 400),
        (1, 300),
        (2, 200),
        (3, 100),
        (4, 10),
        (5, 10),
        (6, 10),
        (7, 10),
    ]);
    let mut contacts = vec![HashMap::new(); lengths.len()];
    for (u, v, links) in [
        (0, 1, 1_000),
        (1, 2, 1_000),
        (2, 3, 1_000),
        (3, 4, 200),
        (4, 5, 100),
        (5, 6, 100),
        (6, 7, 100),
    ] {
        contacts[u].insert(v, links);
        contacts[v].insert(u, links);
    }
    let problem = OptimizeProblem::from_sparse_contacts(&lengths, &contacts)
        .unwrap()
        .with_objective(OrderObjective::LengthTiered);
    let good = [0, 1, 2, 3, 4, 5, 6, 7];
    let bad_anchor_order = [0, 2, 1, 3, 4, 5, 6, 7];
    let reordered_fragments = [0, 1, 2, 3, 7, 6, 5, 4];

    let good_tiers = problem.evaluate_length_tiers(&good);
    let bad_tiers = problem.evaluate_length_tiers(&bad_anchor_order);
    let fragment_tiers = problem.evaluate_length_tiers(&reordered_fragments);

    assert!(good_tiers[0] < bad_tiers[0]);
    assert_eq!(good_tiers[0].to_bits(), fragment_tiers[0].to_bits());
    assert!(problem.evaluate(&good) < problem.evaluate(&bad_anchor_order));
}

#[test]
fn endpoint_multiscale_objective_rewards_the_signed_anchor_path() {
    let lengths = IndexMap::from([
        (0, 400),
        (1, 300),
        (2, 200),
        (3, 100),
        (4, 10),
        (5, 10),
        (6, 10),
        (7, 10),
    ]);
    let contacts = vec![HashMap::new(); lengths.len()];
    let split_contacts = SplitContacts {
        file: "synthetic.split.contacts".to_string(),
        contigs: (0..8).map(|id| format!("c{id}")).collect(),
        data: HashMap::from([
            (
                ContigPair::new("c0".to_string(), "c1".to_string()),
                vec![0.0, 0.0, 100.0, 0.0],
            ),
            (
                ContigPair::new("c1".to_string(), "c2".to_string()),
                vec![0.0, 0.0, 100.0, 0.0],
            ),
            (
                ContigPair::new("c2".to_string(), "c3".to_string()),
                vec![0.0, 0.0, 100.0, 0.0],
            ),
        ]),
    };
    let contig_to_id = (0..8)
        .map(|id| (format!("c{id}"), id))
        .collect::<HashMap<_, _>>();
    let tour = Tour {
        contigs: (0..8).collect(),
        signs: vec![true; 8],
    };
    let problem = OptimizeProblem::from_sparse_contacts(&lengths, &contacts)
        .unwrap()
        .with_endpoint_multiscale(&split_contacts, &contig_to_id, &tour)
        .unwrap();

    assert!(problem.evaluate(&tour.contigs) < problem.evaluate(&[0, 2, 1, 3, 4, 5, 6, 7]));
}

#[test]
fn every_mutation_preserves_the_permutation() {
    let mut rng = SmallRng::seed_from_u64(42);
    let expected: Vec<usize> = (0..50).collect();
    let mut order = expected.clone();
    for _ in 0..10_000 {
        mutate_allhic(&mut order, &mut rng);
        let mut observed = order.clone();
        observed.sort_unstable();
        assert_eq!(observed, expected);
    }
}

#[test]
fn ga_keeps_a_result_no_worse_than_the_initial_tour() {
    let problem = problem();
    let initial = Tour {
        contigs: vec![2, 0, 3, 1],
        signs: vec![true, false, true, false],
    };
    let config = OptimizeConfig {
        population_size: 20,
        stale_generations: 20,
        max_generations: 200,
        seed: 42,
        phases: 2,
        report_interval: 0,
        ..OptimizeConfig::default()
    };
    let result = optimize_order(&initial, &problem, &config).unwrap();

    assert!(result.final_score <= result.initial_score);
    let mut contigs = result.tour.contigs.clone();
    contigs.sort_unstable();
    assert_eq!(contigs, vec![0, 1, 2, 3]);
}

#[test]
fn singleton_units_decode_with_exact_full_problem_score() {
    let problem = problem();
    let backbone = BackboneConfig {
        min_links: 1_000.0,
        ..BackboneConfig::default()
    };
    let plan = BackbonePlan::discover(&problem, &backbone);

    assert!(!plan.is_useful());
    assert_eq!(plan.block_count(), 0);
    assert_eq!(plan.unit_count(), problem.contig_count());

    let expected = vec![3, 0, 2, 1];
    let genes = expected
        .iter()
        .map(|&unit| UnitGene {
            unit,
            reversed: true,
        })
        .collect::<Vec<_>>();
    let decoded = plan.decode(&genes).unwrap();

    assert_eq!(decoded, expected);
    assert_eq!(
        problem.evaluate(&decoded).to_bits(),
        problem.evaluate(&expected).to_bits()
    );
}

#[test]
fn hierarchical_joins_form_seed_contiguous_blocks_and_cut_weak_edges() {
    let seed = vec![3, 1, 4, 0, 2, 5];
    let joins = vec![
        hierarchical_join(1, 3, 20.0, 2.0, true),
        hierarchical_join(1, 4, 10.0, 1.5, true),
        // Enough confidence but not enough raw support.
        hierarchical_join(0, 4, 4.0, 3.0, true),
        // Enough support but not enough confidence.
        hierarchical_join(0, 2, 20.0, 1.09, true),
        // Both numeric thresholds pass, but the join was not reciprocal.
        hierarchical_join(2, 5, 20.0, 3.0, false),
        // Trusted but non-adjacent evidence must not reorder the seed.
        hierarchical_join(3, 4, 20.0, 3.0, true),
    ];
    let config = BackboneConfig {
        min_reduction_fraction: 0.0,
        ..BackboneConfig::default()
    };

    let plan = BackbonePlan::from_hierarchical_joins(&seed, &joins, &config).unwrap();

    assert_eq!(plan.units(), &[vec![3, 1, 4], vec![0], vec![2], vec![5]]);
    assert_eq!(plan.block_count(), 1);
    assert_eq!(plan.accepted_edges(), 2);
    assert!(plan.is_useful());
    let genes = (0..plan.unit_count())
        .map(|unit| UnitGene {
            unit,
            reversed: false,
        })
        .collect::<Vec<_>>();
    assert_eq!(plan.decode(&genes).unwrap(), seed);
}

#[test]
fn hierarchical_confidence_threshold_is_independent_of_clm_margin() {
    let seed = vec![0, 1, 2];
    let joins = vec![
        hierarchical_join(0, 1, 10.0, 1.1, true),
        hierarchical_join(1, 2, 10.0, 1.099, true),
    ];
    let config = BackboneConfig {
        min_margin: 10.0,
        min_reduction_fraction: 0.0,
        ..BackboneConfig::default()
    };

    let plan = BackbonePlan::from_hierarchical_joins(&seed, &joins, &config).unwrap();

    assert_eq!(plan.units(), &[vec![0, 1], vec![2]]);
    assert_eq!(plan.accepted_edges(), 1);
}

#[test]
fn hierarchical_blocks_respect_maximum_size_and_preserve_the_full_permutation() {
    let seed = vec![6, 5, 4, 3, 2, 1, 0];
    let joins = seed
        .windows(2)
        .map(|pair| hierarchical_join(pair[1], pair[0], 20.0, 3.0, true))
        .collect::<Vec<_>>();
    let config = BackboneConfig {
        min_reduction_fraction: 0.0,
        max_block_size: 3,
        ..BackboneConfig::default()
    };

    let first = BackbonePlan::from_hierarchical_joins(&seed, &joins, &config).unwrap();
    let mut reversed_joins = joins.clone();
    reversed_joins.reverse();
    let second = BackbonePlan::from_hierarchical_joins(&seed, &reversed_joins, &config).unwrap();

    assert_eq!(first.units(), &[vec![6, 5, 4], vec![3, 2, 1], vec![0]]);
    assert_eq!(first.units(), second.units());
    assert_eq!(first.accepted_edges(), 4);
    assert_eq!(first.accepted_edges(), second.accepted_edges());
    assert!(first.units().iter().all(|unit| unit.len() <= 3));
    let flattened = first.units().iter().flatten().copied().collect::<Vec<_>>();
    assert_eq!(flattened, seed);

    let reverse_genes = (0..first.unit_count())
        .rev()
        .map(|unit| UnitGene {
            unit,
            reversed: true,
        })
        .collect::<Vec<_>>();
    assert_eq!(
        first.decode(&reverse_genes).unwrap(),
        seed.into_iter().rev().collect::<Vec<_>>()
    );
}

#[test]
fn refining_backbone_splits_only_within_parent_blocks() {
    let seed = (0..20).collect::<Vec<_>>();
    let joins = seed
        .windows(2)
        .filter(|pair| pair[0] != 9)
        .map(|pair| hierarchical_join(pair[0], pair[1], 20.0, 3.0, true))
        .collect::<Vec<_>>();
    let config = BackboneConfig {
        min_reduction_fraction: 0.0,
        max_block_size: 32,
        ..BackboneConfig::default()
    };
    let coarse = BackbonePlan::from_hierarchical_joins(&seed, &joins, &config).unwrap();
    let refined = coarse.refined(8).unwrap();

    assert_eq!(
        refined.units(),
        &[
            (0..8).collect::<Vec<_>>(),
            vec![8, 9],
            (10..18).collect::<Vec<_>>(),
            vec![18, 19],
        ]
    );
    assert_eq!(refined.accepted_edges(), 16);
    assert!(refined.units().iter().all(|unit| unit.len() <= 8));
    let reverse_genes = (0..refined.unit_count())
        .rev()
        .map(|unit| UnitGene {
            unit,
            reversed: true,
        })
        .collect::<Vec<_>>();
    assert_eq!(
        refined.decode(&reverse_genes).unwrap(),
        seed.into_iter().rev().collect::<Vec<_>>()
    );
    assert!(coarse.refined(0).is_err());
    assert_eq!(refined.refined(8).unwrap().units(), refined.units());
    assert_eq!(refined.refined(32).unwrap().units(), refined.units());
}

#[test]
fn path_block_decode_is_always_a_complete_permutation() {
    let problem = path_block_problem();
    let plan = BackbonePlan::discover(&problem, &BackboneConfig::default());

    assert!(plan.is_useful());
    assert_eq!(plan.block_count(), 2);
    assert_eq!(plan.accepted_edges(), 3);
    assert_eq!(plan.unit_count(), 3);

    let forward_genes = (0..plan.unit_count())
        .map(|unit| UnitGene {
            unit,
            reversed: false,
        })
        .collect::<Vec<_>>();
    let forward = plan.decode(&forward_genes).unwrap();
    let reverse_genes = (0..plan.unit_count())
        .rev()
        .map(|unit| UnitGene {
            unit,
            reversed: true,
        })
        .collect::<Vec<_>>();
    let reverse = plan.decode(&reverse_genes).unwrap();

    let mut observed = forward.clone();
    observed.sort_unstable();
    assert_eq!(observed, (0..problem.contig_count()).collect::<Vec<_>>());
    assert_eq!(
        reverse,
        forward.into_iter().rev().collect::<Vec<_>>(),
        "reversing the unit order and every unit must reverse the full tour"
    );
    assert!(
        plan.decode(&forward_genes[..plan.unit_count() - 1])
            .is_err()
    );
}

#[test]
fn weak_competitors_still_count_against_backbone_confidence() {
    let lengths = (0..13).map(|id| (id, 100)).collect::<IndexMap<_, _>>();
    let mut contacts = vec![HashMap::new(); lengths.len()];
    for (u, v, links) in [(0, 1, 5), (0, 2, 5)]
        .into_iter()
        .chain((3..13).map(|v| (0, v, 4)))
    {
        contacts[u].insert(v, links);
        contacts[v].insert(u, links);
    }
    let problem = OptimizeProblem::from_sparse_contacts(&lengths, &contacts).unwrap();
    let config = BackboneConfig {
        min_links: 5.0,
        min_margin: 1.0,
        min_top_two_fraction: 0.5,
        min_reduction_fraction: 0.0,
        ..BackboneConfig::default()
    };
    let plan = BackbonePlan::discover(&problem, &config);

    assert_eq!(plan.accepted_edges(), 0);
    assert!(!plan.is_useful());
}

#[test]
fn no_reliable_blocks_fall_back_to_the_standard_ga() {
    let problem = problem();
    let initial = initial_tour(problem.contig_count());
    let disabled = compact_config(BackboneConfig {
        enabled: false,
        ..BackboneConfig::default()
    });
    let strict = compact_config(BackboneConfig {
        enabled: true,
        min_links: 1_000.0,
        ..BackboneConfig::default()
    });

    let standard = optimize_order(&initial, &problem, &disabled).unwrap();
    let fallback = optimize_order(&initial, &problem, &strict).unwrap();
    let report = fallback.backbone.as_ref().unwrap();

    assert!(!report.used);
    assert_eq!(report.block_count, 0);
    assert_eq!(report.unit_count, problem.contig_count());
    assert_eq!(fallback.tour.contigs, standard.tour.contigs);
    assert_eq!(fallback.tour.signs, standard.tour.signs);
    assert_eq!(
        fallback.final_score.to_bits(),
        standard.final_score.to_bits()
    );
    assert_eq!(fallback.generations, standard.generations);
}

#[test]
fn unusable_hierarchical_evidence_falls_back_to_existing_clm_discovery() {
    let problem = path_block_problem();
    let initial = initial_tour(problem.contig_count());
    let config = compact_config(BackboneConfig::default());
    let weak_joins = initial
        .contigs
        .windows(2)
        .map(|pair| hierarchical_join(pair[0], pair[1], 1.0, 1.0, false))
        .collect::<Vec<_>>();

    let expected = optimize_order(&initial, &problem, &config).unwrap();
    let observed =
        optimize_order_with_hierarchical_joins(&initial, &problem, &config, &weak_joins).unwrap();

    assert!(observed.backbone.as_ref().unwrap().used);
    assert_same_optimization(&observed, &expected);
}

#[test]
fn hierarchical_wrapper_preserves_disabled_backbone_path_exactly() {
    let problem = path_block_problem();
    let initial = initial_tour(problem.contig_count());
    let joins = initial
        .contigs
        .windows(2)
        .map(|pair| hierarchical_join(pair[0], pair[1], 100.0, 3.0, true))
        .collect::<Vec<_>>();
    let config = compact_config(BackboneConfig {
        enabled: false,
        ..BackboneConfig::default()
    });

    let expected = optimize_order(&initial, &problem, &config).unwrap();
    let observed =
        optimize_order_with_hierarchical_joins(&initial, &problem, &config, &joins).unwrap();

    assert_same_optimization(&observed, &expected);
}

#[test]
fn hierarchical_wrapper_preserves_single_phase_path_exactly() {
    let problem = path_block_problem();
    let initial = initial_tour(problem.contig_count());
    let joins = initial
        .contigs
        .windows(2)
        .map(|pair| hierarchical_join(pair[0], pair[1], 100.0, 3.0, true))
        .collect::<Vec<_>>();
    let mut config = compact_config(BackboneConfig::default());
    config.phases = 1;

    let expected = optimize_order(&initial, &problem, &config).unwrap();
    let observed =
        optimize_order_with_hierarchical_joins(&initial, &problem, &config, &joins).unwrap();

    assert_same_optimization(&observed, &expected);
}

#[test]
fn hierarchical_backbone_optimization_is_deterministic_and_complete() {
    let problem = fragmented_path_problem(12);
    let initial = Tour {
        contigs: (0..12).collect(),
        signs: (0..12).map(|index| index % 2 == 0).collect(),
    };
    let joins = initial
        .contigs
        .windows(2)
        .map(|pair| hierarchical_join(pair[0], pair[1], 100.0, 3.0, true))
        .collect::<Vec<_>>();
    let config = compact_config(BackboneConfig {
        max_block_size: 4,
        ..BackboneConfig::default()
    });

    let first =
        optimize_order_with_hierarchical_joins(&initial, &problem, &config, &joins).unwrap();
    let second =
        optimize_order_with_hierarchical_joins(&initial, &problem, &config, &joins).unwrap();

    assert!(first.backbone.as_ref().unwrap().used);
    assert_eq!(first.backbone.as_ref().unwrap().block_count, 3);
    assert_same_optimization(&first, &second);
    let mut observed = first.tour.contigs.clone();
    observed.sort_unstable();
    assert_eq!(observed, (0..problem.contig_count()).collect::<Vec<_>>());
}

#[test]
fn hierarchical_progressive_stages_keep_two_macro_phase_budgets() {
    let problem = fragmented_path_problem(20);
    let initial = Tour {
        contigs: (0..20).collect(),
        signs: vec![true; 20],
    };
    let joins = initial
        .contigs
        .windows(2)
        .map(|pair| hierarchical_join(pair[0], pair[1], 100.0, 3.0, true))
        .collect::<Vec<_>>();
    let config = OptimizeConfig {
        population_size: 12,
        stale_generations: 100,
        max_generations: 3,
        mutation_probability: 1.0,
        seed: 42,
        phases: 2,
        report_interval: 1,
        backbone: BackboneConfig::default(),
    };

    let result =
        optimize_order_with_hierarchical_joins(&initial, &problem, &config, &joins).unwrap();

    assert!(result.backbone.as_ref().unwrap().used);
    assert_eq!(result.backbone.as_ref().unwrap().unit_count, 1);
    assert_eq!(result.generations.len(), 2);
    assert_eq!(result.generations[0], 6, "block32 and block8 share phase 1");
    assert_eq!(result.generations[1], 3, "block1 remains macro phase 2");
    let phase_one_reports = result
        .reports
        .iter()
        .filter(|report| report.phase == 1)
        .map(|report| report.generation)
        .collect::<Vec<_>>();
    let phase_two_reports = result
        .reports
        .iter()
        .filter(|report| report.phase == 2)
        .map(|report| report.generation)
        .collect::<Vec<_>>();
    assert_eq!(phase_one_reports, (1..=6).collect::<Vec<_>>());
    assert_eq!(phase_two_reports, (1..=3).collect::<Vec<_>>());
    let mut observed = result.tour.contigs.clone();
    observed.sort_unstable();
    assert_eq!(observed, (0..20).collect::<Vec<_>>());
}

#[test]
fn hierarchical_final_phase_can_break_a_wrong_locked_join() {
    let problem = problem();
    let initial = Tour {
        contigs: vec![2, 0, 3, 1],
        signs: vec![true; 4],
    };
    let joins = initial
        .contigs
        .windows(2)
        .map(|pair| hierarchical_join(pair[0], pair[1], 100.0, 3.0, true))
        .collect::<Vec<_>>();
    let config = OptimizeConfig {
        population_size: 20,
        stale_generations: 20,
        max_generations: 200,
        mutation_probability: 1.0,
        seed: 42,
        phases: 2,
        report_interval: 0,
        backbone: BackboneConfig::default(),
    };

    let result =
        optimize_order_with_hierarchical_joins(&initial, &problem, &config, &joins).unwrap();
    let backbone = result.backbone.as_ref().unwrap();

    assert_eq!(backbone.unit_count, 1);
    assert!(result.final_score < backbone.coarse_score.unwrap());
    assert_eq!(
        result.final_score.to_bits(),
        problem.evaluate(&result.tour.contigs).to_bits()
    );
}

#[test]
fn disabling_backbone_uses_the_unreported_standard_path() {
    let problem = path_block_problem();
    let initial = initial_tour(problem.contig_count());
    let config = compact_config(BackboneConfig {
        enabled: false,
        ..BackboneConfig::default()
    });

    let first = optimize_order(&initial, &problem, &config).unwrap();
    let second = optimize_order(&initial, &problem, &config).unwrap();

    assert!(first.backbone.is_none());
    assert_eq!(first.tour.contigs, second.tour.contigs);
    assert_eq!(first.tour.signs, second.tour.signs);
    assert_eq!(first.final_score.to_bits(), second.final_score.to_bits());
    assert_eq!(first.generations, second.generations);
}

#[test]
fn standard_ga_is_bitwise_identical_across_serial_and_parallel_scoring() {
    let problem = fragmented_path_problem(128);
    let initial = initial_tour(problem.contig_count());
    let config = OptimizeConfig {
        population_size: 12,
        stale_generations: 4,
        max_generations: 8,
        mutation_probability: 1.0,
        seed: 42,
        phases: 1,
        report_interval: 2,
        backbone: BackboneConfig {
            enabled: false,
            ..BackboneConfig::default()
        },
    };
    let serial_pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let parallel_pool = ThreadPoolBuilder::new().num_threads(4).build().unwrap();

    let serial = serial_pool
        .install(|| optimize_order(&initial, &problem, &config))
        .unwrap();
    let parallel = parallel_pool
        .install(|| optimize_order(&initial, &problem, &config))
        .unwrap();

    assert_same_optimization(&serial, &parallel);
}

#[test]
fn backbone_ga_is_bitwise_identical_across_serial_and_parallel_scoring() {
    let problem = fragmented_path_problem(2_048);
    let initial = initial_tour(problem.contig_count());
    let config = OptimizeConfig {
        population_size: 14,
        stale_generations: 4,
        max_generations: 8,
        mutation_probability: 1.0,
        seed: 42,
        phases: 2,
        report_interval: 2,
        backbone: BackboneConfig::default(),
    };
    let serial_pool = ThreadPoolBuilder::new().num_threads(1).build().unwrap();
    let parallel_pool = ThreadPoolBuilder::new().num_threads(4).build().unwrap();

    let serial = serial_pool
        .install(|| optimize_order(&initial, &problem, &config))
        .unwrap();
    let parallel = parallel_pool
        .install(|| optimize_order(&initial, &problem, &config))
        .unwrap();

    assert!(serial.backbone.as_ref().unwrap().used);
    assert!(parallel.backbone.as_ref().unwrap().used);
    assert_same_optimization(&serial, &parallel);
}

#[test]
fn partial_tours_are_rejected_before_evaluation() {
    let problem = path_block_problem();
    let partial = Tour {
        contigs: vec![0, 1, 2],
        signs: vec![true; 3],
    };
    let config = compact_config(BackboneConfig {
        enabled: false,
        ..BackboneConfig::default()
    });

    let error = optimize_order(&partial, &problem, &config).unwrap_err();
    assert!(error.contains("does not match problem contig count"));
    assert!(problem.evaluate(&partial.contigs).is_infinite());
}

#[test]
fn unlocked_final_ga_never_loses_the_backbone_stage_seed() {
    let problem = path_block_problem();
    let initial = initial_tour(problem.contig_count());
    let config = compact_config(BackboneConfig::default());
    let result = optimize_order(&initial, &problem, &config).unwrap();
    let report = result.backbone.as_ref().unwrap();

    assert!(report.used);
    assert!(result.final_score <= report.seed_score.unwrap());
    assert!(result.final_score <= report.coarse_score.unwrap());
    assert_eq!(
        result.final_score.to_bits(),
        problem.evaluate(&result.tour.contigs).to_bits()
    );
    let mut observed = result.tour.contigs.clone();
    observed.sort_unstable();
    assert_eq!(observed, (0..problem.contig_count()).collect::<Vec<_>>());
}

fn complete_orientation_pair(
    u: usize,
    v: usize,
    distances: [u64; 4],
    links: usize,
) -> Vec<OrientedDistanceRecord> {
    [(true, true), (true, false), (false, true), (false, false)]
        .into_iter()
        .zip(distances)
        .map(
            |((u_forward, v_forward), distance)| OrientedDistanceRecord {
                u,
                v,
                u_forward,
                v_forward,
                distances: vec![distance; links],
            },
        )
        .collect()
}

fn exact_banded_orientation_problem() -> (OrientationProblem, Vec<bool>) {
    let target = vec![true, false, true, false, true];
    let pairs = [
        (0, 1, 1),
        (1, 2, 4),
        (2, 3, 9),
        (3, 4, 16),
        (0, 2, 5),
        (1, 3, 7),
        (2, 4, 3),
        (0, 3, 2),
        (1, 4, 6),
    ];
    let records = pairs
        .into_iter()
        .enumerate()
        .flat_map(|(pair_index, (u, v, links))| {
            let preferred = (usize::from(!target[u]) << 1) | usize::from(!target[v]);
            let mut distances = [10_000, 30_000, 90_000, 270_000];
            distances[preferred] = 3_000 + pair_index as u64 * 100;
            complete_orientation_pair(u, v, distances, links)
        })
        .collect::<Vec<_>>();
    let problem =
        OrientationProblem::from_oriented_distances(vec![101, 1_009, 307, 701, 211], records)
            .unwrap();
    (problem, target)
}

type SourcePriorPairFixture = (usize, usize, Vec<OrientedDistanceRecord>);

fn source_prior_orientation_fixture() -> (Vec<u64>, Vec<SourcePriorPairFixture>, OrientationProblem)
{
    let lengths = vec![100, 1_000, 300];
    let pairs = vec![
        (
            0,
            1,
            complete_orientation_pair(0, 1, [30_000, 3_000, 90_000, 270_000], 1),
        ),
        (
            1,
            2,
            complete_orientation_pair(1, 2, [90_000, 30_000, 3_000, 270_000], 4),
        ),
        (
            0,
            2,
            complete_orientation_pair(0, 2, [3_000, 30_000, 90_000, 270_000], 2),
        ),
    ];
    let records = pairs
        .iter()
        .flat_map(|(_, _, records)| records.iter().cloned())
        .collect::<Vec<_>>();
    let problem = OrientationProblem::from_oriented_distances(lengths.clone(), records).unwrap();
    (lengths, pairs, problem)
}

#[test]
fn orientation_intervening_gap_excludes_endpoints_and_is_strictly_rc_invariant() {
    let problem = OrientationProblem::from_oriented_distances(
        vec![100, 1_000, 300],
        [OrientedDistanceRecord {
            u: 0,
            v: 2,
            u_forward: true,
            v_forward: true,
            distances: vec![3_000],
        }],
    )
    .unwrap();
    let direct = Tour {
        contigs: vec![0, 1, 2],
        signs: vec![true, true, true],
    };
    let reverse_complement = Tour {
        contigs: vec![2, 1, 0],
        signs: vec![false, false, false],
    };
    let exponent = (3_000.0_f64.ln() / 0.481_211_825_059_668_4).round();
    let representative =
        ((0.481_211_825_059_668_4 * exponent).exp().round() as u64).clamp(2_048, 1_u64 << 32);

    let legacy = problem.evaluate(&direct);
    let explicit_legacy =
        problem.evaluate_with_gap_model(&direct, OrientationGapModel::AllhicLegacy);
    let reverse_legacy =
        problem.evaluate_with_gap_model(&reverse_complement, OrientationGapModel::AllhicLegacy);
    assert_eq!(legacy.to_bits(), explicit_legacy.to_bits());
    assert_eq!(
        legacy.to_bits(),
        (-((representative + 100) as f64).ln()).to_bits()
    );
    assert_eq!(
        reverse_legacy.to_bits(),
        (-((representative + 300) as f64).ln()).to_bits()
    );
    assert_ne!(legacy.to_bits(), reverse_legacy.to_bits());

    let intervening = problem.evaluate_with_gap_model(&direct, OrientationGapModel::Intervening);
    let reverse_intervening =
        problem.evaluate_with_gap_model(&reverse_complement, OrientationGapModel::Intervening);
    let only_middle = -((representative + 1_000) as f64).ln();
    assert_eq!(intervening.to_bits(), only_middle.to_bits());
    assert_eq!(intervening.to_bits(), reverse_intervening.to_bits());
}

#[test]
fn orientation_banded_dp_matches_bruteforce_for_every_pair_weight() {
    let (problem, target) = exact_banded_orientation_problem();
    let initial = Tour {
        contigs: (0..5).collect(),
        signs: vec![false, true, false, true, false],
    };
    let mut initial_scores = Vec::new();
    for pair_weight in [
        OrientationPairWeight::Links,
        OrientationPairWeight::SqrtLinks,
        OrientationPairWeight::EqualPair,
    ] {
        let config = BandedOrientationConfig {
            rank_window: 3,
            gap_model: OrientationGapModel::Intervening,
            pair_weight,
            min_links: 1,
            require_complete: true,
        };
        let mut best_score = f64::NEG_INFINITY;
        let mut fewest_changes = usize::MAX;
        for bits in 0..(1usize << initial.signs.len()) {
            let signs = (0..initial.signs.len())
                .map(|position| bits & (1 << position) != 0)
                .collect::<Vec<_>>();
            let changes = signs
                .iter()
                .zip(&initial.signs)
                .filter(|(candidate, input)| candidate != input)
                .count();
            let mut candidate = Tour {
                contigs: initial.contigs.clone(),
                signs,
            };
            let score = problem
                .optimize_banded(&mut candidate, config)
                .unwrap()
                .initial_score;
            if score > best_score || (score == best_score && changes < fewest_changes) {
                best_score = score;
                fewest_changes = changes;
            }
        }

        let mut optimized = initial.clone();
        let result = problem.optimize_banded(&mut optimized, config).unwrap();
        assert_eq!(optimized.contigs, initial.contigs);
        assert_eq!(optimized.signs, target);
        assert_eq!(result.final_score.to_bits(), best_score.to_bits());
        assert_eq!(result.changed_signs, fewest_changes);
        assert!(result.final_score >= result.initial_score);
        assert_eq!(result.used_pairs, 9);
        assert_eq!(result.skipped_incomplete_pairs, 0);
        initial_scores.push(result.initial_score);
    }
    assert_ne!(initial_scores[0].to_bits(), initial_scores[1].to_bits());
    assert_ne!(initial_scores[1].to_bits(), initial_scores[2].to_bits());
}

#[test]
fn historical_block_refinement_rejects_bp_limits_without_mutation() {
    let (problem, target) = exact_banded_orientation_problem();
    let initial = Tour { contigs: (0..target.len()).collect(), signs: target };
    for fraction in [0.0, 0.05, 1.1, f64::NAN] {
        let mut tour = initial.clone();
        let error = problem.refine_signed_blocks_legacy_objective(
            &mut tour,
            BandedOrientationConfig::default(),
            0.05,
            SignedBlockRefineConfig { max_bp_fraction: fraction, ..SignedBlockRefineConfig::default() },
        ).unwrap_err();
        assert!(error.contains("requires max_bp_fraction = 1"));
        assert_eq!(tour.contigs, initial.contigs);
        assert_eq!(tour.signs, initial.signs);
    }
}

#[test]
fn signed_block_refinement_repairs_a_reverse_complemented_internal_block() {
    let target = [true, false, true, false];
    let records = [(0, 1), (1, 2), (2, 3)]
        .into_iter()
        .flat_map(|(u, v)| {
            let preferred = (usize::from(!target[u]) << 1) | usize::from(!target[v]);
            let mut distances = [30_000, 90_000, 150_000, 270_000];
            distances[preferred] = 3_000;
            complete_orientation_pair(u, v, distances, 10)
        })
        .collect::<Vec<_>>();
    let problem = OrientationProblem::from_oriented_distances(
        vec![200_000_000, 10_000_000, 20_000_000, 150_000_000],
        records,
    )
    .unwrap();
    let orientation_config = BandedOrientationConfig {
        rank_window: 1,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::SqrtLinks,
        min_links: 3,
        require_complete: true,
    };
    let block_config = SignedBlockRefineConfig {
        max_span: 3,
        max_bp_fraction: 0.1,
        max_passes: 4,
        min_relative_gain: 1e-4,
    };
    let mut tour = Tour {
        contigs: vec![0, 2, 1, 3],
        signs: vec![true, false, true, false],
    };
    let source_tour = tour.clone();

    problem
        .optimize_banded_with_source_prior(&mut tour, orientation_config, 0.05)
        .unwrap();
    let result = problem
        .refine_signed_blocks_conservative(
            &mut tour,
            &source_tour,
            orientation_config,
            0.05,
            block_config,
        )
        .unwrap();

    assert_eq!(tour.contigs, vec![0, 1, 2, 3]);
    assert_eq!(tour.signs, target);
    assert_eq!(result.accepted_moves, 1);
    assert!(result.evaluated_moves > 0);
    assert!(result.final_score > result.initial_score);
}

#[test]
fn signed_block_refinement_requires_both_boundaries() {
    let target = [true, false, true, false];
    let records = [(0, 1)]
        .into_iter()
        .flat_map(|(u, v)| {
            let preferred = (usize::from(!target[u]) << 1) | usize::from(!target[v]);
            let mut distances = [30_000, 90_000, 150_000, 270_000];
            distances[preferred] = 3_000;
            complete_orientation_pair(u, v, distances, 10)
        })
        .collect::<Vec<_>>();
    let problem = OrientationProblem::from_oriented_distances(
        vec![200_000_000, 10_000_000, 20_000_000, 150_000_000],
        records,
    )
    .unwrap();
    let orientation_config = BandedOrientationConfig {
        rank_window: 1,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::SqrtLinks,
        min_links: 3,
        require_complete: true,
    };
    let block_config = SignedBlockRefineConfig {
        max_span: 3,
        max_bp_fraction: 0.1,
        max_passes: 4,
        min_relative_gain: 1e-4,
    };
    let source = Tour {
        contigs: vec![0, 2, 1, 3],
        signs: vec![true, false, true, false],
    };
    let mut tour = source.clone();

    let result = problem
        .refine_signed_blocks_conservative(
            &mut tour,
            &source,
            orientation_config,
            0.05,
            block_config,
        )
        .unwrap();

    assert_eq!(tour.contigs, source.contigs);
    assert_eq!(tour.signs, source.signs);
    assert_eq!(result.accepted_moves, 0);
}

#[test]
fn signed_block_refinement_respects_bp_fraction() {
    let target = [true, false, true, false];
    let records = [(0, 1), (1, 2), (2, 3)]
        .into_iter()
        .flat_map(|(u, v)| {
            let preferred = (usize::from(!target[u]) << 1) | usize::from(!target[v]);
            let mut distances = [30_000, 90_000, 150_000, 270_000];
            distances[preferred] = 3_000;
            complete_orientation_pair(u, v, distances, 10)
        })
        .collect::<Vec<_>>();
    let problem = OrientationProblem::from_oriented_distances(
        vec![200_000_000, 10_000_000, 20_000_000, 150_000_000],
        records,
    )
    .unwrap();
    let orientation_config = BandedOrientationConfig {
        rank_window: 1,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::SqrtLinks,
        min_links: 3,
        require_complete: true,
    };
    let block_config = SignedBlockRefineConfig {
        max_span: 3,
        max_bp_fraction: 0.05,
        max_passes: 4,
        min_relative_gain: 1e-4,
    };
    let source = Tour {
        contigs: vec![0, 2, 1, 3],
        signs: vec![true, false, true, false],
    };
    let mut tour = source.clone();

    let result = problem
        .refine_signed_blocks_conservative(
            &mut tour,
            &source,
            orientation_config,
            0.05,
            block_config,
        )
        .unwrap();

    assert_eq!(tour.contigs, source.contigs);
    assert_eq!(tour.signs, source.signs);
    assert_eq!(result.accepted_moves, 0);
}

#[test]
fn signed_block_refinement_is_reverse_complement_covariant() {
    let target = [true, false, true, false];
    let records = [(0, 1), (1, 2), (2, 3)]
        .into_iter()
        .flat_map(|(u, v)| {
            let preferred = (usize::from(!target[u]) << 1) | usize::from(!target[v]);
            let mut distances = [30_000, 90_000, 150_000, 270_000];
            distances[preferred] = 3_000;
            complete_orientation_pair(u, v, distances, 10)
        })
        .collect::<Vec<_>>();
    let problem = OrientationProblem::from_oriented_distances(
        vec![200_000_000, 10_000_000, 20_000_000, 150_000_000],
        records,
    )
    .unwrap();
    let orientation_config = BandedOrientationConfig {
        rank_window: 1,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::SqrtLinks,
        min_links: 3,
        require_complete: true,
    };
    let block_config = SignedBlockRefineConfig {
        max_span: 3,
        max_bp_fraction: 0.1,
        max_passes: 4,
        min_relative_gain: 1e-4,
    };
    let mut direct = Tour {
        contigs: vec![0, 2, 1, 3],
        signs: vec![true, false, true, false],
    };
    let mut reverse_complement = Tour {
        contigs: direct.contigs.iter().rev().copied().collect(),
        signs: direct.signs.iter().rev().map(|sign| !*sign).collect(),
    };

    let direct_source = direct.clone();
    let reverse_source = reverse_complement.clone();
    problem
        .refine_signed_blocks_conservative(
            &mut direct,
            &direct_source,
            orientation_config,
            0.05,
            block_config,
        )
        .unwrap();
    problem
        .refine_signed_blocks_conservative(
            &mut reverse_complement,
            &reverse_source,
            orientation_config,
            0.05,
            block_config,
        )
        .unwrap();

    assert_eq!(
        reverse_complement.contigs,
        direct.contigs.iter().rev().copied().collect::<Vec<_>>()
    );
    assert_eq!(
        reverse_complement.signs,
        direct
            .signs
            .iter()
            .rev()
            .map(|sign| !*sign)
            .collect::<Vec<_>>()
    );
}

#[test]
fn orientation_banded_skips_missing_and_unequal_four_cell_pairs() {
    let mut records = complete_orientation_pair(0, 1, [3_000, 4_000, 5_000, 6_000], 2);
    records.pop();
    let mut unequal = complete_orientation_pair(1, 2, [6_000, 5_000, 4_000, 3_000], 2);
    unequal[3].distances.pop();
    records.extend(unequal);
    let problem =
        OrientationProblem::from_oriented_distances(vec![100, 200, 300], records).unwrap();
    let original_signs = vec![false, true, false];
    let mut tour = Tour {
        contigs: vec![0, 1, 2],
        signs: original_signs.clone(),
    };

    let result = problem
        .optimize_banded(&mut tour, BandedOrientationConfig::default())
        .unwrap();

    assert_eq!(tour.signs, original_signs);
    assert_eq!(result.initial_score.to_bits(), 0.0_f64.to_bits());
    assert_eq!(result.final_score.to_bits(), 0.0_f64.to_bits());
    assert_eq!(result.changed_signs, 0);
    assert_eq!(result.used_pairs, 0);
    assert_eq!(result.skipped_incomplete_pairs, 2);
}

#[test]
fn orientation_banded_solution_is_strictly_reverse_complement_covariant() {
    let (problem, target) = exact_banded_orientation_problem();
    let mut direct = Tour {
        contigs: (0..5).collect(),
        signs: vec![false, true, false, true, false],
    };
    let mut reverse_complement = Tour {
        contigs: direct.contigs.iter().rev().copied().collect(),
        signs: direct.signs.iter().rev().map(|sign| !*sign).collect(),
    };
    let direct_order = direct.contigs.clone();
    let reverse_order = reverse_complement.contigs.clone();
    let config = BandedOrientationConfig {
        rank_window: 3,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::SqrtLinks,
        min_links: 1,
        require_complete: true,
    };

    let direct_result = problem.optimize_banded(&mut direct, config).unwrap();
    let reverse_result = problem
        .optimize_banded(&mut reverse_complement, config)
        .unwrap();

    assert_eq!(direct.contigs, direct_order);
    assert_eq!(reverse_complement.contigs, reverse_order);
    assert_eq!(direct.signs, target);
    assert_eq!(
        reverse_complement.signs,
        direct
            .signs
            .iter()
            .rev()
            .map(|sign| !*sign)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        direct_result.initial_score.to_bits(),
        reverse_result.initial_score.to_bits()
    );
    assert_eq!(
        direct_result.final_score.to_bits(),
        reverse_result.final_score.to_bits()
    );
    assert_eq!(direct_result.changed_signs, reverse_result.changed_signs);
    assert_eq!(direct_result.used_pairs, reverse_result.used_pairs);
    assert_eq!(
        direct_result.skipped_incomplete_pairs,
        reverse_result.skipped_incomplete_pairs
    );
}

#[test]
fn orientation_banded_zero_source_prior_is_bit_exact_and_invalid_priors_are_rejected() {
    let (problem, _) = exact_banded_orientation_problem();
    let initial = Tour {
        contigs: (0..5).collect(),
        signs: vec![false, true, false, true, false],
    };
    let config = BandedOrientationConfig {
        rank_window: 3,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::SqrtLinks,
        min_links: 1,
        require_complete: true,
    };
    let mut existing = initial.clone();
    let mut explicit_zero = initial.clone();

    let existing_result = problem.optimize_banded(&mut existing, config).unwrap();
    let explicit_result = problem
        .optimize_banded_with_source_prior(&mut explicit_zero, config, 0.0)
        .unwrap();

    assert_eq!(existing.contigs, explicit_zero.contigs);
    assert_eq!(existing.signs, explicit_zero.signs);
    assert_eq!(existing_result, explicit_result);
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        let mut rejected = initial.clone();
        let error = problem
            .optimize_banded_with_source_prior(&mut rejected, config, invalid)
            .unwrap_err();
        assert!(error.contains("finite and non-negative"));
        assert_eq!(rejected.contigs, initial.contigs);
        assert_eq!(rejected.signs, initial.signs);
    }
}

#[test]
fn orientation_banded_source_prior_matches_small_bruteforce() {
    let (lengths, pair_fixtures, problem) = source_prior_orientation_fixture();
    let initial = Tour {
        contigs: vec![0, 1, 2],
        signs: vec![false, false, false],
    };
    let config = BandedOrientationConfig {
        rank_window: 2,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::EqualPair,
        min_links: 1,
        require_complete: true,
    };
    let prior_strength = 0.35;
    let mut node_scale = vec![0.0; initial.signs.len()];
    for (u, v, pair_records) in pair_fixtures {
        let pair_problem =
            OrientationProblem::from_oriented_distances(lengths.clone(), pair_records).unwrap();
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        for bits in 0..4usize {
            let mut candidate = initial.clone();
            candidate.signs[u] = bits & 1 != 0;
            candidate.signs[v] = bits & 2 != 0;
            let score = pair_problem
                .optimize_banded(&mut candidate, config)
                .unwrap()
                .initial_score;
            minimum = minimum.min(score);
            maximum = maximum.max(score);
        }
        let edge_range = maximum - minimum;
        node_scale[u] += edge_range;
        node_scale[v] += edge_range;
    }

    let mut best_objective = f64::NEG_INFINITY;
    let mut best_raw_score = f64::NEG_INFINITY;
    let mut best_changes = usize::MAX;
    let mut best_signs = Vec::new();
    for bits in 0..(1usize << initial.signs.len()) {
        let signs = (0..initial.signs.len())
            .map(|position| bits & (1 << position) != 0)
            .collect::<Vec<_>>();
        let changes = signs
            .iter()
            .zip(&initial.signs)
            .filter(|(candidate, input)| candidate != input)
            .count();
        let mut candidate = Tour {
            contigs: initial.contigs.clone(),
            signs: signs.clone(),
        };
        let raw_score = problem
            .optimize_banded(&mut candidate, config)
            .unwrap()
            .initial_score;
        let penalty = signs
            .iter()
            .zip(&initial.signs)
            .enumerate()
            .filter(|(_, (candidate, input))| candidate != input)
            .map(|(rank, _)| node_scale[rank])
            .sum::<f64>();
        let objective = raw_score - prior_strength * penalty;
        if objective > best_objective || (objective == best_objective && changes < best_changes) {
            best_objective = objective;
            best_raw_score = raw_score;
            best_changes = changes;
            best_signs = signs;
        }
    }

    let mut optimized = initial;
    let result = problem
        .optimize_banded_with_source_prior(&mut optimized, config, prior_strength)
        .unwrap();

    assert_eq!(optimized.signs, best_signs);
    assert_eq!(result.final_score.to_bits(), best_raw_score.to_bits());
    assert_eq!(result.changed_signs, best_changes);
    assert!(best_objective >= result.initial_score);
    assert!(result.final_score >= result.initial_score);
}

#[test]
fn orientation_banded_source_prior_controls_changes_and_is_strictly_rc_covariant() {
    let (_, _, problem) = source_prior_orientation_fixture();
    let initial = Tour {
        contigs: vec![0, 1, 2],
        signs: vec![false, false, false],
    };
    let config = BandedOrientationConfig {
        rank_window: 2,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::EqualPair,
        min_links: 1,
        require_complete: true,
    };
    let mut changed_counts = Vec::new();
    for strength in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let mut candidate = initial.clone();
        let result = problem
            .optimize_banded_with_source_prior(&mut candidate, config, strength)
            .unwrap();
        assert!(result.final_score >= result.initial_score);
        changed_counts.push(result.changed_signs);
        if strength == 1.0 {
            assert_eq!(candidate.signs, initial.signs);
            assert_eq!(result.changed_signs, 0);
        }
    }
    assert!(
        changed_counts
            .windows(2)
            .all(|counts| counts[1] <= counts[0]),
        "source-prior changes were not monotonic: {changed_counts:?}"
    );

    let mut direct = initial.clone();
    let mut reverse_complement = Tour {
        contigs: direct.contigs.iter().rev().copied().collect(),
        signs: direct.signs.iter().rev().map(|sign| !*sign).collect(),
    };
    let direct_result = problem
        .optimize_banded_with_source_prior(&mut direct, config, 0.5)
        .unwrap();
    let reverse_result = problem
        .optimize_banded_with_source_prior(&mut reverse_complement, config, 0.5)
        .unwrap();

    assert_eq!(
        reverse_complement.signs,
        direct
            .signs
            .iter()
            .rev()
            .map(|sign| !*sign)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        direct_result.initial_score.to_bits(),
        reverse_result.initial_score.to_bits()
    );
    assert_eq!(
        direct_result.final_score.to_bits(),
        reverse_result.final_score.to_bits()
    );
    assert_eq!(direct_result.changed_signs, reverse_result.changed_signs);
}

#[test]
fn orientation_banded_confidence_does_not_charge_the_source_prior_twice() {
    let problem = OrientationProblem::from_oriented_distances(
        vec![100_000; 2],
        complete_orientation_pair(0, 1, [3_000, 3_000, 270_000, 270_000], 64),
    ).unwrap();
    let mut tour = Tour { contigs: vec![0, 1], signs: vec![false, false] };
    let result = problem.optimize_banded_conservative(
        &mut tour,
        BandedOrientationConfig {
            rank_window: 1, gap_model: OrientationGapModel::Intervening,
            pair_weight: OrientationPairWeight::EqualPair, min_links: 1, require_complete: true,
        },
        0.05, 0.95, 1.0,
    ).unwrap();
    // The contact-only margin is maximal; the 0.05 optimization prior must
    // not cap its reported effect at 0.95 and make the >0.95 gate unreachable.
    assert_eq!(tour.signs, vec![true, false]);
    assert_eq!(result.changed_signs, 1);
    assert_eq!(result.rejected_low_confidence, 0);
}

#[test]
fn orientation_banded_gates_preserve_the_joint_objective_on_physical_contacts() {
    use rand::Rng;
    let mut rng = SmallRng::seed_from_u64(20260905);
    let config = BandedOrientationConfig {
        rank_window: 3, gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::SqrtLinks, min_links: 3, require_complete: true,
    };
    let mut confidence_rollbacks = 0;
    let mut bp_rollbacks = 0;
    // Small fixed-seed cases cover interacting sign changes with unequal bp
    // lengths. Every quartet comes from valid positions on its two contigs.
    for case in 0..256 {
        let lengths = (0..5).map(|_| rng.gen_range(10_000..100_000u64)).collect::<Vec<_>>();
        let input = Tour {
            contigs: (0..5).collect(),
            signs: (0..5).map(|_| rng.gen_bool(0.5)).collect(),
        };
        let mut records = Vec::new();
        let mut scales = [0.0; 5];
        for u in 0..5 {
            for v in u + 1..5 {
                let x = rng.gen_range(0..=lengths[u]);
                let y = rng.gen_range(0..=lengths[v]);
                let pair = complete_orientation_pair(u, v, [
                    lengths[u] - x + y, lengths[u] - x + lengths[v] - y,
                    x + y, x + lengths[v] - y,
                ], rng.gen_range(3..=64));
                let pair_problem = OrientationProblem::from_oriented_distances(lengths.clone(), pair.clone()).unwrap();
                let mut minimum = f64::INFINITY;
                let mut maximum = f64::NEG_INFINITY;
                for bits in 0..4 {
                    let mut candidate = input.clone();
                    candidate.signs[u] = bits & 1 != 0;
                    candidate.signs[v] = bits & 2 != 0;
                    let score = pair_problem.optimize_banded(&mut candidate, config).unwrap().initial_score;
                    minimum = minimum.min(score);
                    maximum = maximum.max(score);
                }
                scales[u] += maximum - minimum;
                scales[v] += maximum - minimum;
                records.extend(pair);
            }
        }
        let problem = OrientationProblem::from_oriented_distances(lengths, records).unwrap();
        for (confidence, cap) in [(0.1, 1.0), (0.25, 1.0), (0.4, 1.0), (0.0, 0.25), (0.2, 0.25)] {
            let mut candidate = input.clone();
            let result = problem.optimize_banded_conservative(&mut candidate, config, 0.05, confidence, cap).unwrap();
            let mut measured = candidate.clone();
            let final_score = problem.optimize_banded(&mut measured, config).unwrap().initial_score;
            let penalty = candidate.signs.iter().zip(&input.signs).zip(&scales)
                .filter(|((candidate, input), _)| candidate != input)
                .map(|(_, scale)| 0.05 * scale).sum::<f64>();
            assert!((result.final_score - final_score).abs() < 1e-10);
            assert!(final_score - penalty >= result.initial_score - 1e-10,
                "case {case}, confidence {confidence}, cap {cap}: {} -> {}", result.initial_score, final_score - penalty);
            assert_eq!(result.changed_signs, candidate.signs.iter().zip(&input.signs).filter(|(a, b)| a != b).count());
            if result.rejected_joint_changes > 0 {
                assert_eq!(candidate.signs, input.signs);
                assert_eq!(result.initial_score.to_bits(), result.final_score.to_bits());
                if cap == 1.0 { confidence_rollbacks += 1; }
                if confidence == 0.0 { bp_rollbacks += 1; }
            }
            let mut reverse = Tour {
                contigs: input.contigs.iter().rev().copied().collect(),
                signs: input.signs.iter().rev().map(|sign| !*sign).collect(),
            };
            let reverse_result = problem.optimize_banded_contact_evidence(&mut reverse, config, 0.05, confidence, cap).unwrap();
            assert_eq!(reverse.signs, candidate.signs.iter().rev().map(|sign| !*sign).collect::<Vec<_>>());
            assert_eq!(reverse_result, result);
        }
    }
    assert!(confidence_rollbacks > 0, "fixture must exercise joint regressions caused by confidence filtering");
    assert!(bp_rollbacks > 0, "fixture must exercise joint regressions caused by bp filtering");
}

#[test]
fn orientation_banded_max_marginal_gate_rejects_uncertain_changes() {
    let (_, _, problem) = source_prior_orientation_fixture();
    let initial = Tour {
        contigs: vec![0, 1, 2],
        signs: vec![false, false, false],
    };
    let config = BandedOrientationConfig {
        rank_window: 2,
        gap_model: OrientationGapModel::Intervening,
        pair_weight: OrientationPairWeight::EqualPair,
        min_links: 1,
        require_complete: true,
    };
    let mut ungated = initial.clone();
    let mut zero_gate = initial.clone();
    let mut fully_guarded = initial.clone();
    let mut bp_guarded = initial.clone();

    let ungated_result = problem
        .optimize_banded_with_source_prior(&mut ungated, config, 0.05)
        .unwrap();
    let zero_gate_result = problem
        .optimize_banded_with_source_prior_and_confidence(&mut zero_gate, config, 0.05, 0.0)
        .unwrap();
    let guarded_result = problem
        .optimize_banded_with_source_prior_and_confidence(&mut fully_guarded, config, 0.05, 1.0)
        .unwrap();
    let bp_guarded_result = problem
        .optimize_banded_conservative(&mut bp_guarded, config, 0.05, 0.0, 0.01)
        .unwrap();

    assert_eq!(ungated.signs, zero_gate.signs);
    assert_eq!(ungated_result, zero_gate_result);
    assert_eq!(fully_guarded.signs, initial.signs);
    assert_eq!(guarded_result.changed_signs, 0);
    assert_eq!(
        guarded_result.rejected_low_confidence,
        ungated_result.changed_signs
    );
    assert_eq!(bp_guarded.signs, initial.signs);
    assert_eq!(bp_guarded_result.changed_signs, 0);
    assert_eq!(
        bp_guarded_result.rejected_oversized_changes,
        ungated_result.changed_signs
    );

    for invalid in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        let mut rejected = initial.clone();
        let error = problem
            .optimize_banded_with_source_prior_and_confidence(&mut rejected, config, 0.05, invalid)
            .unwrap_err();
        assert!(error.contains("between 0 and 1"));
        assert_eq!(rejected.signs, initial.signs);
    }
    for invalid in [-0.1, 0.0, 1.1, f64::NAN, f64::INFINITY] {
        let mut rejected = initial.clone();
        let error = problem
            .optimize_banded_conservative(&mut rejected, config, 0.05, 0.0, invalid)
            .unwrap_err();
        assert!(error.contains("in (0, 1]"));
        assert_eq!(rejected.signs, initial.signs);
    }
}

#[test]
fn orientation_refine_with_legacy_gap_is_bit_exact_with_existing_refine() {
    let (problem, _) = exact_banded_orientation_problem();
    let initial = Tour {
        contigs: vec![2, 4, 0, 3, 1],
        signs: vec![false, true, true, false, true],
    };
    let mut existing = initial.clone();
    let mut explicit = initial;

    let existing_result = problem.refine(&mut existing);
    let explicit_result =
        problem.refine_with_gap_model(&mut explicit, OrientationGapModel::AllhicLegacy);

    assert_eq!(existing.contigs, explicit.contigs);
    assert_eq!(existing.signs, explicit.signs);
    assert_eq!(
        existing_result.initial_score.to_bits(),
        explicit_result.initial_score.to_bits()
    );
    assert_eq!(
        existing_result.final_score.to_bits(),
        explicit_result.final_score.to_bits()
    );
    assert_eq!(existing_result.phases, explicit_result.phases);
}

#[test]
fn orientation_refine_with_intervening_gap_is_strictly_rc_covariant_and_monotonic() {
    let (problem, target) = exact_banded_orientation_problem();
    let mut direct = Tour {
        contigs: (0..5).collect(),
        signs: vec![false, true, false, true, false],
    };
    let mut reverse_complement = Tour {
        contigs: direct.contigs.iter().rev().copied().collect(),
        signs: direct.signs.iter().rev().map(|sign| !*sign).collect(),
    };
    let direct_order = direct.contigs.clone();
    let reverse_order = reverse_complement.contigs.clone();

    let direct_result =
        problem.refine_with_gap_model(&mut direct, OrientationGapModel::Intervening);
    let reverse_result =
        problem.refine_with_gap_model(&mut reverse_complement, OrientationGapModel::Intervening);

    assert_eq!(direct.contigs, direct_order);
    assert_eq!(reverse_complement.contigs, reverse_order);
    assert_eq!(direct.signs, target);
    assert_eq!(
        reverse_complement.signs,
        direct
            .signs
            .iter()
            .rev()
            .map(|sign| !*sign)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        direct_result.initial_score.to_bits(),
        reverse_result.initial_score.to_bits()
    );
    assert_eq!(
        direct_result.final_score.to_bits(),
        reverse_result.final_score.to_bits()
    );
    assert_eq!(direct_result.phases, reverse_result.phases);
    assert!(direct_result.final_score >= direct_result.initial_score);
    assert!(reverse_result.final_score >= reverse_result.initial_score);
    assert_eq!(
        direct_result.final_score.to_bits(),
        problem
            .evaluate_with_gap_model(&direct, OrientationGapModel::Intervening)
            .to_bits()
    );
    assert_eq!(
        reverse_result.final_score.to_bits(),
        problem
            .evaluate_with_gap_model(&reverse_complement, OrientationGapModel::Intervening)
            .to_bits()
    );
}

#[test]
fn orientation_prefers_the_shorter_oriented_contact_distances() {
    let records = vec![
        OrientedDistanceRecord {
            u: 0,
            v: 1,
            u_forward: true,
            v_forward: true,
            distances: vec![20_000; 10],
        },
        OrientedDistanceRecord {
            u: 0,
            v: 1,
            u_forward: true,
            v_forward: false,
            distances: vec![3_000; 10],
        },
        OrientedDistanceRecord {
            u: 0,
            v: 1,
            u_forward: false,
            v_forward: true,
            distances: vec![3_000; 10],
        },
        OrientedDistanceRecord {
            u: 0,
            v: 1,
            u_forward: false,
            v_forward: false,
            distances: vec![20_000; 10],
        },
    ];
    let problem =
        OrientationProblem::from_oriented_distances(vec![100_000, 100_000], records).unwrap();
    let mut tour = Tour {
        contigs: vec![0, 1],
        signs: vec![true, true],
    };
    let before = problem.evaluate(&tour);
    let result = problem.optimize(&mut tour);

    assert_ne!(tour.signs[0], tour.signs[1]);
    assert!(result.final_score >= before);
}

#[test]
fn orientation_aggregates_repeated_project_clm_records() {
    let records = (0..2).map(|_| OrientedDistanceRecord {
        u: 0,
        v: 1,
        u_forward: true,
        v_forward: true,
        distances: vec![3_000],
    });
    let problem =
        OrientationProblem::from_oriented_distances(vec![100_000, 100_000], records).unwrap();
    let tour = Tour {
        contigs: vec![0, 1],
        signs: vec![true, true],
    };
    let exponent = (3_000.0_f64.ln() / 0.481_211_825_059_668_4).round();
    let representative = (0.481_211_825_059_668_4 * exponent).exp().round();

    assert!((problem.evaluate(&tour) + 2.0 * representative.ln()).abs() < 1e-10);
}

#[test]
fn clmb_builds_ordering_and_orientation_from_the_same_records() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("contacts.clmb");
    let names = vec!["a".to_string(), "b".to_string()];
    let mut writer =
        ClmbWriter::create_synchronous(&path, &names, 1024, Some(4), Some(40)).unwrap();
    for (a_orientation, b_orientation, distance) in
        [(0, 0, 20_000), (0, 1, 3_000), (1, 0, 3_000), (1, 1, 20_000)]
    {
        writer
            .write_record(
                encode_endpoint(0, a_orientation).unwrap(),
                encode_endpoint(1, b_orientation).unwrap(),
                &vec![distance; 10],
            )
            .unwrap();
    }
    writer.finish().unwrap();

    let problem = AllhicProblem::from_clmb(&path, &names, &[100_000, 100_000]).unwrap();
    assert!((problem.ordering.evaluate(&[0, 1]) + 10.0 / 100_000.0).abs() < 1e-12);

    let mut tour = Tour {
        contigs: vec![0, 1],
        signs: vec![true, true],
    };
    problem.orientation.optimize(&mut tour);
    assert_ne!(tour.signs[0], tour.signs[1]);
}
