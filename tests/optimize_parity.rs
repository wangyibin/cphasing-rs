use cphasing::clm::{ClmbWriter, encode_endpoint};
use cphasing::core::ContigPair;
use cphasing::optimize::{
    AllhicProblem, BackboneConfig, BackbonePlan, OptimizeConfig, OptimizeProblem, OrderObjective,
    OrientationProblem, OrientedDistanceRecord, UnitGene, mutate_allhic, optimize_order,
};
use cphasing::order::{Tour, run_beam_end_initializer, run_hierarchical_end_initializer};
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

    let result = run_hierarchical_end_initializer(&initial, &lengths, &contacts, &names);
    let adjacencies = result
        .contigs
        .windows(2)
        .map(|pair| (pair[0].min(pair[1]), pair[0].max(pair[1])))
        .collect::<HashSet<_>>();

    assert_eq!(adjacencies, HashSet::from([(0, 1), (1, 2)]));
    assert_eq!(result.contigs.len(), result.signs.len());
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
    let problem = fragmented_path_problem(128);
    let initial = initial_tour(problem.contig_count());
    let config = OptimizeConfig {
        population_size: 12,
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
