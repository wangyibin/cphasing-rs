use cphasing::clm::{encode_endpoint, ClmbWriter};
use cphasing::optimize::{
    mutate_allhic, optimize_order, AllhicProblem, BackboneConfig, BackbonePlan, OptimizeConfig,
    OptimizeProblem, OrderObjective, OrientationProblem, OrientedDistanceRecord, UnitGene,
};
use cphasing::order::Tour;
use hashbrown::HashMap;
use indexmap::IndexMap;
use rand::rngs::SmallRng;
use rand::SeedableRng;

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
    assert!(plan
        .decode(&forward_genes[..plan.unit_count() - 1])
        .is_err());
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
