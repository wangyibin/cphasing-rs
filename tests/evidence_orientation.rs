use cphasing::optimize::{EvidenceOrientationConfig, OrientationProblem, OrientedDistanceRecord};
use cphasing::order::Tour;

fn config() -> EvidenceOrientationConfig {
    EvidenceOrientationConfig {
        rank_window: 3,
        min_links: 3,
        min_effect: 0.2,
        input_prior: 0.0,
        block_span: 0,
        passes: 4,
    }
}

fn records(
    lengths: &[u64],
    target: &[bool],
    links: usize,
    centred: bool,
) -> Vec<OrientedDistanceRecord> {
    let mut records = Vec::new();
    for u in 0..lengths.len() {
        for v in u + 1..lengths.len() {
            let (lu, lv) = (lengths[u], lengths[v]);
            let x = if centred {
                lu / 2
            } else if target[u] {
                lu * 9 / 10
            } else {
                lu / 10
            };
            let y = if centred {
                lv / 2
            } else if target[v] {
                lv / 10
            } else {
                lv * 9 / 10
            };
            for ((u_forward, v_forward), distance) in
                [(true, true), (true, false), (false, true), (false, false)]
                    .into_iter()
                    .zip([lu - x + y, lu - x + lv - y, x + y, x + lv - y])
            {
                records.push(OrientedDistanceRecord {
                    u,
                    v,
                    u_forward,
                    v_forward,
                    distances: vec![distance; links],
                });
            }
        }
    }
    records
}

fn fixture(links: usize, centred: bool) -> (OrientationProblem, Tour<usize>) {
    let lengths = vec![100_000, 200_000, 30_000_000, 100_000, 100_000];
    let signs = vec![true, false, true, false, true];
    let data = records(&lengths, &signs, links, centred);
    (
        OrientationProblem::from_oriented_distances(lengths, data).unwrap(),
        Tour {
            contigs: (0..5).collect(),
            signs,
        },
    )
}

#[test]
fn evidence_repairs_arbitrary_signs_including_a_large_contig() {
    let (problem, truth) = fixture(64, false);
    for mask in 0..32 {
        let mut input = Tour {
            contigs: truth.contigs.clone(),
            signs: (0..5).map(|i| mask & (1 << i) != 0).collect(),
        };
        let result = problem.optimize_evidence(&mut input, config()).unwrap();
        assert_eq!(input.signs, truth.signs, "mask {mask}");
        assert_eq!(result.uncertain_contigs, 0);
        assert_eq!(result.valid_pairs, 10);
    }
}

#[test]
fn evidence_is_covariant_under_whole_scaffold_reverse_complement() {
    let (problem, truth) = fixture(64, false);
    let mut direct = Tour {
        contigs: truth.contigs.clone(),
        signs: vec![false; 5],
    };
    let mut reverse = Tour {
        contigs: direct.contigs.iter().rev().copied().collect(),
        signs: vec![true; 5],
    };
    let first = problem.optimize_evidence(&mut direct, config()).unwrap();
    let second = problem.optimize_evidence(&mut reverse, config()).unwrap();
    assert_eq!(
        direct.contigs,
        reverse.contigs.iter().rev().copied().collect::<Vec<_>>()
    );
    assert_eq!(
        direct.signs,
        reverse.signs.iter().rev().map(|s| !s).collect::<Vec<_>>()
    );
    assert_eq!(first.changed_signs, second.changed_signs);
}

#[test]
fn centred_or_missing_evidence_preserves_input_and_reports_uncertainty() {
    for links in [0, 2, 64] {
        let (problem, mut input) = fixture(links, true);
        let initial = input.clone();
        let result = problem.optimize_evidence(&mut input, config()).unwrap();
        assert_eq!(input.signs, initial.signs);
        assert_eq!(result.uncertain_contigs, 5);
        assert_eq!(result.changed_signs, 0);
    }
}

#[test]
fn incomplete_and_inconsistent_quartets_are_not_used() {
    let lengths = vec![100_000, 100_000];
    let base = records(&lengths, &[true, true], 64, false);
    for partial in [false, true] {
        let mut data = base.clone();
        if partial {
            data.pop();
        } else {
            data[0].distances.fill(1);
        }
        let problem = OrientationProblem::from_oriented_distances(lengths.clone(), data).unwrap();
        let mut input = Tour {
            contigs: vec![0, 1],
            signs: vec![false, false],
        };
        let result = problem.optimize_evidence(&mut input, config()).unwrap();
        assert_eq!(result.invalid_pairs, 1);
        assert_eq!(result.valid_pairs, 0);
        assert_eq!(input.signs, vec![false, false]);
    }
}

#[test]
fn duplicating_links_above_the_cap_cannot_increase_endpoint_effect() {
    let (first, mut a) = fixture(64, false);
    let (second, mut b) = fixture(640, false);
    let x = first.optimize_evidence(&mut a, config()).unwrap();
    let y = second.optimize_evidence(&mut b, config()).unwrap();
    for (x, y) in x.decisions.iter().zip(y.decisions.iter()) {
        assert_eq!(x.local.to_bits(), y.local.to_bits());
        assert_eq!(x.supported, y.supported);
    }
}

#[test]
fn explicit_trusted_prior_is_optional() {
    let (problem, _) = fixture(64, false);
    let mut input = Tour {
        contigs: (0..5).collect(),
        signs: vec![false; 5],
    };
    let initial = input.clone();
    let guarded = EvidenceOrientationConfig {
        input_prior: 1.0,
        ..config()
    };
    let result = problem.optimize_evidence(&mut input, guarded).unwrap();
    assert_eq!(input.signs, initial.signs);
    assert_eq!(result.changed_signs, 0);
}

#[test]
fn evidence_repairs_internal_and_terminal_signed_blocks_atomically() {
    let lengths = vec![100_000; 8];
    let mut data = records(&lengths, &[true; 8], 64, true);
    // Adjacent endpoints have a strong geometric signal; non-adjacent
    // contacts are uninformative. All quartets use valid mapped coordinates.
    for record in &mut data {
        if record.v == record.u + 1 {
            let distance = match (record.u_forward, record.v_forward) {
                (true, true) => 20_000,
                (true, false) | (false, true) => 100_000,
                (false, false) => 180_000,
            };
            record.distances.fill(distance);
        }
    }
    let problem = OrientationProblem::from_oriented_distances(lengths, data).unwrap();
    for (start, end) in [(2, 5), (0, 3), (4, 7)] {
        let mut tour = Tour {
            contigs: (0..8).collect(),
            signs: vec![true; 8],
        };
        tour.contigs[start..=end].reverse();
        tour.signs[start..=end].fill(false);
        let result = problem
            .optimize_evidence(
                &mut tour,
                EvidenceOrientationConfig {
                    block_span: 7,
                    ..config()
                },
            )
            .unwrap();
        assert_eq!(tour.contigs, (0..8).collect::<Vec<_>>(), "{start}..{end}");
        assert_eq!(tour.signs, vec![true; 8]);
        assert_eq!(result.accepted_blocks, 1);
    }
}

#[test]
fn empty_singleton_and_invalid_config_are_handled_without_mutation() {
    for n in 0..=1 {
        let problem = OrientationProblem::from_oriented_distances(vec![100; n], []).unwrap();
        let mut tour = Tour {
            contigs: (0..n).collect(),
            signs: vec![true; n],
        };
        assert_eq!(
            problem
                .optimize_evidence(&mut tour, config())
                .unwrap()
                .changed_signs,
            0
        );
    }
    let (problem, mut tour) = fixture(64, false);
    let initial = tour.clone();
    for bad in [f64::NAN, -0.1, 1.1] {
        assert!(
            problem
                .optimize_evidence(
                    &mut tour,
                    EvidenceOrientationConfig {
                        min_effect: bad,
                        ..config()
                    }
                )
                .is_err()
        );
        assert_eq!(tour.signs, initial.signs);
    }
}
