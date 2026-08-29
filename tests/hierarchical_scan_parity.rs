#[path = "../src/order/hierarchical_scan.rs"]
mod hierarchical_scan;

use hierarchical_scan::{
    HierarchicalEndCandidate, OrientedPathHalves, sparse_hierarchical_end_candidates,
};
use rand::prelude::*;
use rand::rngs::SmallRng;

fn oriented_halves(path: &[(usize, bool)], lengths: &[usize]) -> OrientedPathHalves {
    let mut nodes = Vec::with_capacity(path.len() * 2);
    let mut half_lengths = Vec::with_capacity(path.len() * 2);
    for &(contig, forward) in path {
        let half_length = lengths[contig] as f64 * 0.5;
        if forward {
            nodes.extend([2 * contig, 2 * contig + 1]);
        } else {
            nodes.extend([2 * contig + 1, 2 * contig]);
        }
        half_lengths.extend([half_length, half_length]);
    }
    let target = half_lengths.iter().sum::<f64>() * 0.5;
    let mut cumulative = 0.0;
    let mut best_split = 1;
    let mut best_difference = f64::INFINITY;
    for split in 1..nodes.len() {
        cumulative += half_lengths[split - 1];
        let difference = (cumulative - target).abs();
        if difference < best_difference {
            best_difference = difference;
            best_split = split;
        }
    }
    [
        (
            nodes[..best_split].to_vec(),
            half_lengths[..best_split].iter().sum(),
        ),
        (
            nodes[best_split..].to_vec(),
            half_lengths[best_split..].iter().sum(),
        ),
    ]
}

fn all_pair_candidates(
    halves: &[OrientedPathHalves],
    contacts: &[Vec<(usize, f64)>],
) -> Vec<HierarchicalEndCandidate> {
    let contact = |left: usize, right: usize| {
        contacts[left]
            .binary_search_by_key(&right, |&(neighbor, _)| neighbor)
            .ok()
            .map_or(0.0, |index| contacts[left][index].1)
    };
    let mut candidates = Vec::new();
    for left_path in 0..halves.len() {
        for right_path in left_path + 1..halves.len() {
            for left_side in 0..2 {
                for right_side in 0..2 {
                    let (left_nodes, left_length) = &halves[left_path][left_side];
                    let (right_nodes, right_length) = &halves[right_path][right_side];
                    let mut count = 0.0;
                    for &left in left_nodes {
                        for &right in right_nodes {
                            count += contact(left, right);
                        }
                    }
                    if count <= 0.0 {
                        continue;
                    }
                    let score = count / (left_length * right_length).max(1.0);
                    candidates.push(HierarchicalEndCandidate {
                        left_path,
                        left_side,
                        right_path,
                        right_side,
                        normalized_score: score,
                        raw_support: count,
                    });
                }
            }
        }
    }
    candidates
}

#[test]
fn sparse_scan_matches_all_pair_reference() {
    let mut rng = SmallRng::seed_from_u64(0x5eed_cafe);
    for _case in 0..256 {
        let n = rng.gen_range(2..=24);
        let lengths = (0..n)
            .map(|_| rng.gen_range(1..=10_000))
            .collect::<Vec<_>>();

        let mut contigs = (0..n).collect::<Vec<_>>();
        contigs.shuffle(&mut rng);
        let mut paths = Vec::new();
        let mut start = 0;
        while start < n {
            let end = (start + rng.gen_range(1..=4)).min(n);
            paths.push(
                contigs[start..end]
                    .iter()
                    .map(|&contig| (contig, rng.gen_bool(0.5)))
                    .collect::<Vec<_>>(),
            );
            start = end;
        }
        let halves = paths
            .iter()
            .map(|path| oriented_halves(path, &lengths))
            .collect::<Vec<_>>();

        let mut rows = vec![Vec::<(usize, f64)>::new(); 2 * n];
        for left in 0..2 * n {
            for right in left + 1..2 * n {
                if rng.gen_bool(0.12) {
                    let count = rng.gen_range(1..=1_000) as f64 / rng.gen_range(1..=17) as f64;
                    rows[left].push((right, count));
                    rows[right].push((left, count));
                }
            }
        }
        for row in &mut rows {
            row.sort_unstable_by_key(|(neighbor, _)| *neighbor);
        }

        let sparse = sparse_hierarchical_end_candidates(&halves, &rows);
        let reference = all_pair_candidates(&halves, &rows);
        assert_eq!(sparse.len(), reference.len());
        for (observed, expected) in sparse.iter().zip(&reference) {
            assert_eq!(observed.left_path, expected.left_path);
            assert_eq!(observed.left_side, expected.left_side);
            assert_eq!(observed.right_path, expected.right_path);
            assert_eq!(observed.right_side, expected.right_side);
            assert_eq!(
                observed.normalized_score.to_bits(),
                expected.normalized_score.to_bits()
            );
            assert_eq!(
                observed.raw_support.to_bits(),
                expected.raw_support.to_bits()
            );
        }
    }
}

#[test]
fn sparse_scan_handles_fifty_thousand_paths() {
    let n = 50_000;
    let halves = (0..n)
        .map(|contig| [(vec![2 * contig], 1.0), (vec![2 * contig + 1], 1.0)])
        .collect::<Vec<_>>();
    let mut rows = vec![Vec::new(); 2 * n];
    for contig in 0..n - 1 {
        let left = 2 * contig + 1;
        let right = 2 * (contig + 1);
        rows[left].push((right, 1.0));
        rows[right].push((left, 1.0));
    }

    let candidates = sparse_hierarchical_end_candidates(&halves, &rows);

    assert_eq!(candidates.len(), n - 1);
    assert_eq!(
        candidates.first(),
        Some(&HierarchicalEndCandidate {
            left_path: 0,
            left_side: 1,
            right_path: 1,
            right_side: 0,
            normalized_score: 1.0,
            raw_support: 1.0,
        })
    );
    assert_eq!(
        candidates.last(),
        Some(&HierarchicalEndCandidate {
            left_path: n - 2,
            left_side: 1,
            right_path: n - 1,
            right_side: 0,
            normalized_score: 1.0,
            raw_support: 1.0,
        })
    );
}
