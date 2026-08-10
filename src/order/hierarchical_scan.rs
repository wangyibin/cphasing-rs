use hashbrown::HashMap;

pub(crate) type HierarchicalEndCandidate = (usize, usize, usize, usize, f64);
pub(crate) type OrientedPathHalves = [(Vec<usize>, f64); 2];

pub(crate) fn sparse_hierarchical_end_candidates(
    halves: &[OrientedPathHalves],
    contacts: &[Vec<(usize, f64)>],
) -> Vec<HierarchicalEndCandidate> {
    let mut location = vec![(usize::MAX, usize::MAX, usize::MAX); contacts.len()];
    for (path, path_halves) in halves.iter().enumerate() {
        for (side, (nodes, _)) in path_halves.iter().enumerate() {
            for (rank, &node) in nodes.iter().enumerate() {
                location[node] = (path, side, rank);
            }
        }
    }

    // Canonical key order is path-left, path-right, side-left, side-right.
    // Sorting these keys below reproduces the previous nested path-pair scan
    // order before downstream confidence ranking and tie handling.
    let mut contributions =
        HashMap::<(usize, usize, usize, usize), Vec<(usize, usize, f64)>>::new();
    for (left_node, row) in contacts.iter().enumerate() {
        for &(right_node, count) in row {
            if left_node >= right_node {
                continue;
            }
            let (left_path, left_side, left_rank) = location[left_node];
            let (right_path, right_side, right_rank) = location[right_node];
            let (key, ranks) = if left_path < right_path {
                (
                    (left_path, right_path, left_side, right_side),
                    (left_rank, right_rank),
                )
            } else if right_path < left_path {
                (
                    (right_path, left_path, right_side, left_side),
                    (right_rank, left_rank),
                )
            } else {
                continue;
            };
            contributions
                .entry(key)
                .or_default()
                .push((ranks.0, ranks.1, count));
        }
    }

    let mut contributions = contributions.into_iter().collect::<Vec<_>>();
    contributions.sort_unstable_by_key(|(key, _)| *key);
    contributions
        .into_iter()
        .filter_map(
            |((left_path, right_path, left_side, right_side), mut values)| {
                // Match the old nested left-node/right-node summation order so
                // fractional f64 contacts retain bitwise-identical scores.
                values.sort_unstable_by_key(|&(left_rank, right_rank, _)| (left_rank, right_rank));
                let count = values.into_iter().map(|(_, _, count)| count).sum::<f64>();
                if count <= 0.0 {
                    return None;
                }
                let left_length = halves[left_path][left_side].1;
                let right_length = halves[right_path][right_side].1;
                let score = count / (left_length * right_length).max(1.0);
                Some((left_path, left_side, right_path, right_side, score))
            },
        )
        .collect()
}
