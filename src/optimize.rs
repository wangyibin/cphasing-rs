use crate::clm::ClmbReader;
use crate::order::Tour;
use hashbrown::HashMap;
use indexmap::IndexMap;
use rand::rngs::SmallRng;
use rand::Rng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::path::Path;

pub const DEFAULT_POPULATION_SIZE: usize = 100;
pub const DEFAULT_STALE_GENERATIONS: usize = 5_000;
pub const DEFAULT_MUTATION_PROBABILITY: f64 = 0.2;
pub const DEFAULT_MAX_GENERATIONS: usize = 1_000_000;
const GOLDEN_LOWER_BOUND: i32 = 16;
const GOLDEN_UPPER_BOUND: i32 = 50;
const GOLDEN_BINS: usize = (GOLDEN_UPPER_BOUND - GOLDEN_LOWER_BOUND + 1) as usize;
const LOG_GOLDEN_RATIO: f64 = 0.481_211_825_059_668_4;
const MAX_ORIENTATION_DISTANCE: u64 = 500_000_000;

pub type GoldenArray = [u32; GOLDEN_BINS];

#[derive(Debug, Clone, Copy)]
struct ContactEdge {
    u: u32,
    v: u32,
    links: f64,
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub enum OrderObjective {
    /// Project ALLHiC's default (`--logDist=false`).
    #[default]
    ReciprocalDistance,
    /// Project ALLHiC's optional `--logDist` objective.
    LogDistance,
}

#[derive(Debug, Clone)]
pub struct OptimizeProblem {
    lengths: Vec<f64>,
    edges: Vec<ContactEdge>,
    objective: OrderObjective,
}

impl OptimizeProblem {
    pub fn from_sparse_contacts(
        lengths: &IndexMap<usize, usize>,
        contacts: &[HashMap<usize, u32>],
    ) -> Result<Self, String> {
        let n = lengths.len();
        if n > u32::MAX as usize {
            return Err("ordering supports at most u32::MAX contigs".into());
        }
        let mut dense_lengths = vec![0.0; n];
        for (&id, &length) in lengths {
            if id >= n {
                return Err(format!(
                    "contig id {id} is outside the required dense range 0..{n}"
                ));
            }
            dense_lengths[id] = length as f64;
        }

        if contacts.len() != n {
            return Err(format!(
                "contact row count ({}) does not match contig count ({n})",
                contacts.len()
            ));
        }

        let mut edges = Vec::new();
        for (u, row) in contacts.iter().enumerate() {
            for (&v, &links) in row {
                if v >= n {
                    return Err(format!("contact endpoint {v} is outside 0..{n}"));
                }
                if u < v && links > 0 {
                    edges.push(ContactEdge {
                        u: u as u32,
                        v: v as u32,
                        links: f64::from(links),
                    });
                }
            }
        }
        edges.sort_unstable_by_key(|edge| (edge.u, edge.v));

        Ok(Self {
            lengths: dense_lengths,
            edges,
            objective: OrderObjective::default(),
        })
    }

    pub fn with_objective(mut self, objective: OrderObjective) -> Self {
        self.objective = objective;
        self
    }

    pub fn contig_count(&self) -> usize {
        self.lengths.len()
    }

    /// Evaluate the configured project-ALLHiC ordering objective. Lower is
    /// better for both reciprocal-distance and log-distance modes.
    pub fn evaluate(&self, order: &[usize]) -> f64 {
        let mut midpoints = vec![0.0; self.lengths.len()];
        self.evaluate_with_buffer(order, &mut midpoints)
    }

    fn evaluate_with_buffer(&self, order: &[usize], midpoints: &mut [f64]) -> f64 {
        if order.len() != self.lengths.len() || midpoints.len() < self.lengths.len() {
            return f64::INFINITY;
        }
        let mut cumulative = 0.0;

        for &id in order {
            if id >= self.lengths.len() {
                return f64::INFINITY;
            }
            let length = self.lengths[id];
            midpoints[id] = cumulative + length * 0.5;
            cumulative += length;
        }

        match self.objective {
            OrderObjective::ReciprocalDistance => {
                let mut score = 0.0;
                for edge in &self.edges {
                    let distance = (midpoints[edge.u as usize] - midpoints[edge.v as usize]).abs();
                    score -= edge.links / distance;
                }
                score
            }
            OrderObjective::LogDistance => {
                let mut score = 0.0;
                for edge in &self.edges {
                    let distance = (midpoints[edge.u as usize] - midpoints[edge.v as usize]).abs();
                    let distance = if distance <= 1.0 { 1.000001 } else { distance };
                    score += edge.links * distance.ln();
                }
                score
            }
        }
    }
}

/// Conservative controls for the path-block search used by the multilevel GA.
/// Normalized weights are used only to discover temporary blocks; every GA
/// individual is still evaluated with the unmodified ALLHiC objective.
#[derive(Debug, Clone)]
pub struct BackboneConfig {
    pub enabled: bool,
    pub min_links: f64,
    pub min_margin: f64,
    pub min_top_two_fraction: f64,
    pub min_reduction_fraction: f64,
    pub max_block_size: usize,
}

impl Default for BackboneConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            min_links: 5.0,
            min_margin: 1.5,
            // Dense Hi-C graphs contain many long-range edges, so requiring
            // the two strongest neighbors to carry half of all normalized
            // signal rejects nearly every real adjacency. A 0.2 gate still
            // removes diffuse hubs while retaining reciprocal top-two paths.
            min_top_two_fraction: 0.2,
            min_reduction_fraction: 0.1,
            max_block_size: 32,
        }
    }
}

impl BackboneConfig {
    fn validate(&self) -> Result<(), String> {
        if !self.min_links.is_finite() || self.min_links < 0.0 {
            return Err("backbone minimum links must be finite and non-negative".into());
        }
        if !self.min_margin.is_finite() || self.min_margin < 1.0 {
            return Err("backbone top-two/top-three margin must be at least 1".into());
        }
        if !(0.0..=1.0).contains(&self.min_top_two_fraction) {
            return Err("backbone top-two fraction must be between 0 and 1".into());
        }
        if !(0.0..=1.0).contains(&self.min_reduction_fraction) {
            return Err("backbone minimum reduction fraction must be between 0 and 1".into());
        }
        if self.max_block_size < 2 {
            return Err("backbone maximum block size must be at least 2".into());
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct UnitGene {
    pub unit: usize,
    pub reversed: bool,
}

/// A complete partition of the contigs into conservative path blocks and
/// singleton residual units. No contig is removed from a backbone plan.
#[derive(Debug, Clone)]
pub struct BackbonePlan {
    units: Vec<Vec<usize>>,
    unit_of_contig: Vec<usize>,
    block_count: usize,
    accepted_edges: usize,
    useful: bool,
}

#[derive(Debug, Clone, Copy)]
struct RankedBackboneEdge {
    edge_index: usize,
    neighbor: usize,
    score: f64,
}

#[derive(Debug, Clone, Copy)]
struct BackboneCandidate {
    u: usize,
    v: usize,
    confidence: f64,
}

impl BackbonePlan {
    pub fn discover(problem: &OptimizeProblem, config: &BackboneConfig) -> Self {
        let n = problem.contig_count();
        if !config.enabled || n < 3 || problem.edges.is_empty() {
            return Self::singletons(n);
        }

        let mut totals = vec![0.0; n];
        for edge in &problem.edges {
            totals[edge.u as usize] += edge.links;
            totals[edge.v as usize] += edge.links;
        }

        let mut ranked = vec![Vec::<RankedBackboneEdge>::new(); n];
        for (edge_index, edge) in problem.edges.iter().enumerate() {
            let u = edge.u as usize;
            let v = edge.v as usize;
            let denominator = (totals[u] * totals[v]).sqrt();
            if denominator <= 0.0 {
                continue;
            }
            let score = edge.links / denominator;
            ranked[u].push(RankedBackboneEdge {
                edge_index,
                neighbor: v,
                score,
            });
            ranked[v].push(RankedBackboneEdge {
                edge_index,
                neighbor: u,
                score,
            });
        }
        for row in &mut ranked {
            row.sort_unstable_by(|a, b| {
                b.score
                    .total_cmp(&a.score)
                    .then_with(|| a.neighbor.cmp(&b.neighbor))
            });
        }

        let mut edge_ranks = vec![(usize::MAX, usize::MAX); problem.edges.len()];
        let mut top_two_fraction = vec![0.0; n];
        let mut third_score = vec![0.0; n];
        for (node, row) in ranked.iter().enumerate() {
            let total_score: f64 = row.iter().map(|edge| edge.score).sum();
            let top_two_score: f64 = row.iter().take(2).map(|edge| edge.score).sum();
            if total_score > 0.0 {
                top_two_fraction[node] = top_two_score / total_score;
            }
            third_score[node] = row.get(2).map_or(0.0, |edge| edge.score);
            for (rank, edge) in row.iter().enumerate() {
                let endpoints = &mut edge_ranks[edge.edge_index];
                if problem.edges[edge.edge_index].u as usize == node {
                    endpoints.0 = rank;
                } else {
                    endpoints.1 = rank;
                }
            }
        }

        let mut candidates = Vec::new();
        let mut supported_edges = 0usize;
        let mut reciprocal_edges = 0usize;
        let mut concentrated_edges = 0usize;
        for (edge_index, edge) in problem.edges.iter().enumerate() {
            if edge.links < config.min_links {
                continue;
            }
            supported_edges += 1;
            let u = edge.u as usize;
            let v = edge.v as usize;
            let (rank_u, rank_v) = edge_ranks[edge_index];
            if rank_u >= 2 || rank_v >= 2 {
                continue;
            }
            reciprocal_edges += 1;
            if top_two_fraction[u] < config.min_top_two_fraction
                || top_two_fraction[v] < config.min_top_two_fraction
            {
                continue;
            }
            concentrated_edges += 1;
            let denominator = (totals[u] * totals[v]).sqrt();
            if denominator <= 0.0 {
                continue;
            }
            let score = edge.links / denominator;
            let margin_u = if third_score[u] > 0.0 {
                score / third_score[u]
            } else {
                f64::INFINITY
            };
            let margin_v = if third_score[v] > 0.0 {
                score / third_score[v]
            } else {
                f64::INFINITY
            };
            if margin_u < config.min_margin || margin_v < config.min_margin {
                continue;
            }
            let margin = margin_u.min(margin_v).min(10.0);
            let concentration = top_two_fraction[u].min(top_two_fraction[v]);
            candidates.push(BackboneCandidate {
                u,
                v,
                confidence: score * concentration * margin * edge.links.ln_1p(),
            });
        }
        log::info!(
            "Backbone discovery: {} supported, {} reciprocal-top2, {} concentrated, {} margin-qualified edges",
            supported_edges,
            reciprocal_edges,
            concentrated_edges,
            candidates.len()
        );
        candidates.sort_unstable_by(|a, b| {
            b.confidence
                .total_cmp(&a.confidence)
                .then_with(|| (a.u, a.v).cmp(&(b.u, b.v)))
        });

        let mut components = DisjointSet::new(n);
        let mut adjacency = vec![Vec::<usize>::new(); n];
        let mut accepted_edges = 0usize;
        for candidate in candidates {
            if adjacency[candidate.u].len() >= 2 || adjacency[candidate.v].len() >= 2 {
                continue;
            }
            let root_u = components.find(candidate.u);
            let root_v = components.find(candidate.v);
            if root_u == root_v
                || components.size(root_u) + components.size(root_v) > config.max_block_size
            {
                continue;
            }
            adjacency[candidate.u].push(candidate.v);
            adjacency[candidate.v].push(candidate.u);
            components.union_roots(root_u, root_v);
            accepted_edges += 1;
        }

        let mut visited = vec![false; n];
        let mut units = Vec::new();
        for start in 0..n {
            if visited[start] || adjacency[start].len() > 1 {
                continue;
            }
            let mut unit = Vec::new();
            let mut previous = usize::MAX;
            let mut current = start;
            loop {
                visited[current] = true;
                unit.push(current);
                let next = adjacency[current]
                    .iter()
                    .copied()
                    .find(|&neighbor| neighbor != previous && !visited[neighbor]);
                let Some(next) = next else {
                    break;
                };
                previous = current;
                current = next;
            }
            units.push(unit);
        }
        // The degree/cycle guards above guarantee paths, but retain every node
        // even if a future candidate policy violates that assumption.
        for contig in 0..n {
            if !visited[contig] {
                units.push(vec![contig]);
            }
        }
        units.sort_unstable_by_key(|unit| unit.iter().copied().min().unwrap_or(usize::MAX));

        let mut unit_of_contig = vec![usize::MAX; n];
        for (unit_index, unit) in units.iter().enumerate() {
            for &contig in unit {
                unit_of_contig[contig] = unit_index;
            }
        }
        let block_count = units.iter().filter(|unit| unit.len() > 1).count();
        let reduction = n.saturating_sub(units.len());
        let useful =
            block_count > 0 && (reduction as f64) >= config.min_reduction_fraction * n as f64;

        Self {
            units,
            unit_of_contig,
            block_count,
            accepted_edges,
            useful,
        }
    }

    fn singletons(n: usize) -> Self {
        Self {
            units: (0..n).map(|contig| vec![contig]).collect(),
            unit_of_contig: (0..n).collect(),
            block_count: 0,
            accepted_edges: 0,
            useful: false,
        }
    }

    pub fn units(&self) -> &[Vec<usize>] {
        &self.units
    }

    pub fn unit_count(&self) -> usize {
        self.units.len()
    }

    pub fn block_count(&self) -> usize {
        self.block_count
    }

    pub fn accepted_edges(&self) -> usize {
        self.accepted_edges
    }

    pub fn is_useful(&self) -> bool {
        self.useful
    }

    pub fn decode(&self, genes: &[UnitGene]) -> Result<Vec<usize>, String> {
        if genes.len() != self.units.len() {
            return Err(format!(
                "unit tour length ({}) does not match backbone unit count ({})",
                genes.len(),
                self.units.len()
            ));
        }
        let mut seen = vec![false; self.units.len()];
        let mut decoded = Vec::with_capacity(self.unit_of_contig.len());
        for gene in genes {
            if gene.unit >= self.units.len() {
                return Err(format!(
                    "backbone unit {} is outside 0..{}",
                    gene.unit,
                    self.units.len()
                ));
            }
            if std::mem::replace(&mut seen[gene.unit], true) {
                return Err(format!("unit tour contains duplicate unit {}", gene.unit));
            }
            self.append_unit(gene, &mut decoded);
        }
        Ok(decoded)
    }

    fn decode_into(&self, genes: &[UnitGene], decoded: &mut Vec<usize>) {
        decoded.clear();
        for gene in genes {
            self.append_unit(gene, decoded);
        }
    }

    fn append_unit(&self, gene: &UnitGene, decoded: &mut Vec<usize>) {
        let unit = &self.units[gene.unit];
        if gene.reversed {
            decoded.extend(unit.iter().rev().copied());
        } else {
            decoded.extend(unit.iter().copied());
        }
    }

    fn initial_genes(&self, seed_order: &[usize]) -> Vec<UnitGene> {
        let mut rank = vec![usize::MAX; self.unit_of_contig.len()];
        for (position, &contig) in seed_order.iter().enumerate() {
            rank[contig] = position;
        }
        let mut keyed = self
            .units
            .iter()
            .enumerate()
            .map(|(unit_index, unit)| {
                let key = unit
                    .iter()
                    .map(|&contig| rank[contig])
                    .min()
                    .unwrap_or(usize::MAX);
                let reversed = unit.len() > 1 && rank[unit[0]] > rank[*unit.last().unwrap()];
                (
                    key,
                    UnitGene {
                        unit: unit_index,
                        reversed,
                    },
                )
            })
            .collect::<Vec<_>>();
        keyed.sort_unstable_by_key(|(key, gene)| (*key, gene.unit));
        keyed.into_iter().map(|(_, gene)| gene).collect()
    }
}

struct DisjointSet {
    parent: Vec<usize>,
    sizes: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            sizes: vec![1; n],
        }
    }

    fn find(&mut self, node: usize) -> usize {
        if self.parent[node] != node {
            self.parent[node] = self.find(self.parent[node]);
        }
        self.parent[node]
    }

    fn size(&self, root: usize) -> usize {
        self.sizes[root]
    }

    fn union_roots(&mut self, first: usize, second: usize) {
        let (large, small) = if self.sizes[first] >= self.sizes[second] {
            (first, second)
        } else {
            (second, first)
        };
        self.parent[small] = large;
        self.sizes[large] += self.sizes[small];
    }
}

#[derive(Debug, Clone)]
pub struct OrientedDistanceRecord {
    pub u: usize,
    pub v: usize,
    pub u_forward: bool,
    pub v_forward: bool,
    pub distances: Vec<u64>,
}

#[derive(Debug, Clone)]
struct OrientationSummary {
    bins: GoldenArray,
    links: usize,
    sum_log_distance: f64,
    present: bool,
}

impl Default for OrientationSummary {
    fn default() -> Self {
        Self {
            bins: [0; GOLDEN_BINS],
            links: 0,
            sum_log_distance: 0.0,
            present: false,
        }
    }
}

#[derive(Debug, Clone)]
struct PairOrientationData {
    u: usize,
    v: usize,
    orientations: [OrientationSummary; 4],
    signed_links: f64,
}

#[derive(Debug, Clone)]
pub struct OrientationProblem {
    lengths: Vec<u64>,
    pairs: Vec<PairOrientationData>,
    incident_pairs: Vec<Vec<usize>>,
    golden_representatives: GoldenArrayRepresentatives,
}

type GoldenArrayRepresentatives = [u64; GOLDEN_BINS];
type OrientationSummaries = HashMap<(usize, usize), [OrientationSummary; 4]>;

#[derive(Debug, Clone)]
pub struct AllhicProblem {
    pub ordering: OptimizeProblem,
    pub orientation: OrientationProblem,
}

impl AllhicProblem {
    /// Build both ordering and orientation inputs in one streaming pass over
    /// the project's existing CLMB reader. No CLM/CLMB format logic is
    /// duplicated here.
    pub fn from_clmb(
        path: impl AsRef<Path>,
        contig_names: &[String],
        lengths: &[u64],
    ) -> Result<Self, String> {
        if contig_names.len() != lengths.len() {
            return Err(format!(
                "contig-name count ({}) does not match length count ({})",
                contig_names.len(),
                lengths.len()
            ));
        }
        if lengths.len() > u32::MAX as usize {
            return Err("ordering supports at most u32::MAX contigs".into());
        }
        let name_to_id: HashMap<&str, usize> = contig_names
            .iter()
            .enumerate()
            .map(|(id, name)| (name.as_str(), id))
            .collect();
        let mut reader = ClmbReader::open(path.as_ref()).map_err(|error| error.to_string())?;
        let dictionary = reader.header.contigs.clone();
        let dictionary_to_id: Vec<Option<usize>> = dictionary
            .iter()
            .map(|name| name_to_id.get(name.as_str()).copied())
            .collect();
        let mut summaries = OrientationSummaries::new();
        while let Some(block) = reader.next_block().map_err(|error| error.to_string())? {
            for record in block {
                let Some(u) = dictionary_to_id[record.contig1() as usize] else {
                    continue;
                };
                let Some(v) = dictionary_to_id[record.contig2() as usize] else {
                    continue;
                };
                OrientationProblem::add_record(
                    &mut summaries,
                    lengths.len(),
                    OrientedDistanceRecord {
                        u,
                        v,
                        u_forward: record.orientation1() == 0,
                        v_forward: record.orientation2() == 0,
                        distances: record.distances,
                    },
                )?;
            }
        }

        let orientation = OrientationProblem::from_summaries(lengths.to_vec(), summaries);
        let edges = orientation
            .pairs
            .iter()
            .filter_map(|pair| {
                let links = pair.signed_links.abs();
                (links > 0.0).then_some(ContactEdge {
                    u: pair.u as u32,
                    v: pair.v as u32,
                    links,
                })
            })
            .collect();
        let ordering = OptimizeProblem {
            lengths: lengths.iter().map(|&length| length as f64).collect(),
            edges,
            objective: OrderObjective::default(),
        };
        Ok(Self {
            ordering,
            orientation,
        })
    }

    pub fn with_objective(mut self, objective: OrderObjective) -> Self {
        self.ordering.objective = objective;
        self
    }
}

impl OrientationProblem {
    pub fn from_oriented_distances(
        lengths: Vec<u64>,
        records: impl IntoIterator<Item = OrientedDistanceRecord>,
    ) -> Result<Self, String> {
        let n = lengths.len();
        let mut summaries = OrientationSummaries::new();

        for record in records {
            Self::add_record(&mut summaries, n, record)?;
        }

        Ok(Self::from_summaries(lengths, summaries))
    }

    fn add_record(
        summaries: &mut OrientationSummaries,
        n: usize,
        record: OrientedDistanceRecord,
    ) -> Result<(), String> {
        if record.u >= n || record.v >= n {
            return Err(format!(
                "orientation contact endpoint ({}, {}) is outside 0..{n}",
                record.u, record.v
            ));
        }
        if record.u == record.v {
            return Ok(());
        }

        let (u, v, u_forward, v_forward) = if record.u < record.v {
            (record.u, record.v, record.u_forward, record.v_forward)
        } else {
            (record.v, record.u, !record.v_forward, !record.u_forward)
        };
        let orientation = orientation_index(u_forward, v_forward);
        let entry = &mut summaries
            .entry((u, v))
            .or_insert_with(|| std::array::from_fn(|_| OrientationSummary::default()))[orientation];
        entry.present = true;
        entry.links += record.distances.len();
        for distance in record.distances {
            if distance == 0 {
                return Err("orientation contact distance must be greater than zero".into());
            }
            entry.sum_log_distance += (distance as f64).ln();
            entry.bins[golden_bin(distance)] += 1;
        }
        Ok(())
    }

    fn from_summaries(lengths: Vec<u64>, summaries: OrientationSummaries) -> Self {
        let n = lengths.len();
        let mut pairs = Vec::with_capacity(summaries.len());
        for ((u, v), orientations) in summaries {
            let best = orientations
                .iter()
                .enumerate()
                .filter(|(_, summary)| summary.present)
                .min_by(|(_, a), (_, b)| a.sum_log_distance.total_cmp(&b.sum_log_distance));
            let signed_links = if let Some((orientation, summary)) = best {
                let same_orientation = orientation == 0 || orientation == 3;
                summary.links as f64 * if same_orientation { 1.0 } else { -1.0 }
            } else {
                0.0
            };
            pairs.push(PairOrientationData {
                u,
                v,
                orientations,
                signed_links,
            });
        }
        pairs.sort_unstable_by_key(|pair| (pair.u, pair.v));

        let mut incident_pairs = vec![Vec::new(); n];
        for (pair_index, pair) in pairs.iter().enumerate() {
            incident_pairs[pair.u].push(pair_index);
            incident_pairs[pair.v].push(pair_index);
        }

        Self {
            lengths,
            pairs,
            incident_pairs,
            golden_representatives: golden_representatives(),
        }
    }

    /// Equivalent to ALLHiC's `EvaluateQ`. Larger values are better.
    pub fn evaluate(&self, tour: &Tour<usize>) -> f64 {
        let Some((positions, signs_by_id, starts)) = self.layout(tour) else {
            return f64::NEG_INFINITY;
        };

        self.pairs
            .iter()
            .map(|pair| self.evaluate_pair(pair, &positions, &signs_by_id, &starts))
            .sum()
    }

    fn layout(&self, tour: &Tour<usize>) -> Option<(Vec<usize>, Vec<bool>, Vec<u64>)> {
        if tour.contigs.len() != tour.signs.len() {
            return None;
        }
        let n = self.lengths.len();
        let mut positions = vec![usize::MAX; n];
        let mut signs_by_id = vec![true; n];
        let mut starts = Vec::with_capacity(tour.contigs.len());
        let mut cumulative = 0u64;
        for (position, (&id, &sign)) in tour.contigs.iter().zip(&tour.signs).enumerate() {
            if id >= n || positions[id] != usize::MAX {
                return None;
            }
            positions[id] = position;
            signs_by_id[id] = sign;
            starts.push(cumulative);
            cumulative = cumulative.saturating_add(self.lengths[id]);
        }
        Some((positions, signs_by_id, starts))
    }

    fn evaluate_pair(
        &self,
        pair: &PairOrientationData,
        positions: &[usize],
        signs_by_id: &[bool],
        starts: &[u64],
    ) -> f64 {
        let pu = positions[pair.u];
        let pv = positions[pair.v];
        if pu == usize::MAX || pv == usize::MAX {
            return 0.0;
        }
        let (left, right, orientation) = if pu < pv {
            (
                pu,
                pv,
                orientation_index(signs_by_id[pair.u], signs_by_id[pair.v]),
            )
        } else {
            // ALLHiC stores the reverse contact with both orientations
            // flipped: (b, a, flip(bo), flip(ao)).
            (
                pv,
                pu,
                orientation_index(!signs_by_id[pair.u], !signs_by_id[pair.v]),
            )
        };
        let summary = &pair.orientations[orientation];
        if !summary.present {
            return 0.0;
        }
        let gap = starts[right - 1].saturating_sub(starts[left]);
        if gap > MAX_ORIENTATION_DISTANCE {
            return 0.0;
        }
        summary
            .bins
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .map(|(bin, &count)| {
                let distance = self.golden_representatives[bin].saturating_add(gap);
                -(count as f64) * (distance as f64).ln()
            })
            .sum()
    }

    /// ALLHiC-compatible spectral initialization followed by iterative
    /// flip-whole and flip-one refinement.
    pub fn optimize(&self, tour: &mut Tour<usize>) -> OrientationResult {
        self.initialize_spectral(tour);
        self.refine(tour)
    }

    /// De-novo orientation initialization corresponding to ALLHiC `flipAll`.
    /// Resume callers can omit this and call [`Self::refine`] directly.
    pub fn initialize_spectral(&self, tour: &mut Tour<usize>) -> bool {
        self.flip_all(tour)
    }

    /// Repeated ALLHiC `flipWhole` + `flipOne` refinement. This is also the
    /// orientation entry point for an already oriented/resumed tour.
    pub fn refine(&self, tour: &mut Tour<usize>) -> OrientationResult {
        let initial_score = self.evaluate(tour);
        let mut phases = 0;
        loop {
            phases += 1;
            let whole_accepted = self.flip_whole(tour);
            let one_accepted = self.flip_one(tour);
            if !whole_accepted && !one_accepted {
                break;
            }
        }
        OrientationResult {
            initial_score,
            final_score: self.evaluate(tour),
            phases,
        }
    }

    fn flip_all(&self, tour: &mut Tour<usize>) -> bool {
        let old_signs = tour.signs.clone();
        let old_score = self.evaluate(tour);
        let spectral = self.spectral_signs();
        for (position, &id) in tour.contigs.iter().enumerate() {
            tour.signs[position] = spectral[id];
        }
        let new_score = self.evaluate(tour);
        if new_score < old_score {
            tour.signs = old_signs;
            false
        } else {
            true
        }
    }

    fn flip_whole(&self, tour: &mut Tour<usize>) -> bool {
        let old_score = self.evaluate(tour);
        for sign in &mut tour.signs {
            *sign = !*sign;
        }
        if self.evaluate(tour) <= old_score {
            for sign in &mut tour.signs {
                *sign = !*sign;
            }
            false
        } else {
            true
        }
    }

    fn flip_one(&self, tour: &mut Tour<usize>) -> bool {
        let Some((positions, mut signs_by_id, starts)) = self.layout(tour) else {
            return false;
        };
        let mut accepted = false;
        let mut score = self.evaluate(tour);
        for position in 0..tour.contigs.len() {
            let id = tour.contigs[position];
            let old_contribution: f64 = self.incident_pairs[id]
                .iter()
                .map(|&pair_index| {
                    self.evaluate_pair(&self.pairs[pair_index], &positions, &signs_by_id, &starts)
                })
                .sum();
            tour.signs[position] = !tour.signs[position];
            signs_by_id[id] = !signs_by_id[id];
            let new_contribution: f64 = self.incident_pairs[id]
                .iter()
                .map(|&pair_index| {
                    self.evaluate_pair(&self.pairs[pair_index], &positions, &signs_by_id, &starts)
                })
                .sum();
            let new_score = score - old_contribution + new_contribution;
            if new_score > score {
                score = new_score;
                accepted = true;
            } else {
                tour.signs[position] = !tour.signs[position];
                signs_by_id[id] = !signs_by_id[id];
            }
        }
        accepted
    }

    fn spectral_signs(&self) -> Vec<bool> {
        let n = self.lengths.len();
        if n == 0 {
            return Vec::new();
        }
        let mut row_sums = vec![0.0; n];
        for pair in &self.pairs {
            let weight = pair.signed_links.abs();
            row_sums[pair.u] += weight;
            row_sums[pair.v] += weight;
        }
        let shift = row_sums.into_iter().fold(0.0_f64, f64::max);

        let mut vector: Vec<f64> = (1..=n).map(|value| value as f64).collect();
        let initial_norm = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
        for value in &mut vector {
            *value /= initial_norm;
        }
        let mut next = vec![0.0; n];
        for _ in 0..256 {
            for (next_value, &value) in next.iter_mut().zip(&vector) {
                *next_value = shift * value;
            }
            // Adding a scalar diagonal shift preserves eigenvectors while
            // allowing power iteration over the sparse signed O matrix.
            for pair in &self.pairs {
                next[pair.u] += pair.signed_links * vector[pair.v];
                next[pair.v] += pair.signed_links * vector[pair.u];
            }
            let norm = next.iter().map(|value| value * value).sum::<f64>().sqrt();
            if norm == 0.0 {
                break;
            }
            for value in &mut next {
                *value /= norm;
            }
            let delta = next
                .iter()
                .zip(&vector)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            std::mem::swap(&mut vector, &mut next);
            if delta < 1e-12 {
                break;
            }
        }
        vector.into_iter().map(|value| value >= 0.0).collect()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct OrientationResult {
    pub initial_score: f64,
    pub final_score: f64,
    pub phases: usize,
}

fn orientation_index(a_forward: bool, b_forward: bool) -> usize {
    (usize::from(!a_forward) << 1) | usize::from(!b_forward)
}

fn golden_bin(distance: u64) -> usize {
    let exponent = ((distance as f64).ln() / LOG_GOLDEN_RATIO).round() as i32;
    exponent.clamp(GOLDEN_LOWER_BOUND, GOLDEN_UPPER_BOUND) as usize - GOLDEN_LOWER_BOUND as usize
}

fn golden_representatives() -> GoldenArrayRepresentatives {
    std::array::from_fn(|bin| {
        let exponent = GOLDEN_LOWER_BOUND + bin as i32;
        ((LOG_GOLDEN_RATIO * exponent as f64).exp().round() as u64).clamp(2_048, 1_u64 << 32)
    })
}

#[derive(Debug, Clone)]
pub struct OptimizeConfig {
    pub population_size: usize,
    /// Stop after strictly more than this many generations without improvement,
    /// matching ALLHiC's `generation - updated > ngen` condition.
    pub stale_generations: usize,
    pub max_generations: usize,
    pub mutation_probability: f64,
    pub seed: u64,
    pub phases: usize,
    pub report_interval: usize,
    pub backbone: BackboneConfig,
}

impl Default for OptimizeConfig {
    fn default() -> Self {
        Self {
            population_size: DEFAULT_POPULATION_SIZE,
            stale_generations: DEFAULT_STALE_GENERATIONS,
            max_generations: DEFAULT_MAX_GENERATIONS,
            mutation_probability: DEFAULT_MUTATION_PROBABILITY,
            seed: 42,
            phases: 2,
            report_interval: 500,
            backbone: BackboneConfig::default(),
        }
    }
}

impl OptimizeConfig {
    fn validate(&self, contig_count: usize) -> Result<(), String> {
        if self.population_size < 4 {
            return Err("population size must be at least 4 for tournament-3 selection".into());
        }
        if !(0.0..=1.0).contains(&self.mutation_probability) {
            return Err("mutation probability must be between 0 and 1".into());
        }
        if self.phases == 0 {
            return Err("at least one ordering phase is required".into());
        }
        if contig_count == 0 {
            return Err("cannot optimize an empty tour".into());
        }
        self.backbone.validate()?;
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct GenerationReport {
    pub phase: usize,
    pub generation: usize,
    pub best_score: f64,
}

#[derive(Debug, Clone)]
pub struct OptimizeResult {
    pub tour: Tour<usize>,
    pub initial_score: f64,
    pub final_score: f64,
    pub generations: Vec<usize>,
    pub reports: Vec<GenerationReport>,
    pub backbone: Option<BackboneReport>,
}

#[derive(Debug, Clone)]
pub struct BackboneReport {
    pub used: bool,
    pub contig_count: usize,
    pub unit_count: usize,
    pub block_count: usize,
    pub accepted_edges: usize,
    pub seed_score: Option<f64>,
    pub coarse_score: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct AllhicOptimizeResult {
    pub tour: Tour<usize>,
    pub ordering: OptimizeResult,
    pub orientation: OrientationResult,
}

/// Run the functional ALLHiC optimize pipeline without imposing ALLHiC's
/// RE/CLM/tour file formats on callers.
pub fn optimize_allhic(
    initial_tour: &Tour<usize>,
    ordering_problem: &OptimizeProblem,
    orientation_problem: &OrientationProblem,
    config: &OptimizeConfig,
) -> Result<AllhicOptimizeResult, String> {
    let ordering = optimize_order(initial_tour, ordering_problem, config)?;
    let mut tour = ordering.tour.clone();
    let orientation = orientation_problem.optimize(&mut tour);
    Ok(AllhicOptimizeResult {
        tour,
        ordering,
        orientation,
    })
}

#[derive(Clone)]
struct Individual {
    order: Vec<usize>,
    score: f64,
    dirty: bool,
}

#[derive(Clone)]
struct UnitIndividual {
    genes: Vec<UnitGene>,
    score: f64,
    dirty: bool,
}

/// Run ordering with a conservative path-block phase followed by the ordinary
/// fully unlocked ALLHiC GA. Set `config.backbone.enabled = false` to recover
/// the mutation-only ALLHiC-compatible path exactly.
pub fn optimize_order(
    initial_tour: &Tour<usize>,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
) -> Result<OptimizeResult, String> {
    config.validate(problem.contig_count())?;
    validate_tour(&initial_tour.contigs, problem.contig_count())?;
    if initial_tour.signs.len() != problem.contig_count() {
        return Err(format!(
            "tour sign count ({}) does not match contig count ({})",
            initial_tour.signs.len(),
            problem.contig_count()
        ));
    }

    if config.backbone.enabled && config.phases >= 2 {
        let plan = BackbonePlan::discover(problem, &config.backbone);
        if plan.is_useful() {
            return optimize_order_multilevel(initial_tour, problem, config, plan);
        }
        log::info!(
            "Backbone GA fallback: {} blocks / {} units for {} contigs do not meet the reduction threshold",
            plan.block_count(),
            plan.unit_count(),
            problem.contig_count()
        );
        return Ok(optimize_order_standard(
            initial_tour,
            problem,
            config,
            Some(BackboneReport {
                used: false,
                contig_count: problem.contig_count(),
                unit_count: plan.unit_count(),
                block_count: plan.block_count(),
                accepted_edges: plan.accepted_edges(),
                seed_score: None,
                coarse_score: None,
            }),
        ));
    }

    Ok(optimize_order_standard(initial_tour, problem, config, None))
}

fn optimize_order_standard(
    initial_tour: &Tour<usize>,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    backbone: Option<BackboneReport>,
) -> OptimizeResult {
    let initial_score = problem.evaluate(&initial_tour.contigs);
    let mut best_order = initial_tour.contigs.clone();
    let mut best_score = initial_score;
    let mut rng = SmallRng::seed_from_u64(config.seed);
    let mut reports = Vec::new();
    let mut generations = Vec::with_capacity(config.phases);

    for phase in 1..=config.phases {
        let outcome = run_phase(&best_order, problem, config, phase, &mut rng, &mut reports);
        generations.push(outcome.generations);
        if outcome.score < best_score {
            best_score = outcome.score;
            best_order = outcome.order;
        }
    }

    build_optimize_result(
        initial_tour,
        problem.contig_count(),
        initial_score,
        best_order,
        best_score,
        generations,
        reports,
        backbone,
    )
}

fn optimize_order_multilevel(
    initial_tour: &Tour<usize>,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    plan: BackbonePlan,
) -> Result<OptimizeResult, String> {
    let initial_score = problem.evaluate(&initial_tour.contigs);
    let seed_genes = plan.initial_genes(&initial_tour.contigs);
    let seed_order = plan.decode(&seed_genes)?;
    let seed_score = problem.evaluate(&seed_order);
    let mut rng = SmallRng::seed_from_u64(config.seed);
    let mut reports = Vec::new();
    let mut generations = Vec::with_capacity(config.phases);

    log::info!(
        "Backbone GA enabled: {} path blocks, {} units, {} accepted edges for {} contigs",
        plan.block_count(),
        plan.unit_count(),
        plan.accepted_edges(),
        problem.contig_count()
    );
    let coarse = run_unit_phase(
        &seed_genes,
        &plan,
        problem,
        config,
        1,
        &mut rng,
        &mut reports,
    );
    generations.push(coarse.generations);

    let coarse_order = plan.decode(&coarse.genes)?;
    let mut best_order = initial_tour.contigs.clone();
    let mut best_score = initial_score;
    if coarse.score < best_score {
        best_order = coarse_order;
        best_score = coarse.score;
    }

    // Every later phase is the ordinary contig-level ALLHiC GA. All temporary
    // path-block constraints are therefore removable by the original four
    // mutation operators.
    for phase in 2..=config.phases {
        let outcome = run_phase(&best_order, problem, config, phase, &mut rng, &mut reports);
        generations.push(outcome.generations);
        if outcome.score < best_score {
            best_score = outcome.score;
            best_order = outcome.order;
        }
    }

    Ok(build_optimize_result(
        initial_tour,
        problem.contig_count(),
        initial_score,
        best_order,
        best_score,
        generations,
        reports,
        Some(BackboneReport {
            used: true,
            contig_count: problem.contig_count(),
            unit_count: plan.unit_count(),
            block_count: plan.block_count(),
            accepted_edges: plan.accepted_edges(),
            seed_score: Some(seed_score),
            coarse_score: Some(coarse.score),
        }),
    ))
}

#[allow(clippy::too_many_arguments)]
fn build_optimize_result(
    initial_tour: &Tour<usize>,
    contig_count: usize,
    initial_score: f64,
    mut best_order: Vec<usize>,
    mut best_score: f64,
    generations: Vec<usize>,
    reports: Vec<GenerationReport>,
    backbone: Option<BackboneReport>,
) -> OptimizeResult {
    // ALLHiC keeps the initial tour when GA fails to improve it.
    if best_score > initial_score {
        best_score = initial_score;
        best_order.clone_from(&initial_tour.contigs);
    }

    let mut sign_by_id = vec![true; contig_count];
    for (&id, &sign) in initial_tour.contigs.iter().zip(&initial_tour.signs) {
        sign_by_id[id] = sign;
    }
    let signs = best_order.iter().map(|&id| sign_by_id[id]).collect();

    OptimizeResult {
        tour: Tour {
            contigs: best_order,
            signs,
        },
        initial_score,
        final_score: best_score,
        generations,
        reports,
        backbone,
    }
}

struct UnitPhaseOutcome {
    genes: Vec<UnitGene>,
    score: f64,
    generations: usize,
}

fn run_unit_phase(
    seed_genes: &[UnitGene],
    plan: &BackbonePlan,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    phase: usize,
    rng: &mut SmallRng,
    reports: &mut Vec<GenerationReport>,
) -> UnitPhaseOutcome {
    let seed_order = plan
        .decode(seed_genes)
        .expect("validated backbone seed must decode");
    let seed_score = problem.evaluate(&seed_order);
    let mut population = vec![
        UnitIndividual {
            genes: seed_genes.to_vec(),
            score: seed_score,
            dirty: false,
        };
        config.population_size
    ];
    let mut offspring = population.clone();
    let mut hall_of_fame = population[0].clone();
    let mut last_improvement = 0usize;
    let mut completed = 0usize;

    for generation in 1..=config.max_generations {
        let mut offspring_index = 0;
        while offspring_index < config.population_size {
            let (first, second) = unit_tournament_pair(&population, rng);
            offspring[offspring_index]
                .genes
                .clone_from(&population[first].genes);
            offspring[offspring_index].score = population[first].score;
            offspring[offspring_index].dirty = false;
            offspring_index += 1;
            if offspring_index < config.population_size {
                offspring[offspring_index]
                    .genes
                    .clone_from(&population[second].genes);
                offspring[offspring_index].score = population[second].score;
                offspring[offspring_index].dirty = false;
                offspring_index += 1;
            }
        }

        for individual in &mut offspring {
            if rng.gen_bool(config.mutation_probability) {
                mutate_units(&mut individual.genes, plan, rng);
                individual.dirty = true;
            }
        }
        offspring
            .par_iter_mut()
            .filter(|individual| individual.dirty)
            .for_each_init(
                || {
                    (
                        Vec::with_capacity(problem.contig_count()),
                        vec![0.0; problem.contig_count()],
                    )
                },
                |(decoded, midpoints), individual| {
                    plan.decode_into(&individual.genes, decoded);
                    individual.score = problem.evaluate_with_buffer(decoded, midpoints);
                    individual.dirty = false;
                },
            );
        offspring.sort_unstable_by(|a, b| a.score.total_cmp(&b.score));
        std::mem::swap(&mut population, &mut offspring);
        completed = generation;

        if population[0].score < hall_of_fame.score {
            hall_of_fame = population[0].clone();
            last_improvement = generation;
        }
        if config.report_interval > 0 && generation % config.report_interval == 0 {
            log::info!(
                "Current iteration BackboneGA{}-{}: max_score={:.5}",
                phase,
                generation,
                -hall_of_fame.score
            );
            reports.push(GenerationReport {
                phase,
                generation,
                best_score: hall_of_fame.score,
            });
        }
        if generation - last_improvement > config.stale_generations {
            break;
        }
    }

    UnitPhaseOutcome {
        genes: hall_of_fame.genes,
        score: hall_of_fame.score,
        generations: completed,
    }
}

struct PhaseOutcome {
    order: Vec<usize>,
    score: f64,
    generations: usize,
}

fn run_phase(
    seed_order: &[usize],
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    phase: usize,
    rng: &mut SmallRng,
    reports: &mut Vec<GenerationReport>,
) -> PhaseOutcome {
    // ALLHiC's MakeTour clones the same seed for every initial individual.
    let seed_score = problem.evaluate(seed_order);
    let mut population = vec![
        Individual {
            order: seed_order.to_vec(),
            score: seed_score,
            dirty: false,
        };
        config.population_size
    ];
    // Reuse both the Individual and chromosome Vec allocations across
    // generations. Fragmented assemblies otherwise allocate and free
    // population_size full tours on every generation.
    let mut offspring = population.clone();
    let mut hall_of_fame = population[0].clone();
    let mut last_improvement = 0usize;
    let mut completed = 0usize;

    for generation in 1..=config.max_generations {
        let mut offspring_index = 0;
        while offspring_index < config.population_size {
            let (first, second) = tournament_pair(&population, rng);
            offspring[offspring_index]
                .order
                .clone_from(&population[first].order);
            offspring[offspring_index].score = population[first].score;
            offspring[offspring_index].dirty = false;
            offspring_index += 1;
            if offspring_index < config.population_size {
                offspring[offspring_index]
                    .order
                    .clone_from(&population[second].order);
                offspring[offspring_index].score = population[second].score;
                offspring[offspring_index].dirty = false;
                offspring_index += 1;
            }
        }

        for individual in &mut offspring {
            if rng.gen_bool(config.mutation_probability) {
                mutate_allhic(&mut individual.order, rng);
                individual.dirty = true;
            }
        }
        offspring
            .par_iter_mut()
            .filter(|individual| individual.dirty)
            .for_each_init(
                || vec![0.0; problem.contig_count()],
                |midpoints, individual| {
                    individual.score = problem.evaluate_with_buffer(&individual.order, midpoints);
                    individual.dirty = false;
                },
            );
        offspring.sort_unstable_by(|a, b| a.score.total_cmp(&b.score));
        std::mem::swap(&mut population, &mut offspring);
        completed = generation;

        if population[0].score < hall_of_fame.score {
            hall_of_fame = population[0].clone();
            last_improvement = generation;
        }
        if config.report_interval > 0 && generation % config.report_interval == 0 {
            log::info!(
                "Current iteration GA{}-{}: max_score={:.5}",
                phase,
                generation,
                -hall_of_fame.score
            );
            reports.push(GenerationReport {
                phase,
                generation,
                best_score: hall_of_fame.score,
            });
        }
        if generation - last_improvement > config.stale_generations {
            break;
        }
    }

    PhaseOutcome {
        order: hall_of_fame.order,
        score: hall_of_fame.score,
        generations: completed,
    }
}

fn tournament_pair(population: &[Individual], rng: &mut SmallRng) -> (usize, usize) {
    let first = tournament_winner(population, rng, None);
    let second = tournament_winner(population, rng, Some(first));
    (first, second)
}

fn tournament_winner(
    population: &[Individual],
    rng: &mut SmallRng,
    excluded: Option<usize>,
) -> usize {
    let mut candidates = [usize::MAX; 3];
    for index in 0..candidates.len() {
        loop {
            let candidate = rng.gen_range(0..population.len());
            if Some(candidate) != excluded && !candidates[..index].contains(&candidate) {
                candidates[index] = candidate;
                break;
            }
        }
    }
    let mut winner = candidates[0];
    for candidate in candidates.into_iter().skip(1) {
        if population[candidate].score < population[winner].score {
            winner = candidate;
        }
    }
    winner
}

fn unit_tournament_pair(population: &[UnitIndividual], rng: &mut SmallRng) -> (usize, usize) {
    let first = unit_tournament_winner(population, rng, None);
    let second = unit_tournament_winner(population, rng, Some(first));
    (first, second)
}

fn unit_tournament_winner(
    population: &[UnitIndividual],
    rng: &mut SmallRng,
    excluded: Option<usize>,
) -> usize {
    let mut candidates = [usize::MAX; 3];
    for index in 0..candidates.len() {
        loop {
            let candidate = rng.gen_range(0..population.len());
            if Some(candidate) != excluded && !candidates[..index].contains(&candidate) {
                candidates[index] = candidate;
                break;
            }
        }
    }
    let mut winner = candidates[0];
    for candidate in candidates.into_iter().skip(1) {
        if population[candidate].score < population[winner].score {
            winner = candidate;
        }
    }
    winner
}

fn mutate_units<R: Rng + ?Sized>(genes: &mut [UnitGene], plan: &BackbonePlan, rng: &mut R) {
    let n = genes.len();
    if n <= 1 {
        if n == 1 && plan.units[genes[0].unit].len() > 1 {
            genes[0].reversed = !genes[0].reversed;
        }
        return;
    }

    let strategy = rng.r#gen::<f64>();
    if strategy < 0.2 {
        let p = rng.gen_range(0..n);
        let q = rng.gen_range(0..n);
        genes.swap(p, q);
    } else if strategy < 0.4 {
        let k = rng.gen_range(1..n);
        genes.rotate_left(k);
    } else if strategy < 0.7 {
        // Prefer moving a residual singleton; path blocks remain atomic.
        let mut p = rng.gen_range(0..n);
        for _ in 0..4 {
            let candidate = rng.gen_range(0..n);
            if plan.units[genes[candidate].unit].len() == 1 {
                p = candidate;
                break;
            }
        }
        let q = rng.gen_range(0..n);
        if p < q {
            genes[p..=q].rotate_left(1);
        } else if q < p {
            genes[q..=p].rotate_right(1);
        }
    } else if strategy < 0.9 {
        let mut p = rng.gen_range(0..n);
        let mut q = rng.gen_range(0..n);
        if p > q {
            std::mem::swap(&mut p, &mut q);
        }
        if p < q {
            genes[p..=q].reverse();
            for gene in &mut genes[p..=q] {
                gene.reversed = !gene.reversed;
            }
        }
    } else {
        for _ in 0..8 {
            let p = rng.gen_range(0..n);
            if plan.units[genes[p].unit].len() > 1 {
                genes[p].reversed = !genes[p].reversed;
                break;
            }
        }
    }
}

/// ALLHiC mutation mixture: 20% swap, 20% splice, 30% insertion,
/// and 30% inversion.
pub fn mutate_allhic<R: Rng + ?Sized>(order: &mut [usize], rng: &mut R) {
    let n = order.len();
    if n <= 1 {
        return;
    }
    let strategy = rng.r#gen::<f64>();
    if strategy < 0.2 {
        let p = rng.gen_range(0..n);
        let q = rng.gen_range(0..n);
        order.swap(p, q);
    } else if strategy < 0.4 {
        let k = rng.gen_range(1..n);
        order.rotate_left(k);
    } else {
        let mut p = rng.gen_range(0..n);
        let mut q = rng.gen_range(0..n);
        if p > q {
            std::mem::swap(&mut p, &mut q);
        }
        if p == q {
            return;
        }
        if strategy < 0.7 {
            if rng.gen_bool(0.5) {
                order[p..=q].rotate_right(1);
            } else {
                order[p..=q].rotate_left(1);
            }
        } else {
            order[p..=q].reverse();
        }
    }
}

fn validate_tour(order: &[usize], contig_count: usize) -> Result<(), String> {
    if order.is_empty() {
        return Err("tour cannot be empty".into());
    }
    if order.len() != contig_count {
        return Err(format!(
            "tour contig count ({}) does not match problem contig count ({contig_count})",
            order.len()
        ));
    }
    let mut seen = vec![false; contig_count];
    for &id in order {
        if id >= contig_count {
            return Err(format!("tour contig id {id} is outside 0..{contig_count}"));
        }
        if std::mem::replace(&mut seen[id], true) {
            return Err(format!("tour contains duplicate contig id {id}"));
        }
    }
    Ok(())
}
