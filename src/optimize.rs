use crate::clm::ClmbReader;
use crate::order::{HierarchicalJoinEvidence, Tour};
use crate::splitcontacts::SplitContacts;
use hashbrown::HashMap;
use indexmap::IndexMap;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rayon::prelude::*;
use std::ffi::OsStr;
use std::path::Path;
use std::sync::{Arc, OnceLock};
use std::time::{Duration, Instant};

pub const DEFAULT_POPULATION_SIZE: usize = 100;
pub const DEFAULT_STALE_GENERATIONS: usize = 5_000;
pub const DEFAULT_MUTATION_PROBABILITY: f64 = 0.2;
pub const DEFAULT_MAX_GENERATIONS: usize = 1_000_000;
const HIERARCHICAL_REFINED_BLOCK_SIZE: usize = 8;
/// Endpoint competition is distinct from the whole-CLM top-two/top-three
/// margin in `BackboneConfig`, so keep its empirically selected gate separate.
const HIERARCHICAL_MIN_CONFIDENCE: f64 = 1.1;
const HIERARCHICAL_MIN_SIGNIFICANT_RELATIVE_GAIN: f64 = 5e-4;
const HIERARCHICAL_BLOCK_PATIENCE: usize = 2_000;
const HIERARCHICAL_FULL_PATIENCE: usize = 5_000;
const ELITE_COUNT: usize = 1;
const GOLDEN_LOWER_BOUND: i32 = 16;
const GOLDEN_UPPER_BOUND: i32 = 50;
const GOLDEN_BINS: usize = (GOLDEN_UPPER_BOUND - GOLDEN_LOWER_BOUND + 1) as usize;
const LOG_GOLDEN_RATIO: f64 = 0.481_211_825_059_668_4;
const MAX_ORIENTATION_DISTANCE: u64 = 500_000_000;
/// Smaller batches run serially to avoid paying a Rayon barrier on every GA
/// generation. Once a batch crosses this threshold, scoring is balanced over
/// all useful workers because each edge evaluation is relatively expensive.
const MIN_PARALLEL_EVALUATION_WORK: usize = 65_536;
/// Bound persistent per-individual scratch on exceptionally fragmented tours
/// or unusually large populations. Above this budget we retain one reusable
/// scratch buffer per parallel chunk instead.
const MAX_PERSISTENT_EVALUATION_SCRATCH_BYTES: usize = 256 * 1024 * 1024;
/// Reuse chromosome allocations without allowing the pool to dominate memory
/// beside the independent evaluation scratch budget.
const MAX_CHROMOSOME_BUFFER_POOL_BYTES: usize = 64 * 1024 * 1024;
const ANCHOR_COVERAGE_FRACTION: f64 = 0.8;
const PROFILE_GA_ENV: &str = "CPHASING_PROFILE_GA";
const DISABLE_INTERLEAVED_FITNESS_ENV: &str = "CPHASING_DISABLE_INTERLEAVED_FITNESS";
const AFFECTED_EDGE_SAMPLE_INTERVAL: u32 = 16;

pub type GoldenArray = [u32; GOLDEN_BINS];

#[derive(Debug, Clone, Copy)]
struct ContactEdge {
    u: u32,
    v: u32,
    links: f64,
}

#[inline(always)]
fn contact_edge_distance(midpoints: &[f64], edge: &ContactEdge) -> f64 {
    let u = edge.u as usize;
    let v = edge.v as usize;
    debug_assert!(u < midpoints.len() && v < midpoints.len());
    // OptimizeProblem validates every endpoint when it builds the immutable
    // edge list. Avoid repeating two bounds checks for every fitness edge.
    unsafe { (*midpoints.get_unchecked(u) - *midpoints.get_unchecked(v)).abs() }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx")]
unsafe fn reciprocal_scores_avx4(edges: &[ContactEdge], midpoints: [&[f64]; 4]) -> [f64; 4] {
    use std::arch::x86_64::{
        __m256d, _mm256_div_pd, _mm256_set_pd, _mm256_set1_pd, _mm256_setzero_pd, _mm256_storeu_pd,
        _mm256_sub_pd,
    };

    let mut scores: __m256d = _mm256_setzero_pd();
    for edge in edges {
        let u = edge.u as usize;
        let v = edge.v as usize;
        debug_assert!(
            midpoints
                .iter()
                .all(|positions| u < positions.len() && v < positions.len())
        );
        let d0 = unsafe { (*midpoints[0].get_unchecked(u) - *midpoints[0].get_unchecked(v)).abs() };
        let d1 = unsafe { (*midpoints[1].get_unchecked(u) - *midpoints[1].get_unchecked(v)).abs() };
        let d2 = unsafe { (*midpoints[2].get_unchecked(u) - *midpoints[2].get_unchecked(v)).abs() };
        let d3 = unsafe { (*midpoints[3].get_unchecked(u) - *midpoints[3].get_unchecked(v)).abs() };
        let distances = _mm256_set_pd(d3, d2, d1, d0);
        let links = _mm256_set1_pd(edge.links);
        scores = _mm256_sub_pd(scores, _mm256_div_pd(links, distances));
    }
    let mut result = [0.0; 4];
    unsafe { _mm256_storeu_pd(result.as_mut_ptr(), scores) };
    result
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx")]
unsafe fn reciprocal_scores_interleaved_avx4(
    edges: &[ContactEdge],
    interleaved_midpoints: &[f64],
) -> [f64; 4] {
    use std::arch::x86_64::{
        __m256d, _mm256_andnot_pd, _mm256_div_pd, _mm256_loadu_pd, _mm256_set1_pd,
        _mm256_setzero_pd, _mm256_storeu_pd, _mm256_sub_pd,
    };

    let mut scores: __m256d = _mm256_setzero_pd();
    let sign_bits = _mm256_set1_pd(-0.0);
    for edge in edges {
        let u_offset = edge.u as usize * 4;
        let v_offset = edge.v as usize * 4;
        debug_assert!(u_offset + 4 <= interleaved_midpoints.len());
        debug_assert!(v_offset + 4 <= interleaved_midpoints.len());
        let u_positions = unsafe { _mm256_loadu_pd(interleaved_midpoints.as_ptr().add(u_offset)) };
        let v_positions = unsafe { _mm256_loadu_pd(interleaved_midpoints.as_ptr().add(v_offset)) };
        let distances = _mm256_andnot_pd(sign_bits, _mm256_sub_pd(u_positions, v_positions));
        let links = _mm256_set1_pd(edge.links);
        scores = _mm256_sub_pd(scores, _mm256_div_pd(links, distances));
    }
    let mut result = [0.0; 4];
    unsafe { _mm256_storeu_pd(result.as_mut_ptr(), scores) };
    result
}

fn avx_fitness_batch_supported() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::arch::is_x86_feature_detected!("avx")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

fn interleaved_fitness_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os(DISABLE_INTERLEAVED_FITNESS_ENV).is_none())
}

#[derive(Debug, Clone, Copy)]
struct EndpointOrderEdge {
    u: u32,
    v: u32,
    u_before_v: f64,
    v_before_u: f64,
}

#[derive(Debug, Clone)]
struct EndpointMultiscaleData {
    edges: Vec<EndpointOrderEdge>,
    total_best_score: f64,
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub enum OrderObjective {
    /// Project ALLHiC's default (`--logDist=false`).
    #[default]
    ReciprocalDistance,
    /// Project ALLHiC's optional `--logDist` objective.
    LogDistance,
    /// Length-aware, multi-scale objective. Long contigs covering 80% of the
    /// assembly contribute 70%, anchor-fragment placement contributes 20%,
    /// and fragment-fragment ordering contributes 10%.
    LengthTiered,
    /// Orientation-aware endpoint support over anchor rank distances
    /// 1/2/4/8, with a small whole-contig term to avoid fragment drift.
    EndpointMultiscale,
}

#[derive(Debug, Clone)]
pub struct OptimizeProblem {
    lengths: Vec<f64>,
    edges: Vec<ContactEdge>,
    objective: OrderObjective,
    anchors: Vec<bool>,
    endpoint_multiscale: Option<EndpointMultiscaleData>,
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
        let anchors = select_length_anchors(&dense_lengths);

        Ok(Self {
            lengths: dense_lengths,
            edges,
            objective: OrderObjective::default(),
            anchors,
            endpoint_multiscale: None,
        })
    }

    pub fn with_objective(mut self, objective: OrderObjective) -> Self {
        self.objective = objective;
        self
    }

    pub fn with_endpoint_multiscale(
        mut self,
        split_contacts: &SplitContacts,
        contig_to_id: &HashMap<String, usize>,
        tour: &Tour<usize>,
    ) -> Result<Self, String> {
        validate_tour(&tour.contigs, self.contig_count())?;
        if tour.signs.len() != self.contig_count() {
            return Err("endpoint objective requires one sign per contig".into());
        }
        let mut signs = vec![true; self.contig_count()];
        for (&id, &sign) in tour.contigs.iter().zip(&tour.signs) {
            signs[id] = sign;
        }
        let mut endpoint_totals = vec![0.0; self.contig_count() * 2];
        for (pair, counts) in &split_contacts.data {
            let (Some(&u), Some(&v)) = (
                contig_to_id.get(&pair.Contig1),
                contig_to_id.get(&pair.Contig2),
            ) else {
                continue;
            };
            if counts.len() < 4 || u == v {
                continue;
            }
            for u_end in 0..2 {
                for v_end in 0..2 {
                    let count = counts[u_end * 2 + v_end].max(0.0);
                    endpoint_totals[2 * u + u_end] += count;
                    endpoint_totals[2 * v + v_end] += count;
                }
            }
        }
        let mut edges = Vec::new();
        let mut total_best_score = 0.0;
        for (pair, counts) in &split_contacts.data {
            let (Some(&u), Some(&v)) = (
                contig_to_id.get(&pair.Contig1),
                contig_to_id.get(&pair.Contig2),
            ) else {
                continue;
            };
            if counts.len() < 4 || u == v || !self.anchors[u] || !self.anchors[v] {
                continue;
            }
            let u_left = usize::from(!signs[u]);
            let u_right = 1 - u_left;
            let v_left = usize::from(!signs[v]);
            let v_right = 1 - v_left;
            let normalized = |u_end: usize, v_end: usize| {
                let count = counts[u_end * 2 + v_end].max(0.0);
                let denominator =
                    (endpoint_totals[2 * u + u_end] * endpoint_totals[2 * v + v_end]).sqrt();
                if denominator > 0.0 {
                    count / denominator
                } else {
                    0.0
                }
            };
            let u_before_v = normalized(u_right, v_left);
            let v_before_u = normalized(u_left, v_right);
            let best = u_before_v.max(v_before_u);
            if best > 0.0 {
                edges.push(EndpointOrderEdge {
                    u: u as u32,
                    v: v as u32,
                    u_before_v,
                    v_before_u,
                });
                total_best_score += best;
            }
        }
        self.endpoint_multiscale = Some(EndpointMultiscaleData {
            edges,
            total_best_score,
        });
        self.objective = OrderObjective::EndpointMultiscale;
        Ok(self)
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
        self.evaluate_with_buffer_profiled(order, midpoints, None)
    }

    fn evaluate_with_buffer_profiled(
        &self,
        order: &[usize],
        midpoints: &mut [f64],
        mut timings: Option<&mut FitnessDetailTimings>,
    ) -> f64 {
        if order.len() != self.lengths.len() || midpoints.len() < self.lengths.len() {
            return f64::INFINITY;
        }
        let layout_started = timings.is_some().then(Instant::now);
        let mut cumulative = 0.0;

        for &id in order {
            if id >= self.lengths.len() {
                return f64::INFINITY;
            }
            let length = self.lengths[id];
            midpoints[id] = cumulative + length * 0.5;
            cumulative += length;
        }
        if let (Some(timings), Some(layout_started)) = (timings.as_deref_mut(), layout_started) {
            timings.layout += layout_started.elapsed();
        }

        let scoring_started = timings.is_some().then(Instant::now);
        let score = match self.objective {
            OrderObjective::ReciprocalDistance => {
                let mut score = 0.0;
                for edge in &self.edges {
                    let distance = contact_edge_distance(midpoints, edge);
                    score -= edge.links / distance;
                }
                score
            }
            OrderObjective::LogDistance => {
                let mut score = 0.0;
                for edge in &self.edges {
                    let distance = contact_edge_distance(midpoints, edge);
                    let distance = if distance <= 1.0 { 1.000001 } else { distance };
                    score += edge.links * distance.ln();
                }
                score
            }
            OrderObjective::LengthTiered => {
                let components = self.length_tiered_components(order, midpoints);
                encode_tiered_score(components)
            }
            OrderObjective::EndpointMultiscale => {
                let endpoint_cost = self.endpoint_multiscale_cost(order);
                let whole_contig_cost =
                    encode_tiered_score(self.length_tiered_components(order, midpoints));
                0.9 * endpoint_cost + 0.1 * whole_contig_cost
            }
        };
        if let (Some(timings), Some(scoring_started)) = (timings.as_deref_mut(), scoring_started) {
            timings.scoring += scoring_started.elapsed();
            timings.evaluations += 1;
            timings.contact_edges += self.edges.len() as u64;
        }
        score
    }

    fn evaluate_reciprocal_batch4(
        &self,
        orders: [&[usize]; 4],
        midpoints: [&mut [f64]; 4],
        mut timings: Option<&mut FitnessDetailTimings>,
    ) -> Option<[f64; 4]> {
        if self.objective != OrderObjective::ReciprocalDistance || !avx_fitness_batch_supported() {
            return None;
        }
        let layout_started = timings.is_some().then(Instant::now);
        let [positions0, positions1, positions2, positions3] = midpoints;
        for (order, positions) in orders.into_iter().zip([
            &mut *positions0,
            &mut *positions1,
            &mut *positions2,
            &mut *positions3,
        ]) {
            if order.len() != self.lengths.len() || positions.len() < self.lengths.len() {
                return None;
            }
            let mut cumulative = 0.0;
            for &id in order {
                if id >= self.lengths.len() {
                    return None;
                }
                let length = self.lengths[id];
                positions[id] = cumulative + length * 0.5;
                cumulative += length;
            }
        }
        if let (Some(timings), Some(layout_started)) = (timings.as_deref_mut(), layout_started) {
            timings.layout += layout_started.elapsed();
        }

        let scoring_started = timings.is_some().then(Instant::now);
        #[cfg(target_arch = "x86_64")]
        let scores = {
            // Runtime AVX detection above protects this call. Every lane keeps
            // the original edge order and one sequential score accumulator.
            Some(unsafe {
                reciprocal_scores_avx4(
                    &self.edges,
                    [positions0, positions1, positions2, positions3],
                )
            })
        };
        #[cfg(not(target_arch = "x86_64"))]
        let scores = None;
        if let (Some(timings), Some(scoring_started)) = (timings.as_deref_mut(), scoring_started) {
            timings.scoring += scoring_started.elapsed();
            timings.evaluations += 4;
            timings.simd_batches += 1;
            timings.contact_edges += (self.edges.len() as u64).saturating_mul(4);
        }
        scores
    }

    fn evaluate_reciprocal_interleaved_batch4(
        &self,
        orders: [&[usize]; 4],
        interleaved_midpoints: &mut [f64],
        mut timings: Option<&mut FitnessDetailTimings>,
    ) -> Option<[f64; 4]> {
        if self.objective != OrderObjective::ReciprocalDistance || !avx_fitness_batch_supported() {
            return None;
        }
        let required_midpoints = self.lengths.len().checked_mul(4)?;
        if interleaved_midpoints.len() < required_midpoints {
            return None;
        }

        let layout_started = timings.is_some().then(Instant::now);
        for (lane, order) in orders.into_iter().enumerate() {
            if order.len() != self.lengths.len() {
                return None;
            }
            let mut cumulative = 0.0;
            for &id in order {
                if id >= self.lengths.len() {
                    return None;
                }
                let length = self.lengths[id];
                interleaved_midpoints[id * 4 + lane] = cumulative + length * 0.5;
                cumulative += length;
            }
        }
        if let (Some(timings), Some(layout_started)) = (timings.as_deref_mut(), layout_started) {
            timings.layout += layout_started.elapsed();
        }

        let scoring_started = timings.is_some().then(Instant::now);
        #[cfg(target_arch = "x86_64")]
        let scores = {
            // Runtime AVX detection above protects this call. Contig-major
            // storage turns the four midpoint gathers into two vector loads.
            Some(unsafe { reciprocal_scores_interleaved_avx4(&self.edges, interleaved_midpoints) })
        };
        #[cfg(not(target_arch = "x86_64"))]
        let scores = None;
        if let (Some(timings), Some(scoring_started)) = (timings.as_deref_mut(), scoring_started) {
            timings.scoring += scoring_started.elapsed();
            timings.evaluations += 4;
            timings.simd_batches += 1;
            timings.contact_edges += (self.edges.len() as u64).saturating_mul(4);
        }
        scores
    }

    fn endpoint_multiscale_cost(&self, order: &[usize]) -> f64 {
        let Some(endpoint) = &self.endpoint_multiscale else {
            return f64::INFINITY;
        };
        if endpoint.total_best_score <= 0.0 {
            return 1.0;
        }
        let mut anchor_rank = vec![usize::MAX; self.contig_count()];
        let mut rank = 0usize;
        for &id in order {
            if self.anchors[id] {
                anchor_rank[id] = rank;
                rank += 1;
            }
        }
        let mut reward = 0.0;
        for edge in &endpoint.edges {
            let u = edge.u as usize;
            let v = edge.v as usize;
            let (left_rank, right_rank, score) = if anchor_rank[u] < anchor_rank[v] {
                (anchor_rank[u], anchor_rank[v], edge.u_before_v)
            } else {
                (anchor_rank[v], anchor_rank[u], edge.v_before_u)
            };
            let distance = right_rank.saturating_sub(left_rank);
            let scale = match distance {
                1 => 1.0,
                2 => 0.5,
                3..=4 => 0.25,
                5..=8 => 0.125,
                _ => 0.0,
            };
            reward += score * scale;
        }
        (1.0 - reward / endpoint.total_best_score).clamp(0.0, 1.0)
    }

    /// Return normalized primary/secondary/tertiary components for the
    /// length-tiered objective. Each component is in [0, 1], and lower is
    /// better. This is public primarily for diagnostics and benchmark tests.
    pub fn evaluate_length_tiers(&self, order: &[usize]) -> [f64; 3] {
        let mut midpoints = vec![0.0; self.lengths.len()];
        if order.len() != self.lengths.len() {
            return [f64::INFINITY; 3];
        }
        let mut cumulative = 0.0;
        for &id in order {
            if id >= self.lengths.len() {
                return [f64::INFINITY; 3];
            }
            midpoints[id] = cumulative + self.lengths[id] * 0.5;
            cumulative += self.lengths[id];
        }
        self.length_tiered_components(order, &midpoints)
    }

    fn length_tiered_components(&self, order: &[usize], midpoints: &[f64]) -> [f64; 3] {
        let mut anchor_midpoints = vec![0.0; self.lengths.len()];
        let mut cumulative = 0.0;
        for &id in order {
            if self.anchors[id] {
                anchor_midpoints[id] = cumulative + self.lengths[id] * 0.5;
                cumulative += self.lengths[id];
            }
        }

        let mut rewards = [0.0; 3];
        let mut upper_bounds = [0.0; 3];
        for edge in &self.edges {
            let u = edge.u as usize;
            let v = edge.v as usize;
            let (tier, distance) = match (self.anchors[u], self.anchors[v]) {
                (true, true) => (0, (anchor_midpoints[u] - anchor_midpoints[v]).abs()),
                (true, false) | (false, true) => (1, (midpoints[u] - midpoints[v]).abs()),
                (false, false) => (2, (midpoints[u] - midpoints[v]).abs()),
            };
            let minimum_distance = (self.lengths[u] + self.lengths[v]) * 0.5;
            if distance > 0.0 && minimum_distance > 0.0 {
                rewards[tier] += edge.links / distance;
                upper_bounds[tier] += edge.links / minimum_distance;
            }
        }
        let mut costs = [0.0; 3];
        for tier in 0..3 {
            if upper_bounds[tier] > 0.0 {
                costs[tier] = (1.0 - rewards[tier] / upper_bounds[tier]).clamp(0.0, 1.0);
            }
        }
        costs
    }

    /// Build a deterministic, objective-aligned ordering seed by spectral
    /// seriation of the contact graph followed by exact adjacent-swap
    /// polishing. The supplied fallback remains a candidate, so the returned
    /// score can never be worse than the input ordering score.
    pub fn seriation_seed(&self, fallback: &[usize]) -> Result<SeriationResult, String> {
        validate_tour(fallback, self.contig_count())?;
        let n = self.contig_count();
        let fallback_score = self.evaluate(fallback);
        if n < 3 || self.edges.is_empty() {
            return Ok(SeriationResult {
                order: fallback.to_vec(),
                fallback_score,
                spectral_score: fallback_score,
                final_score: fallback_score,
                candidates: 1,
                accepted_swaps: 0,
            });
        }

        // log1p compression limits domination by a few coverage-heavy pairs.
        // Symmetric degree normalization then removes most contig-specific
        // contact-coverage bias.
        let mut degree = vec![0.0; n];
        let mut adjacency = vec![Vec::<(usize, f64)>::new(); n];
        for edge in &self.edges {
            let u = edge.u as usize;
            let v = edge.v as usize;
            let weight = edge.links.ln_1p();
            degree[u] += weight;
            degree[v] += weight;
            adjacency[u].push((v, weight));
            adjacency[v].push((u, weight));
        }

        let mut principal = degree.iter().map(|value| value.sqrt()).collect::<Vec<_>>();
        normalize_vector(&mut principal);
        let mut vector = (0..n)
            .map(|index| {
                let mixed = index.wrapping_mul(2_654_435_761usize) % n;
                mixed as f64 - (n.saturating_sub(1) as f64 * 0.5)
            })
            .collect::<Vec<_>>();
        orthogonalize(&mut vector, &principal);
        normalize_vector(&mut vector);

        let mut best_order = fallback.to_vec();
        let mut best_spectral_order = fallback.to_vec();
        let mut best_score = fallback_score;
        let mut best_spectral_score = fallback_score;
        let mut candidates = 1usize;
        let checkpoints = [2usize, 4, 8, 16, 32, 64, 128, 256];
        let mut checkpoint_index = 0usize;

        for iteration in 1..=checkpoints[checkpoints.len() - 1] {
            // Power iteration on (I + D^-1/2 W D^-1/2) / 2. The lazy
            // transform keeps eigenvalues non-negative; removing the trivial
            // sqrt(degree) vector exposes the Fiedler ordering direction.
            let mut next = vector.iter().map(|value| value * 0.5).collect::<Vec<_>>();
            for u in 0..n {
                if degree[u] <= 0.0 {
                    continue;
                }
                for &(v, weight) in &adjacency[u] {
                    if degree[v] > 0.0 {
                        next[u] += 0.5 * weight * vector[v] / (degree[u] * degree[v]).sqrt();
                    }
                }
            }
            orthogonalize(&mut next, &principal);
            if !normalize_vector(&mut next) {
                break;
            }
            vector = next;

            if checkpoint_index < checkpoints.len() && iteration == checkpoints[checkpoint_index] {
                let mut candidate = (0..n).collect::<Vec<_>>();
                candidate.sort_unstable_by(|&left, &right| {
                    vector[left]
                        .total_cmp(&vector[right])
                        .then_with(|| left.cmp(&right))
                });
                let score = self.evaluate(&candidate);
                candidates += 1;
                if score < best_spectral_score {
                    best_spectral_score = score;
                    best_spectral_order.clone_from(&candidate);
                }
                if score < best_score {
                    best_score = score;
                    best_order = candidate;
                }
                checkpoint_index += 1;
            }
        }

        // Polish the strongest spectral candidate even when the fallback has
        // a better raw score: local corrections can cross that boundary. Keep
        // the fallback as the final safety net.
        let mut polished = best_spectral_order;
        let mut polished_score = best_spectral_score;
        let mut accepted_swaps = 0usize;
        for pass in 0..8 {
            let mut improved = false;
            if pass % 2 == 0 {
                for position in 0..n - 1 {
                    polished.swap(position, position + 1);
                    let score = self.evaluate(&polished);
                    if score < polished_score {
                        polished_score = score;
                        accepted_swaps += 1;
                        improved = true;
                    } else {
                        polished.swap(position, position + 1);
                    }
                }
            } else {
                for position in (0..n - 1).rev() {
                    polished.swap(position, position + 1);
                    let score = self.evaluate(&polished);
                    if score < polished_score {
                        polished_score = score;
                        accepted_swaps += 1;
                        improved = true;
                    } else {
                        polished.swap(position, position + 1);
                    }
                }
            }
            if !improved {
                break;
            }
        }
        if polished_score < best_score {
            best_score = polished_score;
            best_order = polished;
        }

        Ok(SeriationResult {
            order: best_order,
            fallback_score,
            spectral_score: best_spectral_score,
            final_score: best_score,
            candidates,
            accepted_swaps,
        })
    }
}

fn select_length_anchors(lengths: &[f64]) -> Vec<bool> {
    let mut anchors = vec![false; lengths.len()];
    if lengths.is_empty() {
        return anchors;
    }
    let total: f64 = lengths.iter().sum();
    let target = total * ANCHOR_COVERAGE_FRACTION;
    let minimum = lengths.len().min(4);
    let mut ranked = (0..lengths.len()).collect::<Vec<_>>();
    ranked.sort_unstable_by(|&left, &right| {
        lengths[right]
            .total_cmp(&lengths[left])
            .then_with(|| left.cmp(&right))
    });
    let mut covered = 0.0;
    for (rank, id) in ranked.into_iter().enumerate() {
        if rank >= minimum && covered >= target {
            break;
        }
        anchors[id] = true;
        covered += lengths[id];
    }
    anchors
}

fn encode_tiered_score(components: [f64; 3]) -> f64 {
    if components.iter().any(|component| !component.is_finite()) {
        return f64::INFINITY;
    }
    0.7 * components[0] + 0.2 * components[1] + 0.1 * components[2]
}

fn orthogonalize(vector: &mut [f64], basis: &[f64]) {
    let projection = vector
        .iter()
        .zip(basis)
        .map(|(value, base)| value * base)
        .sum::<f64>();
    for (value, base) in vector.iter_mut().zip(basis) {
        *value -= projection * base;
    }
}

fn normalize_vector(vector: &mut [f64]) -> bool {
    let norm = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
    if !norm.is_finite() || norm <= f64::EPSILON {
        return false;
    }
    for value in vector {
        *value /= norm;
    }
    true
}

#[derive(Debug, Clone)]
pub struct SeriationResult {
    pub order: Vec<usize>,
    pub fallback_score: f64,
    pub spectral_score: f64,
    pub final_score: f64,
    pub candidates: usize,
    pub accepted_swaps: usize,
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
    /// Build conservative path blocks directly from the adjacency evidence
    /// produced by the hierarchical end initializer. Only joins that are
    /// adjacent in `seed_order` can become locked edges, so every unit is a
    /// contiguous slice of that seed and its orientation is unambiguous.
    pub fn from_hierarchical_joins(
        seed_order: &[usize],
        joins: &[HierarchicalJoinEvidence],
        config: &BackboneConfig,
    ) -> Result<Self, String> {
        config.validate()?;
        let n = seed_order.len();
        validate_tour(seed_order, n)?;
        if !config.enabled || n < 2 || joins.is_empty() {
            return Ok(Self::seed_singletons(seed_order, config));
        }

        // A hierarchical join connects two cluster endpoints. Once the final
        // seed has been decoded, such a join is useful only when those
        // endpoints are adjacent in the seed. Canonicalizing the pair makes
        // this independent of cluster orientation and also deduplicates
        // repeated evidence deterministically.
        let mut reciprocal_joins = 0usize;
        let mut supported_joins = 0usize;
        let mut confidence_qualified_joins = 0usize;
        let mut trusted_pairs = HashMap::<(usize, usize), ()>::new();
        for join in joins {
            if join.left_contig >= n || join.right_contig >= n {
                return Err(format!(
                    "hierarchical join endpoint ({}, {}) is outside 0..{}",
                    join.left_contig, join.right_contig, n
                ));
            }
            if !join.reciprocal_confident {
                continue;
            }
            reciprocal_joins += 1;
            if !join.raw_support.is_finite() || join.raw_support < config.min_links {
                continue;
            }
            supported_joins += 1;
            if join.confidence.is_nan() || join.confidence < HIERARCHICAL_MIN_CONFIDENCE {
                continue;
            }
            confidence_qualified_joins += 1;
            let pair = if join.left_contig < join.right_contig {
                (join.left_contig, join.right_contig)
            } else {
                (join.right_contig, join.left_contig)
            };
            trusted_pairs.insert(pair, ());
        }

        let mut units = Vec::new();
        let mut current = vec![seed_order[0]];
        let mut accepted_edges = 0usize;
        for adjacent in seed_order.windows(2) {
            let pair = if adjacent[0] < adjacent[1] {
                (adjacent[0], adjacent[1])
            } else {
                (adjacent[1], adjacent[0])
            };
            if current.len() < config.max_block_size && trusted_pairs.contains_key(&pair) {
                current.push(adjacent[1]);
                accepted_edges += 1;
            } else {
                units.push(current);
                current = vec![adjacent[1]];
            }
        }
        units.push(current);

        let plan = Self::from_hierarchical_units(n, units, accepted_edges, config);
        let largest_block = plan.units.iter().map(Vec::len).max().unwrap_or(0);
        log::info!(
            "Hierarchical Backbone evidence: {} joins, {} reciprocal, {} support-qualified, {} confidence-qualified (threshold {:.3}), {} unique trusted pairs, {} accepted seed adjacencies; {} blocks / {} units, largest block {}",
            joins.len(),
            reciprocal_joins,
            supported_joins,
            confidence_qualified_joins,
            HIERARCHICAL_MIN_CONFIDENCE,
            trusted_pairs.len(),
            plan.accepted_edges,
            plan.block_count,
            plan.units.len(),
            largest_block
        );
        Ok(plan)
    }

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

    fn seed_singletons(seed_order: &[usize], config: &BackboneConfig) -> Self {
        let n = seed_order.len();
        Self::from_hierarchical_units(
            n,
            seed_order.iter().map(|&contig| vec![contig]).collect(),
            0,
            config,
        )
    }

    fn from_hierarchical_units(
        n: usize,
        units: Vec<Vec<usize>>,
        accepted_edges: usize,
        config: &BackboneConfig,
    ) -> Self {
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

    /// Split every existing path unit into smaller contiguous units without
    /// joining across a parent boundary. The resulting units can be freely
    /// reordered and reversed by a later unit-GA stage, progressively
    /// releasing constraints while retaining a complete contig partition.
    pub fn refined(&self, max_block_size: usize) -> Result<Self, String> {
        if max_block_size == 0 {
            return Err("refined backbone maximum block size must be at least 1".into());
        }
        if self.units.iter().all(|unit| unit.len() <= max_block_size) {
            return Ok(self.clone());
        }

        let n = self.unit_of_contig.len();
        let units = self
            .units
            .iter()
            .flat_map(|unit| unit.chunks(max_block_size).map(|chunk| chunk.to_vec()))
            .collect::<Vec<_>>();
        let mut unit_of_contig = vec![usize::MAX; n];
        for (unit_index, unit) in units.iter().enumerate() {
            for &contig in unit {
                unit_of_contig[contig] = unit_index;
            }
        }
        let block_count = units.iter().filter(|unit| unit.len() > 1).count();
        let accepted_edges = n.saturating_sub(units.len());
        Ok(Self {
            units,
            unit_of_contig,
            block_count,
            accepted_edges,
            useful: self.useful && block_count > 0,
        })
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
    sum_distance: f64,
    present: bool,
}

impl Default for OrientationSummary {
    fn default() -> Self {
        Self {
            bins: [0; GOLDEN_BINS],
            links: 0,
            sum_log_distance: 0.0,
            sum_distance: 0.0,
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

/// Gap convention used when projecting CLM distances onto a fixed tour.
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub enum OrientationGapModel {
    /// Preserve the historical ALLHiC `EvaluateQ` indexing exactly.
    #[default]
    AllhicLegacy,
    /// Add only complete contigs strictly between the linked contigs.
    Intervening,
}

/// Per-pair support normalization for the banded orientation objective.
#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub enum OrientationPairWeight {
    /// Retain the raw link-weighted CLM log-distance score.
    #[default]
    Links,
    /// Divide each pair score by the square root of its link support.
    SqrtLinks,
    /// Divide each pair score by its link support so every pair has equal mass.
    EqualPair,
}

/// Configuration for exact fixed-order, banded CLM orientation inference.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct BandedOrientationConfig {
    pub rank_window: usize,
    pub gap_model: OrientationGapModel,
    pub pair_weight: OrientationPairWeight,
    pub min_links: usize,
    pub require_complete: bool,
}

impl Default for BandedOrientationConfig {
    fn default() -> Self {
        Self {
            rank_window: 2,
            gap_model: OrientationGapModel::Intervening,
            pair_weight: OrientationPairWeight::Links,
            min_links: 1,
            require_complete: true,
        }
    }
}

/// Audit summary returned by exact banded orientation inference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BandedOrientationResult {
    pub initial_score: f64,
    pub final_score: f64,
    pub changed_signs: usize,
    pub rejected_low_confidence: usize,
    pub rejected_oversized_changes: usize,
    pub used_pairs: usize,
    pub skipped_incomplete_pairs: usize,
}

/// Configuration for post-order signed block refinement. Each accepted move
/// reverse-complements one contiguous block: the order is reversed and every
/// sign in the block is flipped atomically.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SignedBlockRefineConfig {
    /// Largest block, in contigs, considered by one refinement move.
    pub max_span: usize,
    /// Largest block as a fraction of the scaffold length in base pairs.
    pub max_bp_fraction: f64,
    /// Maximum number of accepted best-improvement sweeps.
    pub max_passes: usize,
    /// Required normalized improvement at each of the two block boundaries.
    pub min_relative_gain: f64,
}

impl Default for SignedBlockRefineConfig {
    fn default() -> Self {
        Self {
            max_span: 32,
            max_bp_fraction: 0.05,
            max_passes: 4,
            min_relative_gain: 0.05,
        }
    }
}

/// Audit summary returned by signed block refinement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SignedBlockRefineResult {
    pub initial_score: f64,
    pub final_score: f64,
    pub accepted_moves: usize,
    pub evaluated_moves: usize,
    pub passes: usize,
}

#[derive(Debug, Clone, Copy)]
struct BandedOrientationEdge {
    left_rank: usize,
    potentials: [f64; 4],
}

#[derive(Debug, Clone, Copy)]
struct SignedBoundaryEvidence {
    reward: f64,
    support: usize,
    confidence: f64,
}

#[derive(Debug, Clone, Copy)]
struct SignedBlockMove {
    start: usize,
    end: usize,
    block_bp: u64,
    old_left: SignedBoundaryEvidence,
    new_left: SignedBoundaryEvidence,
    old_right: SignedBoundaryEvidence,
    new_right: SignedBoundaryEvidence,
    left_gain: f64,
    right_gain: f64,
    total_gain: f64,
}

#[derive(Debug, Clone, Copy)]
struct BandedOrientationCell {
    score: f64,
    changed_signs: usize,
}

const MAX_BANDED_ORIENTATION_WINDOW: usize = 16;
const MAX_BANDED_ORIENTATION_DP_CELLS: usize = 32 * 1024 * 1024;

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
            .collect::<Vec<_>>();
        let dense_lengths = lengths
            .iter()
            .map(|&length| length as f64)
            .collect::<Vec<_>>();
        let anchors = select_length_anchors(&dense_lengths);
        let ordering = OptimizeProblem {
            lengths: dense_lengths,
            edges,
            objective: OrderObjective::default(),
            anchors,
            endpoint_multiscale: None,
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
            entry.sum_distance += distance as f64;
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
        self.evaluate_with_gap_model(tour, OrientationGapModel::AllhicLegacy)
    }

    /// Evaluate a tour with an explicit gap convention. The historical
    /// [`Self::evaluate`] entry point deliberately retains ALLHiC's indexing.
    pub fn evaluate_with_gap_model(
        &self,
        tour: &Tour<usize>,
        gap_model: OrientationGapModel,
    ) -> f64 {
        let Some((positions, signs_by_id, starts)) = self.layout(tour) else {
            return f64::NEG_INFINITY;
        };

        self.pairs
            .iter()
            .map(|pair| {
                self.evaluate_pair_with_gap_model(
                    pair,
                    &positions,
                    &signs_by_id,
                    &starts,
                    gap_model,
                )
            })
            .sum()
    }

    /// Exactly maximize a fixed-order CLM orientation objective whose pair
    /// interactions are bounded by `rank_window`.
    pub fn optimize_banded(
        &self,
        tour: &mut Tour<usize>,
        config: BandedOrientationConfig,
    ) -> Result<BandedOrientationResult, String> {
        self.optimize_banded_with_source_prior(tour, config, 0.0)
    }

    /// Exactly maximize the banded CLM objective with a dimensionless prior
    /// favoring the input signs. Returned scores remain raw CLM evidence.
    pub fn optimize_banded_with_source_prior(
        &self,
        tour: &mut Tour<usize>,
        config: BandedOrientationConfig,
        prior_strength: f64,
    ) -> Result<BandedOrientationResult, String> {
        self.optimize_banded_with_source_prior_and_confidence(tour, config, prior_strength, 0.0)
    }

    /// Exactly maximize the banded objective, then retain only sign changes
    /// whose global max-marginal advantage is sufficiently decisive. The
    /// margin is normalized by all retained evidence incident on that contig,
    /// so one threshold is comparable across contact protocols and depths.
    pub fn optimize_banded_with_source_prior_and_confidence(
        &self,
        tour: &mut Tour<usize>,
        config: BandedOrientationConfig,
        prior_strength: f64,
        min_confidence: f64,
    ) -> Result<BandedOrientationResult, String> {
        self.optimize_banded_conservative(tour, config, prior_strength, min_confidence, 1.0)
    }

    /// Apply max-marginal confidence and physical-span safety gates after the
    /// exact banded optimization. Consecutive changed signs whose aggregate
    /// length exceeds `max_flip_bp_fraction` of the scaffold are restored to
    /// their input signs.
    pub fn optimize_banded_conservative(
        &self,
        tour: &mut Tour<usize>,
        config: BandedOrientationConfig,
        prior_strength: f64,
        min_confidence: f64,
        max_flip_bp_fraction: f64,
    ) -> Result<BandedOrientationResult, String> {
        self.optimize_banded_margin(tour, config, prior_strength, min_confidence, max_flip_bp_fraction, false)
    }

    /// Experimental contact-margin variant. Remove the queried node's input
    /// prior from its max-marginal contrast; other nodes retain their priors.
    /// The resulting normalized effect is not a probability of correctness.
    pub fn optimize_banded_contact_evidence(
        &self,
        tour: &mut Tour<usize>,
        config: BandedOrientationConfig,
        prior_strength: f64,
        min_confidence: f64,
        max_flip_bp_fraction: f64,
    ) -> Result<BandedOrientationResult, String> {
        self.optimize_banded_margin(tour, config, prior_strength, min_confidence, max_flip_bp_fraction, true)
    }

    fn optimize_banded_margin(
        &self,
        tour: &mut Tour<usize>,
        config: BandedOrientationConfig,
        prior_strength: f64,
        min_confidence: f64,
        max_flip_bp_fraction: f64,
        contact_margin: bool,
    ) -> Result<BandedOrientationResult, String> {
        if !prior_strength.is_finite() || prior_strength < 0.0 {
            return Err(format!(
                "banded orientation prior_strength must be finite and non-negative, got {prior_strength}"
            ));
        }
        if !min_confidence.is_finite() || !(0.0..=1.0).contains(&min_confidence) {
            return Err(format!(
                "banded orientation min_confidence must be between 0 and 1, got {min_confidence}"
            ));
        }
        if !max_flip_bp_fraction.is_finite()
            || max_flip_bp_fraction <= 0.0
            || max_flip_bp_fraction > 1.0
        {
            return Err(format!(
                "banded orientation max_flip_bp_fraction must be in (0, 1], got {max_flip_bp_fraction}"
            ));
        }
        if config.rank_window == 0 || config.rank_window > MAX_BANDED_ORIENTATION_WINDOW {
            return Err(format!(
                "banded orientation rank_window must be in 1..={MAX_BANDED_ORIENTATION_WINDOW}, got {}",
                config.rank_window
            ));
        }
        let n = self.lengths.len();
        if tour.contigs.len() != n || tour.signs.len() != n {
            return Err(format!(
                "banded orientation requires a complete tour of {n} contigs, got {} contigs and {} signs",
                tour.contigs.len(),
                tour.signs.len()
            ));
        }
        // Intervening gaps are reverse-complement invariant. Canonicalizing the
        // axis also makes floating-point accumulation and tertiary DP ties
        // identical for a tour and its reverse complement.
        let reverse_axis = config.gap_model == OrientationGapModel::Intervening
            && tour.contigs.iter().rev().cmp(tour.contigs.iter()) == std::cmp::Ordering::Less;
        let mut working_tour = tour.clone();
        if reverse_axis {
            working_tour.contigs.reverse();
            working_tour.signs.reverse();
            for sign in &mut working_tour.signs {
                *sign = !*sign;
            }
        }
        let Some((positions, _, starts)) = self.layout(&working_tour) else {
            return Err(
                "banded orientation tour contains an out-of-range or duplicate contig".into(),
            );
        };
        if positions.iter().any(|&position| position == usize::MAX) {
            return Err("banded orientation tour does not contain every problem contig".into());
        }

        let mut edges_by_right = vec![Vec::<BandedOrientationEdge>::new(); n];
        let mut node_scale = vec![0.0; n];
        let mut used_pairs = 0usize;
        let mut skipped_incomplete_pairs = 0usize;
        let mut bandwidth = 0usize;
        for pair in &self.pairs {
            let pu = positions[pair.u];
            let pv = positions[pair.v];
            let (left_rank, right_rank) = if pu < pv { (pu, pv) } else { (pv, pu) };
            let span = right_rank - left_rank;
            if span > config.rank_window {
                continue;
            }

            let first_links = pair.orientations[0].links;
            let complete = pair
                .orientations
                .iter()
                .all(|summary| summary.present && summary.links == first_links);
            let support = if complete {
                first_links
            } else {
                pair.orientations
                    .iter()
                    .map(|summary| summary.links)
                    .max()
                    .unwrap_or(0)
            };
            if config.require_complete && (!complete || support < config.min_links) {
                skipped_incomplete_pairs += 1;
                continue;
            }
            if support == 0 || support < config.min_links {
                continue;
            }

            let gap = orientation_gap(&starts, left_rank, right_rank, config.gap_model);
            if gap > MAX_ORIENTATION_DISTANCE {
                continue;
            }
            let divisor = match config.pair_weight {
                OrientationPairWeight::Links => 1.0,
                OrientationPairWeight::SqrtLinks => (support as f64).sqrt(),
                OrientationPairWeight::EqualPair => support as f64,
            };
            let pair_scores: [f64; 4] = std::array::from_fn(|orientation| {
                self.evaluate_summary(&pair.orientations[orientation], gap) / divisor
            });
            let mut potentials = [0.0; 4];
            for left_bit in 0..2 {
                for right_bit in 0..2 {
                    let left_forward = left_bit != 0;
                    let right_forward = right_bit != 0;
                    let orientation = if pu < pv {
                        orientation_index(left_forward, right_forward)
                    } else {
                        orientation_index(!right_forward, !left_forward)
                    };
                    potentials[(left_bit << 1) | right_bit] = pair_scores[orientation];
                }
            }
            let pair_max = potentials.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            for potential in &mut potentials {
                *potential -= pair_max;
            }
            let edge_min = potentials.iter().copied().fold(f64::INFINITY, f64::min);
            let edge_max = potentials.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let edge_range = edge_max - edge_min;
            node_scale[left_rank] += edge_range;
            node_scale[right_rank] += edge_range;
            edges_by_right[right_rank].push(BandedOrientationEdge {
                left_rank,
                potentials,
            });
            bandwidth = bandwidth.max(span);
            used_pairs += 1;
        }
        for edges in &mut edges_by_right {
            edges.sort_unstable_by_key(|edge| edge.left_rank);
        }
        if used_pairs == 0 {
            return Ok(BandedOrientationResult {
                initial_score: 0.0,
                final_score: 0.0,
                changed_signs: 0,
                rejected_low_confidence: 0,
                rejected_oversized_changes: 0,
                used_pairs,
                skipped_incomplete_pairs,
            });
        }

        let input_signs = working_tour.signs.clone();
        let score_signs = |signs: &[bool]| -> f64 {
            edges_by_right
                .iter()
                .enumerate()
                .map(|(right_rank, edges)| {
                    edges
                        .iter()
                        .map(|edge| {
                            let left_bit = usize::from(signs[edge.left_rank]);
                            let right_bit = usize::from(signs[right_rank]);
                            edge.potentials[(left_bit << 1) | right_bit]
                        })
                        .sum::<f64>()
                })
                .sum()
        };
        let initial_score = score_signs(&input_signs);

        let state_count = 1usize << bandwidth;
        let dp_cells = n.checked_mul(state_count).ok_or_else(|| {
            "banded orientation DP dimensions overflow addressable memory".to_string()
        })?;
        if dp_cells > MAX_BANDED_ORIENTATION_DP_CELLS {
            return Err(format!(
                "banded orientation DP needs {dp_cells} traceback cells for {n} contigs and bandwidth {bandwidth}; reduce rank_window (limit {MAX_BANDED_ORIENTATION_DP_CELLS})"
            ));
        }

        let unreachable = BandedOrientationCell {
            score: f64::NEG_INFINITY,
            changed_signs: usize::MAX,
        };
        let mut current = vec![unreachable; state_count];
        let mut next = vec![unreachable; state_count];
        current[0] = BandedOrientationCell {
            score: 0.0,
            changed_signs: 0,
        };
        let transition_score = |right_rank: usize, state: usize, right_bit: usize| -> f64 {
            let added_score = edges_by_right[right_rank]
                .iter()
                .map(|edge| {
                    let distance = right_rank - edge.left_rank;
                    let left_bit = (state >> (distance - 1)) & 1;
                    edge.potentials[(left_bit << 1) | right_bit]
                })
                .sum::<f64>();
            let sign_changed = (right_bit != 0) != input_signs[right_rank];
            if prior_strength == 0.0 || !sign_changed {
                added_score
            } else {
                added_score - prior_strength * node_scale[right_rank]
            }
        };
        let mut parents = vec![vec![u32::MAX; state_count]; n];
        let mut forward_scores = (min_confidence > 0.0).then(|| {
            let mut scores = Vec::with_capacity(n + 1);
            scores.push(current.iter().map(|cell| cell.score).collect::<Vec<_>>());
            scores
        });
        let state_mask = state_count - 1;
        for right_rank in 0..n {
            next.fill(unreachable);
            for (state, cell) in current.iter().copied().enumerate() {
                if cell.changed_signs == usize::MAX {
                    continue;
                }
                for right_bit in 0..2 {
                    let next_state = ((state << 1) | right_bit) & state_mask;
                    let sign_changed = (right_bit != 0) != input_signs[right_rank];
                    let candidate = BandedOrientationCell {
                        score: cell.score + transition_score(right_rank, state, right_bit),
                        changed_signs: cell.changed_signs + usize::from(sign_changed),
                    };
                    let incumbent_parent = parents[right_rank][next_state] as usize;
                    if banded_orientation_candidate_is_better(
                        candidate,
                        state,
                        next[next_state],
                        incumbent_parent,
                    ) {
                        next[next_state] = candidate;
                        parents[right_rank][next_state] = state as u32;
                    }
                }
            }
            std::mem::swap(&mut current, &mut next);
            if let Some(scores) = &mut forward_scores {
                scores.push(current.iter().map(|cell| cell.score).collect());
            }
        }

        let mut best_state = 0usize;
        let mut best = unreachable;
        for (state, cell) in current.iter().copied().enumerate() {
            if banded_orientation_candidate_is_better(cell, state, best, best_state) {
                best = cell;
                best_state = state;
            }
        }
        let mut optimized_signs = vec![false; n];
        let mut state = best_state;
        for right_rank in (0..n).rev() {
            optimized_signs[right_rank] = state & 1 != 0;
            let parent = parents[right_rank][state];
            debug_assert_ne!(parent, u32::MAX);
            state = parent as usize;
        }
        let unconstrained_changed_signs = optimized_signs
            .iter()
            .zip(&input_signs)
            .filter(|(optimized, input)| optimized != input)
            .count();
        debug_assert_eq!(unconstrained_changed_signs, best.changed_signs);

        let mut rejected_low_confidence = 0usize;
        if let Some(forward_scores) = forward_scores {
            let mut backward_next = vec![0.0_f64; state_count];
            let mut backward_current = vec![f64::NEG_INFINITY; state_count];
            for right_rank in (0..n).rev() {
                let mut marginal = [f64::NEG_INFINITY; 2];
                backward_current.fill(f64::NEG_INFINITY);
                for state in 0..state_count {
                    let forward = forward_scores[right_rank][state];
                    if !forward.is_finite() {
                        continue;
                    }
                    for right_bit in 0..2 {
                        let next_state = ((state << 1) | right_bit) & state_mask;
                        let suffix = backward_next[next_state];
                        if !suffix.is_finite() {
                            continue;
                        }
                        let transition = transition_score(right_rank, state, right_bit);
                        marginal[right_bit] =
                            marginal[right_bit].max(forward + transition + suffix);
                        backward_current[state] = backward_current[state].max(transition + suffix);
                    }
                }

                if optimized_signs[right_rank] != input_signs[right_rank] {
                    let optimized_bit = usize::from(optimized_signs[right_rank]);
                    let input_bit = usize::from(input_signs[right_rank]);
                    let scale = node_scale[right_rank];
                    let confidence = if scale > f64::EPSILON {
                        // Remove this node's own input-sign penalty before
                        // measuring evidence. Otherwise the maximum margin is
                        // 1 - prior_strength, making a 0.95 gate with a 0.05
                        // prior impossible to pass in exact arithmetic.
                        let raw = (marginal[optimized_bit] - marginal[input_bit]) / scale;
                        if contact_margin { (raw + prior_strength).clamp(0.0, 1.0) }
                        else { raw.clamp(0.0, 1.0) }
                    } else {
                        0.0
                    };
                    if confidence <= min_confidence {
                        optimized_signs[right_rank] = input_signs[right_rank];
                        rejected_low_confidence += 1;
                    }
                }
                std::mem::swap(&mut backward_current, &mut backward_next);
            }
        }

        let total_bp = working_tour
            .contigs
            .iter()
            .fold(0u64, |total, &id| total.saturating_add(self.lengths[id]));
        let max_changed_bp = ((total_bp as f64) * max_flip_bp_fraction).ceil() as u64;
        let mut rejected_oversized_changes = 0usize;
        let mut run_start = 0usize;
        while run_start < n {
            if optimized_signs[run_start] == input_signs[run_start] {
                run_start += 1;
                continue;
            }
            let mut run_end = run_start + 1;
            let mut run_bp = self.lengths[working_tour.contigs[run_start]];
            while run_end < n && optimized_signs[run_end] != input_signs[run_end] {
                run_bp = run_bp.saturating_add(self.lengths[working_tour.contigs[run_end]]);
                run_end += 1;
            }
            if run_bp > max_changed_bp {
                optimized_signs[run_start..run_end]
                    .copy_from_slice(&input_signs[run_start..run_end]);
                rejected_oversized_changes += run_end - run_start;
            }
            run_start = run_end;
        }
        let final_score = score_signs(&optimized_signs);
        let changed_signs = optimized_signs
            .iter()
            .zip(&input_signs)
            .filter(|(optimized, input)| optimized != input)
            .count();
        if reverse_axis {
            optimized_signs.reverse();
            for sign in &mut optimized_signs {
                *sign = !*sign;
            }
        }
        tour.signs = optimized_signs;
        Ok(BandedOrientationResult {
            initial_score,
            final_score,
            changed_signs,
            rejected_low_confidence,
            rejected_oversized_changes,
            used_pairs,
            skipped_incomplete_pairs,
        })
    }

    /// Historical local-band block objective used by `optimize --orientation-method
    /// banded-legacy`. Reverse-complement blocks are scored using all affected
    /// banded pairs; each accepted move is followed by ungated banded orientation
    /// with a prior on the current signs. Terminal blocks are included.
    ///
    /// `max_bp_fraction` must be 1: the historical method has no physical-span
    /// limit. Use [`Self::refine_signed_blocks_conservative`] for that policy.
    pub fn refine_signed_blocks_legacy_objective(
        &self,
        tour: &mut Tour<usize>,
        orientation_config: BandedOrientationConfig,
        prior_strength: f64,
        config: SignedBlockRefineConfig,
    ) -> Result<SignedBlockRefineResult, String> {
        if orientation_config.gap_model != OrientationGapModel::Intervening {
            return Err(
                "signed block refinement requires the reverse-complement-invariant intervening gap model"
                    .into(),
            );
        }
        if orientation_config.rank_window == 0
            || orientation_config.rank_window > MAX_BANDED_ORIENTATION_WINDOW
        {
            return Err(format!(
                "signed block refinement rank_window must be in 1..={MAX_BANDED_ORIENTATION_WINDOW}, got {}",
                orientation_config.rank_window
            ));
        }
        if config.max_span < 2 {
            return Err("signed block refinement max_span must be at least 2".into());
        }
        if config.max_bp_fraction != 1.0 {
            return Err(
                "historical signed block refinement requires max_bp_fraction = 1; use conservative refinement for a physical-span limit".into(),
            );
        }
        if config.max_passes == 0 {
            return Err("signed block refinement max_passes must be greater than zero".into());
        }
        if !config.min_relative_gain.is_finite() || config.min_relative_gain < 0.0 {
            return Err(format!(
                "signed block refinement min_relative_gain must be finite and non-negative, got {}",
                config.min_relative_gain
            ));
        }
        if !prior_strength.is_finite() || prior_strength < 0.0 {
            return Err(format!(
                "signed block refinement prior_strength must be finite and non-negative, got {prior_strength}"
            ));
        }

        let Some((mut positions, mut signs_by_id, mut starts)) = self.layout(tour) else {
            return Err(
                "signed block refinement requires a complete tour without duplicate contigs".into(),
            );
        };
        if positions.iter().any(|&position| position == usize::MAX) {
            return Err(
                "signed block refinement tour does not contain every problem contig".into(),
            );
        }

        let mut pair_rewards = self
            .pairs
            .iter()
            .map(|pair| {
                self.banded_signed_pair_reward(
                    pair,
                    &positions,
                    &signs_by_id,
                    &starts,
                    orientation_config,
                )
            })
            .collect::<Vec<_>>();
        let mut current_score = pair_rewards.iter().sum::<f64>();
        let initial_score = current_score;
        let mut accepted_moves = 0usize;
        let mut evaluated_moves = 0usize;
        let mut passes = 0usize;
        let mut pair_marks = vec![0u32; self.pairs.len()];
        let mut mark = 0u32;
        let mut affected_pairs = Vec::new();

        let max_span = config.max_span.min(tour.contigs.len());
        for _ in 0..config.max_passes {
            passes += 1;
            let mut best_move = None;
            let mut best_score = current_score;

            for span in 2..=max_span {
                for start in 0..=tour.contigs.len() - span {
                    let end = start + span - 1;
                    mark = mark.wrapping_add(1);
                    if mark == 0 {
                        pair_marks.fill(0);
                        mark = 1;
                    }
                    affected_pairs.clear();
                    for &id in &tour.contigs[start..=end] {
                        for &pair_index in &self.incident_pairs[id] {
                            if pair_marks[pair_index] == mark {
                                continue;
                            }
                            pair_marks[pair_index] = mark;
                            let pair = &self.pairs[pair_index];
                            let u_inside = (start..=end).contains(&positions[pair.u]);
                            let v_inside = (start..=end).contains(&positions[pair.v]);
                            if u_inside != v_inside {
                                affected_pairs.push(pair_index);
                            }
                        }
                    }

                    reverse_complement_range(tour, start, end);
                    self.refresh_layout_range(
                        tour,
                        &mut positions,
                        &mut signs_by_id,
                        &mut starts,
                        start,
                        end,
                    );
                    let old_local = affected_pairs
                        .iter()
                        .map(|&pair_index| pair_rewards[pair_index])
                        .sum::<f64>();
                    let new_local = affected_pairs
                        .iter()
                        .map(|&pair_index| {
                            self.banded_signed_pair_reward(
                                &self.pairs[pair_index],
                                &positions,
                                &signs_by_id,
                                &starts,
                                orientation_config,
                            )
                        })
                        .sum::<f64>();
                    let candidate_score = current_score - old_local + new_local;
                    evaluated_moves += 1;
                    reverse_complement_range(tour, start, end);
                    self.refresh_layout_range(
                        tour,
                        &mut positions,
                        &mut signs_by_id,
                        &mut starts,
                        start,
                        end,
                    );

                    let local_scale = old_local.abs().max(new_local.abs()).max(1.0);
                    let sufficient_gain =
                        candidate_score - current_score > config.min_relative_gain * local_scale;
                    let better_than_best = candidate_score > best_score
                        || (candidate_score == best_score
                            && best_move.is_some_and(|(best_start, best_end)| {
                                span < best_end - best_start + 1
                                    || (span == best_end - best_start + 1 && start < best_start)
                            }));
                    if sufficient_gain && better_than_best {
                        best_score = candidate_score;
                        best_move = Some((start, end));
                    }
                }
            }

            let Some((start, end)) = best_move else {
                break;
            };
            reverse_complement_range(tour, start, end);
            self.optimize_banded_with_source_prior(tour, orientation_config, prior_strength)?;
            let Some(layout) = self.layout(tour) else {
                return Err("signed block refinement produced an invalid optimized tour".into());
            };
            (positions, signs_by_id, starts) = layout;
            for (pair_index, pair) in self.pairs.iter().enumerate() {
                pair_rewards[pair_index] = self.banded_signed_pair_reward(
                    pair,
                    &positions,
                    &signs_by_id,
                    &starts,
                    orientation_config,
                );
            }
            current_score = pair_rewards.iter().sum();
            accepted_moves += 1;
        }

        Ok(SignedBlockRefineResult {
            initial_score,
            final_score: current_score,
            accepted_moves,
            evaluated_moves,
            passes,
        })
    }

    fn banded_signed_pair_reward(
        &self,
        pair: &PairOrientationData,
        positions: &[usize],
        signs_by_id: &[bool],
        starts: &[u64],
        config: BandedOrientationConfig,
    ) -> f64 {
        let pu = positions[pair.u];
        let pv = positions[pair.v];
        if pu == usize::MAX || pv == usize::MAX {
            return 0.0;
        }
        let span = pu.abs_diff(pv);
        if span == 0 || span > config.rank_window {
            return 0.0;
        }

        let first_links = pair.orientations[0].links;
        let complete = pair
            .orientations
            .iter()
            .all(|summary| summary.present && summary.links == first_links);
        let support = if complete {
            first_links
        } else {
            pair.orientations
                .iter()
                .map(|summary| summary.links)
                .max()
                .unwrap_or(0)
        };
        if (config.require_complete && !complete) || support == 0 || support < config.min_links {
            return 0.0;
        }

        let (left, right, orientation) = if pu < pv {
            (
                pu,
                pv,
                orientation_index(signs_by_id[pair.u], signs_by_id[pair.v]),
            )
        } else {
            (
                pv,
                pu,
                orientation_index(!signs_by_id[pair.u], !signs_by_id[pair.v]),
            )
        };
        let gap = orientation_gap(starts, left, right, config.gap_model);
        if gap > MAX_ORIENTATION_DISTANCE {
            return 0.0;
        }
        let divisor = match config.pair_weight {
            OrientationPairWeight::Links => 1.0,
            OrientationPairWeight::SqrtLinks => (support as f64).sqrt(),
            OrientationPairWeight::EqualPair => support as f64,
        };
        let pair_scores = pair
            .orientations
            .each_ref()
            .map(|summary| self.evaluate_summary(summary, gap) / divisor);
        let pair_min = pair_scores.iter().copied().fold(f64::INFINITY, f64::min);
        (pair_scores[orientation] - pair_min).max(0.0)
    }

    /// Conservative signed-block refinement shared by all contact protocols.
    /// A block is accepted only when both of its replacement adjacencies improve
    /// independently, its physical span is bounded, and the gain overcomes a
    /// sign prior anchored to the immutable source tour.
    pub fn refine_signed_blocks_conservative(
        &self,
        tour: &mut Tour<usize>,
        source_tour: &Tour<usize>,
        orientation_config: BandedOrientationConfig,
        prior_strength: f64,
        config: SignedBlockRefineConfig,
    ) -> Result<SignedBlockRefineResult, String> {
        if orientation_config.gap_model != OrientationGapModel::Intervening {
            return Err(
                "conservative signed block refinement requires the reverse-complement-invariant intervening gap model"
                    .into(),
            );
        }
        if config.max_span < 2 {
            return Err("signed block refinement max_span must be at least 2".into());
        }
        if !config.max_bp_fraction.is_finite()
            || config.max_bp_fraction <= 0.0
            || config.max_bp_fraction > 1.0
        {
            return Err(format!(
                "signed block refinement max_bp_fraction must be in (0, 1], got {}",
                config.max_bp_fraction
            ));
        }
        if config.max_passes == 0 {
            return Err("signed block refinement max_passes must be greater than zero".into());
        }
        if !config.min_relative_gain.is_finite() || config.min_relative_gain < 0.0 {
            return Err(format!(
                "signed block refinement min_relative_gain must be finite and non-negative, got {}",
                config.min_relative_gain
            ));
        }
        if !prior_strength.is_finite() || prior_strength < 0.0 {
            return Err(format!(
                "signed block refinement prior_strength must be finite and non-negative, got {prior_strength}"
            ));
        }

        let Some((source_positions, source_signs_by_id, _)) = self.layout(source_tour) else {
            return Err(
                "signed block refinement source tour must be complete without duplicate contigs"
                    .into(),
            );
        };
        if source_positions
            .iter()
            .any(|&position| position == usize::MAX)
        {
            return Err("signed block refinement source tour is incomplete".into());
        }
        let Some((positions, _, _)) = self.layout(tour) else {
            return Err(
                "signed block refinement requires a complete tour without duplicate contigs".into(),
            );
        };
        if positions.iter().any(|&position| position == usize::MAX) {
            return Err(
                "signed block refinement tour does not contain every problem contig".into(),
            );
        }

        let mismatch_count = |candidate: &Tour<usize>| -> usize {
            candidate
                .contigs
                .iter()
                .zip(&candidate.signs)
                .filter(|(id, sign)| source_signs_by_id[**id] != **sign)
                .count()
        };
        let objective = |candidate: &Tour<usize>| -> f64 {
            self.signed_adjacency_score(candidate, orientation_config)
                - prior_strength * mismatch_count(candidate) as f64
        };
        let total_bp = tour.contigs.iter().map(|&id| self.lengths[id]).sum::<u64>();
        let max_block_bp = ((total_bp as f64) * config.max_bp_fraction).ceil() as u64;
        let mut current_score = objective(tour);
        let initial_score = current_score;
        let mut accepted_moves = 0usize;
        let mut evaluated_moves = 0usize;
        let mut passes = 0usize;

        // Terminal blocks expose only one changed boundary, so they cannot pass
        // the same two-sided evidence rule as internal blocks.
        let max_span = config.max_span.min(tour.contigs.len().saturating_sub(2));
        for _ in 0..config.max_passes {
            passes += 1;
            let mut best_move: Option<SignedBlockMove> = None;
            let mut best_score = current_score;

            for span in 2..=max_span {
                for start in 1..tour.contigs.len() - span {
                    let end = start + span - 1;
                    let block_bp = tour.contigs[start..=end]
                        .iter()
                        .map(|&id| self.lengths[id])
                        .sum::<u64>();
                    if block_bp > max_block_bp {
                        continue;
                    }

                    let old_left = self
                        .signed_boundary_evidence(
                            tour.contigs[start - 1],
                            tour.signs[start - 1],
                            tour.contigs[start],
                            tour.signs[start],
                            orientation_config,
                        )
                        .unwrap_or(SignedBoundaryEvidence {
                            reward: 0.0,
                            support: 0,
                            confidence: 0.0,
                        });
                    let old_right = self
                        .signed_boundary_evidence(
                            tour.contigs[end],
                            tour.signs[end],
                            tour.contigs[end + 1],
                            tour.signs[end + 1],
                            orientation_config,
                        )
                        .unwrap_or(SignedBoundaryEvidence {
                            reward: 0.0,
                            support: 0,
                            confidence: 0.0,
                        });
                    let old_mismatches = tour.contigs[start..=end]
                        .iter()
                        .zip(&tour.signs[start..=end])
                        .filter(|(id, sign)| source_signs_by_id[**id] != **sign)
                        .count();

                    reverse_complement_range(tour, start, end);
                    let new_left = self.signed_boundary_evidence(
                        tour.contigs[start - 1],
                        tour.signs[start - 1],
                        tour.contigs[start],
                        tour.signs[start],
                        orientation_config,
                    );
                    let new_right = self.signed_boundary_evidence(
                        tour.contigs[end],
                        tour.signs[end],
                        tour.contigs[end + 1],
                        tour.signs[end + 1],
                        orientation_config,
                    );
                    let new_mismatches = tour.contigs[start..=end]
                        .iter()
                        .zip(&tour.signs[start..=end])
                        .filter(|(id, sign)| source_signs_by_id[**id] != **sign)
                        .count();
                    let candidate = match (new_left, new_right) {
                        (Some(new_left), Some(new_right)) => {
                            let support_not_weaker = new_left.support >= old_left.support
                                && new_right.support >= old_right.support;
                            let confidence_improves = new_left.confidence - old_left.confidence
                                > config.min_relative_gain
                                && new_right.confidence - old_right.confidence
                                    > config.min_relative_gain;
                            if !support_not_weaker || !confidence_improves {
                                reverse_complement_range(tour, start, end);
                                evaluated_moves += 1;
                                continue;
                            }
                            let left_gain = new_left.reward - old_left.reward;
                            let right_gain = new_right.reward - old_right.reward;
                            let prior_delta =
                                prior_strength * (new_mismatches as f64 - old_mismatches as f64);
                            let total_gain = left_gain + right_gain - prior_delta;
                            Some((left_gain, right_gain, total_gain))
                        }
                        _ => None,
                    };
                    reverse_complement_range(tour, start, end);
                    evaluated_moves += 1;

                    let Some((left_gain, right_gain, total_gain)) = candidate else {
                        continue;
                    };
                    let candidate_score = current_score + total_gain;
                    let both_boundaries_improve = left_gain > config.min_relative_gain
                        && right_gain > config.min_relative_gain;
                    let sufficient_total_gain = total_gain > config.min_relative_gain;
                    let better_than_best = candidate_score > best_score
                        || (candidate_score == best_score
                            && best_move.is_some_and(|best| {
                                span < best.end - best.start + 1
                                    || (span == best.end - best.start + 1 && start < best.start)
                            }));
                    if both_boundaries_improve && sufficient_total_gain && better_than_best {
                        best_score = candidate_score;
                        best_move = Some(SignedBlockMove {
                            start,
                            end,
                            block_bp,
                            old_left,
                            new_left: new_left.unwrap(),
                            old_right,
                            new_right: new_right.unwrap(),
                            left_gain,
                            right_gain,
                            total_gain,
                        });
                    }
                }
            }

            let Some(best) = best_move else {
                break;
            };
            reverse_complement_range(tour, best.start, best.end);
            current_score = objective(tour);
            debug_assert!((current_score - best_score).abs() <= 1e-9);
            accepted_moves += 1;
            log::info!(
                "Accepted signed block: start={}, end={}, contigs={}, bp={}, left_support={}->{}, right_support={}->{}, left_confidence={:.6}->{:.6}, right_confidence={:.6}->{:.6}, left_gain={:.6}, right_gain={:.6}, total_gain={:.6}",
                best.start + 1,
                best.end + 1,
                best.end - best.start + 1,
                best.block_bp,
                best.old_left.support,
                best.new_left.support,
                best.old_right.support,
                best.new_right.support,
                best.old_left.confidence,
                best.new_left.confidence,
                best.old_right.confidence,
                best.new_right.confidence,
                best.left_gain,
                best.right_gain,
                best.total_gain,
            );
        }

        Ok(SignedBlockRefineResult {
            initial_score,
            final_score: current_score,
            accepted_moves,
            evaluated_moves,
            passes,
        })
    }

    fn signed_adjacency_score(&self, tour: &Tour<usize>, config: BandedOrientationConfig) -> f64 {
        tour.contigs
            .windows(2)
            .zip(tour.signs.windows(2))
            .filter_map(|(ids, signs)| {
                self.signed_boundary_evidence(ids[0], signs[0], ids[1], signs[1], config)
            })
            .map(|evidence| evidence.reward)
            .sum()
    }

    /// Return a bounded, support-saturated orientation compatibility score for
    /// one proposed signed adjacency. Normalizing each pair to [0, 1] avoids
    /// comparing arbitrary pair-specific likelihood offsets when a block move
    /// replaces one boundary pair with another.
    fn signed_boundary_evidence(
        &self,
        left_id: usize,
        left_sign: bool,
        right_id: usize,
        right_sign: bool,
        config: BandedOrientationConfig,
    ) -> Option<SignedBoundaryEvidence> {
        let (u, v) = if left_id < right_id {
            (left_id, right_id)
        } else {
            (right_id, left_id)
        };
        let pair_index = self
            .pairs
            .binary_search_by_key(&(u, v), |pair| (pair.u, pair.v))
            .ok()?;
        let pair = &self.pairs[pair_index];
        let first_links = pair.orientations[0].links;
        let complete = pair
            .orientations
            .iter()
            .all(|summary| summary.present && summary.links == first_links);
        let support = if complete {
            first_links
        } else {
            pair.orientations
                .iter()
                .map(|summary| summary.links)
                .max()
                .unwrap_or(0)
        };
        if (config.require_complete && !complete) || support < config.min_links {
            return None;
        }

        let divisor = match config.pair_weight {
            OrientationPairWeight::Links => 1.0,
            OrientationPairWeight::SqrtLinks => (support as f64).sqrt(),
            OrientationPairWeight::EqualPair => support as f64,
        };
        let scores = pair
            .orientations
            .each_ref()
            .map(|summary| self.evaluate_summary(summary, 0) / divisor);
        let minimum = scores.iter().copied().fold(f64::INFINITY, f64::min);
        let maximum = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let range = maximum - minimum;
        if !range.is_finite() || range <= f64::EPSILON {
            return None;
        }
        let orientation = if pair.u == left_id {
            orientation_index(left_sign, right_sign)
        } else {
            orientation_index(!right_sign, !left_sign)
        };
        let runner_up = scores
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != orientation)
            .map(|(_, score)| *score)
            .fold(f64::NEG_INFINITY, f64::max);
        let confidence = ((scores[orientation] - runner_up) / range).clamp(0.0, 1.0);
        Some(SignedBoundaryEvidence {
            reward: confidence * (support as f64).ln_1p(),
            support,
            confidence,
        })
    }

    fn refresh_layout_range(
        &self,
        tour: &Tour<usize>,
        positions: &mut [usize],
        signs_by_id: &mut [bool],
        starts: &mut [u64],
        start: usize,
        end: usize,
    ) {
        let mut cumulative = starts[start];
        for position in start..=end {
            let id = tour.contigs[position];
            positions[id] = position;
            signs_by_id[id] = tour.signs[position];
            starts[position] = cumulative;
            cumulative = cumulative.saturating_add(self.lengths[id]);
        }
        debug_assert!(end + 1 == starts.len() || cumulative == starts[end + 1]);
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
        self.evaluate_pair_with_gap_model(
            pair,
            positions,
            signs_by_id,
            starts,
            OrientationGapModel::AllhicLegacy,
        )
    }

    fn evaluate_pair_with_gap_model(
        &self,
        pair: &PairOrientationData,
        positions: &[usize],
        signs_by_id: &[bool],
        starts: &[u64],
        gap_model: OrientationGapModel,
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
        let gap = orientation_gap(starts, left, right, gap_model);
        if gap > MAX_ORIENTATION_DISTANCE {
            return 0.0;
        }
        self.evaluate_summary(summary, gap)
    }

    fn evaluate_summary(&self, summary: &OrientationSummary, gap: u64) -> f64 {
        if !summary.present {
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

    /// Spectral orientation initialization evaluated with an explicit gap
    /// convention. This avoids accepting a legacy-asymmetric initialization
    /// before intervening-gap refinement.
    pub fn initialize_spectral_with_gap_model(
        &self,
        tour: &mut Tour<usize>,
        gap_model: OrientationGapModel,
    ) -> bool {
        let old_signs = tour.signs.clone();
        let old_score = self.evaluate_with_gap_model(tour, gap_model);
        let spectral = self.spectral_signs();
        for (position, &id) in tour.contigs.iter().enumerate() {
            tour.signs[position] = spectral[id];
        }
        let new_score = self.evaluate_with_gap_model(tour, gap_model);
        if new_score < old_score {
            tour.signs = old_signs;
            false
        } else {
            true
        }
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

    /// Refine signs with an explicit gap convention. The intervening-gap
    /// path operates on a canonical tour axis so scan order and ties are
    /// strictly reverse-complement covariant.
    pub fn refine_with_gap_model(
        &self,
        tour: &mut Tour<usize>,
        gap_model: OrientationGapModel,
    ) -> OrientationResult {
        if gap_model == OrientationGapModel::AllhicLegacy {
            return self.refine(tour);
        }

        let reverse_axis = tour.contigs.len() == tour.signs.len()
            && tour.contigs.iter().rev().cmp(tour.contigs.iter()) == std::cmp::Ordering::Less;
        let mut working_tour = tour.clone();
        if reverse_axis {
            working_tour.contigs.reverse();
            working_tour.signs.reverse();
            for sign in &mut working_tour.signs {
                *sign = !*sign;
            }
        }

        let initial_score = self.evaluate_with_gap_model(&working_tour, gap_model);
        let mut phases = 0;
        loop {
            phases += 1;
            let whole_accepted = self.flip_whole_with_gap_model(&mut working_tour, gap_model);
            let one_accepted = self.flip_one_with_gap_model(&mut working_tour, gap_model);
            if !whole_accepted && !one_accepted {
                break;
            }
        }
        let result = OrientationResult {
            initial_score,
            final_score: self.evaluate_with_gap_model(&working_tour, gap_model),
            phases,
        };

        if reverse_axis {
            working_tour.signs.reverse();
            for sign in &mut working_tour.signs {
                *sign = !*sign;
            }
        }
        tour.signs = working_tour.signs;
        result
    }

    fn flip_all(&self, tour: &mut Tour<usize>) -> bool {
        self.initialize_spectral_with_gap_model(tour, OrientationGapModel::AllhicLegacy)
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
            if orientation_local_improves(old_contribution, new_contribution) {
                accepted = true;
            } else {
                tour.signs[position] = !tour.signs[position];
                signs_by_id[id] = !signs_by_id[id];
            }
        }
        accepted
    }

    fn flip_whole_with_gap_model(
        &self,
        tour: &mut Tour<usize>,
        gap_model: OrientationGapModel,
    ) -> bool {
        let old_score = self.evaluate_with_gap_model(tour, gap_model);
        for sign in &mut tour.signs {
            *sign = !*sign;
        }
        if self.evaluate_with_gap_model(tour, gap_model) <= old_score {
            for sign in &mut tour.signs {
                *sign = !*sign;
            }
            false
        } else {
            true
        }
    }

    fn flip_one_with_gap_model(
        &self,
        tour: &mut Tour<usize>,
        gap_model: OrientationGapModel,
    ) -> bool {
        let Some((positions, mut signs_by_id, starts)) = self.layout(tour) else {
            return false;
        };
        let mut accepted = false;
        for position in 0..tour.contigs.len() {
            let id = tour.contigs[position];
            let old_contribution: f64 = self.incident_pairs[id]
                .iter()
                .map(|&pair_index| {
                    self.evaluate_pair_with_gap_model(
                        &self.pairs[pair_index],
                        &positions,
                        &signs_by_id,
                        &starts,
                        gap_model,
                    )
                })
                .sum();
            tour.signs[position] = !tour.signs[position];
            signs_by_id[id] = !signs_by_id[id];
            let new_contribution: f64 = self.incident_pairs[id]
                .iter()
                .map(|&pair_index| {
                    self.evaluate_pair_with_gap_model(
                        &self.pairs[pair_index],
                        &positions,
                        &signs_by_id,
                        &starts,
                        gap_model,
                    )
                })
                .sum();
            if orientation_local_improves(old_contribution, new_contribution) {
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

/// Protocol-independent endpoint evidence recovered from complete CLM quartets.
/// The four mean distances identify each endpoint's mean mapped coordinate;
/// unlike log-distance potentials this is independent of the other contig's
/// length and of the chosen orientation of that contig. Link counts are capped
/// when weighting evidence: CLM does not retain independent molecule IDs.

#[derive(Clone, Copy, Debug)]
pub struct EvidenceOrientationConfig {
    pub rank_window: usize,
    pub min_links: usize,
    /// Minimum endpoint effect, not a probability of correctness.
    pub min_effect: f64,
    /// Only supplied for an explicitly trusted input tour.
    pub input_prior: f64,
    pub block_span: usize,
    pub passes: usize,
}

#[derive(Clone, Debug, Default)]
pub struct EvidenceOrientationResult {
    pub changed_signs: usize,
    pub uncertain_contigs: usize,
    pub valid_pairs: usize,
    pub invalid_pairs: usize,
    pub accepted_blocks: usize,
    pub evaluated_blocks: usize,
    pub decisions: Vec<EvidenceOrientationDecision>,
}

#[derive(Clone, Debug)]
pub struct EvidenceOrientationDecision {
    pub contig: usize,
    pub local: f64,
    pub context: f64,
    pub left: f64,
    pub right: f64,
    pub neighbours: usize,
    pub supported: bool,
}

#[derive(Clone, Copy, Debug)]
struct EndpointPair {
    u: usize,
    v: usize,
    /// Mean coordinate scaled to [-1, 1]; +1 is the forward tail.
    u_bias: f64,
    v_bias: f64,
    links: usize,
    weight: f64,
}

#[derive(Clone, Copy, Default)]
struct Vote {
    score: f64,
    weight: f64,
}

impl Vote {
    fn add(&mut self, effect: f64, weight: f64) {
        self.score += effect * weight;
        self.weight += weight;
    }

    fn mean(self) -> f64 {
        if self.weight == 0.0 {
            0.0
        } else {
            self.score / self.weight
        }
    }
}

impl OrientationProblem {
    /// Infer signs and signed blocks with the same evidence policy for every
    /// protocol and for resumed and de-novo orders. Unknown signs keep their
    /// input value, but this fallback is not treated as evidence.
    pub fn optimize_evidence(
        &self,
        tour: &mut Tour<usize>,
        config: EvidenceOrientationConfig,
    ) -> Result<EvidenceOrientationResult, String> {
        if config.rank_window == 0
            || config.rank_window > 16
            || config.min_links == 0
            || config.passes == 0
            || !config.min_effect.is_finite()
            || !(0.0..=1.0).contains(&config.min_effect)
            || !config.input_prior.is_finite()
            || config.input_prior < 0.0
        {
            return Err("invalid endpoint-evidence configuration".into());
        }
        let n = self.lengths.len();
        if tour.contigs.len() != n || self.layout(tour).is_none() {
            return Err("endpoint evidence requires a complete, unique signed tour".into());
        }
        let reverse_axis =
            tour.contigs.iter().rev().cmp(tour.contigs.iter()) == std::cmp::Ordering::Less;
        if reverse_axis {
            reverse_complement_range(tour, 0, n - 1);
        }
        let initial = tour.clone();
        let mut original_signs = vec![false; n];
        for (&id, &sign) in initial.contigs.iter().zip(&initial.signs) {
            original_signs[id] = sign;
        }
        let mut result = EvidenceOrientationResult::default();
        let mut pairs = Vec::new();
        for pair in &self.pairs {
            let links = pair.orientations[0].links;
            if links < config.min_links {
                continue;
            }
            if !pair
                .orientations
                .iter()
                .all(|s| s.present && s.links == links)
            {
                result.invalid_pairs += 1;
                continue;
            }
            let means = pair
                .orientations
                .each_ref()
                .map(|s| s.sum_distance / links as f64);
            let lu = self.lengths[pair.u] as f64;
            let lv = self.lengths[pair.v] as f64;
            // D++ = Lu-x+y, D+- = Lu-x+Lv-y, D-+ = x+y,
            // D-- = x+Lv-y. Reject quartets from incompatible coordinates.
            let tolerance = (lu + lv) * 1e-8 + 2.0;
            if lu == 0.0
                || lv == 0.0
                || (means[0] + means[3] - lu - lv).abs() > tolerance
                || (means[1] + means[2] - lu - lv).abs() > tolerance
            {
                result.invalid_pairs += 1;
                continue;
            }
            let x = (means[2] + means[3] - lv) * 0.5;
            let y = (means[0] + means[2] - lu) * 0.5;
            if x < -tolerance || x > lu + tolerance || y < -tolerance || y > lv + tolerance {
                result.invalid_pairs += 1;
                continue;
            }
            // A deliberately bounded support shrinkage, not an independence
            // assumption or a molecule-level confidence interval.
            let allowance = 1.0 / (links.min(64) as f64).sqrt();
            let shrink = |bias: f64| bias.signum() * (bias.abs() - allowance).max(0.0);
            pairs.push(EndpointPair {
                u: pair.u,
                v: pair.v,
                u_bias: shrink((2.0 * x / lu - 1.0).clamp(-1.0, 1.0)),
                v_bias: shrink((2.0 * y / lv - 1.0).clamp(-1.0, 1.0)),
                links,
                weight: (links.min(64) as f64).ln_1p(),
            });
        }
        result.valid_pairs = pairs.len();
        let mut sorted_lengths = self.lengths.clone();
        sorted_lengths.sort_unstable();
        let physical_window = sorted_lengths
            .get(n / 2)
            .copied()
            .unwrap_or(0)
            .saturating_mul(config.rank_window as u64);

        for pass in 0..=config.passes {
            let (positions, _, starts) = self.layout(tour).unwrap();
            let mut near = vec![[Vote::default(); 2]; n];
            let mut wide = vec![[Vote::default(); 2]; n];
            let mut individual = vec![Vec::new(); n];
            for pair in &pairs {
                let pu = positions[pair.u];
                let pv = positions[pair.v];
                let (lo, hi) = if pu < pv { (pu, pv) } else { (pv, pu) };
                let gap = starts[hi] - starts[lo] - self.lengths[tour.contigs[lo]];
                let radius = physical_window
                    .max(self.lengths[pair.u])
                    .max(self.lengths[pair.v]);
                if hi - lo > config.rank_window * 4 || gap > radius {
                    continue;
                }
                for (id, bias, right) in [
                    (pair.u, pair.u_bias, pu < pv),
                    (pair.v, pair.v_bias, pv < pu),
                ] {
                    let side = usize::from(right);
                    let effect = if right { bias } else { -bias };
                    wide[id][side].add(effect, pair.weight);
                    if hi - lo <= config.rank_window {
                        near[id][side].add(effect, pair.weight);
                        individual[id].push((effect, pair.weight));
                    }
                }
            }
            let mut uncertain = 0;
            let mut proposed_signs = tour.signs.clone();
            result.decisions.clear();
            for rank in 0..n {
                let id = tour.contigs[rank];
                let local = combine_sides(near[id]);
                let context = combine_sides(wide[id]);
                let proposed = local > 0.0;
                let direction = if proposed { 1.0 } else { -1.0 };
                let threshold = config.min_effect
                    + if proposed != original_signs[id] {
                        config.input_prior
                    } else {
                        0.0
                    };
                let sides_agree = near[id]
                    .iter()
                    .chain(wide[id].iter())
                    .all(|v| v.weight == 0.0 || direction * v.mean() >= 0.0);
                let total =
                    individual[id]
                        .iter()
                        .fold(Vote::default(), |mut v, &(effect, weight)| {
                            v.add(effect, weight);
                            v
                        });
                let stable = individual[id].len() < 3
                    || individual[id].iter().all(|&(effect, weight)| {
                        let remaining = total.weight - weight;
                        remaining <= f64::EPSILON
                            || direction * (total.score - effect * weight) / remaining > 0.0
                    });
                let supported = local.abs() > threshold
                    && direction * context > config.min_effect * 0.5
                    && sides_agree
                    && stable;
                result.decisions.push(EvidenceOrientationDecision {
                    contig: id,
                    local,
                    context,
                    left: near[id][0].mean(),
                    right: near[id][1].mean(),
                    neighbours: individual[id].len(),
                    supported,
                });
                if supported {
                    proposed_signs[rank] = proposed;
                } else {
                    uncertain += 1;
                }
            }
            result.uncertain_contigs = uncertain;
            if config.block_span < 2 || n < 3 || pass == config.passes {
                tour.signs = proposed_signs;
                break;
            }
            let mut best = None;
            let mut best_gain = 0.0;
            for span in 2..=config.block_span.min(n - 1) {
                for start in 0..=n - span {
                    let end = start + span - 1;
                    let cuts = [start, end + 1];
                    let old = cuts.map(|cut| boundary_evidence(tour, cut, &pairs));
                    let old_context =
                        cuts.map(|cut| boundary_context(tour, cut, &pairs, config.rank_window));
                    reverse_complement_range(tour, start, end);
                    let new = cuts.map(|cut| boundary_evidence(tour, cut, &pairs));
                    let new_context =
                        cuts.map(|cut| boundary_context(tour, cut, &pairs, config.rank_window));
                    let mut gain = 0.0;
                    let mut valid = true;
                    for boundary in 0..2 {
                        if cuts[boundary] == 0 || cuts[boundary] == n {
                            continue;
                        }
                        match new[boundary] {
                            Some((reward, links)) => {
                                let (old_reward, old_links) = old[boundary].unwrap_or((-1.0, 0));
                                let delta = reward - old_reward;
                                valid &= reward > config.min_effect
                                    && delta > config.min_effect
                                    && links.saturating_mul(2) >= old_links
                                    && new_context[boundary]
                                        > old_context[boundary] + config.min_effect * 0.25;
                                gain += delta;
                            }
                            None => valid = false,
                        }
                    }
                    reverse_complement_range(tour, start, end);
                    result.evaluated_blocks += 1;
                    if valid && gain > best_gain + 1e-12 {
                        best_gain = gain;
                        best = Some((start, end));
                    }
                }
            }
            if let Some((start, end)) = best {
                // Keep a signed block atomic: independently fixing its signs
                // before testing the reversal can destroy the candidate.
                reverse_complement_range(tour, start, end);
                result.accepted_blocks += 1;
                log::info!(
                    "Evidence block accepted: start={}, end={}, gain={:.6}",
                    start + 1,
                    end + 1,
                    best_gain
                );
            } else {
                tour.signs = proposed_signs;
                break;
            }
        }
        result.changed_signs = tour
            .contigs
            .iter()
            .zip(&tour.signs)
            .filter(|(id, sign)| original_signs[**id] != **sign)
            .count();
        if reverse_axis {
            reverse_complement_range(tour, 0, n - 1);
            for decision in &mut result.decisions {
                decision.local = -decision.local;
                decision.context = -decision.context;
                let old_left = decision.left;
                decision.left = -decision.right;
                decision.right = -old_left;
            }
            result.decisions.reverse();
        }
        Ok(result)
    }
}

fn combine_sides(sides: [Vote; 2]) -> f64 {
    let count = sides.iter().filter(|v| v.weight > 0.0).count();
    if count == 0 {
        0.0
    } else {
        sides.iter().map(|v| v.mean()).sum::<f64>() / count as f64
    }
}

fn pair_evidence(
    tour: &Tour<usize>,
    left: usize,
    right: usize,
    pairs: &[EndpointPair],
) -> Option<(f64, usize)> {
    let a = tour.contigs[left];
    let b = tour.contigs[right];
    let key = (a.min(b), a.max(b));
    let pair = &pairs[pairs.binary_search_by_key(&key, |p| (p.u, p.v)).ok()?];
    let (a_bias, b_bias) = if a == pair.u {
        (pair.u_bias, pair.v_bias)
    } else {
        (pair.v_bias, pair.u_bias)
    };
    let left_effect = if tour.signs[left] { a_bias } else { -a_bias };
    let right_effect = if tour.signs[right] { -b_bias } else { b_bias };
    Some(((left_effect + right_effect) * 0.5, pair.links))
}

fn boundary_evidence(
    tour: &Tour<usize>,
    cut: usize,
    pairs: &[EndpointPair],
) -> Option<(f64, usize)> {
    if cut == 0 || cut == tour.contigs.len() {
        return None;
    }
    pair_evidence(tour, cut - 1, cut, pairs)
}

fn boundary_context(tour: &Tour<usize>, cut: usize, pairs: &[EndpointPair], window: usize) -> f64 {
    let mut vote = Vote::default();
    for left in cut.saturating_sub(window)..cut {
        for right in cut..(cut + window).min(tour.contigs.len()) {
            if let Some((effect, links)) = pair_evidence(tour, left, right, pairs) {
                vote.add(
                    effect,
                    (links.min(64) as f64).ln_1p() / (right - left) as f64,
                );
            }
        }
    }
    vote.mean()
}

#[derive(Debug, Clone, Copy)]
pub struct OrientationResult {
    pub initial_score: f64,
    pub final_score: f64,
    pub phases: usize,
}

#[inline]
fn orientation_local_improves(old_contribution: f64, new_contribution: f64) -> bool {
    // Comparing the local terms directly avoids cancellation in
    // `global_score - old + new`. A rounded false improvement can otherwise
    // flip back on the next refinement phase and keep the outer loop alive.
    new_contribution > old_contribution
}

fn orientation_index(a_forward: bool, b_forward: bool) -> usize {
    (usize::from(!a_forward) << 1) | usize::from(!b_forward)
}

fn reverse_complement_range(tour: &mut Tour<usize>, start: usize, end: usize) {
    tour.contigs[start..=end].reverse();
    tour.signs[start..=end].reverse();
    for sign in &mut tour.signs[start..=end] {
        *sign = !*sign;
    }
}

#[inline]
fn orientation_gap(
    starts: &[u64],
    left: usize,
    right: usize,
    gap_model: OrientationGapModel,
) -> u64 {
    debug_assert!(left < right);
    match gap_model {
        OrientationGapModel::AllhicLegacy => starts[right - 1].saturating_sub(starts[left]),
        OrientationGapModel::Intervening => starts[right].saturating_sub(starts[left + 1]),
    }
}

#[inline]
fn banded_orientation_candidate_is_better(
    candidate: BandedOrientationCell,
    candidate_parent: usize,
    incumbent: BandedOrientationCell,
    incumbent_parent: usize,
) -> bool {
    candidate.score > incumbent.score
        || (candidate.score == incumbent.score
            && (candidate.changed_signs < incumbent.changed_signs
                || (candidate.changed_signs == incumbent.changed_signs
                    && candidate_parent < incumbent_parent)))
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
    /// matching ALLHiC's `generation - updated > ngen` condition. The
    /// hierarchical fast path additionally caps this patience per resolution
    /// and resets it only after a cumulative meaningful improvement.
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
    order: Arc<Vec<usize>>,
    score: f64,
    dirty: bool,
}

struct ChromosomeBufferPool<T> {
    buffers: Vec<Vec<T>>,
    chromosome_len: usize,
    max_buffers: usize,
}

impl<T: Clone> ChromosomeBufferPool<T> {
    fn new(chromosome_len: usize, max_buffers: usize) -> Self {
        let buffers = (0..max_buffers)
            .map(|_| Vec::with_capacity(chromosome_len))
            .collect();
        Self {
            buffers,
            chromosome_len,
            max_buffers,
        }
    }

    fn copy_from(&mut self, source: &[T]) -> Vec<T> {
        debug_assert_eq!(source.len(), self.chromosome_len);
        let mut buffer = self
            .buffers
            .pop()
            .unwrap_or_else(|| Vec::with_capacity(self.chromosome_len));
        buffer.clear();
        buffer.extend_from_slice(source);
        buffer
    }

    fn recycle(&mut self, chromosome: Arc<Vec<T>>) {
        if self.buffers.len() >= self.max_buffers {
            return;
        }
        if let Ok(mut buffer) = Arc::try_unwrap(chromosome) {
            buffer.clear();
            if buffer.capacity() >= self.chromosome_len {
                self.buffers.push(buffer);
            }
        }
    }

    fn replace(&mut self, target: &mut Arc<Vec<T>>, replacement: Arc<Vec<T>>) {
        let previous = std::mem::replace(target, replacement);
        self.recycle(previous);
    }
}

#[derive(Clone)]
struct UnitIndividual {
    genes: Arc<Vec<UnitGene>>,
    score: f64,
    dirty: bool,
}

#[derive(Clone, Copy, Default)]
struct FitnessDetailTimings {
    layout: Duration,
    scoring: Duration,
    evaluations: u64,
    simd_batches: u64,
    contact_edges: u64,
}

impl FitnessDetailTimings {
    fn add_assign(&mut self, other: Self) {
        self.layout += other.layout;
        self.scoring += other.scoring;
        self.evaluations += other.evaluations;
        self.simd_batches += other.simd_batches;
        self.contact_edges += other.contact_edges;
    }
}

struct OrderEvaluationSlot {
    midpoints: Vec<f64>,
    timings: FitnessDetailTimings,
}

impl OrderEvaluationSlot {
    fn new(contig_count: usize) -> Self {
        Self {
            midpoints: vec![0.0; contig_count],
            timings: FitnessDetailTimings::default(),
        }
    }
}

struct OrderEvaluationWorkspace {
    slots: Vec<OrderEvaluationSlot>,
}

impl OrderEvaluationWorkspace {
    fn new(contig_count: usize) -> Self {
        Self {
            slots: vec![OrderEvaluationSlot::new(contig_count)],
        }
    }

    fn ensure_slots(&mut self, slot_count: usize, contig_count: usize) {
        let required = slot_count.max(1);
        if self.slots.len() < required {
            self.slots
                .resize_with(required, || OrderEvaluationSlot::new(contig_count));
        }
    }

    fn timings(&self) -> FitnessDetailTimings {
        let mut total = FitnessDetailTimings::default();
        for slot in &self.slots {
            total.add_assign(slot.timings);
        }
        total
    }
}

struct UnitEvaluationScratch {
    decoded: [Vec<usize>; 4],
    midpoints: Vec<f64>,
    timings: FitnessDetailTimings,
}

impl UnitEvaluationScratch {
    fn new(contig_count: usize) -> Self {
        Self {
            // The first lane also serves scalar objectives. Allocate the
            // other three lanes lazily when a reciprocal-distance batch is
            // encountered so non-SIMD runs retain their previous footprint.
            decoded: std::array::from_fn(|lane| {
                Vec::with_capacity(if lane == 0 { contig_count } else { 0 })
            }),
            midpoints: vec![0.0; contig_count],
            timings: FitnessDetailTimings::default(),
        }
    }
}

struct UnitEvaluationWorkspace {
    scratch: Vec<UnitEvaluationScratch>,
}

impl UnitEvaluationWorkspace {
    fn new(contig_count: usize) -> Self {
        Self {
            scratch: vec![UnitEvaluationScratch::new(contig_count)],
        }
    }

    fn ensure_slots(&mut self, slot_count: usize, contig_count: usize) {
        let required = slot_count.max(1);
        if self.scratch.len() < required {
            self.scratch
                .resize_with(required, || UnitEvaluationScratch::new(contig_count));
        }
    }

    fn timings(&self) -> FitnessDetailTimings {
        let mut total = FitnessDetailTimings::default();
        for scratch in &self.scratch {
            total.add_assign(scratch.timings);
        }
        total
    }
}

fn evaluation_worker_count(
    dirty_count: usize,
    work_per_evaluation: usize,
    available_threads: usize,
) -> usize {
    if dirty_count == 0 {
        return 0;
    }
    if dirty_count == 1 || available_threads <= 1 {
        return 1;
    }

    let total_work = dirty_count.saturating_mul(work_per_evaluation);
    if total_work < MIN_PARALLEL_EVALUATION_WORK {
        1
    } else {
        available_threads.min(dirty_count)
    }
}

fn scratch_slot_budget(contig_count: usize, bytes_per_contig: usize) -> usize {
    match contig_count.checked_mul(bytes_per_contig) {
        Some(0) => usize::MAX,
        Some(bytes_per_slot) => (MAX_PERSISTENT_EVALUATION_SCRATCH_BYTES / bytes_per_slot).max(1),
        None => 1,
    }
}

fn chromosome_buffer_count(
    population_size: usize,
    chromosome_len: usize,
    bytes_per_gene: usize,
) -> usize {
    match chromosome_len.checked_mul(bytes_per_gene) {
        Some(0) => 0,
        Some(bytes_per_buffer) => {
            population_size.min((MAX_CHROMOSOME_BUFFER_POOL_BYTES / bytes_per_buffer).max(1))
        }
        None => 1,
    }
}

fn balanced_dirty_chunks_mut<'a, T>(
    items: &'a mut [T],
    dirty_count: usize,
    chunk_count: usize,
    is_dirty: impl Fn(&T) -> bool,
) -> Vec<&'a mut [T]> {
    debug_assert!(dirty_count > 0);
    let chunk_count = chunk_count.clamp(1, dirty_count);
    let mut split_ends = Vec::with_capacity(chunk_count);
    let mut dirty_seen = 0usize;
    let mut next_split = 1usize;
    let mut next_target = dirty_count / chunk_count;

    for (index, item) in items.iter().enumerate() {
        if !is_dirty(item) {
            continue;
        }
        dirty_seen += 1;
        if next_split < chunk_count && dirty_seen >= next_target {
            split_ends.push(index + 1);
            next_split += 1;
            next_target =
                ((dirty_count as u128 * next_split as u128) / chunk_count as u128) as usize;
        }
    }
    split_ends.push(items.len());

    let mut chunks = Vec::with_capacity(chunk_count);
    let mut remaining = items;
    let mut consumed = 0usize;
    for end in split_ends {
        let (chunk, tail) = remaining.split_at_mut(end - consumed);
        chunks.push(chunk);
        remaining = tail;
        consumed = end;
    }
    debug_assert!(remaining.is_empty());
    chunks
}

fn score_dirty_order_scalar(
    individual: &mut Individual,
    problem: &OptimizeProblem,
    slot: &mut OrderEvaluationSlot,
    timing_enabled: bool,
) {
    individual.score = if timing_enabled {
        problem.evaluate_with_buffer_profiled(
            individual.order.as_ref(),
            &mut slot.midpoints,
            Some(&mut slot.timings),
        )
    } else {
        problem.evaluate_with_buffer(individual.order.as_ref(), &mut slot.midpoints)
    };
    individual.dirty = false;
}

fn score_dirty_order_chunk(
    individuals: &mut [Individual],
    problem: &OptimizeProblem,
    slot: &mut OrderEvaluationSlot,
    timing_enabled: bool,
) {
    if problem.objective != OrderObjective::ReciprocalDistance || !avx_fitness_batch_supported() {
        for individual in individuals.iter_mut().filter(|individual| individual.dirty) {
            score_dirty_order_scalar(individual, problem, slot, timing_enabled);
        }
        return;
    }

    let mut dirty = individuals.iter_mut().filter(|individual| individual.dirty);
    loop {
        let Some(first) = dirty.next() else {
            break;
        };
        let Some(second) = dirty.next() else {
            score_dirty_order_scalar(first, problem, slot, timing_enabled);
            break;
        };
        let Some(third) = dirty.next() else {
            score_dirty_order_scalar(first, problem, slot, timing_enabled);
            score_dirty_order_scalar(second, problem, slot, timing_enabled);
            break;
        };
        let Some(fourth) = dirty.next() else {
            score_dirty_order_scalar(first, problem, slot, timing_enabled);
            score_dirty_order_scalar(second, problem, slot, timing_enabled);
            score_dirty_order_scalar(third, problem, slot, timing_enabled);
            break;
        };

        let orders = [
            first.order.as_slice(),
            second.order.as_slice(),
            third.order.as_slice(),
            fourth.order.as_slice(),
        ];
        let scores = if let Some(required_midpoints) = problem.contig_count().checked_mul(4) {
            if slot.midpoints.len() < required_midpoints {
                slot.midpoints.resize(required_midpoints, 0.0);
            }
            if interleaved_fitness_enabled() {
                problem.evaluate_reciprocal_interleaved_batch4(
                    orders,
                    &mut slot.midpoints,
                    timing_enabled.then_some(&mut slot.timings),
                )
            } else {
                let contig_count = problem.contig_count();
                let (positions0, remaining) = slot.midpoints.split_at_mut(contig_count);
                let (positions1, remaining) = remaining.split_at_mut(contig_count);
                let (positions2, positions3) = remaining.split_at_mut(contig_count);
                problem.evaluate_reciprocal_batch4(
                    orders,
                    [positions0, positions1, positions2, positions3],
                    timing_enabled.then_some(&mut slot.timings),
                )
            }
        } else {
            None
        };
        if let Some(scores) = scores {
            for (individual, score) in [first, second, third, fourth].into_iter().zip(scores) {
                individual.score = score;
                individual.dirty = false;
            }
        } else {
            score_dirty_order_scalar(first, problem, slot, timing_enabled);
            score_dirty_order_scalar(second, problem, slot, timing_enabled);
            score_dirty_order_scalar(third, problem, slot, timing_enabled);
            score_dirty_order_scalar(fourth, problem, slot, timing_enabled);
        }
    }
}

fn score_dirty_orders(
    individuals: &mut [Individual],
    problem: &OptimizeProblem,
    workspace: &mut OrderEvaluationWorkspace,
    available_threads: usize,
    timing_enabled: bool,
) -> usize {
    let dirty_count = individuals
        .iter()
        .filter(|individual| individual.dirty)
        .count();
    let work_per_evaluation = problem.contig_count().saturating_add(problem.edges.len());
    let worker_count = evaluation_worker_count(dirty_count, work_per_evaluation, available_threads);
    if worker_count == 0 {
        return 0;
    }

    let slot_count = worker_count
        .min(scratch_slot_budget(
            problem.contig_count(),
            std::mem::size_of::<f64>() * 4,
        ))
        .min(individuals.len());
    workspace.ensure_slots(slot_count, problem.contig_count());
    if slot_count == 1 {
        score_dirty_order_chunk(
            individuals,
            problem,
            &mut workspace.slots[0],
            timing_enabled,
        );
    } else {
        let chunks =
            balanced_dirty_chunks_mut(individuals, dirty_count, slot_count, |individual| {
                individual.dirty
            });
        chunks
            .into_par_iter()
            .zip(workspace.slots[..slot_count].par_iter_mut())
            .for_each(|(chunk, slot)| {
                score_dirty_order_chunk(chunk, problem, slot, timing_enabled);
            });
    }
    worker_count
}

fn score_dirty_unit_scalar(
    individual: &mut UnitIndividual,
    plan: &BackbonePlan,
    problem: &OptimizeProblem,
    scratch: &mut UnitEvaluationScratch,
    timing_enabled: bool,
) {
    plan.decode_into(individual.genes.as_ref(), &mut scratch.decoded[0]);
    individual.score = if timing_enabled {
        problem.evaluate_with_buffer_profiled(
            &scratch.decoded[0],
            &mut scratch.midpoints,
            Some(&mut scratch.timings),
        )
    } else {
        problem.evaluate_with_buffer(&scratch.decoded[0], &mut scratch.midpoints)
    };
    individual.dirty = false;
}

fn score_dirty_unit_chunk(
    individuals: &mut [UnitIndividual],
    plan: &BackbonePlan,
    problem: &OptimizeProblem,
    scratch: &mut UnitEvaluationScratch,
    timing_enabled: bool,
) {
    if problem.objective != OrderObjective::ReciprocalDistance || !avx_fitness_batch_supported() {
        for individual in individuals.iter_mut().filter(|individual| individual.dirty) {
            score_dirty_unit_scalar(individual, plan, problem, scratch, timing_enabled);
        }
        return;
    }

    let mut dirty = individuals.iter_mut().filter(|individual| individual.dirty);
    loop {
        let Some(first) = dirty.next() else {
            break;
        };
        let Some(second) = dirty.next() else {
            score_dirty_unit_scalar(first, plan, problem, scratch, timing_enabled);
            break;
        };
        let Some(third) = dirty.next() else {
            score_dirty_unit_scalar(first, plan, problem, scratch, timing_enabled);
            score_dirty_unit_scalar(second, plan, problem, scratch, timing_enabled);
            break;
        };
        let Some(fourth) = dirty.next() else {
            score_dirty_unit_scalar(first, plan, problem, scratch, timing_enabled);
            score_dirty_unit_scalar(second, plan, problem, scratch, timing_enabled);
            score_dirty_unit_scalar(third, plan, problem, scratch, timing_enabled);
            break;
        };

        for (decoded, individual) in scratch
            .decoded
            .iter_mut()
            .zip([&*first, &*second, &*third, &*fourth])
        {
            plan.decode_into(individual.genes.as_ref(), decoded);
        }
        let orders = scratch.decoded.each_ref().map(Vec::as_slice);
        let scores = if let Some(required_midpoints) = problem.contig_count().checked_mul(4) {
            if scratch.midpoints.len() < required_midpoints {
                scratch.midpoints.resize(required_midpoints, 0.0);
            }
            if interleaved_fitness_enabled() {
                problem.evaluate_reciprocal_interleaved_batch4(
                    orders,
                    &mut scratch.midpoints,
                    timing_enabled.then_some(&mut scratch.timings),
                )
            } else {
                let contig_count = problem.contig_count();
                let (positions0, remaining) = scratch.midpoints.split_at_mut(contig_count);
                let (positions1, remaining) = remaining.split_at_mut(contig_count);
                let (positions2, positions3) = remaining.split_at_mut(contig_count);
                problem.evaluate_reciprocal_batch4(
                    orders,
                    [positions0, positions1, positions2, positions3],
                    timing_enabled.then_some(&mut scratch.timings),
                )
            }
        } else {
            None
        };
        if let Some(scores) = scores {
            for (individual, score) in [first, second, third, fourth].into_iter().zip(scores) {
                individual.score = score;
                individual.dirty = false;
            }
        } else {
            score_dirty_unit_scalar(first, plan, problem, scratch, timing_enabled);
            score_dirty_unit_scalar(second, plan, problem, scratch, timing_enabled);
            score_dirty_unit_scalar(third, plan, problem, scratch, timing_enabled);
            score_dirty_unit_scalar(fourth, plan, problem, scratch, timing_enabled);
        }
    }
}

fn score_dirty_units(
    individuals: &mut [UnitIndividual],
    plan: &BackbonePlan,
    problem: &OptimizeProblem,
    workspace: &mut UnitEvaluationWorkspace,
    available_threads: usize,
    timing_enabled: bool,
) -> usize {
    let dirty_count = individuals
        .iter()
        .filter(|individual| individual.dirty)
        .count();
    let work_per_evaluation = problem
        .contig_count()
        .saturating_mul(2)
        .saturating_add(problem.edges.len());
    let worker_count = evaluation_worker_count(dirty_count, work_per_evaluation, available_threads);
    if worker_count == 0 {
        return 0;
    }

    let simd_batching = problem.objective == OrderObjective::ReciprocalDistance
        && avx_fitness_batch_supported()
        && dirty_count >= 4;
    let lanes_per_scratch = if simd_batching { 4 } else { 1 };
    // Keep at least four dirty individuals in every SIMD worker chunk. On
    // high-core-count hosts, splitting one individual per worker would
    // silently defeat batch4 and recreate the scalar unit-GA bottleneck.
    let batch_slot_limit = if simd_batching {
        (dirty_count / 4).max(1)
    } else {
        individuals.len()
    };
    let bytes_per_contig = (std::mem::size_of::<usize>() + std::mem::size_of::<f64>())
        .saturating_mul(lanes_per_scratch);
    let slot_count = worker_count
        .min(batch_slot_limit)
        .min(scratch_slot_budget(
            problem.contig_count(),
            bytes_per_contig,
        ))
        .min(individuals.len());
    workspace.ensure_slots(slot_count, problem.contig_count());
    if slot_count == 1 {
        score_dirty_unit_chunk(
            individuals,
            plan,
            problem,
            &mut workspace.scratch[0],
            timing_enabled,
        );
    } else {
        let chunks =
            balanced_dirty_chunks_mut(individuals, dirty_count, slot_count, |individual| {
                individual.dirty
            });
        chunks
            .into_par_iter()
            .zip(workspace.scratch[..slot_count].par_iter_mut())
            .for_each(|(chunk, scratch)| {
                score_dirty_unit_chunk(chunk, plan, problem, scratch, timing_enabled);
            });
    }
    worker_count
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

/// Run ordering with hierarchical initializer evidence as the preferred
/// coarse partition. If that evidence does not provide the configured amount
/// of dimensionality reduction, the ordinary CLM backbone discovery (and its
/// standard-GA fallback) remains unchanged.
pub fn optimize_order_with_hierarchical_joins(
    initial_tour: &Tour<usize>,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    joins: &[HierarchicalJoinEvidence],
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
        let plan =
            BackbonePlan::from_hierarchical_joins(&initial_tour.contigs, joins, &config.backbone)?;
        if plan.is_useful() {
            log::info!(
                "Using hierarchical initializer evidence for Backbone GA: {} blocks / {} units for {} contigs",
                plan.block_count(),
                plan.unit_count(),
                problem.contig_count()
            );
            return optimize_order_hierarchical_multilevel(initial_tour, problem, config, plan);
        }
        log::info!(
            "Hierarchical Backbone fallback: {} blocks / {} units for {} contigs do not meet the reduction threshold; discovering CLM backbone",
            plan.block_count(),
            plan.unit_count(),
            problem.contig_count()
        );
    }

    optimize_order(initial_tour, problem, config)
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
        let outcome = run_phase(
            &best_order,
            problem,
            config,
            phase,
            &mut rng,
            &mut reports,
            0.0,
        );
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
        0,
        0.0,
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
        let outcome = run_phase(
            &best_order,
            problem,
            config,
            phase,
            &mut rng,
            &mut reports,
            0.0,
        );
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

fn hierarchical_stage_config(
    config: &OptimizeConfig,
    search_dimension: usize,
    generation_multiplier: usize,
    minimum_generations: usize,
    maximum_generations: usize,
    maximum_patience: usize,
) -> OptimizeConfig {
    let mut stage = config.clone();
    let adaptive_limit = search_dimension
        .saturating_mul(generation_multiplier)
        .clamp(minimum_generations, maximum_generations);
    stage.max_generations = stage.max_generations.min(adaptive_limit);
    stage.stale_generations = stage.stale_generations.min(maximum_patience);
    stage
}

fn optimize_order_hierarchical_multilevel(
    initial_tour: &Tour<usize>,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    coarse_plan: BackbonePlan,
) -> Result<OptimizeResult, String> {
    let initial_score = problem.evaluate(&initial_tour.contigs);
    let seed_genes = coarse_plan.initial_genes(&initial_tour.contigs);
    let seed_order = coarse_plan.decode(&seed_genes)?;
    let seed_score = problem.evaluate(&seed_order);
    let mut rng = SmallRng::seed_from_u64(config.seed);
    let mut reports = Vec::new();
    let mut generations = Vec::with_capacity(config.phases);

    let coarse_config = hierarchical_stage_config(
        config,
        coarse_plan.unit_count(),
        10,
        5_000,
        20_000,
        HIERARCHICAL_BLOCK_PATIENCE,
    );
    log::info!(
        "Hierarchical Backbone stage block{}: {} blocks / {} units, max {} generations, patience {}",
        config.backbone.max_block_size,
        coarse_plan.block_count(),
        coarse_plan.unit_count(),
        coarse_config.max_generations,
        coarse_config.stale_generations
    );
    let coarse = run_unit_phase(
        &seed_genes,
        &coarse_plan,
        problem,
        &coarse_config,
        1,
        &mut rng,
        &mut reports,
        0,
        HIERARCHICAL_MIN_SIGNIFICANT_RELATIVE_GAIN,
    );
    log::info!(
        "Hierarchical Backbone block{} completed {} generations with fitness {:.6}",
        config.backbone.max_block_size,
        coarse.generations,
        coarse.score
    );
    let coarse_order = coarse_plan.decode(&coarse.genes)?;
    let mut phase_one_generations = coarse.generations;
    let mut phase_one_score = coarse.score;
    let mut phase_one_order = coarse_order;

    let refined_plan = coarse_plan.refined(HIERARCHICAL_REFINED_BLOCK_SIZE)?;
    if refined_plan.unit_count() > coarse_plan.unit_count() {
        let refined_seed_genes = refined_plan.initial_genes(&phase_one_order);
        let refined_config = hierarchical_stage_config(
            config,
            refined_plan.unit_count(),
            10,
            5_000,
            20_000,
            HIERARCHICAL_BLOCK_PATIENCE,
        );
        log::info!(
            "Hierarchical Backbone stage block{}: {} blocks / {} units, max {} generations, patience {}",
            HIERARCHICAL_REFINED_BLOCK_SIZE,
            refined_plan.block_count(),
            refined_plan.unit_count(),
            refined_config.max_generations,
            refined_config.stale_generations
        );
        let refined = run_unit_phase(
            &refined_seed_genes,
            &refined_plan,
            problem,
            &refined_config,
            1,
            &mut rng,
            &mut reports,
            phase_one_generations,
            HIERARCHICAL_MIN_SIGNIFICANT_RELATIVE_GAIN,
        );
        log::info!(
            "Hierarchical Backbone block{} completed {} generations with fitness {:.6}",
            HIERARCHICAL_REFINED_BLOCK_SIZE,
            refined.generations,
            refined.score
        );
        phase_one_generations = phase_one_generations.saturating_add(refined.generations);
        if refined.score < phase_one_score {
            phase_one_score = refined.score;
            phase_one_order = refined_plan.decode(&refined.genes)?;
        }
    } else {
        log::info!(
            "Hierarchical Backbone block{} refinement skipped: all coarse blocks are already at most {} contigs",
            HIERARCHICAL_REFINED_BLOCK_SIZE,
            HIERARCHICAL_REFINED_BLOCK_SIZE
        );
    }
    // Preserve the public two-phase shape: callers see block32 + block8 as
    // the combined first macro phase, followed by the unlocked phase(s).
    generations.push(phase_one_generations);

    let mut best_order = initial_tour.contigs.clone();
    let mut best_score = initial_score;
    if phase_one_score < best_score {
        best_score = phase_one_score;
        best_order = phase_one_order;
    }

    // The final macro phase works on individual contigs, so every temporary
    // hierarchical constraint can be broken if the exact ALLHiC objective
    // disagrees with an initializer join.
    let full_config = hierarchical_stage_config(
        config,
        problem.contig_count(),
        50,
        10_000,
        100_000,
        HIERARCHICAL_FULL_PATIENCE,
    );
    for phase in 2..=config.phases {
        log::info!(
            "Hierarchical Backbone stage block1 (GA{}): {} contigs, max {} generations, patience {}",
            phase,
            problem.contig_count(),
            full_config.max_generations,
            full_config.stale_generations
        );
        let outcome = run_phase(
            &best_order,
            problem,
            &full_config,
            phase,
            &mut rng,
            &mut reports,
            HIERARCHICAL_MIN_SIGNIFICANT_RELATIVE_GAIN,
        );
        log::info!(
            "Hierarchical Backbone block1 GA{} completed {} generations with fitness {:.6}",
            phase,
            outcome.generations,
            outcome.score
        );
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
            unit_count: coarse_plan.unit_count(),
            block_count: coarse_plan.block_count(),
            accepted_edges: coarse_plan.accepted_edges(),
            seed_score: Some(seed_score),
            coarse_score: Some(phase_one_score),
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

#[derive(Default)]
struct GaPhaseTimings {
    selection: Duration,
    mutation: Duration,
    fitness: Duration,
    sorting: Duration,
    bookkeeping: Duration,
}

#[derive(Clone, Copy)]
enum MutationKind {
    Swap,
    Splice,
    Insertion,
    Inversion,
}

impl MutationKind {
    const ALL: [Self; 4] = [Self::Swap, Self::Splice, Self::Insertion, Self::Inversion];

    fn index(self) -> usize {
        match self {
            Self::Swap => 0,
            Self::Splice => 1,
            Self::Insertion => 2,
            Self::Inversion => 3,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Swap => "swap",
            Self::Splice => "splice",
            Self::Insertion => "insertion",
            Self::Inversion => "inversion",
        }
    }
}

#[derive(Clone, Copy)]
struct MutationRecord {
    kind: MutationKind,
    first: usize,
    second: usize,
    span: usize,
    insertion_rotates_right: Option<bool>,
    changed: bool,
}

impl MutationRecord {
    fn span(self) -> usize {
        self.span
    }
}

struct AffectedEdgeProfiler {
    incident_edges: Vec<Vec<u32>>,
    edge_seen: Vec<u32>,
    vertex_class: Vec<u8>,
    vertex_seen: Vec<u32>,
    generation: u32,
    histograms: [Vec<u64>; 4],
    totals: [u64; 4],
    attempts: [u64; 4],
    samples: [u64; 4],
    edge_count: usize,
}

impl AffectedEdgeProfiler {
    fn new(problem: &OptimizeProblem) -> Self {
        let mut incident_edges = vec![Vec::new(); problem.contig_count()];
        for (edge_index, edge) in problem.edges.iter().enumerate() {
            incident_edges[edge.u as usize].push(edge_index as u32);
            incident_edges[edge.v as usize].push(edge_index as u32);
        }
        Self {
            incident_edges,
            edge_seen: vec![0; problem.edges.len()],
            vertex_class: vec![0; problem.contig_count()],
            vertex_seen: vec![0; problem.contig_count()],
            generation: 0,
            histograms: std::array::from_fn(|_| vec![0; problem.edges.len() + 1]),
            totals: [0; 4],
            attempts: [0; 4],
            samples: [0; 4],
            edge_count: problem.edges.len(),
        }
    }

    fn record(&mut self, order: &[usize], mutation: MutationRecord, edges: &[ContactEdge]) {
        self.begin_mutation();
        let kind_index = mutation.kind.index();
        self.attempts[kind_index] += 1;
        if (self.generation - 1) % AFFECTED_EDGE_SAMPLE_INTERVAL != 0 {
            return;
        }
        self.samples[kind_index] += 1;
        let affected = if !mutation.changed {
            0
        } else {
            match mutation.kind {
                MutationKind::Swap => self.record_swap(order, mutation, edges),
                MutationKind::Splice => self.record_splice(order, mutation, edges),
                MutationKind::Insertion => self.record_insertion(order, mutation, edges),
                MutationKind::Inversion => self.record_inversion(order, mutation, edges),
            }
        };
        self.totals[kind_index] += affected as u64;
        self.histograms[kind_index][affected] += 1;
    }

    fn begin_mutation(&mut self) {
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.edge_seen.fill(0);
            self.vertex_seen.fill(0);
            self.generation = 1;
        }
    }

    fn set_class(&mut self, contig: usize, class: u8) {
        self.vertex_seen[contig] = self.generation;
        self.vertex_class[contig] = class;
    }

    fn class(&self, contig: usize) -> u8 {
        if self.vertex_seen[contig] == self.generation {
            self.vertex_class[contig]
        } else {
            0
        }
    }

    fn count_cross_class_edges(
        &mut self,
        vertices: impl Iterator<Item = usize>,
        edges: &[ContactEdge],
        unchanged_class_pair: Option<(u8, u8)>,
    ) -> usize {
        let mut affected = 0usize;
        for contig in vertices {
            for &edge_index in &self.incident_edges[contig] {
                let edge_index = edge_index as usize;
                if self.edge_seen[edge_index] == self.generation {
                    continue;
                }
                self.edge_seen[edge_index] = self.generation;
                let edge = edges[edge_index];
                let u_class = self.class(edge.u as usize);
                let v_class = self.class(edge.v as usize);
                let unchanged_pair = unchanged_class_pair.is_some_and(|(first, second)| {
                    (u_class == first && v_class == second)
                        || (u_class == second && v_class == first)
                });
                if u_class != v_class && !unchanged_pair {
                    affected += 1;
                }
            }
        }
        affected
    }

    fn record_swap(
        &mut self,
        order: &[usize],
        mutation: MutationRecord,
        edges: &[ContactEdge],
    ) -> usize {
        let (start, end) = ordered_endpoints(mutation.first, mutation.second);
        self.set_class(order[start], 1);
        for &contig in &order[start + 1..end] {
            self.set_class(contig, 2);
        }
        self.set_class(order[end], 3);
        self.count_cross_class_edges(order[start..=end].iter().copied(), edges, Some((1, 3)))
    }

    fn record_splice(
        &mut self,
        order: &[usize],
        mutation: MutationRecord,
        edges: &[ContactEdge],
    ) -> usize {
        let cut = mutation.second;
        for &contig in &order[..cut] {
            self.set_class(contig, 1);
        }
        if cut <= order.len() - cut {
            self.count_cross_class_edges(order[..cut].iter().copied(), edges, None)
        } else {
            self.count_cross_class_edges(order[cut..].iter().copied(), edges, None)
        }
    }

    fn record_insertion(
        &mut self,
        order: &[usize],
        mutation: MutationRecord,
        edges: &[ContactEdge],
    ) -> usize {
        let (start, end) = ordered_endpoints(mutation.first, mutation.second);
        let moved_position = if mutation.insertion_rotates_right == Some(true) {
            end
        } else {
            start
        };
        for (position, &contig) in order[start..=end].iter().enumerate() {
            let position = start + position;
            self.set_class(contig, if position == moved_position { 1 } else { 2 });
        }
        self.count_cross_class_edges(order[start..=end].iter().copied(), edges, None)
    }

    fn record_inversion(
        &mut self,
        order: &[usize],
        mutation: MutationRecord,
        edges: &[ContactEdge],
    ) -> usize {
        let (start, end) = ordered_endpoints(mutation.first, mutation.second);
        for &contig in &order[start..=end] {
            self.set_class(contig, 1);
        }
        let inside_len = end - start + 1;
        let outside_len = order.len() - inside_len;
        if inside_len <= outside_len {
            self.count_cross_class_edges(order[start..=end].iter().copied(), edges, None)
        } else {
            self.count_cross_class_edges(
                order[..start].iter().chain(&order[end + 1..]).copied(),
                edges,
                None,
            )
        }
    }

    fn log(&self, label: &str, phase: usize) {
        for kind in MutationKind::ALL {
            let index = kind.index();
            let attempts = self.attempts[index];
            let samples = self.samples[index];
            let mean = if samples > 0 {
                self.totals[index] as f64 / samples as f64
            } else {
                0.0
            };
            let mean_percent = if self.edge_count > 0 {
                100.0 * mean / self.edge_count as f64
            } else {
                0.0
            };
            log::info!(
                "{}{} affected edges {}: attempts={} sampled={} sample_interval={} mean={:.2} mean_pct={:.2} p50={} p90={} p99={}",
                label,
                phase,
                kind.label(),
                attempts,
                samples,
                AFFECTED_EDGE_SAMPLE_INTERVAL,
                mean,
                mean_percent,
                histogram_quantile(&self.histograms[index], samples, 50, 100),
                histogram_quantile(&self.histograms[index], samples, 90, 100),
                histogram_quantile(&self.histograms[index], samples, 99, 100),
            );
        }
        let attempts: u64 = self.attempts.iter().sum();
        let samples: u64 = self.samples.iter().sum();
        let total: u64 = self.totals.iter().sum();
        let mean = if samples > 0 {
            total as f64 / samples as f64
        } else {
            0.0
        };
        let mean_percent = if self.edge_count > 0 {
            100.0 * mean / self.edge_count as f64
        } else {
            0.0
        };
        log::info!(
            "{}{} affected edges all: attempts={} sampled={} sample_interval={} mean={:.2} mean_pct={:.2} p50={} p90={} p99={}",
            label,
            phase,
            attempts,
            samples,
            AFFECTED_EDGE_SAMPLE_INTERVAL,
            mean,
            mean_percent,
            combined_histogram_quantile(&self.histograms, samples, 50, 100),
            combined_histogram_quantile(&self.histograms, samples, 90, 100),
            combined_histogram_quantile(&self.histograms, samples, 99, 100),
        );
    }
}

fn ordered_endpoints(first: usize, second: usize) -> (usize, usize) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn histogram_quantile(histogram: &[u64], count: u64, numerator: u64, denominator: u64) -> usize {
    if count == 0 {
        return 0;
    }
    let target =
        ((count as u128 * numerator as u128) + denominator as u128 - 1) / denominator as u128;
    let mut cumulative = 0u128;
    for (value, &frequency) in histogram.iter().enumerate() {
        cumulative += frequency as u128;
        if cumulative >= target {
            return value;
        }
    }
    histogram.len().saturating_sub(1)
}

fn combined_histogram_quantile(
    histograms: &[Vec<u64>; 4],
    count: u64,
    numerator: u64,
    denominator: u64,
) -> usize {
    if count == 0 {
        return 0;
    }
    let target =
        ((count as u128 * numerator as u128) + denominator as u128 - 1) / denominator as u128;
    let mut cumulative = 0u128;
    for value in 0..histograms[0].len() {
        cumulative += histograms
            .iter()
            .map(|histogram| histogram[value] as u128)
            .sum::<u128>();
        if cumulative >= target {
            return value;
        }
    }
    histograms[0].len().saturating_sub(1)
}

#[derive(Default)]
struct MutationStats {
    swap: u64,
    splice: u64,
    insertion: u64,
    insertion_left: u64,
    insertion_right: u64,
    inversion: u64,
    no_op: u64,
    changed_span_total: u64,
    changed_endpoint_distance_total: u64,
    changed_span_max: usize,
}

impl MutationStats {
    fn record(&mut self, mutation: MutationRecord) {
        match mutation.kind {
            MutationKind::Swap => self.swap += 1,
            MutationKind::Splice => self.splice += 1,
            MutationKind::Insertion => self.insertion += 1,
            MutationKind::Inversion => self.inversion += 1,
        }
        if let Some(rotates_right) = mutation.insertion_rotates_right {
            if rotates_right {
                self.insertion_right += 1;
            } else {
                self.insertion_left += 1;
            }
        }
        if mutation.changed {
            let span = mutation.span();
            self.changed_span_total += span as u64;
            self.changed_endpoint_distance_total += mutation.first.abs_diff(mutation.second) as u64;
            self.changed_span_max = self.changed_span_max.max(span);
        } else {
            self.no_op += 1;
        }
    }

    fn log(&self, label: &str, phase: usize) {
        let attempts = self.swap + self.splice + self.insertion + self.inversion;
        let changed = attempts.saturating_sub(self.no_op);
        let mean_span = if changed > 0 {
            self.changed_span_total as f64 / changed as f64
        } else {
            0.0
        };
        let mean_endpoint_distance = if changed > 0 {
            self.changed_endpoint_distance_total as f64 / changed as f64
        } else {
            0.0
        };
        log::info!(
            "{}{} mutations: attempts={} swap={} splice={} insertion={} insertion_left={} insertion_right={} inversion={} no_op={} mean_endpoint_distance={:.2} mean_changed_span={:.2} max_changed_span={}",
            label,
            phase,
            attempts,
            self.swap,
            self.splice,
            self.insertion,
            self.insertion_left,
            self.insertion_right,
            self.inversion,
            self.no_op,
            mean_endpoint_distance,
            mean_span,
            self.changed_span_max,
        );
    }
}

impl FitnessDetailTimings {
    fn log(&self, label: &str, phase: usize) {
        let measured = self.layout + self.scoring;
        let measured_seconds = measured.as_secs_f64();
        let scoring_percent = if measured_seconds > 0.0 {
            100.0 * self.scoring.as_secs_f64() / measured_seconds
        } else {
            0.0
        };
        let simd_evaluations = self.simd_batches.saturating_mul(4);
        let scalar_evaluations = self.evaluations.saturating_sub(simd_evaluations);
        log::info!(
            "{}{} fitness detail: evaluations={} simd_batches={} simd_evaluations={} scalar_evaluations={} contact_edges={} layout_cpu_s={:.6} scoring_cpu_s={:.6} scoring_pct={:.2}",
            label,
            phase,
            self.evaluations,
            self.simd_batches,
            simd_evaluations,
            scalar_evaluations,
            self.contact_edges,
            self.layout.as_secs_f64(),
            self.scoring.as_secs_f64(),
            scoring_percent,
        );
    }
}

impl GaPhaseTimings {
    fn log(&self, label: &str, phase: usize, generations: usize, total: Duration) {
        let accounted =
            self.selection + self.mutation + self.fitness + self.sorting + self.bookkeeping;
        let other = total.saturating_sub(accounted);
        let total_seconds = total.as_secs_f64();
        let fitness_percent = if total_seconds > 0.0 {
            100.0 * self.fitness.as_secs_f64() / total_seconds
        } else {
            0.0
        };
        log::info!(
            "{}{} timing: generations={} total_s={:.6} selection_s={:.6} mutation_s={:.6} fitness_s={:.6} fitness_pct={:.2} sort_s={:.6} bookkeeping_s={:.6} other_s={:.6}",
            label,
            phase,
            generations,
            total_seconds,
            self.selection.as_secs_f64(),
            self.mutation.as_secs_f64(),
            self.fitness.as_secs_f64(),
            fitness_percent,
            self.sorting.as_secs_f64(),
            self.bookkeeping.as_secs_f64(),
            other.as_secs_f64(),
        );
    }
}

fn profile_ga_value_enabled(value: Option<&OsStr>) -> bool {
    let Some(value) = value else {
        return false;
    };
    !matches!(
        value.to_string_lossy().trim().to_ascii_lowercase().as_str(),
        "" | "0" | "false" | "no" | "off"
    )
}

fn ga_timing_enabled() -> bool {
    profile_ga_value_enabled(std::env::var_os(PROFILE_GA_ENV).as_deref())
}

struct UnitPhaseOutcome {
    genes: Vec<UnitGene>,
    score: f64,
    generations: usize,
}

fn is_significant_improvement(reference: f64, candidate: f64, minimum_relative_gain: f64) -> bool {
    candidate < reference
        && (minimum_relative_gain <= 0.0
            || reference - candidate >= minimum_relative_gain * reference.abs().max(f64::EPSILON))
}

fn run_unit_phase(
    seed_genes: &[UnitGene],
    plan: &BackbonePlan,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    phase: usize,
    rng: &mut SmallRng,
    reports: &mut Vec<GenerationReport>,
    report_generation_offset: usize,
    minimum_significant_relative_gain: f64,
) -> UnitPhaseOutcome {
    let timing_enabled = ga_timing_enabled();
    let phase_started = timing_enabled.then(Instant::now);
    let mut timings = GaPhaseTimings::default();
    let seed_order = plan
        .decode(seed_genes)
        .expect("validated backbone seed must decode");
    let seed_score = problem.evaluate(&seed_order);
    let mut population = vec![
        UnitIndividual {
            genes: Arc::new(seed_genes.to_vec()),
            score: seed_score,
            dirty: false,
        };
        config.population_size
    ];
    let mut offspring = population.clone();
    let mut hall_of_fame = population[0].clone();
    let mut significant_score = hall_of_fame.score;
    let mut last_improvement = 0usize;
    let mut completed = 0usize;
    let available_threads = rayon::current_num_threads().max(1);
    let mut evaluation_workspace = UnitEvaluationWorkspace::new(problem.contig_count());
    let buffer_count = chromosome_buffer_count(
        config.population_size,
        plan.unit_count(),
        std::mem::size_of::<UnitGene>(),
    );
    let mut chromosome_buffers = ChromosomeBufferPool::new(plan.unit_count(), buffer_count);

    for generation in 1..=config.max_generations {
        let stage_started = timing_enabled.then(Instant::now);
        let mut offspring_index = 0;
        while offspring_index < config.population_size {
            let (first, second) = unit_tournament_pair(&population, rng);
            chromosome_buffers.replace(
                &mut offspring[offspring_index].genes,
                population[first].genes.clone(),
            );
            offspring[offspring_index].score = population[first].score;
            offspring[offspring_index].dirty = false;
            offspring_index += 1;
            if offspring_index < config.population_size {
                chromosome_buffers.replace(
                    &mut offspring[offspring_index].genes,
                    population[second].genes.clone(),
                );
                offspring[offspring_index].score = population[second].score;
                offspring[offspring_index].dirty = false;
                offspring_index += 1;
            }
        }
        if let Some(stage_started) = stage_started {
            timings.selection += stage_started.elapsed();
        }

        let stage_started = timing_enabled.then(Instant::now);
        for individual in &mut offspring {
            if rng.gen_bool(config.mutation_probability) {
                let mut genes = chromosome_buffers.copy_from(individual.genes.as_slice());
                mutate_units(&mut genes, plan, rng);
                chromosome_buffers.replace(&mut individual.genes, Arc::new(genes));
                individual.dirty = true;
            }
        }
        if let Some(stage_started) = stage_started {
            timings.mutation += stage_started.elapsed();
        }
        let stage_started = timing_enabled.then(Instant::now);
        score_dirty_units(
            &mut offspring,
            plan,
            problem,
            &mut evaluation_workspace,
            available_threads,
            timing_enabled,
        );
        if let Some(stage_started) = stage_started {
            timings.fitness += stage_started.elapsed();
        }
        for elite in offspring.iter_mut().take(ELITE_COUNT) {
            chromosome_buffers.replace(&mut elite.genes, hall_of_fame.genes.clone());
            elite.score = hall_of_fame.score;
            elite.dirty = false;
        }
        let stage_started = timing_enabled.then(Instant::now);
        offspring.sort_unstable_by(|a, b| a.score.total_cmp(&b.score));
        std::mem::swap(&mut population, &mut offspring);
        if let Some(stage_started) = stage_started {
            timings.sorting += stage_started.elapsed();
        }
        let stage_started = timing_enabled.then(Instant::now);
        completed = generation;

        if population[0].score < hall_of_fame.score {
            hall_of_fame = population[0].clone();
            if is_significant_improvement(
                significant_score,
                hall_of_fame.score,
                minimum_significant_relative_gain,
            ) {
                significant_score = hall_of_fame.score;
                last_improvement = generation;
            }
        }
        let reported_generation = report_generation_offset.saturating_add(generation);
        if config.report_interval > 0 && reported_generation % config.report_interval == 0 {
            log::info!(
                "Current iteration BackboneGA{}-{}: max_score={:.5}",
                phase,
                reported_generation,
                -hall_of_fame.score
            );
            reports.push(GenerationReport {
                phase,
                generation: reported_generation,
                best_score: hall_of_fame.score,
            });
        }
        if generation - last_improvement > config.stale_generations {
            if let Some(stage_started) = stage_started {
                timings.bookkeeping += stage_started.elapsed();
            }
            break;
        }
        if let Some(stage_started) = stage_started {
            timings.bookkeeping += stage_started.elapsed();
        }
    }

    if let Some(phase_started) = phase_started {
        timings.log("BackboneGA", phase, completed, phase_started.elapsed());
        evaluation_workspace.timings().log("BackboneGA", phase);
    }

    UnitPhaseOutcome {
        genes: hall_of_fame.genes.as_ref().to_vec(),
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
    minimum_significant_relative_gain: f64,
) -> PhaseOutcome {
    let timing_enabled = ga_timing_enabled();
    let phase_started = timing_enabled.then(Instant::now);
    let mut timings = GaPhaseTimings::default();
    let mut mutation_stats = MutationStats::default();
    let mut affected_edge_profiler = timing_enabled.then(|| AffectedEdgeProfiler::new(problem));
    // ALLHiC's MakeTour clones the same seed for every initial individual.
    let seed_score = problem.evaluate(seed_order);
    let mut population = vec![
        Individual {
            order: Arc::new(seed_order.to_vec()),
            score: seed_score,
            dirty: false,
        };
        config.population_size
    ];
    // Clean offspring share their selected parent's immutable chromosome.
    // Only offspring selected for mutation copy the underlying tour.
    let mut offspring = population.clone();
    let mut hall_of_fame = population[0].clone();
    let mut significant_score = hall_of_fame.score;
    let mut last_improvement = 0usize;
    let mut completed = 0usize;
    let available_threads = rayon::current_num_threads().max(1);
    let mut evaluation_workspace = OrderEvaluationWorkspace::new(problem.contig_count());
    let buffer_count = chromosome_buffer_count(
        config.population_size,
        problem.contig_count(),
        std::mem::size_of::<usize>(),
    );
    let mut chromosome_buffers = ChromosomeBufferPool::new(problem.contig_count(), buffer_count);

    for generation in 1..=config.max_generations {
        let stage_started = timing_enabled.then(Instant::now);
        let mut offspring_index = 0;
        while offspring_index < config.population_size {
            let (first, second) = tournament_pair(&population, rng);
            chromosome_buffers.replace(
                &mut offspring[offspring_index].order,
                population[first].order.clone(),
            );
            offspring[offspring_index].score = population[first].score;
            offspring[offspring_index].dirty = false;
            offspring_index += 1;
            if offspring_index < config.population_size {
                chromosome_buffers.replace(
                    &mut offspring[offspring_index].order,
                    population[second].order.clone(),
                );
                offspring[offspring_index].score = population[second].score;
                offspring[offspring_index].dirty = false;
                offspring_index += 1;
            }
        }
        if let Some(stage_started) = stage_started {
            timings.selection += stage_started.elapsed();
        }

        let stage_started = timing_enabled.then(Instant::now);
        for individual in &mut offspring {
            if rng.gen_bool(config.mutation_probability) {
                let mut order = chromosome_buffers.copy_from(individual.order.as_slice());
                let mutation = mutate_allhic_recorded(&mut order, rng);
                if timing_enabled {
                    mutation_stats.record(mutation);
                }
                if let Some(profiler) = affected_edge_profiler.as_mut() {
                    profiler.record(individual.order.as_slice(), mutation, &problem.edges);
                }
                chromosome_buffers.replace(&mut individual.order, Arc::new(order));
                individual.dirty = true;
            }
        }
        if let Some(stage_started) = stage_started {
            timings.mutation += stage_started.elapsed();
        }
        let stage_started = timing_enabled.then(Instant::now);
        score_dirty_orders(
            &mut offspring,
            problem,
            &mut evaluation_workspace,
            available_threads,
            timing_enabled,
        );
        if let Some(stage_started) = stage_started {
            timings.fitness += stage_started.elapsed();
        }
        for elite in offspring.iter_mut().take(ELITE_COUNT) {
            chromosome_buffers.replace(&mut elite.order, hall_of_fame.order.clone());
            elite.score = hall_of_fame.score;
            elite.dirty = false;
        }
        let stage_started = timing_enabled.then(Instant::now);
        offspring.sort_unstable_by(|a, b| a.score.total_cmp(&b.score));
        std::mem::swap(&mut population, &mut offspring);
        if let Some(stage_started) = stage_started {
            timings.sorting += stage_started.elapsed();
        }
        let stage_started = timing_enabled.then(Instant::now);
        completed = generation;

        if population[0].score < hall_of_fame.score {
            hall_of_fame = population[0].clone();
            if is_significant_improvement(
                significant_score,
                hall_of_fame.score,
                minimum_significant_relative_gain,
            ) {
                significant_score = hall_of_fame.score;
                last_improvement = generation;
            }
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
            if let Some(stage_started) = stage_started {
                timings.bookkeeping += stage_started.elapsed();
            }
            break;
        }
        if let Some(stage_started) = stage_started {
            timings.bookkeeping += stage_started.elapsed();
        }
    }

    if let Some(phase_started) = phase_started {
        timings.log("GA", phase, completed, phase_started.elapsed());
        evaluation_workspace.timings().log("GA", phase);
        mutation_stats.log("GA", phase);
        if let Some(profiler) = &affected_edge_profiler {
            profiler.log("GA", phase);
        }
    }

    PhaseOutcome {
        order: hall_of_fame.order.as_ref().to_vec(),
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
    let _ = mutate_allhic_recorded(order, rng);
}

fn mutate_allhic_recorded<R: Rng + ?Sized>(order: &mut [usize], rng: &mut R) -> MutationRecord {
    let n = order.len();
    if n <= 1 {
        return MutationRecord {
            kind: MutationKind::Swap,
            first: 0,
            second: 0,
            span: 1,
            insertion_rotates_right: None,
            changed: false,
        };
    }
    let strategy = rng.r#gen::<f64>();
    if strategy < 0.2 {
        let p = rng.gen_range(0..n);
        let q = rng.gen_range(0..n);
        order.swap(p, q);
        MutationRecord {
            kind: MutationKind::Swap,
            first: p,
            second: q,
            span: p.abs_diff(q).saturating_add(1),
            insertion_rotates_right: None,
            changed: p != q,
        }
    } else if strategy < 0.4 {
        let k = rng.gen_range(1..n);
        order.rotate_left(k);
        MutationRecord {
            kind: MutationKind::Splice,
            first: 0,
            second: k,
            span: n,
            insertion_rotates_right: None,
            changed: true,
        }
    } else {
        let mut p = rng.gen_range(0..n);
        let mut q = rng.gen_range(0..n);
        if p > q {
            std::mem::swap(&mut p, &mut q);
        }
        if p == q {
            return MutationRecord {
                kind: if strategy < 0.7 {
                    MutationKind::Insertion
                } else {
                    MutationKind::Inversion
                },
                first: p,
                second: q,
                span: 1,
                insertion_rotates_right: None,
                changed: false,
            };
        }
        if strategy < 0.7 {
            let rotates_right = rng.gen_bool(0.5);
            if rotates_right {
                order[p..=q].rotate_right(1);
            } else {
                order[p..=q].rotate_left(1);
            }
            MutationRecord {
                kind: MutationKind::Insertion,
                first: p,
                second: q,
                span: q - p + 1,
                insertion_rotates_right: Some(rotates_right),
                changed: true,
            }
        } else {
            order[p..=q].reverse();
            MutationRecord {
                kind: MutationKind::Inversion,
                first: p,
                second: q,
                span: q - p + 1,
                insertion_rotates_right: None,
                changed: true,
            }
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

#[cfg(test)]
mod orientation_precision_tests {
    use super::orientation_local_improves;

    #[test]
    fn equal_local_contributions_are_not_rounded_into_an_improvement() {
        let global_score = f64::from_bits(0xc0df_5ce4_ac5f_529d);
        let contribution = f64::from_bits(0xc0c0_8f20_d41e_e71d);

        // This is the former acceptance expression. Although the local
        // contribution is unchanged, cancellation rounds it one ULP upward.
        assert!(global_score - contribution + contribution > global_score);
        assert!(!orientation_local_improves(contribution, contribution));
    }
}

#[cfg(test)]
mod timing_tests {
    use super::{
        AffectedEdgeProfiler, BackbonePlan, ContactEdge, FitnessDetailTimings, MutationKind,
        MutationRecord, OptimizeProblem, OrderEvaluationWorkspace, OrderObjective,
        UnitEvaluationWorkspace, UnitGene, UnitIndividual, avx_fitness_batch_supported,
        combined_histogram_quantile, histogram_quantile, profile_ga_value_enabled,
        score_dirty_units,
    };
    use std::ffi::OsStr;
    use std::sync::Arc;

    #[test]
    fn ga_timing_environment_values_are_parsed_explicitly() {
        for disabled in [
            None,
            Some(""),
            Some("0"),
            Some("false"),
            Some("NO"),
            Some("off"),
        ] {
            assert!(!profile_ga_value_enabled(disabled.map(OsStr::new)));
        }
        for enabled in ["1", "true", "yes", "profile"] {
            assert!(profile_ga_value_enabled(Some(OsStr::new(enabled))));
        }
    }

    #[test]
    fn affected_edge_histogram_quantiles_use_nearest_rank() {
        let histogram = vec![0, 2, 1, 0, 1];
        assert_eq!(histogram_quantile(&histogram, 4, 50, 100), 1);
        assert_eq!(histogram_quantile(&histogram, 4, 90, 100), 4);
        let histograms = [histogram, vec![1, 0, 0, 0, 0], vec![0; 5], vec![0; 5]];
        assert_eq!(combined_histogram_quantile(&histograms, 5, 50, 100), 1);
        assert_eq!(combined_histogram_quantile(&histograms, 5, 99, 100), 4);
    }

    #[test]
    fn affected_edge_profiler_counts_cross_block_edges() {
        let mut edges = Vec::new();
        for u in 0..5u32 {
            for v in u + 1..5u32 {
                edges.push(ContactEdge { u, v, links: 1.0 });
            }
        }
        let problem = OptimizeProblem {
            lengths: vec![2.0, 3.0, 5.0, 7.0, 11.0],
            edges,
            objective: OrderObjective::ReciprocalDistance,
            anchors: vec![false; 5],
            endpoint_multiscale: None,
        };
        let order = [0, 1, 2, 3, 4];
        let cases = [
            (MutationKind::Swap, 1, 3, None, 8),
            (MutationKind::Splice, 0, 2, None, 6),
            (MutationKind::Insertion, 1, 3, Some(false), 8),
            (MutationKind::Inversion, 1, 3, None, 6),
        ];
        for (kind, first, second, insertion_rotates_right, expected) in cases {
            let mut profiler = AffectedEdgeProfiler::new(&problem);
            profiler.record(
                &order,
                MutationRecord {
                    kind,
                    first,
                    second,
                    span: first.abs_diff(second) + 1,
                    insertion_rotates_right,
                    changed: true,
                },
                &problem.edges,
            );
            assert_eq!(profiler.totals[kind.index()], expected);
        }
    }

    #[test]
    fn interleaved_avx4_matches_scalar_and_lane_major_bits() {
        if !avx_fitness_batch_supported() {
            return;
        }
        let problem = OptimizeProblem {
            lengths: vec![2.0, 3.0, 5.0, 7.0, 11.0],
            edges: vec![
                ContactEdge {
                    u: 0,
                    v: 1,
                    links: 13.0,
                },
                ContactEdge {
                    u: 3,
                    v: 4,
                    links: 17.0,
                },
                ContactEdge {
                    u: 0,
                    v: 4,
                    links: 19.0,
                },
                ContactEdge {
                    u: 1,
                    v: 3,
                    links: 23.0,
                },
            ],
            objective: OrderObjective::ReciprocalDistance,
            anchors: vec![false; 5],
            endpoint_multiscale: None,
        };
        let order_storage = [
            vec![0, 1, 2, 3, 4],
            vec![4, 3, 2, 1, 0],
            vec![1, 3, 0, 4, 2],
            vec![2, 0, 4, 1, 3],
        ];
        let orders = order_storage.each_ref().map(Vec::as_slice);

        let mut lane0 = vec![0.0; 5];
        let mut lane1 = vec![0.0; 5];
        let mut lane2 = vec![0.0; 5];
        let mut lane3 = vec![0.0; 5];
        let lane_major = problem
            .evaluate_reciprocal_batch4(
                orders,
                [&mut lane0, &mut lane1, &mut lane2, &mut lane3],
                None,
            )
            .unwrap();

        let mut interleaved = vec![0.0; 20];
        let mut timings = FitnessDetailTimings::default();
        let contig_major = problem
            .evaluate_reciprocal_interleaved_batch4(orders, &mut interleaved, Some(&mut timings))
            .unwrap();

        for lane in 0..4 {
            let mut scalar_midpoints = vec![0.0; 5];
            let scalar = problem.evaluate_with_buffer(orders[lane], &mut scalar_midpoints);
            assert_eq!(contig_major[lane].to_bits(), lane_major[lane].to_bits());
            assert_eq!(contig_major[lane].to_bits(), scalar.to_bits());
        }
        assert_eq!(timings.evaluations, 4);
        assert_eq!(timings.simd_batches, 1);
        assert_eq!(timings.contact_edges, 16);

        let mut short = vec![0.0; 19];
        assert!(
            problem
                .evaluate_reciprocal_interleaved_batch4(orders, &mut short, None)
                .is_none()
        );
        let invalid = [0, 1, 2, 3, 9];
        assert!(
            problem
                .evaluate_reciprocal_interleaved_batch4(
                    [invalid.as_slice(); 4],
                    &mut interleaved,
                    None,
                )
                .is_none()
        );
    }

    #[test]
    fn unit_batch4_matches_scalar_bits_across_parallel_chunks_and_tail() {
        if !avx_fitness_batch_supported() {
            return;
        }
        let edge_pattern = [
            ContactEdge {
                u: 0,
                v: 1,
                links: 13.0,
            },
            ContactEdge {
                u: 3,
                v: 4,
                links: 17.0,
            },
            ContactEdge {
                u: 0,
                v: 4,
                links: 19.0,
            },
            ContactEdge {
                u: 1,
                v: 3,
                links: 23.0,
            },
        ];
        let problem = OptimizeProblem {
            lengths: vec![2.0, 3.0, 5.0, 7.0, 11.0],
            // Enough work to select the parallel chunk path. Repeating an
            // edge is valid for this scorer and keeps the fixture compact.
            edges: (0..8_192)
                .map(|index| edge_pattern[index % edge_pattern.len()])
                .collect(),
            objective: OrderObjective::ReciprocalDistance,
            anchors: vec![false; 5],
            endpoint_multiscale: None,
        };
        let plan = BackbonePlan::singletons(5);
        let permutations = [
            [0, 1, 2, 3, 4],
            [4, 3, 2, 1, 0],
            [1, 3, 0, 4, 2],
            [2, 0, 4, 1, 3],
            [3, 1, 4, 0, 2],
            [4, 0, 2, 3, 1],
            [0, 2, 4, 3, 1],
            [1, 4, 2, 0, 3],
            [3, 0, 1, 4, 2],
        ];
        let mut individuals = permutations
            .into_iter()
            .map(|permutation| UnitIndividual {
                genes: Arc::new(
                    permutation
                        .into_iter()
                        .map(|unit| UnitGene {
                            unit,
                            reversed: false,
                        })
                        .collect(),
                ),
                score: f64::INFINITY,
                dirty: true,
            })
            .collect::<Vec<_>>();
        let expected = individuals
            .iter()
            .map(|individual| {
                problem
                    .evaluate(&plan.decode(individual.genes.as_ref()).unwrap())
                    .to_bits()
            })
            .collect::<Vec<_>>();

        let mut workspace = UnitEvaluationWorkspace::new(problem.contig_count());
        let workers = score_dirty_units(&mut individuals, &plan, &problem, &mut workspace, 4, true);

        assert_eq!(workers, 4);
        assert_eq!(workspace.scratch.len(), 2);
        for (individual, expected_bits) in individuals.iter().zip(expected) {
            assert!(!individual.dirty);
            assert_eq!(individual.score.to_bits(), expected_bits);
        }
        let timings = workspace.timings();
        assert_eq!(timings.evaluations, 9);
        assert_eq!(timings.simd_batches, 2);
        assert_eq!(timings.contact_edges, 8_192 * 9);
    }

    #[test]
    fn unit_non_reciprocal_objective_remains_scalar() {
        let problem = OptimizeProblem {
            lengths: vec![2.0, 3.0, 5.0, 7.0, 11.0],
            edges: vec![
                ContactEdge {
                    u: 0,
                    v: 1,
                    links: 13.0,
                },
                ContactEdge {
                    u: 0,
                    v: 4,
                    links: 19.0,
                },
            ],
            objective: OrderObjective::LogDistance,
            anchors: vec![false; 5],
            endpoint_multiscale: None,
        };
        let plan = BackbonePlan::singletons(5);
        let mut individuals = [
            [0, 1, 2, 3, 4],
            [4, 2, 0, 1, 3],
            [2, 3, 1, 4, 0],
            [1, 0, 4, 3, 2],
        ]
        .into_iter()
        .map(|permutation| UnitIndividual {
            genes: Arc::new(
                permutation
                    .into_iter()
                    .map(|unit| UnitGene {
                        unit,
                        reversed: false,
                    })
                    .collect(),
            ),
            score: f64::INFINITY,
            dirty: true,
        })
        .collect::<Vec<_>>();
        let expected = individuals
            .iter()
            .map(|individual| {
                problem
                    .evaluate(&plan.decode(individual.genes.as_ref()).unwrap())
                    .to_bits()
            })
            .collect::<Vec<_>>();

        let mut workspace = UnitEvaluationWorkspace::new(problem.contig_count());
        score_dirty_units(&mut individuals, &plan, &problem, &mut workspace, 1, true);

        for (individual, expected_bits) in individuals.iter().zip(expected) {
            assert!(!individual.dirty);
            assert_eq!(individual.score.to_bits(), expected_bits);
        }
        let timings = workspace.timings();
        assert_eq!(timings.evaluations, 4);
        assert_eq!(timings.simd_batches, 0);
        assert!(workspace.scratch[0].decoded[1..].iter().all(Vec::is_empty));
    }

    #[test]
    fn profiling_workspaces_preserve_accumulated_slots_when_work_shrinks() {
        let mut orders = OrderEvaluationWorkspace::new(5);
        orders.ensure_slots(3, 5);
        orders.slots[2].timings.evaluations = 7;
        orders.ensure_slots(1, 5);
        assert_eq!(orders.slots.len(), 3);
        assert_eq!(orders.timings().evaluations, 7);

        let mut units = UnitEvaluationWorkspace::new(5);
        units.ensure_slots(3, 5);
        units.scratch[2].timings.evaluations = 11;
        units.ensure_slots(1, 5);
        assert_eq!(units.scratch.len(), 3);
        assert_eq!(units.timings().evaluations, 11);
    }
}
