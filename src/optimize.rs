use crate::clm::ClmbReader;
use crate::order::Tour;
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
    decoded: Vec<usize>,
    midpoints: Vec<f64>,
    timings: FitnessDetailTimings,
}

impl UnitEvaluationScratch {
    fn new(contig_count: usize) -> Self {
        Self {
            decoded: Vec::with_capacity(contig_count),
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

fn per_individual_scratch_fits(
    population_size: usize,
    contig_count: usize,
    bytes_per_contig: usize,
) -> bool {
    population_size
        .checked_mul(contig_count)
        .and_then(|elements| elements.checked_mul(bytes_per_contig))
        .is_some_and(|bytes| bytes <= MAX_PERSISTENT_EVALUATION_SCRATCH_BYTES)
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

    if worker_count == 1 {
        let scratch = &mut workspace.scratch[0];
        for individual in individuals.iter_mut().filter(|individual| individual.dirty) {
            plan.decode_into(individual.genes.as_ref(), &mut scratch.decoded);
            individual.score = if timing_enabled {
                problem.evaluate_with_buffer_profiled(
                    &scratch.decoded,
                    &mut scratch.midpoints,
                    Some(&mut scratch.timings),
                )
            } else {
                problem.evaluate_with_buffer(&scratch.decoded, &mut scratch.midpoints)
            };
            individual.dirty = false;
        }
    } else {
        let bytes_per_contig = std::mem::size_of::<usize>() + std::mem::size_of::<f64>();
        if per_individual_scratch_fits(individuals.len(), problem.contig_count(), bytes_per_contig)
        {
            workspace.ensure_slots(individuals.len(), problem.contig_count());
            individuals
                .par_iter_mut()
                .zip(workspace.scratch.par_iter_mut())
                .filter(|(individual, _)| individual.dirty)
                .for_each(|(individual, scratch)| {
                    plan.decode_into(individual.genes.as_ref(), &mut scratch.decoded);
                    individual.score = if timing_enabled {
                        problem.evaluate_with_buffer_profiled(
                            &scratch.decoded,
                            &mut scratch.midpoints,
                            Some(&mut scratch.timings),
                        )
                    } else {
                        problem.evaluate_with_buffer(&scratch.decoded, &mut scratch.midpoints)
                    };
                    individual.dirty = false;
                });
        } else {
            let slot_count = worker_count
                .min(scratch_slot_budget(
                    problem.contig_count(),
                    bytes_per_contig,
                ))
                .min(individuals.len());
            workspace.ensure_slots(slot_count, problem.contig_count());
            if slot_count == 1 {
                let scratch = &mut workspace.scratch[0];
                for individual in individuals.iter_mut().filter(|individual| individual.dirty) {
                    plan.decode_into(individual.genes.as_ref(), &mut scratch.decoded);
                    individual.score = if timing_enabled {
                        problem.evaluate_with_buffer_profiled(
                            &scratch.decoded,
                            &mut scratch.midpoints,
                            Some(&mut scratch.timings),
                        )
                    } else {
                        problem.evaluate_with_buffer(&scratch.decoded, &mut scratch.midpoints)
                    };
                    individual.dirty = false;
                }
            } else {
                let chunks =
                    balanced_dirty_chunks_mut(individuals, dirty_count, slot_count, |individual| {
                        individual.dirty
                    });
                chunks
                    .into_par_iter()
                    .zip(workspace.scratch[..slot_count].par_iter_mut())
                    .for_each(|(chunk, scratch)| {
                        for individual in chunk.iter_mut().filter(|individual| individual.dirty) {
                            plan.decode_into(individual.genes.as_ref(), &mut scratch.decoded);
                            individual.score = if timing_enabled {
                                problem.evaluate_with_buffer_profiled(
                                    &scratch.decoded,
                                    &mut scratch.midpoints,
                                    Some(&mut scratch.timings),
                                )
                            } else {
                                problem
                                    .evaluate_with_buffer(&scratch.decoded, &mut scratch.midpoints)
                            };
                            individual.dirty = false;
                        }
                    });
            }
        }
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

fn run_unit_phase(
    seed_genes: &[UnitGene],
    plan: &BackbonePlan,
    problem: &OptimizeProblem,
    config: &OptimizeConfig,
    phase: usize,
    rng: &mut SmallRng,
    reports: &mut Vec<GenerationReport>,
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
mod timing_tests {
    use super::{
        AffectedEdgeProfiler, ContactEdge, FitnessDetailTimings, MutationKind, MutationRecord,
        OptimizeProblem, OrderEvaluationWorkspace, OrderObjective, UnitEvaluationWorkspace,
        avx_fitness_batch_supported, combined_histogram_quantile, histogram_quantile,
        profile_ga_value_enabled,
    };
    use std::ffi::OsStr;

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
