use cphasing::optimize::mutate_allhic;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

const CONTIG_COUNT: usize = 50_000;
const POPULATION_SIZE: usize = 100;
const DIRTY_PER_GENERATION: usize = 20;
const GENERATIONS: usize = 250;
const SEED: u64 = 0x5eed_b0ff;

fn parent_population() -> Vec<Arc<[usize]>> {
    let base = (0..CONTIG_COUNT).collect::<Vec<_>>();
    (0..POPULATION_SIZE)
        .map(|parent| {
            let mut order = base.clone();
            order.rotate_left(parent * 499 % CONTIG_COUNT);
            Arc::from(order)
        })
        .collect()
}

fn observe(order: &[usize], checksum: &mut u64) {
    black_box(order);
    for position in [0, order.len() / 3, 2 * order.len() / 3, order.len() - 1] {
        *checksum = checksum.rotate_left(7).wrapping_add(order[position] as u64);
    }
}

fn run_arc_cow(parents: &[Arc<[usize]>], generations: usize) -> (Duration, u64) {
    let mut rng = SmallRng::seed_from_u64(SEED);
    let mut checksum = 0u64;
    let start = Instant::now();
    for _ in 0..generations {
        let mut children = Vec::with_capacity(DIRTY_PER_GENERATION);
        for _ in 0..DIRTY_PER_GENERATION {
            let parent = rng.gen_range(0..parents.len());
            let mut child = parents[parent].clone();
            mutate_allhic(Arc::make_mut(&mut child), &mut rng);
            observe(&child, &mut checksum);
            children.push(child);
        }
        black_box(&children);
    }
    (start.elapsed(), checksum)
}

fn run_reused_buffers(parents: &[Arc<[usize]>], generations: usize) -> (Duration, u64) {
    let mut rng = SmallRng::seed_from_u64(SEED);
    let mut checksum = 0u64;
    let mut pool = (0..DIRTY_PER_GENERATION)
        .map(|_| Vec::with_capacity(CONTIG_COUNT))
        .collect::<Vec<_>>();
    let start = Instant::now();
    for _ in 0..generations {
        let mut children = Vec::with_capacity(DIRTY_PER_GENERATION);
        for _ in 0..DIRTY_PER_GENERATION {
            let parent = rng.gen_range(0..parents.len());
            let mut buffer = pool.pop().unwrap();
            buffer.extend_from_slice(&parents[parent]);
            mutate_allhic(&mut buffer, &mut rng);
            observe(&buffer, &mut checksum);
            children.push(Arc::new(buffer));
        }
        black_box(&children);
        for child in children {
            let mut buffer = Arc::try_unwrap(child).unwrap();
            buffer.clear();
            pool.push(buffer);
        }
    }
    (start.elapsed(), checksum)
}

fn median(mut values: Vec<Duration>) -> Duration {
    values.sort_unstable();
    values[values.len() / 2]
}

#[test]
fn reused_buffers_preserve_every_mutated_gene() {
    let parents = parent_population();
    let mut arc_rng = SmallRng::seed_from_u64(SEED);
    let mut pooled_rng = SmallRng::seed_from_u64(SEED);
    let mut buffer = Vec::with_capacity(CONTIG_COUNT);

    for _ in 0..100 {
        let arc_parent = arc_rng.gen_range(0..parents.len());
        let pooled_parent = pooled_rng.gen_range(0..parents.len());
        assert_eq!(arc_parent, pooled_parent);

        let mut child = parents[arc_parent].clone();
        mutate_allhic(Arc::make_mut(&mut child), &mut arc_rng);
        buffer.clear();
        buffer.extend_from_slice(&parents[pooled_parent]);
        mutate_allhic(&mut buffer, &mut pooled_rng);

        assert_eq!(child.as_ref(), buffer.as_slice());
    }
}

#[test]
#[ignore = "50k chromosome allocation microbenchmark"]
fn benchmark_reused_chromosome_buffers() {
    let parents = parent_population();

    // Warm allocator pages and mutation code before measuring either path.
    let (_, arc_warmup) = run_arc_cow(&parents, 5);
    let (_, pooled_warmup) = run_reused_buffers(&parents, 5);
    assert_eq!(arc_warmup, pooled_warmup);

    let mut arc_times = Vec::new();
    let mut pooled_times = Vec::new();
    for repetition in 0..3 {
        let ((arc_time, arc_checksum), (pooled_time, pooled_checksum)) = if repetition % 2 == 0 {
            (
                run_arc_cow(&parents, GENERATIONS),
                run_reused_buffers(&parents, GENERATIONS),
            )
        } else {
            let pooled = run_reused_buffers(&parents, GENERATIONS);
            let arc = run_arc_cow(&parents, GENERATIONS);
            (arc, pooled)
        };
        assert_eq!(arc_checksum, pooled_checksum);
        arc_times.push(arc_time);
        pooled_times.push(pooled_time);
    }

    let arc = median(arc_times);
    let pooled = median(pooled_times);
    let speedup = arc.as_secs_f64() / pooled.as_secs_f64();
    eprintln!(
        "50k chromosome, {} dirty offspring: Arc COW {:?}, reused buffers {:?}, speedup {:.3}x",
        DIRTY_PER_GENERATION * GENERATIONS,
        arc,
        pooled,
        speedup
    );
}
