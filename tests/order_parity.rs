use cphasing::order::{ContactMatrix, calculate_allhic_fitness, calculate_fitness};
use hashbrown::HashMap;
use indexmap::IndexMap;

fn add_symmetric_contact(contacts: &mut [HashMap<usize, u32>], u: usize, v: usize, weight: u32) {
    contacts[u].insert(v, weight);
    contacts[v].insert(u, weight);
}

fn fixture() -> (IndexMap<usize, usize>, Vec<HashMap<usize, u32>>) {
    let lengths = IndexMap::from([(0, 100), (1, 200), (2, 300), (3, 400)]);
    let mut contacts = vec![HashMap::new(); lengths.len()];
    add_symmetric_contact(&mut contacts, 0, 2, 5);
    add_symmetric_contact(&mut contacts, 2, 3, 7);
    add_symmetric_contact(&mut contacts, 0, 1, 11);
    add_symmetric_contact(&mut contacts, 1, 3, 13);
    (lengths, contacts)
}

#[test]
fn allhic_fitness_matches_fixed_dense_reference() {
    let (lengths, contacts) = fixture();
    let matrix = ContactMatrix::new(&lengths, contacts);
    let genes = [2, 0, 3, 1];

    // Midpoints in this tour are 150, 350, 600 and 900. This expression
    // mirrors ALLHiC's pairwise nlinks * ln(midpoint distance) objective.
    let expected =
        5.0 * 200.0_f64.ln() + 7.0 * 450.0_f64.ln() + 11.0 * 550.0_f64.ln() + 13.0 * 300.0_f64.ln();
    let observed = calculate_allhic_fitness(&genes, &matrix);

    assert!((observed - expected).abs() < 1e-10);
    assert_eq!(
        calculate_fitness(&genes, &matrix),
        (expected * 1_000_000.0) as isize
    );
}

#[test]
fn allhic_fitness_is_invariant_to_reversing_the_tour() {
    let (lengths, contacts) = fixture();
    let matrix = ContactMatrix::new(&lengths, contacts);
    let forward = calculate_allhic_fitness(&[2, 0, 3, 1], &matrix);
    let reverse = calculate_allhic_fitness(&[1, 3, 0, 2], &matrix);

    assert!((forward - reverse).abs() < 1e-10);
}

#[test]
fn subset_fitness_ignores_edges_to_absent_contigs_between_calls() {
    let (lengths, contacts) = fixture();
    let matrix = ContactMatrix::new(&lengths, contacts);

    let first_subset = calculate_allhic_fitness(&[0, 1], &matrix);
    let second_subset = calculate_allhic_fitness(&[2, 3], &matrix);
    let first_subset_again = calculate_allhic_fitness(&[0, 1], &matrix);

    assert!((first_subset - 11.0 * 150.0_f64.ln()).abs() < 1e-10);
    assert!((second_subset - 7.0 * 350.0_f64.ln()).abs() < 1e-10);
    assert!((first_subset - first_subset_again).abs() < 1e-10);
}
