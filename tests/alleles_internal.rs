// Exercise the private allele/sketch invariants without compiling unrelated
// legacy unit tests in lib.rs. These are the production source modules.
pub use cphasing::core;
#[path = "../src/sketch.rs"]
mod sketch;
#[path = "../src/alleles.rs"]
mod alleles;
