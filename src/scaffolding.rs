// a contig scaffolder, including order and orientation
// of contigs in a scaffold

use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;


struct SplitContigUnit {
    pub contig: String,
    pub cis1: u32,
    pub cis2: u32,
    pub orientation: u8,
    pub order: u32,
}

impl SplitContigUnit {
    pub fn new() -> Self {
        Self {
            contig: String::new(),
            cis1: 0,
            cis2: 0,
            orientation: 0,
            order: 0,
        }
    }
}

