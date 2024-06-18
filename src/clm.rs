use anyhow::Result as anyResult;
use std::collections::{HashMap, HashSet};
use std::io::{ Write, BufReader, BufRead };

// developing


use crate::core::{ common_reader, common_writer };


#[derive(Debug, Clone)]
pub struct Clm {
    file: String
}

impl BaseTable for Clm {
    fn new(name: &String) -> Self {
        Clm {file: name.clone()}
    }

    fn file_name(&self) -> Cow<'_, str> {
        let path = Path::new(&self.file);
        path.file_name().expect("REASON").to_string_lossy()
    }

    fn prefix(&self) -> String {
        let binding = self.file_name().to_string();
        let file_path = Path::new(&binding);
        let file_prefix = file_path.file_stem().unwrap().to_str().unwrap();

        (*file_prefix).to_string()
    }

}

impl Clm {
    pub fn split_clm(&self) {
        
    }
}