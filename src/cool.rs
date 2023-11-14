use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use std::collections::HashMap;


use hdf5::File as H5File;
use hdf5::Group as H5Group;
use hdf5::Dataset as H5Dataset;
use hdf5::types::VarLenArray;
use hdf5::types::VarLenArray as H5VarLenArray;


pub struct Cool {
    file_name: String,
    chromsizes: HashMap<String, usize>,
    binsize: usize,
    bincount: usize,
    bin2chrom: Vec<String>,
    bin2pos: Vec<usize>,
    bin2index: HashMap<(String, usize), usize>,
    index2bin: HashMap<usize, (String, usize)>,
    matrix: Vec<Vec<f64>>,
}

impl Cool {
    pub fn new(file_name: &String) -> Cool {
        let mut cool = Cool {
            file_name: file_name.clone(),
            chromsizes: HashMap::new(),
            binsize: 0,
            bincount: 0,
            bin2chrom: Vec::new(),
            bin2pos: Vec::new(),
            bin2index: HashMap::new(),
            index2bin: HashMap::new(),
            matrix: Vec::new(),
        };
        cool.read_cool();
        cool
    }

    pub fn file_name(&self) -> &String {
        &self.file_name
    }

    pub fn chromsizes(&self) -> &HashMap<String, usize> {
        &self.chromsizes
    }

    pub fn binsize(&self) -> &usize {
        &self.binsize
    }

    pub fn bincount(&self) -> &usize {
        &self.bincount
    }

    pub fn bin2chrom(&self) -> &Vec<String> {
        &self.bin2chrom
    }

    pub fn bin2pos(&self) -> &Vec<usize> {
        &self.bin2pos
    }

    pub fn bin2index(&self) -> &HashMap<(String, usize), usize> {
        &self.bin2index
    }

    pub fn index2bin(&self) -> &HashMap<usize, (String, usize)> {
        &self.index2bin
    }

    pub fn matrix(&self) -> &Vec<Vec<f64>> {
        &self.matrix
    }

    pub fn read_cool(&mut self) {
        let file = H5File::open(&self.file_name).unwrap();
        let group = file.group("resolutions/10000").unwrap();
        let dataset = group.dataset("bins").unwrap();
        let bins: Vec<Vec<String>> = dataset.read_2d().unwrap();
        let dataset = group.dataset("bin1_offset").unwrap();
        let bin1_offset: Vec<usize> = dataset.read_1d().unwrap();
        let dataset = group.dataset("bin2_offset").unwrap();
        let bin2_offset: Vec<usize> = dataset.read_1d().unwrap();
        let dataset = group.dataset("count").unwrap();
        let count: Vec<f64> = dataset.read_1d().unwrap();
        let dataset = group.dataset("chromsize").unwrap();
        let chromsize: Vec<usize> = dataset.read_1d().unwrap();

        for i in 0..bins.len() {
            let chrom = bins[i][0].clone();
            let pos = bins[i][1].parse::<usize>().unwrap();
            self.bin2chrom.push(chrom.clone());
            self.bin2pos.push(pos);
            self.bin2index.insert((chrom.clone(), pos), i);
            self.index2bin.insert(i, (chrom.clone(), pos));
        }

        self.bincount = bins.len();
        self.binsize = 10000;
        
        for i in 0..chromsize.len() {
            let chrom = bins[i][0].clone();
            let size = chromsize[i];
            self.chromsizes.insert(chrom, size);
        }

        self.matrix = vec![vec![0.0; self.bincount]; self.bincount];
        for i in 0..self.bincount {
            for j in bin1_offset[i]..bin2_offset[i] {
                self.matrix[i][bins[j][1].parse::<usize>().unwrap()] = count[j];
            }
        }
    }

    pub fn fetch(&self, chrom1: &String, pos1: usize, chrom2: &String, pos2: usize) -> f64 {
        let bin1 = self.bin2index.get(&(chrom1.clone(), pos1)).unwrap();
        let bin2 = self.bin2index.get(&(chrom2.clone(), pos2)).unwrap();
        self.matrix[*bin1][*bin2]
    }

    pub fn fetch2(&self, bin1: usize, bin2: usize) -> f64 {
        self.matrix[bin1][bin2]
    }

    pub fn write_cool(&self, file_name: &String) {
        // write a hdf5 file which is the same as the input cool file
        
    }

}

