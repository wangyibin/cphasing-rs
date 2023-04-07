#![allow(unused_variables, unused_imports)]
pub mod cli;
pub mod cutsite;
pub mod core;
pub mod paf;
pub mod pairs;
pub mod porec;

#[cfg(test)]

mod tests {
    use super::*;
    use crate::core::BaseTable;
    #[test]
    fn test_reader() {
        let fastq = String::from("test/test.fq");
        core::common_reader(&fastq);
        
    }

    #[test]
    fn test_cut_site() {
        let fastq = String::from("test/test.fq.gz");
        let pattern = String::from("GATCGATC");
        println!("{}", fastq);
        cutsite::cut_site(fastq, pattern, String::from("test/test.out.fq.gz"));
    }

    #[test]
    fn test_chrom_size() {
        let file = String::from("test/test.contigsizes");
        
        let chromsizes: core::ChromSize = core::ChromSize::new(&file);
        // println!("{:?}", chromsizes.to_vec());
        // println!("{:?}", chromsizes.data());
    }


    #[test]
    fn test_pairs() {
        let file = String::from("test/test.pairs");
        let mut ph:pairs::PairHeader = pairs::PairHeader::new();

        ph.from_pairs(&file);
        // println!("{:?}", ph);
        
    }

    #[test]
    fn test_paftable() {
        let pt = paf::PAFTable::new(&String::from("test/test.paf"));
        println!("{}", pt.file_name());
        assert_eq!(pt.file_name(), "test.paf");

        // pt.paf2table(&"test/test.out.csv.gz".to_string(), &1, &0.75, &10).unwrap();
        
    }

    #[test]
    fn test_pore_c_table() {
        let pct = porec::PoreCTable::new(&String::from("test/test.out.csv.gz"));
        let output = String::from("test/test.output.pairs.gz");
        let chromsizes = String::from("test/test.contigsizes");

        pct.to_pairs(&chromsizes, &output).unwrap();
    }
}