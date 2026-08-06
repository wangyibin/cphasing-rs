use cphasing::alleles::*;
use cphasing::core::BaseTable;

#[cfg(test)]
mod tests {
    use super::*;
    use std::borrow::Cow;
    use std::collections::HashMap;
    use std::error::Error;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    // #[test]
    // fn test_allele_table() {
    //     let file = "test//ploidy-2.2.100k.allele.table";
    //     let mut allele_table = AlleleTable::new(&file.to_string());
    //     let allele_records = allele_table.allele_records().unwrap();
    //     println!("{:?}", allele_table.header.to_unique_minimizer_density());
    //     // println!("{:?}", allele_records);
    // }

    // #[test]
    // fn test_allelic_contigs_by_cliques() {
    //     let file = "test/ploidy-2.2.100k.allele.table";
    //     let mut allele_table = AlleleTable::new(&file.to_string());

    //     let mut contigs = allele_table.header.contigs.clone();
    //     let contig_cliques = allele_table.get_allelic_contig_groups_by_cliques();

    // }

    // #[test]
    // fn test_allele_strand_table() {
    //     let file = "/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-4.2/50k/test_wfmash_allele/ploidy-4.2.50k.allele.strand.table";
    //     let mut allele_table = AlleleStrandTable::new(&file.to_string());
    //     let allele_records = allele_table.allele_strand_records().unwrap();
    //     let info = allele_table.get_info();
    //     println!("{:?}", info);
    // }

    // #[test]
    // fn test_allele_table() {
    //     let file = "/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-4.2/50k/test_wfmash_allele/ploidy-4.2.50k.allele.table";
    //     let mut allele_table = AlleleTable::new(&file.to_string());
    //     let allele_records = allele_table.allele_records().unwrap();
    //     println!("{:?}", allele_records);
    // }
    #[test]
    fn test_alleles() {
        let fasta = String::from(
            "/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-2.2/2m/ploidy-2.2.2m.fasta",
        );
        let mut alleles = AllelesFasta::new(&fasta);
        let mut allele = alleles.seqs();
        alleles.run(19, 19, 0.85, &String::from("output.allele.table"));
    }
}
