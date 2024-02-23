use cphasing::core::{ BaseTable };
use cphasing::alleles::{ AlleleRecord, AlleleTable };

#[cfg(test)]    
mod tests {
    use super::*;
    use std::path::Path;
    use std::fs::File;
    use std::io::{ BufReader, BufRead };
    use std::collections::HashMap;
    use std::error::Error;
    use std::borrow::Cow;

    // #[test]
    // fn test_allele_table() {
    //     let file = "test//ploidy-2.2.100k.allele.table";
    //     let mut allele_table = AlleleTable::new(&file.to_string());
    //     let allele_records = allele_table.allele_records().unwrap();
    //     println!("{:?}", allele_table.header.to_unique_minimizer_density());
    //     // println!("{:?}", allele_records);
    // }
    
    #[test]
    fn test_allelic_contigs_by_cliques() {
        let file = "test/ploidy-2.2.100k.allele.table";
        let mut allele_table = AlleleTable::new(&file.to_string());
       
        let mut contigs = allele_table.header.contigs.clone();
        let contig_cliques = allele_table.get_allelic_contig_groups_by_cliques();
        
      
    }
    
}