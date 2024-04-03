use std::collections::HashSet;
use cphasing::core::*;
use cphasing::prune::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prune() {
        let alleletable = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-4.2/50k/test_wfmash_allele/ploidy-4.2.50k.allele.table");
        let allele_strand_table = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-4.2/50k/test_wfmash_allele/ploidy-4.2.50k.allele.strand.table");
        let contacts = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-4.2/50k/test_wfmash_allele/ploidy-4.2.contacts");
        let prunetable = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-4.2/50k/test_wfmash_allele/prune.contig.table2");
        let normalization_method = String::from("precise");
        let whitehash: HashSet<&String> = HashSet::new();

        let mut pruner = Pruner::new(&alleletable, &allele_strand_table, &contacts,
                                    &prunetable, &normalization_method);
        
        pruner.prune();
    }
}