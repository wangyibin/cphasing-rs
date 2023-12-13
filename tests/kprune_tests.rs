use cphasing::core::*;
use cphasing::kprune::*;

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_prune_table() {
    //     let file = String::from("test/prune.contig.table");
    //     let pt = PruneTable::new(&file);
    //     println!("{}", pt.file_name());
    //     pt.contig_pairs().unwrap();
    // }

    // #[test]
    // fn test_maximum_biparbite_matching() {
    //     maximum_bipartite_matching();
    // }

    #[test]
    fn test_kpruner() {
        let alleletable = String::from("test/ploidy-2.3.n500k.allele.table");
        let pixels = String::from("test/test.contacts");
        let count_re = String::from("test/test_counts_GATC.txt");
        let prunetable = String::from("test/prune.contig.table");
        let method: &str = "greedy";
        let mut kpruner = KPruner::new(&alleletable, &pixels, &count_re, &prunetable);
        kpruner.prune(method);
        kpruner.prunetable.write(&String::from("test/prune.contig.table"));
    }
}