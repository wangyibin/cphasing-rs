use cphasing::core::*;
use cphasing::count_re::*;
use cphasing::fastx::Fastx;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_re() {
        let file = String::from("test/test_counts_GATC.txt");
        let mut cr = CountRE::new(&file);
        cr.parse();
        println!("{:?}", cr.file_name());
        println!("{:?}", cr.to_data());
        
    }

    #[test]
    fn test_count_re_from_fasta() {
        let fasta = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT/2.3/n50_500k/ploidy-2.3.n500k.fasta");
        let motif = String::from("GATC");
        let mut cr = CountRE::new(&fasta);
        let fasta = Fastx::new(&fasta);
        let contigsizes = fasta.get_chrom_size().unwrap();
        let counts = fasta.count_re(&motif).unwrap();
        cr.from_hashmap(counts, contigsizes);
        
        cr.write(&String::from("test/test_counts_GATC.txt"));
        
    }
}