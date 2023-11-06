use cphasing::fastx::*;
use cphasing::core::BaseTable;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_get_seq() {
        let fasta = String::from("test/test.fa");
        let fa = Fastx::new(&fasta);
        let seqs = fa.get_chrom_seqs().unwrap();
        println!("{:?}", seqs);
    }

    #[test]
    fn test_split_fastq() {
        let fastq = String::from("test/test.fq");
        split_fastq(&fastq, &String::from("test/out.split"), 1).unwrap();
    }
}