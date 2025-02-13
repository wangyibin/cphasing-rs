#![cfg_attr(debug_assertions, allow(dead_code, unused_imports))]
// use std::io::Write;
// use std::collections::HashMap;

use cphasing::fastx::*;
use cphasing::core::BaseTable;

#[cfg(test)]
mod tests {
    use super::*;
    
    // #[test]
    // fn test_get_seq() {
    //     let fasta = String::from("test/test.fa");
    //     let fa = Fastx::new(&fasta);
    //     let seqs = fa.get_chrom_seqs().unwrap();
    //     println!("{:?}", seqs);
    // }

    // #[test]
    // fn test_split_fastq() {
    //     let fastq = String::from("test/test.fq");
    //     split_fastq(&fastq, &String::from("test/out.split"), 1).unwrap();
    // }

    // #[test]
    // fn test_count_re() {
    //     // let fasta = String::from("test/test.fa");
    //     let fasta = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-2.2/500k/test.fa");
    //     let fa = Fastx::new(&fasta);
    //     let motif = String::from("GATC");
    //     let (counts, lengths) = fa.count_re(&motif).unwrap();
    //     println!("{:?}", lengths);
    // }

    // #[test] 
    // fn test_digest() {
    //     let fasta = String::from("test/test.fa");
    //     let fa = Fastx::new(&fasta);
    //     let motif = String::from("GATC");
    //     let pos = fa.digest(&motif, 50).unwrap();
    //     println!("{:?}", pos);
        
    // }

    #[test]
    fn test_slidefasta() {
        // let fasta = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-2.2/500k/test.fa");
        let fasta = String::from("/data3/wangyb/0.CPhasing/82-114/82-114_hifi.bp.p_utg.fasta");

        let fa = Fastx::new(&fasta);
        let window: u64 = 10000;
        let step: u64 = 0;
        // let output = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-2.2/500k/10k.fa");
        let output = String::from("/data3/wangyb/0.CPhasing/82-114/10k.fa");
        let _ = fa.slidefasta(&output, window, step);
    }

    // #[test]
    // fn test_slide() {
    //     let fastq = String::from("test/test.fq");
    //     let fa = Fastx::new(&fastq);
    //     let window: u64 = 10;
    //     let step: u64 = 0;
    //     let min_length: u64 = 4;
    //     let output = String::from("test/out.slide");
    //     fa.slide(&output, window, step, min_length);

    // }

    // #[test]
    // fn test_kmer_count() {
    //     let fasta = String::from("test/test.fa");
    //     let k = 3 as usize;
    //     let output = String::from("test/test.kmer.count");
    //     let fa = Fastx::new(&fasta);
    //     let kmer_count: HashMap<String, u64> = fa.kmer_count(k).unwrap();
        
    //     let writer = std::fs::File::create(&output).unwrap();
    //     let mut writer = std::io::BufWriter::new(writer);
    //     for (kmer, count) in kmer_count.iter() {
    //         writer.write_all(format!("{}\t{}\n", kmer, count).as_bytes()).unwrap();
    //     }
    // }

    // #[test]
    // fn test_mask_high_frequency_kmer() {
    //     let fasta = String::from("test/test.fa");
    //     let k = 3 as usize;
    //     let output = String::from("test/test.masked.fa");
    //     let fa = Fastx::new(&fasta);
        
    //     let _ = fa.mask_high_frequency_kmer(k, 10, &output);
    // }
    // #[test]
    // fn test_kmer_position() {
    //     let fasta = String::from("test/test.fa");
    //     let k = 3 as usize;
    //     let kmer_list = String::from("test/kmer.list");
    //     let output = String::from("test/test.kmer.position.bed");

    //     let fa = Fastx::new(&fasta);
        
    //     fa.kmer_positions(k, &kmer_list, &output);
    // }
}