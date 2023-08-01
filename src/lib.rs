#![allow(unused_variables, unused_imports)]
pub mod bed;
pub mod cli;
pub mod cutsite;
pub mod core;
pub mod fastx;
pub mod methy;
pub mod mnd;
pub mod optimize;
pub mod paf;
pub mod pairs;
pub mod porec;
pub mod prune;
pub mod sketch;

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
        let pattern = b"GATCGATC";
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

    #[test]
    fn test_sketch() {
        let seq = String::from("AAAAACAAAAATAAGCGGGTTGACTTTTTTATATTCCCCCCCGAACCGGAACCGGGGGGGGATACGA");
        let rid: u64 = 0;
        let k: usize = 3;
        let w: usize = 3;

        let sketch = sketch::sketch(&seq, rid, k, w);
        println!("{:?}", sketch);

    }

    // #[test]
    // fn test_seq2kminmers() {
    //     use rust_seq2kminmers::KminmersIterator;
    //     let seq = b"AAAAAAAAAATAAGCTTGACTTTTTTATATTCCCCCCCGAACCGGGGGGGGGGATACGA";
    //     let iter = KminmersIterator::new(seq, 3, 5, 1.0, false).unwrap();
    //     for kmer in iter {
    //         println!("{:?}", kmer);
    //     }
    // }

    #[test]
    fn test_nthash() {
        use nthash::ntc64;
        let seq = b"AAAAAAAAAATAAGCTTGACTTTTTTATATTCCCCCCCGAACCGGGGGGGGGGATACGA";
        let n = 5;
        let hashes = ntc64(b"GCTT", 0, 4);
        println!("{:?}", hashes);
    }

    #[test]
    fn test_contig_score_table() {
        use crate::optimize::ContigScoreTable;
        let cst = ContigScoreTable::new(&String::from("test/so.score.txt"));
        let mut co = cst.read();
        // println!("{:?}", co.contig_units[2].orderidx);
        // println!("{:?}", co.contig_units[3].orderidx);
        // co.swap(2 as usize, 3 as usize);
        // println!("{:?}", co.contig_units[2].orderidx);
        // println!("{:?}", co.contig_units[3].orderidx);

        // println!("{:?}", co.contig_units[0].orderidx);
        // println!("{:?}", co.contig_units[1].orderidx);
        // println!("{:?}", co.contig_units[2].orderidx);
        // println!("{:?}", co.contig_units[3].orderidx);
        // co.reverse(0 as usize, 3 as usize);
        // println!("{:?}", co.contig_units[0].orderidx);
        // println!("{:?}", co.contig_units[1].orderidx);
        // println!("{:?}", co.contig_units[2].orderidx);
        // println!("{:?}", co.contig_units[3].orderidx);

        // println!("{:?}", co.contig_units[2].orientation);
        co.rotate(2 as usize);
        // println!("{:?}", co.contig_units[2].orientation);
        
        // println!("{}", co.cost());
    }

    #[test]
    fn test_simulated_annealing() {
        use crate::optimize::ContigScoreTable;
        let cst = ContigScoreTable::new(&String::from("test/test.so.score.txt"));
        let co = cst.read();
        let mut sa = optimize::SimulatedAnnealing::new(co, 1000.0, 0.999, 0.01, 1000000);
        let best = sa.run();
        
        println!("{:.?}", best.contigs());

        // write contig order to file
        use std::io::Write;
        let mut file = std::fs::File::create("test/test.so.order").unwrap();
        for i in 0..best.contigs().len() {
            let line = format!("{}\n", best.contigs()[i]);
            file.write_all(line.as_bytes()).unwrap();
        }
    }

    // #[test]
    // fn test_genetic_algorithm() {
    //     use crate::optimize::ContigScoreTable;
    //     let mut cst = ContigScoreTable::new(&String::from("test/so.score.txt"));
    //     let mut co = cst.read();
    //     let mut ga = optimize::GeneticAlgorithm::new(co, 100, 100, 0.01, 0.01);
    //     let best = ga.run();
        
    //     println!("{:.?}", best.contigs());
    // }

    #[test]
    fn test_fastx() {
        use crate::fastx::Fastx;
        let fastq = String::from("test/test.fq.gz");
        let fastq = Fastx::new(&fastq);
   
    }

    #[test]
    fn test_prunetable() {
        use crate::prune::PruneTable;
        let pt = PruneTable::new(&String::from("test/prune.contig.table"));

        println!("{:?}", pt.contig_pairs());
    
    }

    #[test]
    fn test_remove_by_contig_pairs() {
        use std::collections::HashSet;
        use crate::core::ContigPair;
        use crate::prune::PruneTable;
        
        let pt = PruneTable::new(&String::from("test/prune.contig.table"));
        let pairs = String::from("test/test.pairs");
        let mut pairs = pairs::Pairs::new(&pairs);
        let output = String::from("test/test.pairs.prune");
        let contigs: HashSet<ContigPair> = pt.contig_pairs().unwrap().into_iter().collect();
        pairs.remove_by_contig_pairs(contigs, &output);


    }

    #[test]
    fn test_methy_bam2fastq() {
        use crate::methy::modbam2fastq;
        let bam = String::from("test/output.hifi.call_mods.modbam.bam");
        let fastq = String::from("test/output.hifi.call_mods.modbam.fastq.gz");
        let min_prob = 0.1f32;
        let _ = modbam2fastq(&bam, min_prob, &fastq);
        
    }

    #[test]
    fn test_modify_fasta() {
        use crate::methy::modify_fasta;
        let fasta = String::from("test/hprc.chr20.100k.fasta");
        let bed = String::from("test/hprc.output.hifi.call_mods.modbam.pbmm2.freq.aggregate.all.bed");
        let output = String::from("test/hprc.chr20.100k.modify.fasta");
        let min_score = 3u32;
        let min_frac = 0.8f32;
        let _ = modify_fasta(&fasta, &bed, min_score, min_frac, &output);
    }

    #[test]
    fn test_bed() {
        use crate::bed::Bed3;
        let bed = String::from("test/output.hifi.call_mods.modbam.pbmm2.freq.aggregate.all.bed");
        let bed = Bed3::new(&bed);
        for b in bed {
            println!("{:?}", b);
        }
    }
}
