use cphasing::aligner::*;
use cphasing::fastx::*;
use cphasing::core::BaseTable;
use std::io::{Cursor, Write};

#[cfg(test)]
mod tests {
    use super::*;
 

    #[test]
    fn test_input_bam() {
        let input_bam = String::from("test/test.10000.align.sort.bam");
        let fasta = String::from("/data3/wangyb/DATA/GS/UL/20230411-UNL230107-P4-PAK71591-sup/mod.nog.fa");
        let fa = Fastx::new(&fasta);
        let seqs = fa.get_chrom_seqs().unwrap();
        let min_quality: u8 = 10;
        let output = String::from("test/test.10000.align.sort.realign.bam");
        read_bam(&input_bam, &seqs, min_quality, &output);
    }

    use minimap2::*;

    // #[test]
    // fn test_minimap2() {
    // let seq = b"CGGCAmmAGGTTAAAATCTmAGTGCTGCAATAGGCGATTACAGTACAGCACCCAGCCTCCC";
    // let aligner = Aligner::builder()
    //                         .map_ont()
    //                         .with_cigar()
    //                         .with_sam_out()
    //                         .with_seq(seq)
    //                         .expect("Unable to build index");
    // let query = b"CGGCAmmAGGTTAAAATCTmAGTGCTGCAATAGGCGATTACAGTACAGCACCCAGCCTCCC";
    // let hits = aligner.map(query, true, true, None, None);
    // println!("hits: {:?}", hits);
    // }

    // #[test]
    // fn test_pesudo_ref() {
    //     let mut buffer = Curso::new(Vec::new());
    //     buffer.write_all(b">ctg1\n").unwrap();
    //     buffer.write_all(b"CGGCAmmAGGTTAAAATCTmAGTGCTGCAATAGGCGATTACAGTACAGCACCCAGCCTCCC\n").unwrap();
    //     buffer.write_all(b">ctg2\b").unwrap();
    //     buffer.write_all(b"CGGCAmmAGGTTAAAATCTmAGTGCTGCAATAGGCAAGCTTTTTACAGCACCCAGCCTCCC\n").unwrap();
    //     let contents = buffer.into_inner();

    //     let aligner = Aligner::builder()
    //                             .map_ont()
    //                             .with_cigar()
    //                             .with_index(&contents[..])
    //                             .expect("Unable to build index");
    //     let query = b"CGGCAmmAGGTTAAAATCTmAGTGCTGCAATAGGCGATTACAGTACAGCACCCAGCCTCCC";
    //     let hits = aligner.map(query, true, true, None, None);
    //     println!("hits: {:?}", hits);
    // }
    // #[test]
    // fn test_mapping() {
    //     let index_file = String::from("test/utg000111l.mmi");
    //     let input_bam = String::from("test/test.10000.bam");
    //     mapping(&index_file, &input_bam);
    // }
}