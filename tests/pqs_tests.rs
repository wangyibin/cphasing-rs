use cphasing::pqs;
use cphasing::core::BaseTable;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parquet_clm() {
        // let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/");
        // let output = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/output.clm.gz");
        let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/82-114.merge.pairs.pqs/");
        let output = String::from("/data3/wangyb/0.CPhasing/pqs/82-114.merge.pairs.pqs/output.clm.gz");
        // let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/82-114/hic/cphasing_output_hcr_v0.2.6_chimeric_precision/0_2.correct/test2.pqs");
        // let output = String::from("/data3/wangyb/0.CPhasing/82-114/hic/cphasing_output_hcr_v0.2.6_chimeric_precision/0_2.correct/output.clm.gz");
        let min_contacts: u32 = 1;
        let min_quality: u8 = 0;
        let output_split_contacts = true;
        let output_depth = false;
        let binsize = 10000;
        let threads = 10;
        let p = pqs::PQS::new(&input_pq_dir);

        let _ = p.to_clm(min_contacts, min_quality, &output, output_split_contacts, 
                            output_depth, binsize, threads);
    }

    // #[test]
    // fn test_parquet_mnd() {
    //     let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/");
    //     let output = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/output.mnd.txt");
    //     // let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/82-114.merge.pairs.pqs/");
    //     // let output = String::from("/data3/wangyb/0.CPhasing/pqs/82-114.merge.pairs.pqs/output.clm.gz");
        
    //     let min_quality: u8 = 0;
        
    //     let p = pqs::PQS::new(&input_pq_dir);

    //     let _ = p.to_mnd(min_quality, &output);
    // }

    // #[test]
    // fn test_parquet_intersect() {
    //     let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/");
    //     let bed = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.10000.hcr.bed");

    //     let output = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/output.intersect.pqs");
    //     let min_quality: u8 = 0;
    //     let p = pqs::PQS::new(&input_pq_dir);

    //     let _ = p.intersect(&bed, min_quality, &output);

    // }

    // #[test]
    // fn test_break_contig() {
    //     let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/");
    //     let bed = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.break.contigs.bed");
    //     let output = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/output.break.contigs.pqs");
        
    //     let mut p = pqs::PQS::new(&input_pq_dir);
    //     let _ = p.break_contigs(&bed, &output);
    // }
}