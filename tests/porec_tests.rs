use cphasing::porec::*;
use cphasing::core::BaseTable;
#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_merge() {
    //     let porec_files: Vec<String> = vec![
    //         String::from("/data3/wangyb/0.CPhasing/0.simulation/AT/2/ploidy-2.porec.gz"),
    //         String::from("/data3/wangyb/0.CPhasing/0.simulation/AT/2/ploidy-2.porec.gz"),
    //         String::from("/data3/wangyb/0.CPhasing/0.simulation/AT/2/ploidy-2.porec.gz"),
    //     ];

    //     merge_porec_tables(&porec_files, &String::from("test/test.merge.out.csv.gz"));
    // }

    #[test]
    fn test_to_pairs_pqs() {
        // let porec_file = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-2.2/500k/ploidy-2.2.porec.gz");
        // let output = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-2.2/500k/test.pairs.pqs");
        // let contigsizes = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-2.2/500k/ploidy-2.2.500k.contigsizes");
        let porec_file = String::from("/data3/wangyb/0.CPhasing/82-114/porec/porec.merge.porec.gz");
        let output = String::from("/data3/wangyb/0.CPhasing/82-114/porec/test.2.pqs");
        let contigsizes = String::from("/data3/wangyb/0.CPhasing/82-114/porec/82-114_hifi.bp.p_utg.contigsizes");
        let chunksize = 1000000;
        let chunksize = 1000000;
        
        let min_quality = 0 as u8;
        let min_order = 2 as usize;
        let max_order = 50 as usize;
        let mut p = PoreCTable::new(&porec_file);
        let _ = p.to_pairs_pqs(&contigsizes, &output, chunksize, min_quality, min_order, max_order);
        


    }
}