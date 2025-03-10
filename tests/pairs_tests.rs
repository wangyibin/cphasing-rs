use cphasing::pairs::*;
use cphasing::core::BaseTable;
#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_pairs_to_bam() {
    //     let mut pairs = Pairs::new(&String::from("test/test.pairs"));
        
    //     let bam = String::from("test/tests.pairs2bam.bam");
    //     pairs.to_bam(&bam);
        
    // }

    // #[test]
    // fn test_pairs_to_contacts() {
    //     let mut pairs = Pairs::new(&String::from("test/test.pairs"));
    //     let max_counts = 3;
    //     let contacts = pairs.to_contacts(max_counts).unwrap();
    //     contacts.write(&String::from("test/test.pairs2contacts.contacts"));
    // }

    // #[test]
    // fn test_pairs_to_clm() {
    //     let mut pairs = Pairs::new(&String::from("test/test.pairs"));
    //     let clm = pairs.to_clm(3, &String::from("test/test.clm"), );
    // }
    
    // #[test]
    // fn test_intersect() {
    //     let mut pairs = Pairs::new(&String::from("test/test.pairs"));
    //     let hcr_bed = String::from("test/test.hcr.bed");
    //     pairs.intersect(&hcr_bed, &String::from("test/test.hcr.pairs"));
    // }

    // #[test]
    // fn test_split_contacts() {
    //     let mut pairs = Pairs::new(&String::from("test/test.pairs"));
        
    //     let _ = pairs.to_split_contacts(3, 2);
        
    // }

    // #[test]
    // fn test_collapse() {
    //     let mut pairs = Pairs::new(&String::from("/data3/wangyb/0.CPhasing/5.collapsed/B9/B9_all_porec_reads.corrected.pairs.gz"));
    //     let collapsed_list = String::from("/data3/wangyb/0.CPhasing/5.collapsed/B9/tour/asm.collapsed.contig.list");

    //     println!("{}", pairs.file_name());
    //     pairs.dup(&collapsed_list, 123,& String::from("-"));

    // }

    #[test]
    fn test_to_pqs() {
        let mut pairs = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.gz");
        let mut output = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.rust.pqs");
        let chunksize = 1000000 as usize;
        let mut pairs = Pairs::new(&pairs);
        
        let _ = pairs.to_pqs(chunksize, &output);
    }

}

