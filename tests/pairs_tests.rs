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

    #[test]
    fn test_split_contacts() {
        let mut pairs = Pairs::new(&String::from("test/test.pairs"));
        
        let _ = pairs.to_split_contacts(3, 2);
        
    }

}

