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

    #[test]
    fn test_pairs_to_clm() {
        let mut pairs = Pairs::new(&String::from("test/test.pairs"));
        let clm = pairs.to_clm(&String::from("test/test.clm"));
    }

}

