use cphasing::pairs::*;
use cphasing::core::BaseTable;
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pairs_to_bam() {
        let mut pairs = Pairs::new(&String::from("test/test.pairs"));
        
        let bam = String::from("test/tests.pairs2bam.bam");
        pairs.to_bam(&bam);
        
    }
}

