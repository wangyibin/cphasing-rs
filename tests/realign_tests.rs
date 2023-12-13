use cphasing::realign::*;
use cphasing::core::BaseTable;


#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_input_bam() {
    //     let bam = String::from("test/test.secondary.bam");
    //     let output = String::from("test/test.secondary.realign.bam");
    //     read_bam(&bam, &output);

    // }

    #[test]
    fn test_input_paf() {
        let paf = String::from("test/test.secondary.paf.gz");
        let output = String::from("test/test.secondary.realign.paf.gz");
        read_paf(&paf, 1, &output);
    }

}
