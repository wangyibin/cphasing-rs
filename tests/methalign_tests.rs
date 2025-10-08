use cphasing::methalign;

#[cfg(test)]
mod tests {
    use statrs::function::factorial;

    use super::*;

    #[test]
    fn test_parse_bam() {
        let input_bam = String::from("/data3/wangyb/5.POJ/2.CPhasing_ul200k/mappability/methalign/test/cphasing_rs_test/contigs.reads.align.bam");
        let output_bam = String::from("/data3/wangyb/5.POJ/2.CPhasing_ul200k/mappability/methalign/test/cphasing_rs_test/output.bam");
        // let fa: Option<String> = None;
        // let bg: Option<String> = None;
        let fa: Option<String> = Some(String::from("/data3/wangyb/5.POJ/2.CPhasing_ul200k/mappability/methalign/test/cphasing_rs_test/contigs.fa"));
        let bg : Option<String> = Some(String::from("/data3/wangyb/5.POJ/2.CPhasing_ul200k/mappability/methalign/test/cphasing_rs_test/contigs.bg"));
        let threads = 30;
        methalign::parse_bam(&input_bam, 
                            &fa, &bg,
                            0,
                            2,
                            2,
                            0.5,
                            128,
                            2,
                            true, 
                            true, 
                            &output_bam,
                            threads
                            );
    }
}