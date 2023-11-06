use cphasing::bam;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_bam() {
        let input_bam = String::from("test/test.10000.align.sort.bam");
        let output_prefix = String::from("test/test.10000.align.sort.split");
        let record_num = 1000;
        bam::split_bam(&input_bam, &output_prefix, record_num);
    }
}