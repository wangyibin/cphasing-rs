use cphasing::simulation::*;


#[cfg(test)]
mod tests {
    use super::*;
 

    // #[test]
    // fn test_simulation_from_split_read() {
    //     let bam = String::from("test/test.10000.align.sort.bam");
    //     let out = String::from("test/test.10000.align.sort.split.bam");
    //     let min_quality: u8 = 40;
        
    //     simulation_from_split_read(&bam, &out, min_quality);
    // }

    #[test]
    fn test_simulate_porec() {
        let fasta = String::from("/data3/wangyb/DATA/AT/Simulates_tetraploid/data/tair10.fasta");
        let vcf = String::from("/data3/wangyb/DATA/AT/Simulates_tetraploid/data/Cvi-0.vcf.gz");
        let bed = String::from("/data3/wangyb/DATA/AT/Simulates_tetraploid/data/dorado_hac_9.4.q10.DpnII.bed");
        let output = String::from("/data3/wangyb/DATA/AT/Simulates_tetraploid/data/Cvi-0.porec.fasta");
        simulate_porec(&fasta, &vcf, &bed, &output);
    }
}