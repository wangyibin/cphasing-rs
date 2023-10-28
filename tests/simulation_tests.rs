use cphasing::simulation::*;


#[cfg(test)]
mod tests {
    use super::*;
 

    #[test]
    fn test_simulation_from_split_read() {
        let bam = String::from("test/test.10000.align.sort.bam");
        let out = String::from("test/test.10000.align.sort.split.bam");
        let min_quality: u8 = 40;
        
        simulation_from_split_read(&bam, &out, min_quality);
    }
}