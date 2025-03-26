use cphasing::paf::*;
use cphasing::core::BaseTable;



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paf2depth() {
        let paf = PAFTable::new(&String::from("/data3/wangyb/0.CPhasing/alfalfa/1.default_porec/M1-2.paf.gz"));
        let winsize = 5000;
        let stepsize = 1000;
        let chromsize = String::from("/data3/wangyb/0.CPhasing/alfalfa/1.default_porec/M1_hifi.p_utg.contigsizes");
        let output = String::from("/data3/wangyb/0.CPhasing/alfalfa/1.default_porec/M1-2.depth");
        let min_mapq = 0 as u8;
        let _ = paf.to_depth(&chromsize, winsize, stepsize, min_mapq, &output);
    }
}