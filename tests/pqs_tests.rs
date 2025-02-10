use cphasing::pqs;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parquet_clm() {
        let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/ploidy-2.2.pairs.pqs/q0");
        // let input_pq_dir = String::from("/data3/wangyb/0.CPhasing/pqs/82-114.output.pairs.pqs/q0");
        let _ = pqs::parquets2clm(&input_pq_dir);
    }
}