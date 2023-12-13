use cphasing::core::*;
use cphasing::bed::*;
use cphasing::pairs::*;
use std::collections::HashMap;
use rust_lapper::{Interval, Lapper};

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_lapper() {
        use rust_lapper::{Interval, Lapper};
        type Iv = Interval<usize, u32>;
        let data: Vec<Iv> = vec![
            Iv {
                start: 1,
                stop: 10,
                val: 1,
            },
            Iv {
                start: 5,
                stop: 15,
                val: 2,
            },
            Iv {
                start: 12,
                stop: 20,
                val: 3,
            },
        ];
        let mut lapper = Lapper::new(data);
        lapper.find(11, 15).collect::<Vec<&Iv>>();

    }

    #[test]
    fn test_bed3() {
        type Iv = Interval<usize, u8>;
        let bed = String::from("test/hcr.bed");
        let mut bed3 = Bed3::new(&bed);
        let mut hcr: HashMap<String, Vec<Iv>> = HashMap::new();
        for i in bed3 {
            if hcr.contains_key(&i.chrom) {
                hcr.get_mut(&i.chrom).unwrap().push(Iv {
                    start: i.start,
                    stop: i.end,
                    val: 0,
                });
            } else {
                hcr.insert(i.chrom, vec![Iv {
                    start: i.start,
                    stop: i.end,
                    val: 0,
                }]);
            }
        }
        let mut lapper = Lapper::new(hcr.get(&"ctg1".to_string()).unwrap().to_vec());
        let ivs = lapper.find(100, 20000).collect::<Vec<&Iv>>();
        println!("{:?}", ivs);
    }

    #[test]
    fn test_get_interval_hash() {
        type Iv_u8 = Interval<usize, u8>;
        let bed = String::from("test/hcr.bed");
        let bed3 = Bed3::new(&bed);
        let hcr = bed3.to_interval_hash();
        let mut lapper = hcr.get(&"ctg1".to_string()).unwrap();
        let ivs = lapper.count(10000, 200000);
        println!("{:?}", ivs);
    }
}