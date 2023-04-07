pub mod core {

use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::{ HashMap, HashSet };
use std::fmt::Display;
use std::fs::File;
use std::io::{ Write, BufReader, BufRead, Error };
use std::io;
use std::ffi::OsStr;
use std::path::Path; 
use std::hash::Hash;
use std::error::Error as OtherError;
use serde::{ Deserialize, Serialize };
use itertools::Itertools;


pub trait BaseTable {
    fn new(name: &String) -> Self;

    fn file_name(&self) -> Cow<'_, str>;
    
}

#[derive(Debug)]
pub struct CountReRecord {
    Contig: String, 
    RECounts: u32,
    Length: u32,
}

pub struct CountRE {
    file: String,
}

impl BaseTable for CountRE {
    fn new(name: &String) -> CountRE {
        CountRE { file: name.clone() }
    }

    fn file_name(&self) -> Cow<'_, str> {
        let path = Path::new(&self.file);
        path.file_name().expect("REASON").to_string_lossy()
    }
}

impl CountRE {
    pub fn parse(&self) -> Result<HashMap<String, CountReRecord>, Error> {
        let input = File::open(&self.file)?;
        let buffered = BufReader::new(input);
        let mut data: HashMap<String, CountReRecord> = HashMap::new();

        for line in buffered.lines() {
            let line_inner = line.unwrap();
            if line_inner.starts_with("#") {
                continue
            }
            let values: Vec<&str> = line_inner.split("\t").collect();
            if values.len() == 3 {
                let res = CountReRecord { 
                    Contig: values[0].to_string(),
                    RECounts: values[1].parse::<u32>().unwrap_or(0),
                    Length: values[2].parse::<u32>().unwrap_or(0),
                    };

                data.insert(values[0].to_string(), res);
                
            }    
        }
        Ok(data)
    }
}

#[derive(Debug)]
pub struct AlleleRecord {
    Chrom: String, 
    Pos: u32,
    Alleles: Vec<String>,
}

#[derive(Debug)]
pub struct AlleleTable {
    file: String, // with path
}


impl BaseTable for AlleleTable {
    fn new(name: &String) -> AlleleTable {
        AlleleTable { file: name.clone() }
    }

    fn file_name(&self) -> Cow<'_, str> {
        let path = Path::new(&self.file);
        path.file_name().expect("REASON").to_string_lossy()
    }

}

impl AlleleTable {
    pub fn allele_num(&self) -> Result<u32, Error> {
        let input = File::open(&self.file)?;
        let buffered = BufReader::new(input);
        
        let mut maximum_allele_num: u32 = 0;
        for line in buffered.lines() {
            let line_inner = line.unwrap();
            if line_inner.starts_with("#") {
                continue
            }
            
            let values: Vec<&str> = line_inner.trim_end().split("\t").collect();
            if values.len() - 2 > maximum_allele_num.try_into().unwrap() {
                let length = values.len();
                maximum_allele_num = length as u32 - 2;
            }
        }
        Ok(maximum_allele_num)
    }

    pub fn parse(&self) -> Result<Vec<AlleleRecord>, Error> {
        let input = File::open(&self.file)?;
        let buffered = BufReader::new(input);

        let mut data: Vec<AlleleRecord> = Vec::new();

        for line in buffered.lines() {
            let line_inner = line.unwrap();
            if line_inner.starts_with("#") {
                continue
            }
            
            let values: Vec<String> = line_inner
                                        .trim_end()
                                        .split("\t")
                                        .map(| s| s.to_string())
                                        .collect();
            
            data.push(AlleleRecord { Chrom: values[0].clone(), 
                                     Pos: values[1].parse::<u32>().unwrap_or(0), 
                                     Alleles: values[2..].to_vec()});
        }
        Ok(data)
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct PairTableRecord {
    X: u32,
    Y: u32,
    Contig1: String,
    Contig2: String,
    RE1: u32,
    RE2: u32,
    ObservedLinks: u32,
    ExpectedLinksIfAdjacent: f32,
    Label: String,
}

#[derive(PartialEq, Eq, Hash, Debug, Clone)]
pub struct ContigPair {
    pub Contig1: String,
    pub Contig2: String,
}

impl ContigPair {
    fn from_vec(vec: Vec<&String>) -> ContigPair {
        ContigPair { Contig1: (*vec[0].clone()).to_string(), 
                    Contig2: (*vec[1].clone()).to_string()}
    }

    fn reverse(&self) -> ContigPair {
        ContigPair {
            Contig1: self.Contig2.clone(),
            Contig2: self.Contig1.clone()
        }
    }

    fn reverse_by_order(&self) -> ContigPair {
        match self.Contig1.cmp(&self.Contig2) {
            Ordering::Greater => return self.reverse(),
            Ordering::Less => return self.clone(),
            _ => return self.clone()
        }
    }
}

#[derive(Debug, Clone)]
pub struct Contacts {
    // X: u32,
    // Y: u32,
    // RE1: u32,
    // RE2: u32,
    ObservedLinks: u32,
    // ExpectedLinksIfAdjacent: f32,
    // Label: String,
}

// impl PairRecord {
//     fn to_pair_table_record(&self, contig_pair: &ContigPair) -> PairTableRecord {
//         PairTableRecord {
//             X: self.X,
//             Y: self.Y,
//             Contig1: contig_pair.Contig1.clone(), 
//             Contig2: contig_pair.Contig2.clone(),
//             RE1: self.RE1,
//             RE2: self.RE2,
//             ObservedLinks: self.ObservedLinks,
//             ExpectedLinksIfAdjacent: self.ExpectedLinksIfAdjacent,
//             Label: self.Label.clone(),
//         }
//     }
// }

type PairHash = HashMap<ContigPair, Contacts>;


pub struct PairTable {
    file: String,
}   

impl BaseTable for PairTable {
    fn new(name: &String) -> PairTable {
        PairTable { file: name.clone() }
    }

    fn file_name(&self) -> Cow<'_, str> {
        let path = Path::new(&self.file);
        path.file_name().expect("REASON").to_string_lossy()
    }

}

impl PairTable {
    pub fn parse(&self, symmetric: Option<bool>) -> Result<PairHash, Box<dyn OtherError>> {
        
        let input = File::open(&self.file)?;
        let mut rdr = csv::ReaderBuilder::new()
                            .has_headers(false)
                            .comment(Some(b'#'))
                            .delimiter(b'\t')
                            .from_reader(input);

        let mut data: PairHash = HashMap::new();

        for result in rdr.deserialize() {
            let record: PairTableRecord = result?;
            let contig_pair = ContigPair { Contig1: record.Contig1.clone(), 
                                            Contig2: record.Contig2.clone()};
            let pair_record = Contacts {
                // X: record.X,
                // Y: record.Y,
                // RE1: record.RE1,
                // RE2: record.RE2,
                ObservedLinks: record.ObservedLinks,
                // ExpectedLinksIfAdjacent: record.ExpectedLinksIfAdjacent,
                // Label: record.Label,
            };
            data.insert(contig_pair, pair_record.clone());
            if symmetric.unwrap_or(false) {
                let contig_pair = ContigPair { Contig1: record.Contig2, 
                                                Contig2: record.Contig1};
                data.insert(contig_pair, pair_record);
            }
        }

        Ok(data)
    } 
}


pub fn pruner(at: AlleleTable, pt: PairTable) -> Result<(), Error> {
    
    let at_data = at.parse().unwrap();
    let pt_data = pt.parse(Some(true)).unwrap();
    let mut remove_data = HashSet::new();
    let mut retain_data = HashSet::new(); 
    
    let mut contig_db: HashMap<&String, HashSet<&String>> = HashMap::new();
    
    for (k, _) in pt_data.iter() {
        if !contig_db.contains_key(&k.Contig1) {
            let mut set: HashSet<&String> = HashSet::new();
            set.insert(&k.Contig2);
            contig_db.insert(&k.Contig1, set);
        } else {
            contig_db.get_mut(&k.Contig1).map(|val| val.insert(&k.Contig2));
        }
        
    }

    'outer: for ar in at_data.iter() {
        let alleles = &ar.Alleles;

        if alleles.len() == 1 {
            continue 'outer;
            
        } else if alleles.len() == 2 {
            let contig_pair = ContigPair { Contig1: alleles[0].clone(), 
                                            Contig2: alleles[1].clone() };
            
            if pt_data.contains_key(&contig_pair) {
                // remove_data.insert(contig_pair.reverse());
                remove_data.insert(contig_pair); 
            }

        } else {
            for allele_pair in alleles.iter().combinations(2) {
                let contig_pair = ContigPair::from_vec(allele_pair);
                
                if pt_data.contains_key(&contig_pair) {
                    // remove_data.insert(contig_pair.reverse());
                    remove_data.insert(contig_pair); 
                }
            }
        }
        let mut set: HashSet<&String> = HashSet::new();
       
        for allele1 in alleles.iter() {
            match contig_db.get(allele1) {
                Some(v) => {
                    for s in v.iter() {
                        set.insert(&s);
                    }
                },
                None => {
                   continue 'outer;
                }
            }
        }
        
        'inner: for allele2 in set.iter() {
            
            let mut values: HashMap<usize, u32> = HashMap::new();
            
            'inner2: for (i, allele1) in alleles.iter().enumerate() {
                let contig_pair = ContigPair { Contig1: allele1.clone().clone(), 
                                                Contig2: allele2.clone().clone() };
                
                
                match pt_data.get(&contig_pair) {
                    Some(v) => {
                        values.insert(i, v.ObservedLinks);
                    },
                    None => {
                        // println!("Not found {:?}", &contig_pair);
                        continue 'inner2;
                    }
                }
                
            }

            if values.len() > 1 {
                let mut max_value: &u32 = &0;
                let mut max_i: &usize = &(0 as usize);

                for (i, value) in values.iter() {
                    if value > max_value {
                        max_value = value;
                        max_i = i;
                    }    
                }
                let contig_pair = ContigPair { 
                    Contig1: alleles[*max_i].clone(), 
                    Contig2: allele2.clone().clone()
                };
                // retain_data.insert(contig_pair.reverse());
                retain_data.insert(contig_pair);
            }
            
        }
    }

    let retain_data: HashSet<_> = retain_data
                                    .difference(&remove_data)
                                    .map(|c| c.reverse_by_order())
                                    .collect();

    let mut output = File::create(&pt.file.replace(".pairs.txt", ".pairs.prune.txt"))?;

    // temp
    output.write(b"#");

    let mut rdr = csv::ReaderBuilder::new()
                        .delimiter(b'\t')
                        .has_headers(false)
                        .comment(Some(b'#'))
                        .flexible(true)
                        .from_path(&pt.file)?;

    let mut wtr = csv::WriterBuilder::new()
                        .delimiter(b'\t')
                        .from_writer(output);
   
    for result in rdr.deserialize() {
        let record: PairTableRecord = result?;
        
        let contig_pair = ContigPair { Contig1: record.Contig1.clone(), 
                                        Contig2: record.Contig2.clone()};
        
        if retain_data.contains(&contig_pair.reverse_by_order()) {
            wtr.serialize(&record)?;
        }
    }
    // wtr.write_record(&["#X", "Y", "Contig1", "Contig2", 
    //                     "RE1", "RE2", "ObservedLinks", 
    //                     "ExpectedLinksIfAdjacent", "Label"]);
    

    // for retain_pair in retain_data.iter() {
    //     let retain_pair = retain_pair.reverse_by_order();
    //     match pt_data.get(&retain_pair) {
    //         Some(v) => {
    //             wtr.serialize(v.to_pair_table_record(&retain_pair))?;
    //         },
    //         None => {
    //             continue;
    //         }
    //     }
    // }
    Ok(())
}
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::core::BaseTable;
    #[test]
    fn alleletable_filename() {
        let at = core::AlleleTable::new(&String::from("tests/Chr01.allele.table"));
        assert_eq!(at.file_name(), "Chr01.allele.table");
    }

    #[test]
    fn count_re_filename() {
        let cr = core::CountRE::new(&String::from("tests/Chr01.counts_AAGCTT.txt"));
        assert_eq!(cr.file_name(), "Chr01.counts_AAGCTT.txt");
    }

    #[test]
    fn count_re_hash() {
        let cr = core::CountRE::new(&String::from("tests/Chr01.counts_AAGCTT.txt"));
        let data = cr.parse();
        println!("{:?}", data.expect("REASON").get("utg002278l|740000|1189480"));
    }

    #[test]
    fn allele_number() {
        let at = core::AlleleTable::new(&String::from("tests/Chr01.allele.table"));
        let data = at.allele_num();
        println!("{:?}", data);
    }

    #[test]
    fn allele_data() {
        let at = core::AlleleTable::new(&String::from("tests/Chr01.allele.table"));
        let data = at.parse().unwrap();
        println!("{:?}", data[0]);
    }

    #[test]
    fn pairtable_data() {
        let pt = core::PairTable::new(&String::from("tests/Chr01.pairs.txt"));
        let data = pt.parse(Some(false)).unwrap();
        let contig_pair = core::ContigPair { Contig1: String::from("utg000023l"), 
                                        Contig2: String::from("utg001036l|0|600000") };
        println!("{:?}", &contig_pair);
        println!("{:?}", &data[&contig_pair]);
    }

    #[test]
    fn prune_test() {
        let at = core::AlleleTable::new(&String::from("tests/Chr01.allele.table"));
        let pt = core::PairTable::new(&String::from("tests/Chr01.pairs.txt"));
        match core::pruner(at, pt) {
            Ok(()) => println!("Successful..."),
            Err(e) => println!("Failaled to write data: {}", {e}),
        }
    }
}
