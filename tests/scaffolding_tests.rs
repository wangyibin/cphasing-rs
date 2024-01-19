use std::collections::{ HashMap, HashSet };
use cphasing::scaffolding::*;
use cphasing::pairs::*;
use cphasing::core::BaseTable;
use cphasing::core::ContigPair;

use rand::{thread_rng, Rng};
use rand::seq::SliceRandom;

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_split_tour_from_contacts() {
    //     let mut pairs = Pairs::new(&String::from("test/test.pairs"));
    //     let contacts = pairs.to_split_contacts(3, 2).unwrap();
    //     let mut st = SplitTour::new();
    //     let unique_min: HashMap<String, f64> = HashMap::new();
    //     let normalization_method = String::from("cis");
    //     st.from_contacts(contacts.to_data(&unique_min, &normalization_method));
    //     println!("{:?}", st);
    // }

    // #[test]
    // fn test_split_tour_get_split_contact_hash() {
    //     let mut pairs = Pairs::new(&String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_f0.9/12/n50_500k/test_pipeline/test_test_test/ploidy-12.pairs.gz"));
    //     let contacts = pairs.to_split_contacts(3, 2).unwrap();
    //     let mut st = SplitTour::new();
    //     let unique_min: HashMap<String, f64> = HashMap::new();
    //     let normalization_method = String::from("cis");
    //     let contact_hash = contacts.to_data(&unique_min, &normalization_method);

    //     st.from_contacts(&contact_hash);
    //     st.get_split_contact_hash(&contact_hash);
    //     println!("{:?}", st.split_contact_hash);
    // }

    #[test]
    fn test_split_tour_cost() {
        let mut pairs = Pairs::new(&String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_f0.9/12/n50_500k/ploidy-12.q1.pairs.gz"));
        let contacts = pairs.to_split_contacts(3, 2).unwrap();
        let mut st = SplitTour::new();
        let unique_min: HashMap<String, f64> = HashMap::new();
        let normalization_method = String::from("cis");
        let contact_hash = contacts.to_data(&unique_min, &normalization_method);
        let mut contig_set: HashSet<String> = HashSet::new();
        contig_set.insert(String::from("1A.ctg1_0"));
        contig_set.insert(String::from("1A.ctg1_1"));
        contig_set.insert(String::from("1A.ctg2_0"));
        contig_set.insert(String::from("1A.ctg2_1"));
        contig_set.insert(String::from("1A.ctg3_0"));
        contig_set.insert(String::from("1A.ctg3_1"));
        contig_set.insert(String::from("1A.ctg4_0"));
        contig_set.insert(String::from("1A.ctg4_1"));
        contig_set.insert(String::from("1A.ctg5_0"));
        contig_set.insert(String::from("1A.ctg5_1"));
        contig_set.insert(String::from("1A.ctg6_0"));
        contig_set.insert(String::from("1A.ctg6_1"));
        contig_set.insert(String::from("1A.ctg7_0"));
        contig_set.insert(String::from("1A.ctg7_1"));
        contig_set.insert(String::from("1A.ctg8_0"));
        contig_set.insert(String::from("1A.ctg8_1"));
        contig_set.insert(String::from("1A.ctg9_0"));
        contig_set.insert(String::from("1A.ctg9_1"));
        contig_set.insert(String::from("1A.ctg10_0"));
        contig_set.insert(String::from("1A.ctg10_1"));
        contig_set.insert(String::from("1A.ctg11_0"));
        contig_set.insert(String::from("1A.ctg11_1"));
        contig_set.insert(String::from("1A.ctg12_0"));
        contig_set.insert(String::from("1A.ct12_1"));
        contig_set.insert(String::from("1A.ctg13_0"));
        contig_set.insert(String::from("1A.ctg13_1"));
        contig_set.insert(String::from("1A.ctg14_0"));
        contig_set.insert(String::from("1A.ctg14_1"));


        let mut contig_list: Vec<String> = Vec::new();
        contig_list.push(String::from("1A.ctg1"));
        contig_list.push(String::from("1A.ctg2"));
        contig_list.push(String::from("1A.ctg3"));
        contig_list.push(String::from("1A.ctg4"));
        contig_list.push(String::from("1A.ctg5"));
        contig_list.push(String::from("1A.ctg6"));
        contig_list.push(String::from("1A.ctg7"));
        contig_list.push(String::from("1A.ctg8"));
        contig_list.push(String::from("1A.ctg9"));
        contig_list.push(String::from("1A.ctg10"));
        contig_list.push(String::from("1A.ctg11"));
        contig_list.push(String::from("1A.ctg12"));
        contig_list.push(String::from("1A.ctg13"));
        contig_list.push(String::from("1A.ctg14"));
        
        // filter contact_hash which Contig1 or Contig2 not in contig_list
        let mut filtered_contact_hash: HashMap<ContigPair, f64> = HashMap::new();
        for (contig_pair, contact) in contact_hash.iter() {

            if contig_set.contains(&contig_pair.Contig1) && contig_set.contains(&contig_pair.Contig2) {
                filtered_contact_hash.insert(contig_pair.clone(), *contact);
            }
        }


        st.from_contacts(&filtered_contact_hash, &contig_list);
        st.get_split_contact_hash(&filtered_contact_hash);
        println!("{:?}", st.contigs);
        println!("{}", st.get_cost());
        // random 
        let mut rng = thread_rng();
        st.contigs.shuffle(&mut rng);
        println!("{:?}", st.contigs);
        st.get_split_contact_hash(&filtered_contact_hash);

        let mut stp = SplitTourPopulation::new();
        stp.contact_hash = filtered_contact_hash;
        stp.init_population(st); 
        stp.best();
        stp.selection();
        println!("{:?}", stp.best_cost);
        println!("{:?}", stp.best_tour.contigs);
    }
}