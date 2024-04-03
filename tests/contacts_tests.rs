use cphasing::core::*;
use cphasing::count_re::*;
use cphasing::contacts::*;


#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_contacts() {
    //     let file = String::from("test/test.contacts");
    //     let mut pt = Contacts::new(&file);
    //     let mut cr = CountRE::new(&String::from("test/test_counts_GATC.txt"));
    //     // cr.parse();
    //     pt.parse();
    //     println!("{}", pt.file_name());
    //     pt.to_data();
    // }
    
    #[test]
    fn test_contacts() {
        let file = String::from("/data3/wangyb/0.CPhasing/0.simulation/AT_remove_inter_raw/align_data/ploidy-4.2/1m/C-Phasing/ploidy-4.2.contacts");
        let mut pt = Contacts2::new(&file);
        pt.parse();
        let contact_data = pt.to_data(&String::from("none") );
        
        for (contig_pair, count) in contact_data.iter() {
            println!("{}\t{}\t{}", contig_pair.Contig1, contig_pair.Contig2, count);
        }
    }

    // #[test]
    // fn test_contig_swap() {
    //     let contig1 = String::from("1A.ctg2");
    //     let contig2 = String::from("1A.ctg13");
    //     let mut contig1_ = &contig1; 
    //     let mut contig2_ = &contig2;
    //     let mut contig_pair = ContigPair2::new(&contig1, &contig2);
    //     contig_pair.order();
    //     println!("{:?} {:?}", contig1, contig2);
    //     println!("{:?}", contig_pair);

    //     println!("{} {}", contig1_, contig2_);
    //     if contig1_ > contig2_ {
    //         std::mem::swap(&mut contig1_, &mut contig2_);
    //     }
    //     println!("{} {}", contig1_, contig2_);
    //     println!("{:?} {:?}", contig1, contig2);

        
    // }

    // #[test]
    // fn test_from_clm() {
    //     let clm = String::from("test/test.clm");
    //     let contacts = Contacts::from_clm(&clm);
    //     contacts.write(&String::from("test/test.contacts"));
    // }

}