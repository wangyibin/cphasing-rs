use cphasing::core::*;
use cphasing::count_re::*;
use cphasing::contacts::*;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contacts() {
        let file = String::from("test/test.contacts");
        let mut pt = Contacts::new(&file);
        let mut cr = CountRE::new(&String::from("test/test_counts_GATC.txt"));
        // cr.parse();
        pt.parse();
        println!("{}", pt.file_name());
        pt.to_data();
    }

    // #[test]
    // fn test_from_clm() {
    //     let clm = String::from("test/test.clm");
    //     let contacts = Contacts::from_clm(&clm);
    //     contacts.write(&String::from("test/test.contacts"));
    // }

}