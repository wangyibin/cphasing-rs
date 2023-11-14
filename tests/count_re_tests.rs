use cphasing::core::*;
use cphasing::count_re::*;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_re() {
        let file = String::from("test/test_counts_GATC.txt");
        let mut cr = CountRE::new(&file);
        cr.parse();
        println!("{:?}", cr.file_name());
        println!("{:?}", cr.to_data());
        
    }
}