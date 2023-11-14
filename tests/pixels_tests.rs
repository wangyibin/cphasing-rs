use cphasing::core::*;
use cphasing::count_re::*;
use cphasing::pixels::*;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pixels() {
        let file = String::from("test/test.pixels");
        let mut pt = Pixels::new(&file);
        let mut cr = CountRE::new(&String::from("test/test_counts_GATC.txt"));
        // cr.parse();
        pt.parse();
        println!("{}", pt.file_name());
        pt.to_data(cr.to_data(), true);
    }
}