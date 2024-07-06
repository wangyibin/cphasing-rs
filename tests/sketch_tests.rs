use::cphasing::sketch;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_sketch() {
        let seq = String::from("AAAAACAAAAATAAGCGGGTTGACTTTTTTATATTCCCCCCCGAACCGGAACCGGGGGGGGATACGA");
        let rid: u64 = 0;
        let k: usize = 3;
        let w: usize = 3;

        let sketch = sketch::sketch(&seq, rid, k, w);
        println!("{:?}", sketch);

    }
}