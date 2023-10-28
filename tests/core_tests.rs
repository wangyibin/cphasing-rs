use cphasing::core::*;


#[cfg(test)]
mod tests {
    use supper::*;
    #[test]
    fn test_get_tag() {
        let mut record = bam::Record::new();
        record.push_aux(b"mm", &Aux::String(b"test".to_vec())).unwrap();
        let tag_keys = ["mm", "ml"];
        let tag = get_tag(&record, &tag_keys, &parse_mm_tag);
        assert_eq!(tag, Some(Ok(("test".to_string(), "mm"))));
    }
}
