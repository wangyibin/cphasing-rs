
use rust_htslib::bam::{ 
    self,
    record::Aux, record::CigarStringView, 
    record::Cigar, record::CigarString,
    Read, Reader, Record, HeaderView, 
    Header, header::HeaderRecord,
    Writer};

// split bam by record number and write to different files 
pub fn split_bam(input_bam: &String, output_prefix: &String, 
             record_num: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut bam = Reader::from_path(input_bam).unwrap();
    let header = Header::from_template(bam.header());
    let mut i = 0;
    let mut j = 0;
    let mut wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam).unwrap();
    log::info!("write {} records to {}", i, format!("{}_{}.bam", output_prefix, j));
    for r in bam.records() {
        let record = r?;
        wtr.write(&record)?;
        i += 1;
        if i == record_num {
            j += 1;
            wtr = Writer::from_path(format!("{}_{}.bam", output_prefix, j), &header, bam::Format::Bam)?;
            log::info!("write {} records to {}", i, format!("{}_{}.bam", output_prefix, j));
            
            i = 0;
        }
    }
    
    Ok(())
}


