
use std::collections::HashMap;
use rust_htslib::bam::{IndexedReader, Read, Record};
use std::error::Error;
use std::io::ErrorKind;

pub fn read_pair_generator<'a>(
    bam: &'a mut IndexedReader,
    refname: &str,
    min_site: i64,
    max_site: i64,
) -> impl Iterator<Item = Result<(Option<Record>, Option<Record>), Box<dyn Error + 'a>>> + 'a {

    let tid = match bam.header().tid(refname.as_bytes())
            .ok_or_else(|| std::io::Error::new(
                ErrorKind::NotFound, 
                format!("Reference '{}' not found", refname)
            )) {
        Ok(t) => t,
        Err(e) => return std::iter::once(Err(Box::new(e) as Box<dyn Error>)),
    };

    // First pass: identify paired reads
    //let mut is_paired = HashMap::new();

    bam.fetch((tid, min_site, max_site + 1));
    for record in bam.records() {
        println!("{:?}", record.unwrap().cigar()); // Debugging line to check the query name
    }

    std::iter::once(Ok((None, None))) // Placeholder for the iterator
}


