
use std::collections::HashMap;
use std::error::Error;
use std::io::ErrorKind;

use rust_htslib::bam::{IndexedReader, Read, Record};

pub fn read_pair_generator(
    bam: & mut IndexedReader,
    refname: &str,
    min_site: i64,
    max_site: i64,
) -> Result<Vec<(Option<Record>, Option<Record>)>, Box<dyn Error>> {

    let tid = match bam.header().tid(refname.as_bytes())
            .ok_or_else(|| std::io::Error::new(
                ErrorKind::NotFound, 
                format!("Reference '{}' not found", refname)
            )) {
        Ok(t) => t,
        Err(e) => return Err(Box::new(e) as Box<dyn Error>),
    };

    let _ = bam.fetch((tid, min_site as i64, max_site as i64 + 1))
        .map_err(|e| Box::new(e) as Box<dyn Error>);

    // Collect read pairs by query name
    let mut read_pairs: HashMap<Vec<u8>, (Option<Record>, Option<Record>)> = HashMap::new();
    for record_result in bam.records() {
        let record = match record_result {
            Ok(r) => r,
            Err(e) => return Err(Box::new(e) as Box<dyn Error>),
        };

        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }

        let query_name = record.qname().to_owned();
        
        let entry = read_pairs.entry(query_name).or_insert((None, None));
        
        if record.is_first_in_template() {
            entry.0 = Some(record);
        } else if record.is_last_in_template() {
            entry.1 = Some(record);
        } else {
            // Handle singleton reads if needed
            if entry.0.is_none() {
                entry.0 = Some(record);
            } else {
                entry.1 = Some(record);
            }
        }
    }

    // Convert the HashMap into a Vec and return
    Ok(read_pairs.into_iter().map(|(_, pair)| pair).collect())
}


