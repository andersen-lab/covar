use std::collections::HashMap;
use std::error::Error;
use std::io::ErrorKind;
use std::path::PathBuf;

use bio::io::{fasta, gff};
use bio::io::fasta::FastaRead;

use rust_htslib::bam::{IndexedReader, Read, Record};


pub fn read_reference(path: &PathBuf) -> Result<fasta::Record, Box<dyn Error>> {
    let mut reader = fasta::Reader::from_file(path)?;
    let mut reference = fasta::Record::new();
    reader.read(&mut reference)?;
    Ok(reference)
}

pub fn read_annotation(path: &PathBuf) -> Result<HashMap<(u32, u32), String>, Box<dyn Error>> {
    let mut reader = gff::Reader::from_file(path, gff::GffType::GFF3)?;
    let mut gene_regions: HashMap<(u32, u32), String> = HashMap::new();

    for record in reader.records() {
        let rec = record?;
        if rec.feature_type() == "CDS" {
            if let Some(gene) = rec.attributes().get("gene") {
                gene_regions.insert((*rec.start() as u32, *rec.end() as u32), gene.to_string());
            }
        }
    }
    Ok(gene_regions)
}

pub fn read_pair_generator(
    bam: & mut IndexedReader,
    refname: &str,
    min_site: u32,
    max_site: u32,) -> Vec<(Option<Record>, Option<Record>)> {

    let tid = match bam.header().tid(refname.as_bytes())
            .ok_or_else(|| std::io::Error::new(
                ErrorKind::NotFound, 
                format!("Reference '{}' not found", refname)
            )) {
        Ok(t) => t,
        Err(e) => panic!("Error: {}", e),
    };

    let _ = bam.fetch((tid, min_site , max_site + 1))
        .map_err(|e| Box::new(e) as Box<dyn Error>);

    // Get read pairs by query name
    let mut read_pairs: HashMap<Vec<u8>, (Option<Record>, Option<Record>)> = HashMap::new();
    for record_result in bam.records() {
        let record = match record_result {
            Ok(r) => r,
            Err(e) => panic!("Error reading BAM record: {}", e),
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
            // Handle singleton reads
            if entry.0.is_none() {
                entry.0 = Some(record);
            } else {
                entry.1 = Some(record);
            }
        }
    }

    read_pairs.into_values().collect()
}

pub fn get_coverage_map(read_pairs: &[(Option<Record>, Option<Record>)]) -> Vec<(u32, u32)> {
    let mut coverage_map = Vec::new();

    for (read1, read2) in read_pairs.iter() {
        let coverage_range = match (read1, read2) {
            (Some(r1), Some(r2)) => {
                let positions = [r1.pos(), r1.cigar().end_pos(), r2.pos(), r2.cigar().end_pos()];
                let min_pos = *positions.iter().min().unwrap();
                let max_pos = *positions.iter().max().unwrap();
                (Some(min_pos), Some(max_pos))
            },
            (Some(r1), None) => {
                let positions = [r1.pos(), r1.cigar().end_pos()];
                let min_pos = *positions.iter().min().unwrap();
                let max_pos = *positions.iter().max().unwrap();
                (Some(min_pos), Some(max_pos))
            },
            (None, Some(r2)) => {
                let positions = [r2.pos(), r2.cigar().end_pos()];
                let min_pos = *positions.iter().min().unwrap();
                let max_pos = *positions.iter().max().unwrap();
                (Some(min_pos), Some(max_pos))
            },
            (None, None) => continue,
        };
        // unwrap options and convert to u32
        let coverage_range = match (coverage_range.0, coverage_range.1) {
            (Some(start), Some(end)) => (start as u32, end as u32),
            _ => continue, // Skip if either start or end is None
        };
        coverage_map.push(coverage_range);
    }
    // Sort the coverage map by start position
    coverage_map.sort_by_key(|&(start, _)| start);

    coverage_map
}
