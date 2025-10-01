use std::collections::HashMap;
use std::error::Error;
use std::io::ErrorKind;
use std::path::PathBuf;

use bio::io::{fasta, gff};
use bio::io::fasta::FastaRead;

use rust_htslib::bam::{IndexedReader, Read, Record};


pub fn read_reference<P: AsRef<std::path::Path>>(path: P) -> Result<fasta::Record, Box<dyn Error>> {
    let mut reader = fasta::Reader::from_file(path.as_ref())?;
    let mut reference = fasta::Record::new();
    reader.read(&mut reference)?;
    Ok(reference)
}

pub fn read_annotation(path: &PathBuf) -> Result<HashMap<(u32, u32), String>, Box<dyn Error>> {
    let mut reader = gff::Reader::from_file(path, gff::GffType::GFF3)?;
    let mut gene_regions: HashMap<(u32, u32), String> = HashMap::new();

    for record in reader.records() {
        let rec = record?;
        if rec.feature_type() == "gene" {
            if let Some(gene) = rec.attributes().get("Name") {
                if gene == "ORF1ab" { // SARS-CoV-2 specific
                    gene_regions.insert((266, 13468), "ORF1a".to_string());
                    gene_regions.insert((13468, 21555), "ORF1b".to_string());
                    continue;
                }
                gene_regions.insert((*rec.start() as u32, *rec.end() as u32), gene.to_string());
            }
        }
    }
    Ok(gene_regions)
}

pub fn read_pair_generator(
    bam: &mut IndexedReader,
    refname: &str,
    min_site: u32,
    max_site: u32,
) -> Vec<(Option<Record>, Option<Record>)> {
    let tid = match bam.header().tid(refname.as_bytes())
        .ok_or_else(|| std::io::Error::new(
            ErrorKind::NotFound,
            format!("Reference '{}' not found", refname)
        )) {
            Ok(t) => t,
            Err(e) => panic!("Error: {}", e),
        };

    let _ = bam.fetch((tid, min_site, max_site + 1))
        .map_err(|e| Box::new(e) as Box<dyn Error>);

    // Keep only unmatched reads in memory; when a mate arrives emit the pair immediately.
    let mut waiting: HashMap<Vec<u8>, Record> = HashMap::new();
    let mut results: Vec<(Option<Record>, Option<Record>)> = Vec::new();

    for record_result in bam.records() {
        let record = match record_result {
            Ok(r) => r,
            Err(e) => panic!("Error reading BAM record: {}", e),
        };

        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }

        let qname = record.qname().to_vec();

        if let Some(prev) = waiting.remove(&qname) {
            // We have a mate; decide ordering based on flags if possible.
            let pair = if prev.is_first_in_template() {
                (Some(prev), Some(record))
            } else if prev.is_last_in_template() {
                (Some(record), Some(prev))
            } else {
                // Fallback: treat the previously seen as first.
                (Some(prev), Some(record))
            };
            results.push(pair);
        } else {
            // No mate seen yet: store this read and wait for its mate.
            waiting.insert(qname, record);
        }
    }

    // Any remaining reads in `waiting` are singletons (mate not seen in the region).
    for (_k, rec) in waiting.into_iter() {
        if rec.is_first_in_template() {
            results.push((Some(rec), None));
        } else if rec.is_last_in_template() {
            results.push((None, Some(rec)));
        } else {
            // Unknown orientation: put into the first slot.
            results.push((Some(rec), None));
        }
    }

    results
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
    // Sort the coverage map by start position, then end position
    coverage_map.sort_by_key(|&(start, end)| (start, end));
    coverage_map
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_read_reference() {
        let reference = read_reference(&PathBuf::from("tests/data/NC_045512_Hu-1.fasta")).unwrap();

        assert_eq!(reference.id(), "NC_045512.2");
        assert_eq!(reference.seq().len(), 29903);
    }

    #[test]
    fn test_read_annotation() {
        let gene_regions = read_annotation(&PathBuf::from("tests/data/NC_045512_Hu-1.gff")).unwrap();

        assert_eq!(gene_regions.len(), 12);
        assert_eq!(gene_regions.get(&(21563, 25384)).unwrap(), "S");
        assert_eq!(gene_regions.get(&(266, 13468)).unwrap(), "ORF1a");
        assert_eq!(gene_regions.get(&(13468, 21555)).unwrap(), "ORF1b");
    }

    #[test]
    fn test_read_bam_and_coverage() {
        let mock_bam_path = PathBuf::from("tests/data/test.bam");
        let mut bam = IndexedReader::from_path(&mock_bam_path).unwrap();
        let read_pairs = read_pair_generator(&mut bam, "NC_045512.2", 23000, 24000);
        let coverage_map = get_coverage_map(&read_pairs);

        assert!(!read_pairs.is_empty());
        assert!(!coverage_map.is_empty());
    }
}   