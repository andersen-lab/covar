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
        let mock_reference_str = b">chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let mock_reference_file = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(mock_reference_file.path(), mock_reference_str).unwrap();

        let reference = read_reference(&PathBuf::from(mock_reference_file.path())).unwrap();

        assert_eq!(reference.id(), "chr1");
        assert_eq!(reference.seq().len(), 80);
    }

    #[test]
    fn test_read_annotation() {
        let mock_gff_str = b"##gff-version 3\nchr1\t.\tCDS\t1\t1000\t.\t+\t0\tID=cds1;gene=gene1\nchr1\t.\tCDS\t2000\t3000\t.\t+\t0\tID=cds2;gene=gene2\nchr1\t.\tCDS\t4000\t5000\t.\t+\t0\tID=cds3;gene=ORF1ab\n";
        let mock_gff_file = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(mock_gff_file.path(), mock_gff_str).unwrap();

        let gene_regions = read_annotation(&PathBuf::from(mock_gff_file.path())).unwrap();

        assert_eq!(gene_regions.len(), 4);
        assert_eq!(gene_regions.get(&(1, 1000)).unwrap(), "gene1");
        assert_eq!(gene_regions.get(&(2000, 3000)).unwrap(), "gene2");
        assert_eq!(gene_regions.get(&(266, 13468)).unwrap(), "ORF1a");
        assert_eq!(gene_regions.get(&(13468, 21555)).unwrap(), "ORF1b");
    }

    #[test]
    fn test_read_bam_and_coverage() {
        let mock_bam_path = PathBuf::from("tests/data/mock.bam");
        let mut bam = IndexedReader::from_path(&mock_bam_path).unwrap();
        let read_pairs = read_pair_generator(&mut bam, "NC_045512.2", 23000, 24000);
        let coverage_map = get_coverage_map(&read_pairs);

        assert!(!read_pairs.is_empty());
        assert!(!coverage_map.is_empty());
    }
}   