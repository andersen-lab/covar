use std::collections::HashMap;
use std::error::Error;
use std::path::{self, PathBuf};

use bio::io::{fasta, gff};
use bio::io::fasta::FastaRead;
use rust_htslib::{bam, bam::Read};
use clap::Parser;

#[derive(Parser)]
pub struct Cli {
    #[arg(short = 'i', long = "input")]
    pub input_bam: std::path::PathBuf,

    #[arg(short = 'r', long = "reference")]
    pub reference_fasta: std::path::PathBuf,

    #[arg(short = 'a', long = "annotation")]
    pub annotation_gff: std::path::PathBuf,
}

fn read_reference(path: &PathBuf) -> Result<fasta::Record, Box<dyn Error>> {
    let mut reader = fasta::Reader::from_file(path)?;
    let mut reference = fasta::Record::new();
    reader.read(&mut reference)?;
    Ok(reference)
}

fn read_annotation(path: &PathBuf) -> Result<HashMap<String, (u64, u64)>, Box<dyn Error>> {
    let mut reader = gff::Reader::from_file(path, gff::GffType::GFF3)?;
    let mut gene_regions: HashMap<String, (u64, u64)> = HashMap::new();

    for record in reader.records() {
        let rec = record.ok().expect("Error reading GFF record");
        if rec.feature_type() == "CDS" {
            if let Some(gene) = rec.attributes().get("gene") {
                gene_regions.insert(gene.to_string(), (*rec.start(), *rec.end()));
            }
        }
    }

    Ok(gene_regions)
}

pub fn run(args: Cli) -> Result<(), Box<dyn Error>> {
    let reference = read_reference(&args.reference_fasta)?;
    let annotation = read_annotation(&args.annotation_gff)?;

    println!("{:?}", annotation.get("S"));

    let input_bam = bam::Reader::from_path(&args.input_bam)?;
    let _bam_header = bam::Header::from_template(input_bam.header());

    
    // Process the files here...
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_reference() {
        let path = PathBuf::from("src/assets/sars-cov-2/NC_045512_Hu-1.fasta");
        let reference = read_reference(&path).unwrap();
        assert_eq!(reference.seq().len(), 29903);
    }

    #[test]
    fn test_read_annotation() {
        let path = PathBuf::from("src/assets/sars-cov-2/NC_045512_Hu-1.gff");
        let annot = read_annotation(&path).unwrap();
        assert!(!annot.is_empty());
        assert_eq!(*annot.get("S").unwrap(), (21563, 25384));
    }
}

