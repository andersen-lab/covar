use std::collections::HashMap;
use std::error::Error;
use std::path::PathBuf;

use bio::io::{fasta, gff};
use bio::io::fasta::FastaRead;
use clap::Parser;
use rust_htslib::bam;

mod bam_utils;
use bam_utils::read_pair_generator;

mod mutation;
pub use mutation::{call_variants, translate_variants};


#[derive(Parser)]
pub struct Cli {
    #[arg(short = 'i', long = "input")]
    /// Input BAM file (must be sorted and indexed).
    pub input_bam: std::path::PathBuf,

    #[arg(short = 'r', long = "reference")]
    /// Reference genome FASTA file.
    pub reference_fasta: std::path::PathBuf,

    #[clap(short = 'a', long = "annotation")]
    /// Optional annotation GFF file. If provided, will output amino acid translations in addition to nucleotide variant calls.
    pub annotation_gff: std::path::PathBuf,
}

pub fn run(args: Cli) -> Result<(), Box<dyn Error>> {
    let reference = read_reference(&args.reference_fasta)?;
    let annotation = read_annotation(&args.annotation_gff)?;
    let mut bam = bam::IndexedReader::from_path(&args.input_bam)?;

    let read_pairs = read_pair_generator(
        &mut bam,
        reference.id(),
        0,
        reference.seq().len().try_into()? // Whole genome for now
    )?;

    for pair in read_pairs {
        let variants = call_variants(pair, &reference)
            .map_err(|e| {
                eprintln!("Error calling variants: {}", e);
                e
            })?;

        let aa_variants = translate_variants(variants, &reference, &annotation);

        for variant in aa_variants {
            println!("{:?}", variant);
        }
    }

    Ok(())
}

fn read_reference(path: &PathBuf) -> Result<fasta::Record, Box<dyn Error>> {
    let mut reader = fasta::Reader::from_file(path)?;
    let mut reference = fasta::Record::new();
    reader.read(&mut reference)?;
    Ok(reference)
}

fn read_annotation(path: &PathBuf) -> Result<HashMap<(u32, u32), String>, Box<dyn Error>> {
    let mut reader = gff::Reader::from_file(path, gff::GffType::GFF3)?;
    let mut gene_regions: HashMap<(u32, u32), String> = HashMap::new();

    for record in reader.records() {
        let rec = record.ok().expect("Error reading GFF record");
        if rec.feature_type() == "CDS" {
            if let Some(gene) = rec.attributes().get("gene") {
                gene_regions.insert((*rec.start() as u32, *rec.end() as u32), gene.to_string());
            }
        }
    }
    Ok(gene_regions)
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
        assert_eq!(*annot.get(&(21563, 25384)).unwrap(), "S");
    }
}

