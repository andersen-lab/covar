use std::error::Error;
use clap::Parser;
use rust_htslib::bam;

mod mutation;
mod gene;
mod utils;
use utils::{read_reference, read_annotation, read_pair_generator, call_variants};

#[derive(Parser)]
pub struct Cli {
    #[arg(short = 'i', long = "input")] // Add stdin support?
    /// Input BAM file (must be primer trimmed, sorted and indexed).
    pub input_bam: std::path::PathBuf,

    #[arg(short = 'r', long = "reference")]
    /// Reference genome FASTA file.
    pub reference_fasta: std::path::PathBuf,

    #[arg(short = 'a', long = "annotation")]
    /// Annotation GFF3 file. Used for translating mutations to respective amino acid mutation.
    pub annotation_gff: std::path::PathBuf,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();

    println!("input: {:?}\nreference: {:?}\nannotation: {:?}",
    args.input_bam, args.reference_fasta, args.annotation_gff);

    run(args)?;

    Ok(())
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
        let variants = call_variants(pair, &reference, &annotation);

        for variant in variants {
            println!("{:?}", variant);
        }
    }

    Ok(())
}
