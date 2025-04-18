use std::error::Error;
use clap::Parser;
use rust_htslib::bam;

mod cluster;
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

    #[clap(short = 'o', long = "output")]
    /// Output file. If not provided, output will be printed to stdout.
    pub output: Option<std::path::PathBuf>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();

    println!("input: {:?}\nreference: {:?}\nannotation: {:?}",
    args.input_bam, args.reference_fasta, args.annotation_gff);

    run(args)?;

    Ok(())
}


pub fn run(args: Cli) -> Result<(), Box<dyn Error>> {
    let mut bam = bam::IndexedReader::from_path(&args.input_bam)?;
    let reference = read_reference(&args.reference_fasta)?;
    let annotation = read_annotation(&args.annotation_gff)?;

    let read_pairs = read_pair_generator(
        &mut bam,
        reference.id(),
        0,
        // TODO: make this a cli argument
        reference.seq().len().try_into()? // Whole genome for now
    )?;

    for pair in read_pairs {
        let variants = call_variants(pair, &reference, &annotation);

        if variants.len() == 0 {
            continue;
        }
        println!("{:?}", variants.nt_mutations());
    }

    Ok(())
}
