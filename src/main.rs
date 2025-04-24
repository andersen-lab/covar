mod cluster;
mod gene;
mod mutation;
mod utils;

use std::error::Error;
use clap::Parser;
use rust_htslib::bam;

use cluster::{Cluster, call_variants};
use utils::{read_reference, read_annotation, read_pair_generator};

#[derive(Parser, Debug)]
#[command(name = "coVar")]
#[command(version, about)]
struct Cli {
    #[arg(short = 'i', long = "input")] // Add stdin support?
    /// Input BAM file (must be primer trimmed, sorted and indexed).
    pub input_bam: std::path::PathBuf,

    #[arg(short = 'r', long = "reference")]
    /// Reference genome FASTA file.
    pub reference_fasta: std::path::PathBuf,

    #[arg(short = 'a', long = "annotation")]
    /// Annotation GFF3 file. Used for translating mutations to respective amino acid mutation.
    pub annotation_gff: std::path::PathBuf,

    #[arg(short = 'o', long = "output")]
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


fn run(args: Cli) -> Result<(), Box<dyn Error>> {
    let mut bam = bam::IndexedReader::from_path(&args.input_bam)?;
    let reference = read_reference(&args.reference_fasta)?;
    let annotation = read_annotation(&args.annotation_gff)?;

    // Get read pairs
    let read_pairs = read_pair_generator(
        &mut bam,
        reference.id(),
        // 21563,
        // 25384,
        0,
        reference.seq().len().try_into()? // Whole genome
    );

    println!("Done fetching read pairs");

    // Call variants
    let mut clusters = Vec::<Cluster>::new();
    for pair in read_pairs {
        let variants = call_variants(pair, &reference, &annotation);
        //if variants.len() == 0 { continue; }
        clusters.push(variants);
    }

    println!("Done calling variants");

    // Aggregate unique clusters
    let clusters_merged = cluster::merge_clusters(&clusters); //57min

    println!("Done merging clusters");

    let mut output = String::new();

    if let Some(output_path) = args.output {
        std::fs::write(output_path, output)?;
    } else {
        print!("{}", output);
    }
    
    Ok(())
}
