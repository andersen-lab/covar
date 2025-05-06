mod cluster;
mod gene;
mod mutation;
mod utils;

use std::{error::Error, fs::File};
use clap::Parser;
use rust_htslib::bam;
use polars::prelude::*;

use cluster::{Cluster, call_variants};
use utils::{read_reference, read_annotation, read_pair_generator, get_coverage_map};

#[derive(Parser, Debug)]
#[command(name = "coVar", version, about)]
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
    /// Optional output file path. If not provided, output will be printed to stdout.
    pub output: Option<std::path::PathBuf>,

    #[arg(short = 's', long = "start_site", default_value_t = 0)]
    /// Genomic start site for variant calling. Default is 0.
    pub start_site: u32,

    #[arg(short = 'e', long = "end_site")]
    /// Genomic end site for variant calling. Default is the length of the reference genome.
    pub end_site: Option<u32>,

    #[arg(short = 'c', long = "min-count", default_value_t = 1)]
    /// Minimum occurrences to include a cluster in output.
    pub min_count: u32,
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = Cli::parse();

    if args.end_site.is_none() {
        // Set end site to the length of the reference genome if not provided
        let reference = read_reference(&args.reference_fasta)?;
        args.end_site = Some(reference.seq().len() as u32);
    } else if let Some(end_site) = args.end_site {
        if end_site < args.start_site {
            return Err("End site must be greater than or equal to start site.".into());
        }
    }

    eprintln!("input: {:?}\nreference: {:?}\nannotation: {:?}",
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
        args.start_site,
        args.end_site.unwrap() // Whole genome
    );
        
    let coverage_map = get_coverage_map(&read_pairs, reference.seq().len() as u32);


    // Call variants
    eprintln!("Calling variants...");
    let pb = indicatif::ProgressBar::new(read_pairs.len() as u64);
    let mut clusters = Vec::<Cluster>::new();
    for pair in read_pairs {
        let variants = call_variants(pair, &reference, &annotation, &coverage_map);
        clusters.push(variants);
        pb.inc(1);
    }
    pb.finish_and_clear();

    // Aggregate unique clusters
    let mut clusters_merged = cluster::merge_clusters(&clusters, &args)?;

    if let Some(output_path) = args.output { // Write to file if provided
        let file = File::create(output_path)?;
        let mut writer = CsvWriter::new(file);
        writer.finish(&mut clusters_merged)?;
    } else { // Write to stdout
        let mut buffer = Vec::new();
        let mut writer = CsvWriter::new(&mut buffer);
        writer.finish(&mut clusters_merged)?;
        println!("{}", String::from_utf8(buffer)?);
    }
    
    Ok(())
}
