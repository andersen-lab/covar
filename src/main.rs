mod cluster;
mod gene;
mod mutation;
mod utils;

use std::sync::{Arc, Mutex};
use std::{thread, vec};
use std::{error::Error, fs::File};
use clap::Parser;
use indicatif::ProgressBar;
use crossbeam::channel;
use rust_htslib::bam;
use polars::prelude::*;

use cluster::{Cluster, call_variants};
use utils::{read_reference, read_annotation, read_pair_generator, get_coverage_map};

#[derive(Parser, Debug)]
#[command(name = "coVar", version, about)]
struct Cli {
    #[arg(short = 'i', long = "input")] // Add stdin support
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

    #[arg(short = 'c', long = "min_count", default_value_t = 1)]
    /// Minimum occurrences to include a cluster in output. Default is 1.
    pub min_count: u32,

    #[arg(short = 't', long = "threads", default_value_t = 1)]
    /// Number of threads to spawn for variant calling. Default is 1.
    pub threads: u32,
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
    let coverage_map = get_coverage_map(&read_pairs);

    // Call variants
    let pb = Arc::new(indicatif::ProgressBar::new(read_pairs.len() as u64));
    let (sender, receiver) = channel::unbounded();
    for pair in read_pairs {
        sender.send(pair).unwrap();
    }
    drop(sender); // Close the channel

    let mut handles = vec![];
    let reference = Arc::new(reference);
    let annotation = Arc::new(annotation);
    let coverage_map = Arc::new(coverage_map);

    for _ in 0..args.threads {
        let receiver = receiver.clone();
        let reference = Arc::clone(&reference);
        let annotation = Arc::clone(&annotation);
        let coverage_map = Arc::clone(&coverage_map);
        let pb = Arc::clone(&pb);

        let handle = thread::spawn(move || {
            let mut local_clusters = Vec::new();
            while let Ok(pair) = receiver.recv() {
                let variants = call_variants(&pair, &reference, &annotation, &coverage_map);
                local_clusters.push(variants);
                pb.inc(1);
            }
            local_clusters
        });
        handles.push(handle);
    }

    // Collect results from threads
    let mut clusters = Vec::new();
    for handle in handles {
        let mut local_clusters = handle.join().unwrap();
        clusters.append(&mut local_clusters);
    }
    pb.finish_and_clear();

    // Aggregate unique clusters
    let mut clusters_merged = cluster::merge_clusters(&clusters, &args)?;

    if let Some(output_path) = args.output { // Write to file if provided
        let file = File::create(output_path)?;
        let mut writer = CsvWriter::new(file)
            .with_separator(b'\t');
        writer.finish(&mut clusters_merged)?;
    } else { // Write to stdout
        let mut buffer = Vec::new();
        let mut writer = CsvWriter::new(&mut buffer)
            .with_separator(b'\t');
        writer.finish(&mut clusters_merged)?;
        println!("{}", String::from_utf8(buffer)?);
    }
    
    Ok(())
}
