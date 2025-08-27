mod cli;
mod cluster;
mod gene;
mod mutation;
mod utils;

use std::sync::Arc;
use std::{thread, vec, error::Error, fs::File};
use clap::Parser;
use crossbeam::channel;
use polars::prelude::*;
use rust_htslib::bam;

use cli::{Cli, Config};
use cluster::call_variants;
use utils::{read_reference, read_annotation, read_pair_generator, get_coverage_map};


fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    let config = Config::from_cli(args)?;

    eprintln!("input: {:?}\nreference: {:?}\nannotation: {:?}",
        config.input_bam, config.reference_fasta, config.annotation_gff);

    run(config)?;

    Ok(())
}

fn run(config: Config) -> Result<(), Box<dyn Error>> {
    let mut bam = bam::IndexedReader::from_path(&config.input_bam)?;
    let reference = read_reference(&config.reference_fasta)?;
    let annotation = read_annotation(&config.annotation_gff)?;

    // Get read pairs
    // TODO: stream read pairs instead of loading all into memory 
    let read_pairs = read_pair_generator(
        &mut bam,
        reference.id(),
        config.start_site,
        config.end_site // Whole genome
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

    for _ in 0..config.threads {
        let receiver = receiver.clone();
        let reference = Arc::clone(&reference);
        let annotation = Arc::clone(&annotation);
        let coverage_map = Arc::clone(&coverage_map);
        let pb = Arc::clone(&pb);

        let handle = thread::spawn(move || {
            let mut local_clusters = Vec::new();
            while let Ok(pair) = receiver.recv() {
                let variants = call_variants(&pair, &reference, &annotation, &coverage_map, config.min_quality);
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
    let mut clusters_merged = cluster::merge_clusters(&clusters, &config)?;

    if let Some(output_path) = config.output { // Write to file if provided
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
