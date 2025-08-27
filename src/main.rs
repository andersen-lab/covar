mod cli;
mod cluster;
mod gene;
mod mutation;
mod utils;

use std::sync::Arc;
use std::{error::Error, fs::File};
use clap::Parser;
use polars::prelude::*;
use rust_htslib::bam;
use rayon::prelude::*;

use cli::{Cli, Config};
use cluster::call_variants;
use utils::{read_reference, read_annotation, read_pair_generator, get_coverage_map};


fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    let config = Config::from_cli(args)?;

    eprintln!(
        "input: {:?}\nreference: {:?}\nannotation: {:?}",
        config.input_bam,
        config.reference_fasta, 
        config.annotation_gff
    );

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
        config.end_site
    );
    let coverage_map = get_coverage_map(&read_pairs);

    // Call variants in parallel
    rayon::ThreadPoolBuilder::new().num_threads(config.threads).build_global().unwrap();

    let pb = Arc::new(indicatif::ProgressBar::new(read_pairs.len() as u64));
    let reference = Arc::new(reference);
    let annotation = Arc::new(annotation);
    let coverage_map = Arc::new(coverage_map);

    let clusters: Vec<_> = read_pairs
        .into_par_iter()
        .map(|pair| {
            let variants = call_variants(&pair, &reference, &annotation, &coverage_map, config.min_quality);
            pb.inc(1);
            variants
        })
        .collect();

    pb.finish_and_clear();

    // Aggregate unique clusters
    let mut clusters_merged = cluster::merge_clusters(&clusters, &config)?;
    
    // Write to file if provided, else to stdout
    if let Some(output_path) = config.output { 
        let file = File::create(output_path)?;
        let mut writer = CsvWriter::new(file)
            .with_separator(b'\t');
        writer.finish(&mut clusters_merged)?;
    } else {
        let mut buffer = Vec::new();
        let mut writer = CsvWriter::new(&mut buffer)
            .with_separator(b'\t');
        writer.finish(&mut clusters_merged)?;
        println!("{}", String::from_utf8(buffer)?);
    }
    
    Ok(())
}
