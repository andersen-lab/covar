use std::error::Error;
use clap::Parser;
use crate::utils::read_reference;

#[derive(Parser, Debug)]
#[command(name = "coVar", version, about)]
pub struct Cli {
    #[arg(short = 'i', long = "input")]
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
    /// Genomic start site for variant calling.
    pub start_site: u32,

    #[arg(short = 'e', long = "end_site")]
    /// Genomic end site for variant calling. Defaults to the length of the reference genome.
    pub end_site: Option<u32>,

    #[arg(short = 'd', long = "min_depth", default_value_t = 1)]
    /// Minimum coverage depth for a mutation cluster to be considered.
    pub min_depth: u32,

    #[arg(short = 'f', long = "min_frequency", default_value_t = 0.001)]
    /// Minimum frequency (cluster depth / total depth) to include a cluster in output.
    pub min_frequency: f64,

    #[arg(short = 'q', long = "min_quality", default_value_t = 20)]
    /// Minimum base quality for variant calling.
    pub min_quality: u8,

    #[arg(short = 't', long = "threads", default_value_t = 1)]
    /// Number of threads to spawn for variant calling.
    pub threads: usize,
}

pub struct Config {
    pub input_bam: std::path::PathBuf,
    pub reference_fasta: std::path::PathBuf,
    pub annotation_gff: std::path::PathBuf,
    pub output: Option<std::path::PathBuf>,
    pub start_site: u32,
    pub end_site: u32,
    pub min_depth: u32,
    pub min_frequency: f64,
    pub min_quality: u8,
    pub threads: usize,
}

impl Config {
    pub fn from_cli(args: Cli) -> Result<Self, Box<dyn Error>> {
        let end_site = if let Some(end_site) = args.end_site {
            if end_site < args.start_site {
                return Err("End site must be greater than or equal to start site.".into());
            }
            end_site
        } else {
            let reference = read_reference(&args.reference_fasta)?;
            reference.seq().len() as u32
        };

        Ok(Config {
            input_bam: args.input_bam,
            reference_fasta: args.reference_fasta,
            annotation_gff: args.annotation_gff,
            output: args.output,
            start_site: args.start_site,
            end_site,
            min_depth: args.min_depth,
            min_frequency: args.min_frequency,
            min_quality: args.min_quality,
            threads: args.threads,
        })
    }
}