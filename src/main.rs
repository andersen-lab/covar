use std::process;
use clap::Parser;


#[derive(Parser)]
struct Cli {
    #[arg(short = 'i', long = "input")]
    input_bam: std::path::PathBuf,

    #[arg(short = 'r', long = "reference")]
    reference_fasta: std::path::PathBuf,

    #[arg(short = 'a', long = "annotation")]
    annotation_gff: std::path::PathBuf,
}

fn main() {
    let args = Cli::parse();

    println!("input: {:?}, reference: {:?}, annotation: {:?}",
    args.input_bam, args.reference_fasta, args.annotation_gff);

    
}
