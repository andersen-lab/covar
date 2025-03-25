use std::{error::Error, fs, io::BufReader};

use rust_htslib::{bam, bam::Read};
use clap::Parser;

#[derive(Parser)]
pub struct Cli {
    #[arg(short = 'i', long = "input")]
    pub input_bam: std::path::PathBuf,

    #[arg(short = 'r', long = "reference")]
    pub reference_fasta: std::path::PathBuf,

    #[arg(short = 'a', long = "annotation")]
    pub annotation_gff: std::path::PathBuf,
}

pub fn run(args: Cli) -> Result<(), Box<dyn Error>> {
    let input_bam = bam::Reader::from_path(&args.input_bam)
        .map_err(|e| format!("Failed to open input input BAM file: {}", e))?;
    
    let _reference_fasta_reader = BufReader::new(fs::File::open(&args.reference_fasta)
        .map_err(|e| format!("Failed to open reference FASTA file: {}", e))?);
    
    let _annotation_gff_reader = BufReader::new(fs::File::open(&args.annotation_gff)
        .map_err(|e| format!("Failed to open annotation GFF file: {}", e))?);
 
    let header = bam::Header::from_template(input_bam.header());
    
    // print header records to the terminal, akin to samtools
    for (key, records) in header.to_hashmap() {
        for record in records {
                println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
        }
    }
    
    println!("Successfully opened all files");
    
    // Process the files here...
    
    Ok(())
}
