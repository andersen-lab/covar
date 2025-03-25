use std::error::Error;

use bio::io::{fasta, gff};
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
    let input_bam = bam::Reader::from_path(&args.input_bam)?;
    
    let _reference_fasta_reader = fasta::Reader::from_file(&args.reference_fasta)?;
    
    let _annotation_gff_reader = gff::Reader::from_file(&args.annotation_gff, gff::GffType::GFF3)?;
 
    let bam_header = bam::Header::from_template(input_bam.header());
    
    // print header records to the terminal, akin to samtools
    for (key, records) in bam_header.to_hashmap() {
        for record in records {
                println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
        }
    }
    
    println!("Successfully opened all files");
    
    // Process the files here...
    
    Ok(())
}
