use clap::Parser;
use covar::{Cli, run};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();

    println!("input: {:?}, reference: {:?}, annotation: {:?}",
    args.input_bam, args.reference_fasta, args.annotation_gff);


    run(args)?;

    Ok(())
}
