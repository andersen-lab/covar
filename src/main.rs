use clap::Parser;
use covar::{Cli, run};

fn main() {
    let args = Cli::parse();

    println!("input: {:?}, reference: {:?}, annotation: {:?}",
    args.input_bam, args.reference_fasta, args.annotation_gff);


    run(args);
}
