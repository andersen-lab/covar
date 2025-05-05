# coVar
coVar is a tool for detecting physically-linked variants in genomic data. Given a sorted, indexed BAM file, reference genome and gene annotation, coVar identifies and counts reads with unique physically linked variants.

## Installation

Currently, to install coVar, you need to have [cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html) installed. On release we plan to publish a precompiled binary to crates.io and bioconda.

### Local build from source (experimental)
```
git clone https://github.com/andersen-lab/covar.git
cd covar
cargo install --path .
covar --version
```

## Usage
```
covar --input <INPUT_BAM> --reference <REFERENCE_FASTA> --annotation <ANNOTATION_GFF>
```