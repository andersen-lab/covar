# coVar
coVar is a tool for detecting physically-linked mutations in genomic data. Given a sorted, indexed BAM file, reference genome and gene annotation, coVar identifies and counts sequencing reads with unique physically linked mutations.

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

### Required arguments
| Flag                                  | Description                                                              |
| ------------------------------------- | ------------------------------------------------------------------------ |
| `-i`, `--input <INPUT_BAM>`           | Input BAM file (**must** be primer trimmed, sorted, and indexed).        |
| `-r`, `--reference <REFERENCE_FASTA>` | Reference genome in FASTA format.                                        |
| `-a`, `--annotation <ANNOTATION_GFF>` | Annotation GFF3 file for translating nucleotide to amino acid mutations. |

### Optional arguments
| Flag                           | Default            | Description                                                                               |
| ------------------------------ | ------------------ | ----------------------------------------------------------------------------------------- |
| `-o`, `--output <OUTPUT>`      | *stdout*           | Output file path. If not provided, results will be printed to stdout.                     |
| `-s`, `--start_site <START>`   | `0`                | Genomic start position for variant calling.                                               |
| `-e`, `--end_site <END>`       | *reference length* | Genomic end position for variant calling. Defaults to the length of the reference genome. |
| `-d`, `--min_depth <DEPTH>`    | `1`                | Minimum coverage depth for a mutation cluster to be considered.                           |
| `-f`, `--min_frequency <FREQ>` | `0.001`            | Minimum mutation frequency (cluster depth / total depth).                                 |
| `-q`, `--min_quality <QUAL>`   | `20`               | Minimum base quality score for variant calling.                                           |
| `-t`, `--threads <THREADS>`    | `1`                | Number of threads to use for processing.                                                  |

## Example Commands
### Basic run
```bash
coVar \
  -i sample.bam \
  -r reference.fasta \
  -a annotation.gff3
```

### Specify genomic region and output file
```bash
coVar \
  -i sample.bam \
  -r reference.fasta \
  -a annotation.gff3 \
  -s 1000 \
  -e 5000 \
  -o output.tsv
```

### Multi-threaded run with custom depth, quality and frequency thresholds
```bash
coVar \
  -i sample.bam \
  -r reference.fasta \
  -a annotation.gff3 \
  -d 5 \
  -q 30 \
  -f 0.01 \
  -t 4
```

## Output
The output is a tab-delimited file (.tsv) with the following columns:
| Column           | Description                                                 |
| ---------------- | ------------------------------------------------------------|
| `nt_mutations`   | Nucleotide mutations for this cluster                       |
| `aa_mutations`   | Corresponding amino acid translations (where possible*)     |
| `cluster_depth`  | Total number of read pairs with this cluster of mutations   |
| `total_depth`    | Total number of reads spanning this cluster                 |
| `frequency`      | Mutation frequency (cluster depth / total depth)            |
| `coverage_start` | Maximum read start site for which this cluster was detected |
| `coverage_end`   | Minimum read end site for which this cluster was detected   |

\*Note: Not all nucleotide mutations will have a corresponding amino acid mutations. For example, SNPs in codons that span reads or framshift indels will be translated as 'Unknown' and 'NA', respectively.