# coVar
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/covar?color=green)](https://anaconda.org/bioconda/covar)
[![crates.io version](https://img.shields.io/crates/v/covar)](https://crates.io/crates/covar)
[![crates.io downloads](https://img.shields.io/crates/d/covar?color=orange&label=downloads)](https://crates.io/crates/covar)


coVar is a tool for detecting physically-linked mutations in wastewater genomic sequencing data. Given a sorted, indexed BAM file, reference genome and gene annotation, coVar identifies and counts sequencing reads with unique physically linked mutations.

coVar also ships with [Freyja](https://github.com/andersen-lab/freyja), and can be accessed via the `freyja covariants` subcommand.

## Installation

### Bioconda (recommended)
```
conda create -n covar-env
conda activate covar-env
conda install -c bioconda covar
covar --version
```

### crates.io (recommended)
```
cargo install covar
covar --version
```

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
covar \
  -i sample.bam \
  -r reference.fasta \
  -a annotation.gff3
```

### Specify genomic region and output file
```bash
covar \
  -i sample.bam \
  -r reference.fasta \
  -a annotation.gff3 \
  -s 1000 \
  -e 5000 \
  -o output.tsv
```

### Multi-threaded run with custom depth, quality and frequency thresholds
```bash
covar \
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
| `aa_mutations`   | Corresponding amino acid translations (where possible)      |
| `cluster_depth`  | Total number of read pairs with this cluster of mutations   |
| `total_depth`    | Total number of reads spanning this cluster                 |
| `frequency`      | Mutation frequency (cluster depth / total depth)            |
| `coverage_start` | Maximum read start site for which this cluster was detected |
| `coverage_end`   | Minimum read end site for which this cluster was detected   |

**Note**: Not all nucleotide mutations will have a corresponding amino acid mutations. For example, SNPs in codons that span reads or frameshift indels will be translated as 'Unknown' and 'NA', respectively.

Additionally, in the case of multiple nucleotide mutations in a single codon, the amino acid translation will occur multiple times. For example, if a codon has two SNPs, the amino acid translation will be repeated for each SNP. See the following example:

```
nt_mutations:
G22813T T22882G G22895C T22896C G22898A A22948C A22948-CCT
aa_mutations:
S:K417N S:N440K S:V445P S:V445P S:G446S S:K462N S:DEL463
```

Here, the SNPS G22895C and T22896C are in the same codon, which result in Proline (P) at position 445. In order to keep a one-to-one relationship between nucleotide and amino acid mutations, the amino acid translation will be repeated so we see S:V445P twice.
