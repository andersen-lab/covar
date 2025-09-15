use std::collections::HashMap;

use crate::cluster::Cluster;
use crate::mutation::{deletion::Deletion, insertion::Insertion, snp::SNP, Mutation};

use bio::io::fasta;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

pub fn call_variants(
    read_pair: &(Option<Record>, Option<Record>),
    reference: &fasta::Record,
    annotation: &HashMap<(u32, u32), String>,
    coverage_map: &[(u32, u32)],
    min_quality: u8,
) -> Cluster { 
    let (read1, read2) = read_pair;

    // Function to process a single read in the pair
    let process_read = |read: &Record| -> Vec<Mutation> {
        let mut local_variants: Vec<Mutation> = Vec::new();

        let ref_seq = reference.seq();

        let read_seq = String::from_utf8(read.seq().as_bytes().to_vec()).expect("Invalid UTF-8 in read sequence"); // panics
        let read_qual = read.qual();
        let cigar = read.cigar();

        let start_pos: u32 = read.pos() as u32;
        let mut read_pos: u32 = 0;
        let mut ref_pos: u32 = start_pos;

        for c in cigar.iter() {
            match c {
                Cigar::Match(len) => { // Call SNPs

                    for match_idx  in 0..*len {

                        let ref_base = ref_seq[(ref_pos + match_idx) as usize] as char;
                        let read_base = read_seq.chars().nth((read_pos + match_idx) as usize).unwrap();
                        
                        if ref_base != read_base {
                            let mut snp = SNP::new(
                                ref_pos + match_idx,
                                ref_base,
                                 read_base,
                                  read_qual[(read_pos + match_idx) as usize]
                            );
                            if let Some(gene) = snp.get_gene(annotation) {
                                snp.translate(&read_seq, read_pos + match_idx, read_qual, reference, &gene);
                                local_variants.push(Mutation::SNP(snp));        
                            }
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                },
                Cigar::Ins(len) => { // Call insertions
                    let ref_base = ref_seq[(ref_pos - 1) as usize] as char;
                    let ins_seq = read_seq[(read_pos as usize)..(read_pos as usize + *len as usize)].to_string();
                    let mut insertion = Insertion::new(
                        ref_pos - 1,
                        ref_base,
                         ins_seq,
                         *read_qual[read_pos as usize..read_pos as usize + (*len as usize)].iter().min().unwrap_or(&0)
                    );

                    if let Some(gene) = insertion.get_gene(annotation) {
                        insertion.translate(&gene);
                        local_variants.push(Mutation::Insertion(insertion));                                    
                    }
                    
                    read_pos += len;
                },
                Cigar::Del(len) => { // Call deletions
                    let del_qual = if (read_pos as usize + 1) > read_qual.len() {
                        read_qual[read_pos as usize]
                    } else {
                        *read_qual[read_pos as usize..read_pos as usize + 1].iter().min().unwrap_or(&0)
                    };
                    
                    let deletion_site = ref_pos - 1;

                    let ref_base = ref_seq[deletion_site as usize] as char;
                    let del_seq = std::str::from_utf8(&ref_seq[(ref_pos as usize)..(ref_pos as usize + *len as usize)])
                        .expect("Invalid UTF-8 sequence in reference")
                        .to_string();
                    let mut deletion = Deletion::new(
                        deletion_site,
                        ref_base,
                        del_seq,
                        del_qual
                    );
                    if let Some(gene) = deletion.get_gene(annotation) {
                        deletion.translate(&gene);
                        local_variants.push(Mutation::Deletion(deletion));                                    
                    }
                    
                    ref_pos += len;
                },
                Cigar::SoftClip(len) => {
                    read_pos += len;
                },
                _ => {},
            }
        }
        local_variants
    };

    let mut variants: Vec<Mutation> = Vec::new();
    let mut range: (u32, u32) = (0, u32::MAX); // Consider using htslib range type here

    // Process read1 if present
    if let Some(r1) = read1 {
        let r1_variants = process_read(r1);
        variants.extend(r1_variants);

        let start_pos = r1.pos() as u32;
        let end_pos = r1.cigar().end_pos() as u32;
        if range.0 == 0 && range.1 == u32::MAX {
            range.0 = start_pos;
            range.1 = end_pos;
        } else {
            if start_pos < range.0 {
                range.0 = start_pos;
            }
            if end_pos > range.1 {
                range.1 = end_pos;
            }
        }
    }
    
    // Process read2 if present
    if let Some(r2) = read2 {
        let r2_variants = process_read(r2);
        variants.extend(r2_variants);
        let start_pos = r2.pos() as u32;
        let end_pos = r2.cigar().end_pos() as u32;

        if range.0 == 0 && range.1 == u32::MAX {
            range.0 = start_pos;
            range.1 = end_pos;
        } else {
            if start_pos < range.0 {
                range.0 = start_pos;
            }
            if end_pos > range.1 {
                range.1 = end_pos;
            }
        }
    }
    
    // Basic deduplication - keep unique variants only
    let mut unique_variants = Vec::new();
    let mut seen_variants = std::collections::HashSet::new();
    
    for var in variants {
        let var_str = var.to_string();
        if !seen_variants.contains(&var_str) && !var_str.is_empty() {
            seen_variants.insert(var_str);
            unique_variants.push(var);
        }
    }

    // Filter by quality
    unique_variants = filter_quality(unique_variants, min_quality);

    // Sort unique variants by mut position
    unique_variants.sort_by(|a, b| {
        let a_pos = a.get_position();
        let b_pos = b.get_position();
        a_pos.cmp(&b_pos)
    });

    // Separate nt and aa mutations and parse as string
    let mut nt_mutations: Vec<String> = Vec::new();
    let mut aa_mutations: Vec<String> = Vec::new();

    let mut mutations_start = 0;
    let mut mutations_end = u32::MAX;

    for mutation in unique_variants {
        nt_mutations.push(mutation.to_string());
        aa_mutations.push(match mutation {
            Mutation::SNP(_) => mutation.get_aa_mutation().unwrap_or("Unknown".to_string()),
            _ => mutation.get_aa_mutation().unwrap_or("NA".to_string()),
        });
        let pos = mutation.get_position();
        if mutations_start == 0 || pos < mutations_start {
            mutations_start = pos;
        }
        if mutations_end == u32::MAX || pos > mutations_end {
            mutations_end = pos;
        }
    }

    let max_count = get_max_count(
        (mutations_start, mutations_end),
        coverage_map
    );

    Cluster::new(
        nt_mutations.join(" "),
        aa_mutations.join(" "),
        max_count,
        range.0, 
        range.1,
        mutations_start,
        mutations_end,
    )
}

fn get_max_count(mut_range: (u32, u32), coverages: &[(u32, u32)]) -> u32 {
    let (mut_start, mut_end) = mut_range;

    // Binary search for the first coverage where cov_start > mut_start
    let left = coverages.partition_point(|&(cov_start, _)| cov_start <= mut_start);

    coverages[..left]
        .iter()
        .filter(|&&(_, cov_end)| cov_end >= mut_end)
        .count() as u32
}

fn filter_quality(variants: Vec<Mutation>, min_quality: u8) -> Vec<Mutation> {
    variants.into_iter()
        .filter(|var| var.get_quality() >= min_quality)
        .collect()
}