use std::fmt;
use std::collections::HashMap;

use bio::io::fasta;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

use crate::mutation::{Mutation, snp::SNP, insertion::Insertion, deletion::Deletion};

#[derive(Clone)]
pub struct Cluster {
    nt_mutations: Vec<String>,
    aa_mutations: Vec<String>,
    count: u32,
    max_count: u32,
    frequency: f32,
    coverage_start: u32,
    coverage_end: u32,
    mutations_start: u32,
    mutations_end: u32,
}

impl Cluster {
    pub fn new(nt_mutations: Vec<String>, aa_mutations: Vec<String>, coverage_start: u32, coverage_end: u32, mutations_start: u32, mutations_end: u32) -> Self {
        Self {
            nt_mutations,
            aa_mutations,
            count: 1,
            max_count: 1,
            frequency: 1.0,
            coverage_start,
            coverage_end,
            mutations_start,
            mutations_end,
        }
    }

    pub fn nt_mutations(&self) -> String {
        self.nt_mutations.join(",")
    }

    pub fn aa_mutations(&self) -> String {
        self.aa_mutations.join(",")
    }

    // pub fn count(&self) -> u32 {
    //     self.count
    // }

    // pub fn max_count(&self) -> u32 {
    //     self.max_count
    // }

    // pub fn start(&self) -> u32 {
    //     self.coverage_start
    // }

    // pub fn end(&self) -> u32 {
    //     self.coverage_end
    // }
}

impl fmt::Display for Cluster {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        self.nt_mutations(),
        self.aa_mutations(),
        self.count,
        self.max_count,
        self.frequency,
        self.coverage_start,
        self.coverage_end)
    }
}

pub fn call_variants(
    read_pair: (Option<Record>, Option<Record>),
    reference: &fasta::Record,
    annotation: &HashMap<(u32, u32), String>
) -> Cluster { 
    let (read1, read2) = read_pair;

    if read1.is_none() && read2.is_none() {
        panic!("Both reads are missing, cannot proceed with variant calling.")
    }
    
    // Function to process a single read in the pair
    let process_read = |read: &Record| -> Vec<(Mutation, String)> {
        let mut local_variants: Vec<(Mutation, String)> = Vec::new();

        let ref_seq = reference.seq();
        let ref_seq_len = ref_seq.len() as u32;

        let read_seq = String::from_utf8(read.seq().as_bytes().to_vec()).expect("Invalid UTF-8 in read sequence"); // panics
        let read_seq_len = read_seq.len() as u32;
        let cigar = read.cigar();

        let start_pos: u32 = read.pos() as u32;
        let mut read_pos: u32 = 0;
        let mut ref_pos: u32 = start_pos;

        for c in cigar.iter() { // TODO: handle indel offsets
            match c {
                Cigar::Match(len) => { // Call SNPs
                    for i in 0..*len {
                        if ref_pos + i < ref_seq_len && read_pos + i < read_seq_len {
                            let ref_base = ref_seq[(ref_pos + i) as usize] as char;
                            let read_base = read_seq.chars().nth((read_pos + i) as usize).unwrap();
                            
                            if ref_base != read_base && read_base != 'N' {
                                let snp = SNP::new(ref_pos + i, ref_base, read_base);
                                if let Some(gene) = snp.get_gene(annotation) {
                                    let aa_mutation = snp.translate(&read_seq, read_pos, reference, &gene)
                                        .unwrap_or_else(|| "Unknown".to_string());
                                    local_variants.push((Mutation::SNP(snp), aa_mutation));                                    
                                }
                            }
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                },
                Cigar::Ins(len) => { // Call insertions
                    if ref_pos > 0 && read_pos + *len <= read_seq_len {
                        let ref_base = ref_seq[(ref_pos - 1) as usize] as char;
                        let ins_seq = read_seq[(read_pos as usize)..(read_pos as usize + *len as usize)].to_string();
                        let insertion = Insertion::new(ref_pos - 1, ref_base, ins_seq);
                        if let Some(gene) = insertion.get_gene(annotation) {
                            let aa_mutation = insertion.translate(&read_seq, read_pos, reference, &gene)
                                .unwrap_or_else(|| "Unknown".to_string());
                            local_variants.push((Mutation::Insertion(insertion), aa_mutation));                                    
                        }
                    }
                    read_pos += len;
                },
                Cigar::Del(len) => { // Call deletions
                    if ref_pos > 0 && ref_pos + *len <= ref_seq_len {
                        let ref_base = ref_seq[(ref_pos - 1) as usize] as char;
                        let del_seq = std::str::from_utf8(&ref_seq[(ref_pos as usize)..(ref_pos as usize + *len as usize)])
                            .expect("Invalid UTF-8 sequence in reference")
                            .to_string();
                        let deletion = Deletion::new(ref_pos - 1, ref_base, del_seq);
                        if let Some(gene) = deletion.get_gene(annotation) {
                            let aa_mutation = deletion.translate(&read_seq, read_pos, reference, &gene)
                                .unwrap_or_else(|| "Unknown".to_string());
                            local_variants.push((Mutation::Deletion(deletion), aa_mutation));                                    
                        }
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

    let mut variants: Vec<(Mutation, String)> = Vec::new();
    let mut range: (u32, u32) = (0, u32::MAX); // Consider using htslib range type here

    // Process read1 if present
    if let Some(r1) = read1 {
        let r1_variants = process_read(&r1);
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
        let r2_variants = process_read(&r2);
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
    
    for (var, aa_mut) in variants {
        let var_str = var.to_string();
        if !seen_variants.contains(&var_str) {
            seen_variants.insert(var_str);
            unique_variants.push((var, aa_mut));
        }
    }

    // Sort unique variants by mut position
    unique_variants.sort_by(|a, b| {
        let a_pos = a.0.get_position();
        let b_pos = b.0.get_position();
        a_pos.cmp(&b_pos)
    });

    // Separate nt and aa mutations and parse as string
    let mut nt_mutations: Vec<String> = Vec::new();
    let mut aa_mutations: Vec<String> = Vec::new();
    let mut mutations_start = 0;
    let mut mutations_end = u32::MAX;

    for (mutation, aa_mut) in unique_variants {
        nt_mutations.push(mutation.to_string());
        aa_mutations.push(aa_mut);
        let pos = mutation.get_position();
        if mutations_start == 0 || pos < mutations_start {
            mutations_start = pos;
        }
        if mutations_end == u32::MAX || pos > mutations_end {
            mutations_end = pos;
        }
    }

    Cluster::new(
        nt_mutations,
        aa_mutations,
        range.0, 
        range.1,
        mutations_start,
        mutations_end,
    )
}


pub fn merge_clusters(clusters: &Vec<Cluster>) -> Vec<Cluster> {

    let mut merged_map: HashMap<String, Cluster> = HashMap::new();

    for cluster in clusters {
        let key = cluster.nt_mutations();
        if let Some(merged_cluster) = merged_map.get_mut(&key) {
            // Merge existing clusters
            merged_cluster.count += cluster.count;
            merged_cluster.coverage_start = merged_cluster.coverage_start.max(cluster.coverage_start);
            merged_cluster.coverage_end = merged_cluster.coverage_end.min(cluster.coverage_end);
        } else {
            // Add new cluster
            merged_map.insert(key.clone(), cluster.clone());
        }

        // Update max_count for read pairs that span the cluster
        for merged_cluster in merged_map.values_mut() {
            if cluster.coverage_start <= merged_cluster.mutations_start && cluster.coverage_end >= merged_cluster.mutations_end {
                merged_cluster.max_count += cluster.count;
            }
        }
    }

    // Update frequency
    for merged_cluster in merged_map.values_mut() {
        merged_cluster.frequency = (merged_cluster.count as f32) / (merged_cluster.max_count as f32);
    }

    merged_map.into_values().collect()
}
