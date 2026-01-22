use std::fmt;
use std::collections::HashSet;
use bio::io::fasta;

use crate::mutation::Mutation;


#[derive(Clone)]
pub struct Cluster {
    pub mutations: Vec<Mutation>,
    pub nt_mutations: String,
    pub cluster_depth: u32,
    pub total_depth: u32,
    pub coverage_start: u32,
    pub coverage_end: u32,
    pub mutations_start: u32,
    pub mutations_end: u32,
}

impl Cluster {
    pub fn new(mutations: Vec<Mutation>, nt_mutations: String, max_count: u32, coverage_start: u32, coverage_end: u32, mutations_start: u32, mutations_end: u32) -> Self {
        Self {
            mutations,
            nt_mutations,
            cluster_depth: 1,
            total_depth: max_count,
            coverage_start,
            coverage_end,
            mutations_start,
            mutations_end,
        }
    }
    pub fn translate_cluster(&self, reference: &fasta::Record) -> String {

        // create mock read sequence with mutations applied
        let mut read_seq = String::from_utf8(reference.seq().to_vec())
            .expect("Invalid UTF-8 in reference sequence");
        let start = self.coverage_start;
        let end = self.coverage_end;

        for mutation in &self.mutations {
            match mutation {
                Mutation::SNP(snp) => {
                    let pos = snp.get_position() as usize;
                    read_seq.replace_range(pos..pos+1, &snp.get_alternate_base());
                },
                Mutation::Insertion(ins) => {
                    let pos = ins.get_position() as usize;
                    let seq = ins.get_alternate_base();
                    read_seq.insert_str(pos, &seq);
                },
                _ => {},
            }
        }

        let mut aa_mutations: Vec<String> = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();

        for mutation in &self.mutations {
            let aa = match mutation {
                Mutation::SNP(snp) => snp.translate(&read_seq, reference, start, end).unwrap_or("Unknown".to_string()),
                Mutation::Deletion(del) => del.translate().unwrap_or("NA".to_string()),
                Mutation::Insertion(ins) => ins.translate().unwrap_or("NA".to_string()),
            };
            if aa == "NA" || seen.insert(aa.clone()) {
                aa_mutations.push(aa);
            }
        }
        aa_mutations.join(" ")
    }
}

impl fmt::Display for Cluster {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
        "{}\t{}\t{}\t{}",
        self.nt_mutations,
        self.cluster_depth,
        self.coverage_start,
        self.coverage_end)
    }
}