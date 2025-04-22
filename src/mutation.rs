pub mod snp;
pub mod insertion;
pub mod deletion;

use bio::io::fasta;
use std::{collections::HashMap, fmt};
use crate::gene::Gene;

pub use snp::SNP;
pub use insertion::Insertion;
pub use deletion::Deletion;

#[derive(Clone, Debug)]
pub enum Mutation {
    SNP(SNP),
    Insertion(Insertion),
    Deletion(Deletion),
}

impl Mutation {
    pub fn get_position(&self) -> u32 {
        match self {
            Mutation::SNP(snp) => snp.get_position(),
            Mutation::Insertion(ins) => ins.get_position(),
            Mutation::Deletion(del) => del.get_position(),
        }
    }

    pub fn get_reference_base(&self) -> char {
        match self {
            Mutation::SNP(snp) => snp.get_reference_base(),
            Mutation::Insertion(ins) => ins.get_reference_base(),
            Mutation::Deletion(del) => del.get_reference_base(),
        }
    }
    pub fn get_alternate_base(&self) -> String {
        match self {
            Mutation::SNP(snp) => snp.get_alternate_base(),
            Mutation::Insertion(ins) => ins.get_alternate_base(),
            Mutation::Deletion(del) => del.get_alternate_base(),
        }
    }

    pub fn get_gene(&self, annotation: &HashMap<(u32, u32), String>) -> Option<Gene> {
        match self {
            Mutation::SNP(snp) => snp.get_gene(annotation),
            Mutation::Insertion(ins) => ins.get_gene(annotation),
            Mutation::Deletion(del) => del.get_gene(annotation),
        }
    }

    pub fn translate(&self, read: &str, read_pos: u32, reference: &fasta::Record, gene: &Gene) -> Option<String> {
        match self {
            Mutation::SNP(snp) => snp.translate(read, read_pos, reference, gene),
            Mutation::Insertion(ins) => ins.translate(read, read_pos, reference, gene),
            Mutation::Deletion(del) => del.translate(read, read_pos, reference, gene),
        }
    }
}

impl fmt::Display for Mutation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Mutation::SNP(snp) => write!(f, "{}{}{}", snp.get_reference_base(), snp.get_position() + 1, snp.get_alternate_base()),
            Mutation::Insertion(ins) => write!(f, "{}{}+{}", ins.get_reference_base(), ins.get_position() + 1, ins.get_alternate_base()),
            Mutation::Deletion(del) => write!(f, "{}{}-{}", del.get_reference_base(), del.get_position() + 1, del.get_alternate_base()),
        }
    }
}