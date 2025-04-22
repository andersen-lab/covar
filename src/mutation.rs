
pub mod snp;
pub mod insertion;
pub mod deletion;

use bio::io::fasta;
use std::collections::HashMap;

use crate::gene::Gene;

pub trait Mutation {
    fn get_position(&self) -> u32; // Add 1 to position for human-readable format, but still 0-based internally
    fn get_reference_base(&self) -> char;
    fn get_alternate_base(&self) -> String;
    fn translate(&self, read: &str, read_pos: u32, reference: &fasta::Record, gene: &Gene) -> Option<String>;
    fn to_string(&self) -> String;

    fn get_gene(&self, annotation: &HashMap<(u32, u32), String>) -> Option<Gene> {
        for (&(start, end), gene_name) in annotation.iter() {
            if self.get_position() >= start && self.get_position() <= end {
                return Some(Gene {
                    start,
                    end,
                    name: gene_name.clone(),
                });
            }
        }
        None
    }
}



