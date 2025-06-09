use std::collections::HashMap;

use crate::gene::Gene;


#[derive(Clone, Debug)]
pub struct Deletion {
    pos: u32,
    ref_base: char,
    alt_sequence: String,
}

impl Deletion {
    pub fn new(pos: u32, ref_base: char, alt_sequence: String) -> Self {
        Deletion {
            pos,
            ref_base,
            alt_sequence,
        }
    }

    pub fn get_position(&self) -> u32 {
        self.pos
    }

    pub fn get_reference_base(&self) -> char {
        self.ref_base
    }

    pub fn get_alternate_base(&self) -> String {
        self.alt_sequence.clone()
    }

    pub fn get_gene(&self, annotation: &HashMap<(u32, u32), String>) -> Option<Gene> {
        for (&(start, end), gene_name) in annotation.iter() {
            if self.get_position() >= start && self.get_position() <= end {
                return Some(Gene {
                    start,
                    name: gene_name.clone(),
                });
            }
        }
        None
    }
    
    pub fn translate(&self, gene: &Gene) -> Option<String> {
        
        // check if deletion is in frame
        let deletion_seq = self.alt_sequence.as_bytes();
        if deletion_seq.len() % 3 != 0 { return None } // not in frame
        let deletion_aa_len = (deletion_seq.len() / 3) as u32;

        let codon_pos = (self.pos + deletion_seq.len() as u32 - gene.get_start()) / 3;

        let translated = if deletion_aa_len > 1 {
            format!("{}:DEL{}/{}", gene.get_name(), codon_pos, codon_pos + deletion_aa_len - 1)
        } else {
            format!("{}:DEL{}", gene.get_name(), codon_pos)
        };

        Some(translated)
    }
}