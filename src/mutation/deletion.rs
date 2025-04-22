use std::fmt;

use bio::io::fasta;

use super::Mutation;
use crate::gene::Gene;


#[derive(Debug)]
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
}

impl Mutation for Deletion {
    fn get_position(&self) -> u32 {
        self.pos
    }

    fn get_reference_base(&self) -> char {
        self.ref_base
    }

    fn get_alternate_base(&self) -> String {
        self.alt_sequence.clone()
    }

    fn translate(&self, read: &str, read_pos: u32, reference: &fasta::Record, gene: &Gene) -> Option<String> {
        let codon_pos = (self.pos - gene.get_start()) / 3;
        
        // check if deletion is in frame
        let deletion_seq = self.alt_sequence.as_bytes();
        if deletion_seq.len() % 3 != 0 { return None }
        let deletion_aa_len = (deletion_seq.len() / 3) as u32;

        let mut translated = String::new();

        if deletion_aa_len > 1 {
            translated = format!("{}:DEL{}-{}", gene.get_name(), codon_pos + 1, codon_pos + deletion_aa_len + 1);
        } else {
            translated = format!("{}:DEL{}", gene.get_name(), codon_pos + 1);
        }


        Some(translated)
    }

    fn to_string(&self) -> String {
        format!("{}{}-{}", self.ref_base, self.pos + 1, self.alt_sequence)
    }
}
