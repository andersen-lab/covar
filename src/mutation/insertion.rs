use std::fmt;

use bio_seq::prelude::*;
use bio_seq::translation::STANDARD;
use bio_seq::translation::TranslationTable;
use bio::io::fasta;

use super::Mutation;
use crate::gene::Gene;


#[derive(Debug)]
pub struct Insertion {
    pos: u32,
    ref_base: char,
    alt_sequence: String,
}

impl Insertion {
    pub fn new(pos: u32, ref_base: char, alt_sequence: String) -> Self {
        Insertion {
            pos,
            ref_base,
            alt_sequence,
        }
    }
}

impl Mutation for Insertion {

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
        
        // check if insertion is in frame
        let insertion_seq = self.alt_sequence.as_bytes();
        if insertion_seq.len() % 3 != 0 { return None }

        let alt_codon: Seq<Dna> = match insertion_seq.try_into() {
            Ok(codon) => codon,
            Err(_) => return None,
        };

        let aa_insertion = alt_codon
            .windows(3)
            .map(|codon| STANDARD.to_amino(&codon))
            .collect::<Seq<Amino>>()
            .to_string();

        let translated = format!("{}:INS{}{}", gene.get_name(), codon_pos + 1, aa_insertion);

        Some(translated)
    }

    fn to_string(&self) -> String {
        format!("{}{}+{}", self.ref_base, self.pos + 1, self.alt_sequence)
    }
}
