use bio_seq::prelude::*;
use bio_seq::translation::STANDARD;
use bio_seq::translation::TranslationTable;

use crate::gene::Gene;

#[derive(Clone, Debug)]
pub struct Insertion {
    pos: u32,
    ref_base: char,
    alt_sequence: String,
    quality: u8,
    gene: Gene,
}

impl Insertion {
    pub fn new(pos: u32, ref_base: char, alt_sequence: String, quality: u8, gene: Gene) -> Self {
        Insertion {
            pos,
            ref_base,
            alt_sequence,
            quality,
            gene
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

    pub fn get_quality(&self) -> u8 {
        self.quality
    }

    pub fn translate(&self) -> Option<String> {
        
        let insertion_seq = self.alt_sequence.as_bytes();
        if insertion_seq.len() % 3 != 0 { return None }
        
        let codon_pos = ((self.pos + 1 - self.gene.get_start()) / 3) + 2;

        let alt_codon: Seq<Dna> = insertion_seq.try_into().expect("Invalid UTF-8 in reference sequence"); // panics

        let aa_insertion = alt_codon
            .chunks(3)
            .map(|codon| STANDARD.to_amino(codon))
            .collect::<Seq<Amino>>()
            .to_string();

        let translated = format!("{}:INS{}{}", self.gene.get_name(), codon_pos, aa_insertion);

        Some(translated)
    }
}