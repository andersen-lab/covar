use std::collections::HashMap;

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
    aa_mutation: Option<String>,
}

impl Insertion {
    pub fn new(pos: u32, ref_base: char, alt_sequence: String, quality: u8) -> Self {
        Insertion {
            pos,
            ref_base,
            alt_sequence,
            quality,
            aa_mutation: None,
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

    pub fn get_aa_mutation(&self) -> Option<String> {
        self.aa_mutation.clone()
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

    pub fn translate(&mut self, gene: &Gene) {
        
        let insertion_seq = self.alt_sequence.as_bytes();
        if insertion_seq.len() % 3 != 0 { return; }
        
        let codon_pos = ((self.pos + 1 - gene.get_start()) / 3) + 2;

        let alt_codon: Seq<Dna> = insertion_seq.try_into().expect("Invalid UTF-8 in reference sequence"); // panics

        let aa_insertion = alt_codon
            .chunks(3)
            .map(|codon| STANDARD.to_amino(codon))
            .collect::<Seq<Amino>>()
            .to_string();

        let translated = format!("{}:INS{}{}", gene.get_name(), codon_pos, aa_insertion);

        self.aa_mutation = Some(translated);
    }
}