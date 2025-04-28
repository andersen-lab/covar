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
}

impl Insertion {
    pub fn new(pos: u32, ref_base: char, alt_sequence: String) -> Self {
        Insertion {
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
                    end,
                    name: gene_name.clone(),
                });
            }
        }
        None
    }

    pub fn translate(&self, gene: &Gene) -> Option<String> {
        
        let insertion_seq = self.alt_sequence.as_bytes();
        if insertion_seq.len() % 3 != 0 { return None }
        
        let codon_pos = (self.pos + insertion_seq.len() as u32 - gene.get_start()) / 3;

        let alt_codon: Seq<Dna> = match insertion_seq.try_into() {
            Ok(codon) => codon,
            Err(_) => return None,
        };

        let aa_insertion = alt_codon
            .windows(3)
            .map(|codon| STANDARD.to_amino(codon))
            .collect::<Seq<Amino>>()
            .to_string();

        let translated = format!("{}:INS{}{}", gene.get_name(), codon_pos, aa_insertion);

        Some(translated)
    }
}