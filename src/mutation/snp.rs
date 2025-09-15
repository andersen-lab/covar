use bio_seq::prelude::*;
use bio_seq::translation::STANDARD;
use bio_seq::translation::TranslationTable;
use bio::io::fasta;

use crate::gene::Gene;

#[derive(Debug, Clone)]
pub struct SNP {
    pos: u32,
    ref_base: char,
    alt_base: char,
    quality: u8,
    gene: Gene,
}

impl SNP {
    pub fn new(pos: u32, ref_base: char, alt_base: char, quality: u8, gene: Gene) -> Self {
        SNP {
            pos,
            ref_base,
            alt_base,
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
        self.alt_base.to_string()
    }

    pub fn get_quality(&self) -> u8 {
        self.quality
    }


    pub fn translate(&self, read: &str, reference: &fasta::Record) -> Option<String> {
        let mut_pos_one_based = self.pos + 1; // 1-based for gene positions
        
        let codon_phase = (mut_pos_one_based - self.gene.get_start()) % 3;
        let codon_pos = (mut_pos_one_based - self.gene.get_start()) / 3;
        
        let ref_start_pos = (mut_pos_one_based - codon_phase - 1) as usize;
        let ref_end_pos = (mut_pos_one_based - codon_phase + 2) as usize;
        let ref_codon: Seq<Dna> = reference.seq()[ref_start_pos..ref_end_pos].try_into().expect("Invalid UTF-8 in reference sequence"); // panics
        let ref_aa = STANDARD.to_amino(&ref_codon).to_string();
        
        let alt_start_pos = (self.pos - codon_phase) as usize;
        let alt_end_pos = (self.pos - codon_phase + 3) as usize;
        let alt_codon: Seq<Dna> = read[alt_start_pos..alt_end_pos].try_into().expect("Invalid UTF-8 in read sequence"); // panics

        let alt_aa = STANDARD.to_amino(&alt_codon).to_string();
        let translated = format!("{}:{}{}{}", self.gene.get_name(), ref_aa, codon_pos + 1, alt_aa);

        Some(translated)
    }
}