use bio_seq::prelude::*;
use bio_seq::translation::STANDARD;
use bio_seq::translation::TranslationTable;
use bio::io::fasta;

use super::Mutation;
use crate::gene::Gene;


#[derive(Debug)]
pub struct SNP {
    pos: u32,
    ref_base: char,
    alt_base: char,
}

impl SNP {
    pub fn new(pos: u32,ref_base: char, alt_base: char) -> Self {
        SNP {
            pos,
            ref_base,
            alt_base,
        }
    }
}

impl Mutation for SNP {

    fn get_position(&self) -> u32 {
        self.pos
    }


    fn get_reference_base(&self) -> char {
        self.ref_base
    }

    fn get_alternate_base(&self) -> String {
        self.alt_base.to_string()
    }

    fn to_string(&self) -> String {
        format!("{}{}{}", self.ref_base, self.pos + 1, self.alt_base)
    }

    fn translate(&self, read: &str, read_pos: u32, reference: &fasta::Record, gene: &Gene) -> Option<String> {
        let codon_phase = (self.pos - gene.get_start()) % 3;

        let ref_start_pos = (self.pos - codon_phase - 1) as usize;
        let ref_end_pos = (self.pos - codon_phase + 2) as usize;
        let ref_codon: Seq<Dna> = reference.seq()[ref_start_pos..ref_end_pos].try_into().unwrap(); // panics
        let ref_aa = STANDARD.to_amino(&ref_codon).to_string();

        
        let alt_start_pos = (read_pos - codon_phase - 1) as usize;
        let alt_end_pos = (read_pos - codon_phase + 2) as usize;
        let alt_codon: Seq<Dna> = read[alt_start_pos..alt_end_pos].try_into().unwrap(); // panics
        let alt_aa = STANDARD.to_amino(&alt_codon).to_string();

        let translated = format!("{}:{}{}{}", gene.get_name(), ref_aa, self.pos + 1, alt_aa);

        Some(translated)
    }
}