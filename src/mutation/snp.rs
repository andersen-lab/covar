use std::collections::HashMap;

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
    aa_mutation: Option<String>,
}

impl SNP {
    pub fn new(pos: u32,ref_base: char, alt_base: char, quality: u8) -> Self {
        SNP {
            pos,
            ref_base,
            alt_base,
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
        self.alt_base.to_string()
    }

    pub fn get_quality(&self) -> u8 {
        self.quality
    }

    pub fn set_aa_mutation(&mut self, aa_mutation: Option<String>) {
        self.aa_mutation = aa_mutation;
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

    pub fn translate(&mut self, read: &str, read_pos: u32, read_qual: &[u8], reference: &fasta::Record, gene: &Gene) {
        let mut_pos_one_based = self.pos + 1; // 1-based for gene positions
        
        let codon_phase = (mut_pos_one_based - gene.get_start()) % 3;
        let codon_pos = (mut_pos_one_based - gene.get_start()) / 3;
        
        // handle codon spanning reads
        if (read_pos as i32 - codon_phase as i32) < 0 { return }
        if read_pos - codon_phase + 3  > read.len() as u32 { return }
        
        let ref_start_pos = (mut_pos_one_based - codon_phase - 1) as usize;
        let ref_end_pos = (mut_pos_one_based - codon_phase + 2) as usize;
        let ref_codon: Seq<Dna> = reference.seq()[ref_start_pos..ref_end_pos].try_into().expect("Invalid UTF-8 in reference sequence"); // panics
        let ref_aa = STANDARD.to_amino(&ref_codon).to_string();
        
        let alt_start_pos = (read_pos - codon_phase) as usize;
        let alt_end_pos = (read_pos - codon_phase + 3) as usize;
        let alt_codon: Seq<Dna> = read[alt_start_pos..alt_end_pos].try_into().expect("Invalid UTF-8 in read sequence"); // panics
        
        // check min quality of all bases in codon 
        if read_qual[alt_start_pos..alt_end_pos].iter().any(|&qual| qual < self.quality) {
            self.quality = *read_qual[alt_start_pos..alt_end_pos].iter().min().unwrap();
        }

        let alt_aa = STANDARD.to_amino(&alt_codon).to_string();
        let translated = format!("{}:{}{}{}", gene.get_name(), ref_aa, codon_pos + 1, alt_aa);

        self.aa_mutation = Some(translated);
    }
}