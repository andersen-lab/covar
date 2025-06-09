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
}

impl SNP {
    pub fn new(pos: u32,ref_base: char, alt_base: char) -> Self {
        SNP {
            pos,
            ref_base,
            alt_base,
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

    pub fn translate(&self, read: &str, read_pos: u32, reference: &fasta::Record, gene: &Gene) -> Option<String> {
        let mut_pos_one_based = self.pos + 1; // 1-based for gene positions
        
        let codon_phase = (mut_pos_one_based - gene.get_start()) % 3;
        let codon_pos = (mut_pos_one_based - gene.get_start()) / 3;

        // handle codon spanning reads
        if (read_pos as i32 - codon_phase as i32) < 0 { return None }
        if read_pos - codon_phase + 3  > read.len() as u32 { return None }

        let ref_start_pos = (mut_pos_one_based - codon_phase - 1) as usize;
        let ref_end_pos = (mut_pos_one_based - codon_phase + 2) as usize;
        let ref_codon: Seq<Dna> = match reference.seq()[ref_start_pos..ref_end_pos].try_into() {
            Ok(codon) => codon,
            Err(_) => return None,
        }; 
        let ref_aa = STANDARD.to_amino(&ref_codon).to_string();

        let alt_start_pos = (read_pos - codon_phase) as usize;
        let alt_end_pos = (read_pos - codon_phase + 3) as usize;
        let alt_codon: Seq<Dna> = match read[alt_start_pos..alt_end_pos].try_into() {
            Ok(codon) => codon,
            Err(_) => return None,
        };

        let alt_aa = STANDARD.to_amino(&alt_codon).to_string();
        let translated = format!("{}:{}{}{}", gene.get_name(), ref_aa, codon_pos + 1, alt_aa);

        Some(translated)
    }
}