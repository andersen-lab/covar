use super::mutation;
use bio::io::fasta;
use mutation::Gene;

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

impl mutation::Mutation for Deletion {
    fn get_position(&self) -> u32 {
        self.pos
    }

    fn get_reference_base(&self) -> char {
        self.ref_base
    }

    fn get_alternate_base(&self) -> String {
        self.alt_sequence.clone()
    }

    fn to_string(&self) -> String {
        format!("{}{}-{}", self.ref_base, self.pos + 1, self.alt_sequence)
    }

    fn translate(&self, read: &str, read_pos: u32, reference: &fasta::Record, gene: &Gene) -> Option<String> {
        // Placeholder for translation logic
        Some(format!("{} -> {} (deletion)", self.ref_base, self.alt_sequence))
    }
}