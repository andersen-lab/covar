use super::mutation;

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

impl mutation::Mutation for Insertion {

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
        format!("{}{}+{}", self.ref_base, self.pos + 1, self.alt_sequence)
    }

    fn translate(&self, reference: &str, gene_info: (&str, u32)) -> Option<String> {
        // Placeholder for translation logic
        Some(format!("{} -> {} (insertion)", self.ref_base, self.alt_sequence))
    }
}