use super::mutation;

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

impl mutation::Mutation for SNP {
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

    fn translate(&self, ref_codon: &str, gene_info: (&str, u32)) -> Option<String> {
        Some(format!("{} -> {} (SNP)", self.ref_base, self.alt_base))
    }
}