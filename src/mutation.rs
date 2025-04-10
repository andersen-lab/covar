
pub mod snp;
pub mod insertion;
pub mod deletion;

pub mod mutation {
    use std::collections::HashMap;
    pub trait Mutation: std::fmt::Debug {
        fn get_position(&self) -> u32; // Add 1 to position for human-readable format, but still 0-based internally
        fn get_reference_base(&self) -> char;
        fn get_alternate_base(&self) -> String;
        fn to_string(&self) -> String;
        fn translate(&self, ref_codon: &str, gene_info: (&str, u32)) -> Option<String>;
    
        fn get_gene<'a>(&self, annotation: &'a HashMap<(u32, u32), String>) -> Option<(&'a str, u32)> {
            for (&(start, end), gene) in annotation.iter() {
                if self.get_position() >= start && self.get_position() <= end {
                    return Some((gene, start));
                }
            }
            None
        }
    }
}

