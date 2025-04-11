
pub mod snp;
pub mod insertion;
pub mod deletion;

pub mod mutation {
    use bio::io::fasta;
    use std::collections::HashMap;
    pub trait Mutation: std::fmt::Debug {
        fn get_position(&self) -> u32; // Add 1 to position for human-readable format, but still 0-based internally
        fn get_reference_base(&self) -> char;
        fn get_alternate_base(&self) -> String;
        fn to_string(&self) -> String;
        fn translate(&self, read: &str, read_pos: u32, reference: &fasta::Record, gene: &Gene) -> Option<String>;
    
        fn get_gene<'a>(&self, annotation: &'a HashMap<(u32, u32), String>) -> Option<Gene> {
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
    }

    pub struct Gene {
        start: u32,
        end: u32,
        name: String,
    }

    impl Gene {
        pub fn get_name(&self) -> &str {
            &self.name
        }

        pub fn get_start(&self) -> u32 {
            self.start
        }

        pub fn get_end(&self) -> u32 {
            self.end
        }
    }
}

