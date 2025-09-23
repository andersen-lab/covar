pub mod snp;
pub mod insertion;
pub mod deletion;

use std::fmt;

pub use snp::SNP;
pub use insertion::Insertion;
pub use deletion::Deletion;

#[derive(Clone, Debug)]
pub enum Mutation {
    SNP(SNP),
    Insertion(Insertion),
    Deletion(Deletion),
}

impl Mutation {
    pub fn get_position(&self) -> u32 {
        match self {
            Mutation::SNP(snp) => snp.get_position(),
            Mutation::Insertion(ins) => ins.get_position(),
            Mutation::Deletion(del) => del.get_position(),
        }
    }

    pub fn get_quality(&self) -> u8 {
        match self {
            Mutation::SNP(snp) => snp.get_quality(),
            Mutation::Insertion(ins) => ins.get_quality(),
            Mutation::Deletion(del) => del.get_quality(),
        }
    }
}

impl fmt::Display for Mutation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Mutation::SNP(snp) => write!(f, "{}{}{}", snp.get_reference_base(), snp.get_position() + 1, snp.get_alternate_base()),
            Mutation::Insertion(ins) => write!(f, "{}{}+{}", ins.get_reference_base(), ins.get_position() + 1, ins.get_alternate_base()),
            Mutation::Deletion(del) => write!(f, "{}{}-{}", del.get_reference_base(), del.get_position() + 1, del.get_alternate_base()),
        }
    }
}

#[cfg(test)]
mod mutation_tests {
    use super::*;
    use crate::gene::Gene;
    use crate::utils;
    use tempfile::NamedTempFile;

    #[test]
    fn test_snp_translate() {
        // Write the mock reference to a temporary FASTA file
        let mock_reference = NamedTempFile::new().unwrap();
        let mock_reference_str = "AAAAAAAAAAAAAAAAAAAAAAAAAA";
        let fasta_content = format!(">chr1\n{}", mock_reference_str);
        std::fs::write(mock_reference.path(), fasta_content).unwrap();

        let mock_reference = utils::read_reference(mock_reference.path()).unwrap();

        let read = "AAAAAAAAAATAAAAAAAAAAAAAAA";
        let mock_snp = SNP::new(10, 'A', 'T', 30, Gene::new("gene1".to_string(), 1));
        
        let translated_aa = mock_snp.translate(read, &mock_reference, 1, read.len() as u32).unwrap();
        
        assert_eq!(translated_aa.to_string(), "gene1:K4I");
    }

    #[test]
    fn test_single_insertion_translate() {
        let mock_insertion = Insertion::new(10, 'A', "CCC".to_string(), 30, Gene::new("gene1".to_string(), 1));
        
        let translated_aa = mock_insertion.translate().unwrap();
        
        assert_eq!(translated_aa.to_string(), "gene1:INS5P");
    }

    #[test]
    fn test_multiple_insertion_translate() {
        let mock_insertion = Insertion::new(10, 'A', "CCCCCC".to_string(), 30, Gene::new("gene1".to_string(), 1));
        
        let translated_aa = mock_insertion.translate().unwrap();
        
        assert_eq!(translated_aa.to_string(), "gene1:INS5PP");
    }

    #[test]
    fn test_deletion_translate() {
        let mock_deletion = Deletion::new(10, 'A', "CCC".to_string(), 30, Gene::new("gene1".to_string(), 1));
        let translated_aa = mock_deletion.translate().unwrap();

        assert_eq!(translated_aa.to_string(), "gene1:DEL5");
    }

    #[test]
    fn test_multiple_deletion_translate() {
        let mock_deletion = Deletion::new(10, 'A', "CCCCCCCCC".to_string(), 30, Gene::new("gene1".to_string(), 1));
        let translated_aa = mock_deletion.translate().unwrap();

        assert_eq!(translated_aa.to_string(), "gene1:DEL5/7");
    }
}