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