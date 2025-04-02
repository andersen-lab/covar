use rust_htslib::bam::Record;
use bio::io::fasta;
use std::collections::HashMap;
use std::error::Error;

pub trait Mutation: std::fmt::Debug {
    fn position(&self) -> u32; // Add 1 to position for human-readable format, but still 0-based internally
    fn reference_base(&self) -> char; 
    fn alternate_base(&self) -> String;
    fn to_string(&self) -> String; 
    fn translate(&self, ref_codon: &str, read_codon: &str, gene: &str) -> String; // TODO: add read/reference info as parameters
}

#[derive(Debug)]
pub struct SNP {
    pos: u32,
    ref_base: char,
    alt_base: char,
}

impl SNP {
    pub fn new(pos: u32, ref_base: char, alt_base: char) -> Self {
        SNP {
            pos,
            ref_base,
            alt_base,
        }
    }
}

impl Mutation for SNP {
    fn position(&self) -> u32 {
        self.pos
    }

    fn reference_base(&self) -> char {
        self.ref_base
    }

    fn alternate_base(&self) -> String {
        self.alt_base.to_string()
    }

    fn to_string(&self) -> String {
        format!("{}{}{}", self.ref_base, self.pos + 1, self.alt_base)
    }

    fn translate(&self , ref_codon: &str, read_codon: &str, gene: &str) -> String {
        // Placeholder for translation logic
        format!("{} -> {}", self.ref_base, self.alt_base)
    }
}
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

impl Mutation for Insertion {

    fn position(&self) -> u32 {
        self.pos
    }

    fn reference_base(&self) -> char {
        self.ref_base
    }

    fn alternate_base(&self) -> String {
        self.alt_sequence.clone()
    }

    fn to_string(&self) -> String {
        format!("{}{}+{}", self.ref_base, self.pos + 1, self.alt_sequence)
    }

    fn translate(&self, ref_codon: &str, read_codon: &str, gene: &str) -> String {
        // Placeholder for translation logic
        format!("{} -> {} (insertion)", self.ref_base, self.alt_sequence)
    }
}
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

impl Mutation for Deletion {
    fn position(&self) -> u32 {
        self.pos
    }

    fn reference_base(&self) -> char {
        self.ref_base
    }

    fn alternate_base(&self) -> String {
        self.alt_sequence.clone()
    }

    fn to_string(&self) -> String {
        format!("{}{}-{}", self.ref_base, self.pos + 1, self.alt_sequence)
    }

    fn translate(&self, ref_codon: &str, read_codon: &str, gene: &str) -> String {
        // Placeholder for translation logic
        format!("{} -> {} (insertion)", self.ref_base, self.alt_sequence)
    }
}



pub fn call_variants(
    read_pair: (Option<Record>, Option<Record>),
    reference: &fasta::Record,
    annotation: &HashMap<(u32, u32), String>,
) -> Result<Vec<Box<dyn Mutation>>, Box<dyn Error>> {
    let (read1, read2) = read_pair;
    if read1.is_none() && read2.is_none() {
        return Err("Both reads are missing, cannot proceed with variant calling.".into());
    }
    
    let mut variants: Vec<Box<dyn Mutation>> = Vec::new();
    
    // Function to process a single read
    let process_read = |read: &Record| -> Result<Vec<Box<dyn Mutation>>, Box<dyn Error>> {
        let mut local_variants = Vec::new();

        let ref_seq = reference.seq();
        let ref_seq_len = ref_seq.len() as u32;

        let read_seq = String::from_utf8(read.seq().as_bytes().to_vec())
            .map_err(|e| Box::new(e) as Box<dyn Error>)?;
        let read_seq_len = read_seq.len() as u32;
        let cigar = read.cigar();

        let start_pos: u32 = read.pos() as u32;
        let mut read_pos: u32 = 0;
        let mut ref_pos: u32 = start_pos as u32;
        for c in cigar.iter() {
            match c {
                rust_htslib::bam::record::Cigar::Match(len) => {
                    // Check for SNPs in matched regions
                    for i in 0..*len {
                        if ref_pos + i < ref_seq_len && read_pos + i < read_seq_len {
                            let ref_base = ref_seq[(ref_pos + i) as usize] as char;
                            let read_base = read_seq.chars().nth((read_pos + i) as usize).unwrap();
                            
                            if ref_base != read_base && read_base != 'N' {
                                let snp = SNP::new((ref_pos + i) as u32, ref_base, read_base);
                                local_variants.push(Box::new(snp) as Box<dyn Mutation>);
                            }
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                },
                rust_htslib::bam::record::Cigar::Ins(len) => {
                    if ref_pos > 0 && read_pos + *len <= read_seq_len {
                        let ref_base = ref_seq[(ref_pos - 1) as usize] as char;
                        let ins_seq = read_seq[(read_pos as usize)..(read_pos as usize + *len as usize)].to_string();
                        let insertion = Insertion::new((ref_pos - 1) as u32, ref_base, ins_seq);
                        local_variants.push(Box::new(insertion) as Box<dyn Mutation>);
                    }
                    read_pos += len;
                },
                rust_htslib::bam::record::Cigar::Del(len) => {
                    if ref_pos > 0 && ref_pos + *len <= ref_seq_len {
                        let ref_base = ref_seq[(ref_pos - 1) as usize] as char;
                        let del_seq = std::str::from_utf8(&ref_seq[(ref_pos as usize)..(ref_pos as usize + *len as usize)])?.to_string();
                        let deletion = Deletion::new((ref_pos - 1) as u32, ref_base, del_seq);
                        local_variants.push(Box::new(deletion) as Box<dyn Mutation>);
                    }
                    ref_pos += len;
                },
                rust_htslib::bam::record::Cigar::SoftClip(len) => {
                    read_pos += len;
                },
                _ => {},  // Handle other CIGAR operations as needed
            }
        }
        
        Ok(local_variants)
    };
    
    // Process read1 if present
    if let Some(r1) = read1 {
        let r1_variants = process_read(&r1)?;
        variants.extend(r1_variants);
    }
    
    // Process read2 if present
    if let Some(r2) = read2 {
        let r2_variants = process_read(&r2)?;
        variants.extend(r2_variants);
    }
    
    // Basic deduplication - keep unique variants only
    let mut unique_variants = Vec::new();
    let mut seen_variants = std::collections::HashSet::new();
    
    for var in variants {
        let var_str = var.to_string();
        if !seen_variants.contains(&var_str) {
            seen_variants.insert(var_str);
            unique_variants.push(var);
        }
    }
    
    Ok(unique_variants)
}