use std::collections::HashMap;
use std::error::Error;
use std::io::ErrorKind;
use std::path::PathBuf;

use bio::io::{fasta, gff};
use bio::io::fasta::FastaRead;

use rust_htslib::bam::{IndexedReader, Read, Record};
use rust_htslib::bam::record::Cigar;

use crate::mutation::{Mutation, snp::SNP, insertion::Insertion, deletion::Deletion};
use crate::cluster::Cluster;

pub fn read_reference(path: &PathBuf) -> Result<fasta::Record, Box<dyn Error>> {
    let mut reader = fasta::Reader::from_file(path)?;
    let mut reference = fasta::Record::new();
    reader.read(&mut reference)?;
    Ok(reference)
}

pub fn read_annotation(path: &PathBuf) -> Result<HashMap<(u32, u32), String>, Box<dyn Error>> {
    let mut reader = gff::Reader::from_file(path, gff::GffType::GFF3)?;
    let mut gene_regions: HashMap<(u32, u32), String> = HashMap::new();

    for record in reader.records() {
        let rec = record?;
        if rec.feature_type() == "CDS" {
            if let Some(gene) = rec.attributes().get("gene") {
                gene_regions.insert((*rec.start() as u32, *rec.end() as u32), gene.to_string());
            }
        }
    }
    Ok(gene_regions)
}

pub fn read_pair_generator(
    bam: & mut IndexedReader,
    refname: &str,
    min_site: i64,
    max_site: i64,
) -> Result<Vec<(Option<Record>, Option<Record>)>, Box<dyn Error>> {

    let tid = match bam.header().tid(refname.as_bytes())
            .ok_or_else(|| std::io::Error::new(
                ErrorKind::NotFound, 
                format!("Reference '{}' not found", refname)
            )) {
        Ok(t) => t,
        Err(e) => return Err(Box::new(e) as Box<dyn Error>),
    };

    let _ = bam.fetch((tid, min_site , max_site + 1))
        .map_err(|e| Box::new(e) as Box<dyn Error>);

    // Get read pairs by query name
    let mut read_pairs: HashMap<Vec<u8>, (Option<Record>, Option<Record>)> = HashMap::new();
    for record_result in bam.records() {
        let record = match record_result {
            Ok(r) => r,
            Err(e) => return Err(Box::new(e) as Box<dyn Error>),
        };

        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }

        let query_name = record.qname().to_owned();
        
        let entry = read_pairs.entry(query_name).or_insert((None, None));
        
        if record.is_first_in_template() {
            entry.0 = Some(record);
        } else if record.is_last_in_template() {
            entry.1 = Some(record);
        } else {
            // Handle singleton reads
            if entry.0.is_none() {
                entry.0 = Some(record);
            } else {
                entry.1 = Some(record);
            }
        }
    }

    Ok(read_pairs.into_iter().map(|(_, pair)| pair).collect())
}

pub fn call_variants(
    read_pair: (Option<Record>, Option<Record>),
    reference: &fasta::Record,
    annotation: &HashMap<(u32, u32), String>
) -> Cluster { 
    let (read1, read2) = read_pair;

    if read1.is_none() && read2.is_none() {
        panic!("Both reads are missing, cannot proceed with variant calling.")
    }
    
    // Function to process a single read in the pair
    let process_read = |read: &Record| -> Vec<(Box<dyn Mutation>, String)> {
        let mut local_variants = Vec::new();

        let ref_seq = reference.seq();
        let ref_seq_len = ref_seq.len() as u32;

        let read_seq = String::from_utf8(read.seq().as_bytes().to_vec()).expect("Invalid UTF-8 in read sequence"); // panics
        let read_seq_len = read_seq.len() as u32;
        let cigar = read.cigar();

        let start_pos: u32 = read.pos() as u32;
        let mut read_pos: u32 = 0;
        let mut ref_pos: u32 = start_pos;

        for c in cigar.iter() { // TODO: handle indel offsets
            match c {
                Cigar::Match(len) => { // Call SNPs
                    for i in 0..*len {
                        if ref_pos + i < ref_seq_len && read_pos + i < read_seq_len {
                            let ref_base = ref_seq[(ref_pos + i) as usize] as char;
                            let read_base = read_seq.chars().nth((read_pos + i) as usize).unwrap();
                            
                            if ref_base != read_base && read_base != 'N' {
                                let snp = SNP::new(ref_pos + i, ref_base, read_base);
                                if let Some(gene) = snp.get_gene(annotation) {
                                    let aa_mutation = snp.translate(&read_seq, read_pos, reference, &gene)
                                        .unwrap_or_else(|| "Unknown".to_string());
                                    local_variants.push((Box::new(snp) as Box<dyn Mutation>, aa_mutation));                                    
                                }
                            }
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                },
                Cigar::Ins(len) => { // Call insertions
                    if ref_pos > 0 && read_pos + *len <= read_seq_len {
                        let ref_base = ref_seq[(ref_pos - 1) as usize] as char;
                        let ins_seq = read_seq[(read_pos as usize)..(read_pos as usize + *len as usize)].to_string();
                        let insertion = Insertion::new(ref_pos - 1, ref_base, ins_seq);
                        if let Some(gene) = insertion.get_gene(annotation) {
                            let aa_mutation = insertion.translate(&read_seq, read_pos, reference, &gene)
                                .unwrap_or_else(|| "Unknown".to_string());
                            local_variants.push((Box::new(insertion) as Box<dyn Mutation>, aa_mutation));                                    
                        }
                    }
                    read_pos += len;
                },
                Cigar::Del(len) => { // Call deletions
                    if ref_pos > 0 && ref_pos + *len <= ref_seq_len {
                        let ref_base = ref_seq[(ref_pos - 1) as usize] as char;
                        let del_seq = std::str::from_utf8(&ref_seq[(ref_pos as usize)..(ref_pos as usize + *len as usize)])
                            .expect("Invalid UTF-8 sequence in reference")
                            .to_string();
                        let deletion = Deletion::new((ref_pos - 1), ref_base, del_seq);
                        if let Some(gene) = deletion.get_gene(annotation) {
                            let aa_mutation = deletion.translate(&read_seq, read_pos, reference, &gene)
                                .unwrap_or_else(|| "Unknown".to_string());
                            local_variants.push((Box::new(deletion) as Box<dyn Mutation>, aa_mutation));                                    
                        }
                    }
                    ref_pos += len;
                },
                Cigar::SoftClip(len) => {
                    read_pos += len;
                },
                _ => {},
            }
        }
        local_variants
    };

    let mut variants: Vec<(Box<dyn Mutation>, String)> = Vec::new();
    let mut range: (u32, u32) = (0, u32::MAX); // Consider using htslib range type here

    // Process read1 if present
    if let Some(r1) = read1 {
        let r1_variants = process_read(&r1);
        variants.extend(r1_variants);

        let start_pos = r1.pos() as u32;
        let end_pos = r1.cigar().end_pos() as u32;

        if start_pos < range.0 {
            range.0 = start_pos;
        }
        if end_pos > range.1 {
            range.1 = end_pos;
        }
    }
    
    // Process read2 if present
    if let Some(r2) = read2 {
        let r2_variants = process_read(&r2);
        variants.extend(r2_variants);
        let start_pos = r2.pos() as u32;
        let end_pos = r2.cigar().end_pos() as u32;

        if start_pos < range.0 {
            range.0 = start_pos;
        }
        if end_pos > range.1 {
            range.1 = end_pos;
        }
    }
    
    
    // Basic deduplication - keep unique variants only
    let mut unique_variants = Vec::new();
    let mut seen_variants = std::collections::HashSet::new();
    
    for (var, aa_mut) in variants {
        let var_str = var.to_string();
        if !seen_variants.contains(&var_str) {
            seen_variants.insert(var_str);
            unique_variants.push((var, aa_mut));
        }
    }

    // Unpack nt and aa mutations
    let mut nt_mutations = Vec::new();
    let mut aa_mutations = Vec::new();
    for (var, aa_mut) in unique_variants {
        nt_mutations.push(var.to_string());
        aa_mutations.push(aa_mut);
    }
    
    Cluster::new(
        nt_mutations,
        aa_mutations,
        1,
        1,
        range.0, 
        range.1,
    )
}

