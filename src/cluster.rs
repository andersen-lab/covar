use std::fmt;

#[derive(Clone)]
pub struct Cluster {
    pub aa_mutations: String,
    pub nt_mutations: String,
    pub cluster_depth: u32,
    pub total_depth: u32,
    pub coverage_start: u32,
    pub coverage_end: u32,
    pub mutations_start: u32,
    pub mutations_end: u32,
}

impl Cluster {
    pub fn new(nt_mutations: String, aa_mutations: String, max_count: u32, coverage_start: u32, coverage_end: u32, mutations_start: u32, mutations_end: u32) -> Self {
        Self {
            nt_mutations,
            aa_mutations,
            cluster_depth: 1,
            total_depth: max_count,
            coverage_start,
            coverage_end,
            mutations_start,
            mutations_end,
        }
    }
}

impl fmt::Display for Cluster {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
        "{}\t{}\t{}\t{}\t{}",
        self.nt_mutations,
        self.aa_mutations,
        self.cluster_depth,
        self.coverage_start,
        self.coverage_end)
    }
}