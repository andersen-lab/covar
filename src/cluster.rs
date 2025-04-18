pub struct Cluster {
    nt_mutations: Vec<String>,
    aa_mutations: Vec<String>,
    count: u32,
    max_count: u32,
    start: u32,
    end: u32,
}

impl Cluster {
    pub fn new(nt_mutations: Vec<String>, aa_mutations: Vec<String>, count: u32, max_count: u32, start: u32, end: u32) -> Self {
        Self {
            nt_mutations,
            aa_mutations,
            count,
            max_count,
            start,
            end,
        }
    }

    pub fn nt_mutations(&self) -> &Vec<String> {
        &self.nt_mutations
    }

    pub fn count(&self) -> u32 {
        self.count
    }

    pub fn max_count(&self) -> u32 {
        self.max_count
    }

    pub fn start(&self) -> u32 {
        self.start
    }

    pub fn end(&self) -> u32 {
        self.end
    }
}

// pub fn merge_clusters(clusters: &Vec<Cluster>);