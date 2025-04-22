use std::fmt;

use crate::mutation::Mutation;

#[derive(Clone)]
pub struct Cluster {
    nt_mutations: Vec<Mutation>,
    aa_mutations: Vec<String>,
    count: u32,
    max_count: u32,
    start: u32,
    end: u32,
}

impl Cluster {
    pub fn new(nt_mutations: Vec<Mutation>, aa_mutations: Vec<String>, count: u32, max_count: u32, start: u32, end: u32) -> Self {
        Self {
            nt_mutations,
            aa_mutations,
            count,
            max_count,
            start,
            end,
        }
    }
    
    pub fn nt_mutations(&self) -> String {
        let mut sorted_nt_mutations = self.nt_mutations.clone();
        sorted_nt_mutations.sort_by(|a, b| {
            let a_pos = a.get_position();
            let b_pos = b.get_position();
            a_pos.cmp(&b_pos)
        });
        sorted_nt_mutations
            .iter()
            .map(|mutation| mutation.to_string())
            .collect::<Vec<_>>()
            .join(",")
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

    pub fn len(&self) -> usize {
        self.nt_mutations.len()
    }
}

impl fmt::Display for Cluster {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f,
        "{}\t{}\t{}\t{}\t{}\t{}",
        self.nt_mutations(),
        self.aa_mutations.join(","),
        self.count,
        self.max_count,
        self.start,
        self.end)
    }
}

pub fn merge_clusters(clusters: &Vec<Cluster>) -> Vec<Cluster> {
    let mut merged_clusters = Vec::<Cluster>::new();

    for cluster in clusters {
        let mut found = false;
        for merged_cluster in &mut merged_clusters {
            if cluster.nt_mutations() == merged_cluster.nt_mutations() {
                // Merge clusters
                merged_cluster.count += cluster.count;
                merged_cluster.start = merged_cluster.start.max(cluster.start);
                merged_cluster.end = merged_cluster.end.min(cluster.end);
                found = true;
                break;
            }
        }
        if !found {
            // Add new cluster
            merged_clusters.push(cluster.clone());
        }
    }
    merged_clusters
}