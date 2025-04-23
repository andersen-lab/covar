use std::fmt;
use std::collections::HashMap;

#[derive(Clone)]
pub struct Cluster {
    nt_mutations: Vec<String>,
    aa_mutations: Vec<String>,
    count: u32,
    max_count: u32,
    coverage_start: u32,
    coverage_end: u32,
    mutations_start: u32,
    mutations_end: u32,
}

impl Cluster {
    pub fn new(nt_mutations: Vec<String>, aa_mutations: Vec<String>, count: u32, max_count: u32, coverage_start: u32, coverage_end: u32, mutations_start: u32, mutations_end: u32) -> Self {
        Self {
            nt_mutations,
            aa_mutations,
            count,
            max_count,
            coverage_start,
            coverage_end,
            mutations_start,
            mutations_end,
        }
    }

    pub fn nt_mutations(&self) -> String {
        self.nt_mutations.join(",")
    }

    pub fn aa_mutations(&self) -> String {
        self.aa_mutations.join(",")
    }

    pub fn count(&self) -> u32 {
        self.count
    }

    pub fn max_count(&self) -> u32 {
        self.max_count
    }

    pub fn start(&self) -> u32 {
        self.coverage_start
    }

    pub fn end(&self) -> u32 {
        self.coverage_end
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
        self.aa_mutations(),
        self.count,
        self.max_count,
        self.coverage_start,
        self.coverage_end)
    }
}

pub fn merge_clusters(clusters: &Vec<Cluster>) -> Vec<Cluster> {

    let mut merged_map: HashMap<String, Cluster> = HashMap::new();

    for cluster in clusters {
        let key = cluster.nt_mutations();
        if let Some(merged_cluster) = merged_map.get_mut(&key) {
            // Merge clusters
            merged_cluster.count += cluster.count;
            merged_cluster.coverage_start = merged_cluster.coverage_start.max(cluster.coverage_start);
            merged_cluster.coverage_end = merged_cluster.coverage_end.min(cluster.coverage_end);
        } else {
            // Add new cluster
            merged_map.insert(key.clone(), cluster.clone());
        }

        // Update max_count for overlapping clusters
        for merged_cluster in merged_map.values_mut() {
            if cluster.coverage_start <= merged_cluster.mutations_start && cluster.coverage_end >= merged_cluster.mutations_end {
                merged_cluster.max_count += cluster.count;
            }
        }
    }

    merged_map.into_values().collect()
}
