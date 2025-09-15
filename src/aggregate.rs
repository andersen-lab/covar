use std::error::Error;
use std::collections::HashMap;

use bio::io::fasta;
use polars::prelude::*;

use crate::Config;
use crate::cluster::Cluster;

macro_rules! struct_to_dataframe {
    ($input:expr, [$($field:ident),+]) => {
        {
            let len = $input.len().to_owned();

            // Extract the field values into separate vectors
            $(let mut $field = Vec::with_capacity(len);)*

            for e in $input.iter() {
                $($field.push(e.$field.clone());)*
            }
            df! {
                $(stringify!($field) => $field,)*
            }
        }
    };
}


pub fn aggregate_clusters(clusters: &[Cluster], reference: &fasta::Record, config: &Config) -> Result<DataFrame, Box<dyn Error>> {

    let df = match struct_to_dataframe!(clusters,
            [nt_mutations, cluster_depth, total_depth, coverage_start, coverage_end, mutations_start, mutations_end]) {
            Ok(df) => df.lazy(),
            Err(e) => panic!("Error creating DataFrame: {}", e),
        };

    let df = df
        .group_by_stable([col("nt_mutations")])
        .agg([
            col("cluster_depth").sum().alias("cluster_depth"),
            col("total_depth").max().alias("total_depth"),
            col("coverage_start").max().alias("coverage_start"),
            col("coverage_end").min().alias("coverage_end"),
        ])

        // Calculate frequency column
        .with_column(
            (col("cluster_depth").cast(DataType::Float64) / col("total_depth").cast(DataType::Float64))
            .round(6)
            .alias("frequency")
        )

        // Filter by CLI parameters
        .filter(col("total_depth").gt(lit(0)))
        .filter(col("cluster_depth").gt_eq(lit(config.min_depth))) 
        .filter(col("frequency").gt_eq(lit(config.min_frequency)))

        .collect()?;

    // Translate
    let mut nt_to_aa: HashMap<String, String> = HashMap::new();
    for cluster in clusters {
        nt_to_aa
            .entry(cluster.nt_mutations.clone())
            .or_insert_with(|| cluster.translate_cluster(reference));
    }

    let aa_mutations: Vec<String> = df
        .column("nt_mutations")?
        .str()?
        .into_iter()
        .filter_map(|nt_mut| nt_mut.and_then(|nt| nt_to_aa.get(nt).cloned()))
        .collect(); 

    let aa_mutations_series = Series::new(PlSmallStr::from("aa_mutations"), aa_mutations);
    let mut df = df.clone();
    df.with_column(aa_mutations_series)?;


    // Reorder columns and sort
    let df = df.select(["nt_mutations", "aa_mutations", "cluster_depth", "total_depth", "frequency", "coverage_start", "coverage_end"])?
        .sort(["cluster_depth", "nt_mutations"], SortMultipleOptions::default().with_order_descending(true))?;
    Ok(df)
}