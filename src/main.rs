use anyhow::{Context, Result};
use clap::{Arg, Command};
use ndarray::{Array1, Array2};
use phylotree::tree::Tree;
use rayon::prelude::*;
use std::{fs::File, io::{BufRead, BufReader, Write}, path::Path};

fn main() -> Result<()> {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();
    let matches = Command::new("Unweighted_UniFrac")
        .version("0.1.0")
        .about("Fast Unweighted UniFrac")
        .arg(Arg::new("tree")
            .short('t')
            .long("tree")
            .value_name("TREE_FILE")
            .help("Input newick format tree file")
            .required(true))
        .arg(Arg::new("table")
            .short('i')
            .long("input")
            .value_name("TABLE_FILE")
            .help("Input tab-delimited sample-feature table")
            .required(true))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .value_name("OUTPUT_FILE")
            .help("Output file for distance matrix")
            .required(true))
        .get_matches();

    let tree_file = matches.get_one::<String>("tree").unwrap();
    let table_file = matches.get_one::<String>("table").unwrap();
    let output_file = matches.get_one::<String>("output").unwrap();

    // Read the tree
    let tree = Tree::from_file(Path::new(tree_file))?;

    // Read the sample-feature table
    let (taxa_order, sample_names, presence_matrix) = read_sample_table(table_file)?;

    let n_samples = sample_names.len();

    // Compute distance matrix: n_samples x n_samples
    let mut dist_matrix = vec![0.0; n_samples * n_samples];

    for i in 0..n_samples {
        dist_matrix[i * n_samples + i] = 0.0; // distance to itself = 0
        for j in i+1..n_samples {
            let uni = compute_unifrac_for_pair(&tree, &taxa_order, &presence_matrix, i, j)?;
            dist_matrix[i * n_samples + j] = uni;
            dist_matrix[j * n_samples + i] = uni; // symmetric
        }
    }

    // Write output matrix
    write_matrix(&sample_names, &dist_matrix, n_samples, output_file)?;

    Ok(())
}

/// Read the sample-feature table.
/// First line: ignore the first element, subsequent elements are sample names
/// Example:
/// Anything  SampleA  SampleB  SampleC
/// T1        10       0        5
/// T2        0        25       0
/// ...
///
/// Any value > 0 is converted to 1.0, else 0.0.
fn read_sample_table(filename: &str) -> Result<(Vec<String>, Vec<String>, Vec<Vec<f64>>)> {
    let f = File::open(filename)?;
    let mut lines = BufReader::new(f).lines();

    // First line: parse sample names
    let header = lines.next().context("No header in table")??;
    let mut hdr_split = header.split_whitespace();
    hdr_split.next(); // ignore the first element in the header line
    let sample_names: Vec<String> = hdr_split.map(|s| s.to_string()).collect();

    let mut taxa_order = Vec::new();
    let mut presence_matrix = Vec::new();

    for line in lines {
        let line = line?;
        let mut parts = line.split_whitespace();
        let taxon = parts.next().context("Taxon missing in a line")?.to_string();
        taxa_order.push(taxon);
        let values: Vec<f64> = parts.map(|x| {
            let val: f64 = x.parse().unwrap_or(0.0);
            if val > 0.0 {1.0} else {0.0} // ensure binary
        }).collect();
        presence_matrix.push(values);
    }

    Ok((taxa_order, sample_names, presence_matrix))
}

/// Compute UniFrac for a given pair of samples i,j
fn compute_unifrac_for_pair(tree: &Tree, taxa_order: &[String], presence_matrix: &[Vec<f64>], i: usize, j: usize) -> Result<f64> {
    // Determine which taxa are present in either sample i or j
    let mut present_taxa = Vec::new();
    for (t_idx, taxon) in taxa_order.iter().enumerate() {
        let val_i = presence_matrix[t_idx][i];
        let val_j = presence_matrix[t_idx][j];
        if val_i > 0.0 || val_j > 0.0 {
            present_taxa.push(taxon.clone());
        }
    }

    let mut sub_tree = tree.clone();
    // prune taxa not in present_taxa
    {
        let leaves = sub_tree.get_leaves();
        let set: std::collections::HashSet<_> = present_taxa.iter().cloned().collect();
        for l in leaves {
            let name = sub_tree.get(&l).unwrap().name.clone().unwrap();
            if !set.contains(&name) {
                sub_tree.prune(&l).context("Prune failed")?;
            }
        }
    }

    let leaves = sub_tree.get_leaves();
    let mut leaf_order = vec![0; sub_tree.size()];
    let mut leaf_names = Vec::new();
    for (l_ord, l_idx) in leaves.into_iter().enumerate() {
        leaf_order[l_idx] = l_ord;
        leaf_names.push(sub_tree.get(&l_idx).unwrap().name.clone().unwrap());
    }

    let (mat_b, brlens) = construct_b(&sub_tree, &leaf_order)?;

    let p_a = get_sample_vec(&mat_b, &presence_matrix, &taxa_order, &leaf_names, i)?;
    let p_b = get_sample_vec(&mat_b, &presence_matrix, &taxa_order, &leaf_names, j)?;

    let sum_shared = parallel_elementwise_sum(&p_a, &p_b, &brlens);
    let l_total = brlens.sum();
    let unifrac = 1.0 - (sum_shared / l_total);

    Ok(unifrac)
}

/// Construct B and brlens
fn construct_b(tree: &Tree, leaf_order: &[usize]) -> Result<(Array2<u8>, Array1<f64>)> {
    let n_tips = tree.n_leaves();
    let n_branches = tree.size();
    let root = tree.get_root()?;

    let mut mat_b = Array2::<u8>::zeros((n_branches, n_tips));
    let mut brlens = Array1::zeros(n_branches);
    for idx in tree.postorder(&root)? {
        let node = tree.get(&idx)?;
        brlens[idx] = node.parent_edge.unwrap_or_default();
        if node.is_tip() {
            let t_ord = leaf_order[idx];
            mat_b[(idx, t_ord)] = 1;
        } else {
            for c in node.children.iter() {
                let merged = &mat_b.row(idx) + &mat_b.row(*c);
                mat_b.row_mut(idx).assign(&merged);
            }
        }
    }

    Ok((mat_b, brlens))
}

/// Construct p_a (or p_b) for a given sample index
fn get_sample_vec(
    mat: &Array2<u8>,
    presence_matrix: &[Vec<f64>],
    taxa_order: &[String],
    leaf_names: &[String],
    sample_idx: usize
) -> Result<Array1<f64>> {
    let s = mat.shape();
    let mut p: Array1<f64> = Array1::zeros(s[0]);

    // For each leaf_name, find its taxon index in taxa_order, check presence in sample_idx
    for (col, lname) in leaf_names.iter().enumerate() {
        let t_idx = taxa_order.iter().position(|x| x == lname).unwrap();
        let val = presence_matrix[t_idx][sample_idx];
        if val > 0.0 {
            // Convert u8 to f64 before addition
            p = &p + &mat.column(col).mapv(|x| x as f64);
        }
    }

    // clamp to 0,1
    Ok(p.mapv(|v: f64| if v > 0.0 {1.0} else {0.0}))
}


/// Parallelize the element-wise multiply and sum (p_a * p_b * brlens)
fn parallel_elementwise_sum(p_a: &Array1<f64>, p_b: &Array1<f64>, brlens: &Array1<f64>) -> f64 {
    let data = p_a
        .iter()
        .zip(p_b.iter())
        .zip(brlens.iter())
        .map(|((pa, pb), b)| (*pa, *pb, *b))
        .collect::<Vec<_>>();

    data.par_iter()
        .map(|(pa, pb, b)| pa * pb * b)
        .sum()
}

/// Write the resulting matrix to a file
fn write_matrix(sample_names: &[String], dist_matrix: &[f64], n: usize, output_file: &str) -> Result<()> {
    let mut file = File::create(output_file)?;
    // Write header
    // The first column header is sample names as well
    write!(file, "Sample")?;
    for sn in sample_names {
        write!(file, "\t{}", sn)?;
    }
    writeln!(file)?;

    for i in 0..n {
        write!(file, "{}", sample_names[i])?;
        for j in 0..n {
            write!(file, "\t{:.6}", dist_matrix[i * n + j])?;
        }
        writeln!(file)?;
    }

    Ok(())
}
