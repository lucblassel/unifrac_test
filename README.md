# UniFrac implememtation in Rust

This is an example repo to show how to compute the [UniFrac](https://en.wikipedia.org/wiki/UniFrac) distance between a pair of samples containing taxa. 
It uses the [phylotree](https://github.com/lucblassel/phylotree-rs) crate to parse the tree file and feature-sample table (OTU table for example) and then compute the unifrac distance.

## Install
```bash
git clone https://github.com/jianshu93/unifrac-rs-test
cd unifrac-rs-test
cargo build --release
./target/release/unifrac -h
```

## Usage 
```bash
 ************** initializing logger *****************

Fast Unweighted UniFrac

Usage: unifrac --tree <TREE_FILE> --input <TABLE_FILE> --output <OUTPUT_FILE>

Options:
  -t, --tree <TREE_FILE>      Input newick format tree file
  -i, --input <TABLE_FILE>    Input tab-delimited sample-feature table
  -o, --output <OUTPUT_FILE>  Output file for distance matrix
  -h, --help                  Print help
  -V, --version               Print version
```

### example
```bash
### remove bootstrap support first if you have it

### Then run unifrac like this:
unifrac -t data/test_rot_new2.nwk -i data/table.txt -o try.txt
cat try.txt
```

## References
1.Lozupone, C. and Knight, R., 2005. UniFrac: a new phylogenetic method for comparing microbial communities. Applied and environmental microbiology, 71(12), pp.8228-8235.

2.Hamady, M., Lozupone, C. and Knight, R., 2010. Fast UniFrac: facilitating high-throughput phylogenetic analyses of microbial communities including analysis of pyrosequencing and PhyloChip data. The ISME journal, 4(1), pp.17-27.
