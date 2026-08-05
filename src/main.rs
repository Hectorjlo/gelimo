use clap::Parser;

/// gelimo: Generalized Phylogenetic Linear Models CLI
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to Newick phylogenetic tree file
    #[arg(short, long)]
    tree: Option<String>,

    /// Path to CSV trait data file
    #[arg(short, long)]
    data: Option<String>,
}

fn main() {
    let _args = Args::parse();
    println!("gelimo v0.1.0 - Generalized Phylogenetic Linear Models Engine");
}
