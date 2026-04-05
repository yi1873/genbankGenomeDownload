//! GenBank/RefSeq Genome Downloader - Rust implementation
//!
//! Downloads genome FASTA files from NCBI based on a local pre-filtered
//! assembly_summary.txt file (genbank format). Organizes output in NCBI
//! FTP-style directory structure and generates .map files.
//!
//! Author: Xiangzhi Liang (xiang_zhi_@126.com)
//! Date: 2026-04-04
//!
//! This script follows bioinformatics guidelines for strict validation,
//! explicit failure, and reproducibility.

#![allow(dead_code, unused_variables)]
mod assembly;
mod downloader;

use clap::Parser;
use downloader::DownloadConfig;
#[derive(Parser, Debug)]
#[command(name = "genbankGenomeDownload")]
#[command(version, about)]
struct Args {
    /// Input assembly_summary.txt file (pre-filtered genbank assembly format)
    #[arg(short, long, value_name = "FILE")]
    input: String,

    /// Output directory for downloaded genomes
    #[arg(short, long, value_name = "DIR", default_value = ".")]
    output: String,

    /// Number of concurrent download threads
    #[arg(long, default_value_t = 8)]
    threads: usize,

    /// Comma-separated sequence types to download (e.g., genomic,rna,cds_from_genomic)
    #[arg(long, default_value = "genomic")]
    fna: String,

    /// Overwrite existing files instead of skipping them
    #[arg(long)]
    overwrite: bool,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() {
    let args = Args::parse();

    // Parse fna types from comma-separated string
    let fna_types: Vec<String> = args
        .fna
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();

    let config = DownloadConfig {
        output_dir: args.output.clone(),
        threads: args.threads,
        fna_types,
        overwrite: args.overwrite,
        verbose: args.verbose,
    };

    eprintln!("Reading assembly file: {}", args.input);

    let records = match assembly::parse_assembly_file(&args.input) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error reading '{}': {}", args.input, e);
            std::process::exit(1);
        }
    };

    if records.is_empty() {
        eprintln!("No valid assembly records found in '{}'.", args.input);
        eprintln!("Ensure the file is in NCBI assembly_summary.txt tab-delimited format.");
        std::process::exit(1);
    }

    eprintln!(
        "Parsed {} unique assembly records. Output directory: {}\n",
        records.len(),
        config.output_dir
    );

    if config.verbose {
        for rec in &records {
            eprintln!(
                "  {} | {} | {} | level={}",
                rec.assembly_accession,
                rec.organism_name,
                rec.taxid,
                rec.assembly_level
            );
        }
        eprintln!();
    }

    match downloader::run_download(records, &config) {
        Ok((downloaded, skipped)) => {
            // Detailed summary is already printed by run_download
            eprintln!("\n=== Quick Summary ===");
            eprintln!("  Successfully processed genomes: {}", downloaded);
            eprintln!("  Skipped files (already exist):  {}", skipped);
        }
        Err(e) => {
            eprintln!("\nFatal error during download: {}", e);
            std::process::exit(1);
        }
    }
}
