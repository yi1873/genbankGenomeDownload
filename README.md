# GenBank/RefSeq Genome Downloader (Rust)

A high-performance, concurrent downloader for GenBank and RefSeq genome assemblies from NCBI FTP servers. Written in Rust.

## Features

- **Fast concurrent downloads**: Uses multi-threading to download multiple genomes simultaneously.
- **Resumable downloads**: Supports resume on interruption with retry logic.
- **NCBI directory structure**: Organizes downloaded genomes in the same directory layout as NCBI FTP (`GCA/XXX/XXX/XXX/`).
- **Map file generation**: Creates `.map` files mapping sequence accessions to taxonomy IDs and metadata.
- **Flexible filtering**: Accepts a pre-filtered `assembly_summary.txt` file (NCBI assembly summary format).
- **Configurable sequence types**: Download genomic, RNA, CDS, etc. (default: `genomic`).
- **Progress bars**: Real‑time progress indication with `indicatif`.
- **Verbose logging**: Optional detailed output for debugging.

## Installation

### Prerequisites

- Rust toolchain (stable, ≥ 1.70). Install via [rustup](https://rustup.rs/).

### Build from source

```bash
git clone <repository-url>
cd genbankGenomeDownload
cargo build --release
```

The executable will be placed at `target/release/genbankGenomeDownload`.

## Usage

### Basic command

```bash
genbankGenomeDownload -i assembly_summary.txt -o ./
```

### All options

```
genbankGenomeDownload [OPTIONS] -i <FILE>

OPTIONS:
    -i, --input <FILE>      Input assembly_summary.txt file (pre‑filtered genbank assembly format)
    -o, --output <DIR>      Output directory for downloaded genomes [default: .]
    --threads <NUM>         Number of concurrent download threads [default: 8]
    --fna <TYPES>           Comma‑separated sequence types to download
                            (e.g., genomic,rna,cds_from_genomic) [default: genomic]
    --overwrite             Overwrite existing files instead of skipping them
    -v, --verbose           Verbose output
    -h, --help              Print help
    -V, --version           Print version
```

### Example workflow

1. Obtain an `assembly_summary.txt` file from NCBI (e.g., from `ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt`).
2. Filter it to keep only the assemblies you need (by taxid, assembly level, etc.).
3. Run the downloader:

```bash
genbankGenomeDownload -i select_assembly_summary.txt -o ./  --threads 12 --fna genomic
```

The tool will:
- Parse the input file,
- Create the NCBI‑style directory tree under `./`,
- Download the `.fna.gz` files,
- Decompress them to `.fna`,
- Generate a `.map` file for each FASTA with accession‑to‑taxid mappings.

## Project structure

```
genbankGenomeDownload/
├── Cargo.toml                    # Project metadata and dependencies
├── src/
│   ├── main.rs                   # CLI entry point, argument parsing
│   ├── assembly.rs               # Assembly record parsing (assembly_summary.txt)
│   └── downloader.rs             # Core download logic, threading, retry, map generation
├── GCA/                          Example/downloaded genome directories (NCBI layout)
├── test/
    └── test.genbank              Test data
```

### Key dependencies

- `reqwest` (blocking) – HTTP client for downloading
- `clap` – command‑line argument parsing
- `flate2` – gzip decompression
- `regex` – regular expressions for parsing FASTA headers
- `indicatif` – progress bars and spinners

## Performance

- **Concurrency**: Default 8 threads, adjustable with `--threads`.
- **Retry logic**: Up to 10 retries with exponential backoff.
- **Resume support**: If a `.gz` file exists partially, the download resumes from the last byte.

## Compatibility

- Input file must follow the NCBI `assembly_summary.txt` tab‑delimited format (GenBank or RefSeq).
- Output directory structure mirrors the NCBI FTP path after `/all/` (e.g., `GCA/022/175/585/GCA_022175585.2_ASM2217558v2`).

## License

This project is provided under the terms of the [GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.html) license (same as the original Perl script).

## Author

Xiangzhi Liang (xiang_zhi_@126.com)  
Date: 2026‑04‑04

## Acknowledgments

- Based on `centrifuge‑download` by Florian Breitwieser.
- Uses NCBI’s public FTP servers for genome data.
