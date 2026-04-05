#![allow(dead_code, unused_variables)]
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use regex::Regex;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

const MAX_RETRIES: u32 = 10;

/// Download configuration
#[derive(Clone)]
pub struct DownloadConfig {
    pub output_dir: String,
    pub threads: usize,
    pub fna_types: Vec<String>,
    pub overwrite: bool,
    pub verbose: bool,
}

impl Default for DownloadConfig {
    fn default() -> Self {
        DownloadConfig {
            output_dir: ".".to_string(),
            threads: 8,
            fna_types: vec!["genomic".to_string()],
            overwrite: false,
            verbose: false,
        }
    }
}

/// Outcome of processing a single file (per fna_type)
#[derive(Debug, Clone, Copy)]
enum FileOutcome {
    Success { bytes: u64 },
    Skipped,
    Failed,
}

/// Outcome of processing a single assembly record (may include multiple fna_types)
#[derive(Debug, Clone, Copy)]
struct RecordOutcome {
    successes: usize,
    skips: usize,
    failures: usize,
    total_bytes: u64,
}

impl RecordOutcome {
    fn new() -> Self {
        RecordOutcome {
            successes: 0,
            skips: 0,
            failures: 0,
            total_bytes: 0,
        }
    }

    fn add_file_outcome(&mut self, outcome: FileOutcome) {
        match outcome {
            FileOutcome::Success { bytes } => {
                self.successes += 1;
                self.total_bytes += bytes;
            }
            FileOutcome::Skipped => self.skips += 1,
            FileOutcome::Failed => self.failures += 1,
        }
    }
}

/// Build FTP/HTTP download URL for a given assembly and sequence type
fn build_download_url(ftp_path: &str, asm_basename: &str, fna_type: &str) -> String {
    format!("{}/{}_{}.fna.gz", ftp_path, asm_basename, fna_type)
}

/// Download a file with retry logic and resume support.
fn download_file(
    client: &reqwest::blocking::Client,
    url: &str,
    dest_path: &Path,
    config: &DownloadConfig,
) -> std::io::Result<FileOutcome> {
    if dest_path.exists() && !config.overwrite {
        if config.verbose {
            eprintln!("  [SKIP] {} already exists", dest_path.display());
        }
        return Ok(FileOutcome::Skipped);
    }

    if let Some(parent) = dest_path.parent() {
        fs::create_dir_all(parent)?;
    }

    let mut last_err = String::new();

    for attempt in 1..=MAX_RETRIES {
        let resume_from = if dest_path.exists() {
            fs::metadata(dest_path).map(|m| m.len()).unwrap_or(0)
        } else {
            0
        };

        if config.verbose {
            if resume_from > 0 {
                eprintln!(
                    "  [RETRY {}/{}] Resuming {} from {} bytes ...",
                    attempt, MAX_RETRIES, url, resume_from
                );
            } else {
                eprintln!(
                    "  [DOWNLOAD {}/{}] {} -> {}",
                    attempt, MAX_RETRIES, url, dest_path.display()
                );
            }
        }

        let mut req = client.get(url);
        if resume_from > 0 {
            req = req.header("Range", format!("bytes={}-", resume_from));
        }

        match req.send() {
            Ok(resp) => {
                let status = resp.status();

                if resume_from > 0 && status == 200 {
                    if config.verbose {
                        eprintln!("  Server doesn't support Range header, restarting...");
                    }
                    let _ = fs::remove_file(dest_path);
                } else if status.as_u16() >= 400 {
                    last_err = format!("HTTP {}", status);
                    if attempt < MAX_RETRIES {
                        let wait = std::time::Duration::from_secs(2u64.pow(attempt));
                        eprintln!("  Retry after {:?} - error: {}", wait, last_err);
                        thread::sleep(wait);
                    }
                    continue;
                }

                let data = resp
                    .bytes()
                    .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;

                if resume_from > 0 && dest_path.exists() {
                    let mut file = OpenOptions::new().append(true).open(dest_path)?;
                    file.write_all(&data)?;
                } else {
                    let mut file = File::create(dest_path)?;
                    file.write_all(&data)?;
                }

                let bytes_downloaded = data.len() as u64;
                if config.verbose {
                    eprintln!("  OK ({} bytes)", bytes_downloaded);
                }
                return Ok(FileOutcome::Success { bytes: bytes_downloaded });
            }
            Err(e) => {
                last_err = e.to_string();
                if attempt < MAX_RETRIES {
                    let wait = std::time::Duration::from_secs(2u64.pow(attempt));
                    eprintln!(
                        "  Error: {}, retrying after {:?} ({}/{})",
                        last_err, wait, attempt, MAX_RETRIES
                    );
                    thread::sleep(wait);
                }
            }
        }
    }

    Err(std::io::Error::new(
        std::io::ErrorKind::Other,
        format!("Failed after {} attempts: {}", MAX_RETRIES, last_err),
    ))
}

/// Decompress a .gz file to destination, then remove the .gz file
fn decompress_gzip(gz_path: &Path, out_path: &Path, config: &DownloadConfig) -> std::io::Result<()> {
    if out_path.exists() && !config.overwrite {
        if config.verbose {
            eprintln!("  [SKIP] Decompressed file exists: {}", out_path.display());
        }
        if gz_path.exists() {
            let _ = fs::remove_file(gz_path);
        }
        return Ok(());
    }

    if config.verbose {
        eprintln!("  Gunzipping {} -> {}", gz_path.display(), out_path.display());
    }

    let gz_file = File::open(gz_path)?;
    let decoder = flate2::read::GzDecoder::new(BufReader::new(gz_file));
    let mut reader = BufReader::new(decoder);
    let mut out_file = File::create(out_path)?;

    let mut buf = [0u8; 64 * 1024];
    loop {
        match reader.read(&mut buf)? {
            0 => break,
            n => out_file.write_all(&buf[..n])?,
        }
    }

    if gz_path.exists() {
        fs::remove_file(gz_path)?;
    }

    Ok(())
}

/// Generate a .map file mapping each fasta sequence accession to taxid + metadata
pub fn generate_map_file(
    fasta_path: &Path,
    taxid: &str,
    description: &str,
    config: &DownloadConfig,
) -> std::io::Result<()> {
    let map_path_str = format!("{}.map", fasta_path.display());
    let map_path = Path::new(&map_path_str);

    if map_path.exists() && !config.overwrite {
        if config.verbose {
            eprintln!("  [SKIP] Map file exists: {}", map_path.display());
        }
        return Ok(());
    }

    if config.verbose {
        eprintln!("  Generating map: {}", map_path.display());
    }

    let fasta_file = File::open(fasta_path)?;
    let reader = BufReader::new(fasta_file);

    let mut map_file = File::create(&map_path)?;
    let re = Regex::new(r"^>(\S+)").unwrap();

    for line in reader.lines().flatten() {
        if let Some(cap) = re.captures(&line) {
            let accession = cap.get(1).unwrap().as_str();
            writeln!(map_file, "{}\t{}\t{}", accession, taxid, description)?;
        }
    }

    if config.verbose {
        eprintln!("  Map written: {}", map_path.display());
    }
    Ok(())
}

/// Process a single assembly record: download .gz, decompress to .fna, generate .map
fn process_record(
    client: &reqwest::blocking::Client,
    record: &crate::assembly::AssemblyRecord,
    domain_dir: &Path,
    config: &DownloadConfig,
) -> std::io::Result<RecordOutcome> {
    let ftp_path = &record.ftp_path;
    let asm_basename = match record.asm_basename() {
        Some(b) => b,
        None => {
            eprintln!("  Warning: no basename for ftp_path '{}'", ftp_path);
            return Ok(RecordOutcome::new()); // empty outcome
        }
    };

    let subdir = match record.download_subdir() {
        Some(d) => d,
        None => {
            eprintln!("  Warning: cannot extract subdir from ftp_path '{}'", ftp_path);
            return Ok(RecordOutcome::new());
        }
    };

    let download_dir = domain_dir.join(&subdir);
    let taxid = &record.taxid;

    // Handle infraspecific_name for organism naming (mirrors Perl logic)
    let mut organism_name = record.organism_name.clone();
    if !record.infraspecific_name.is_empty() {
        let strain = record.infraspecific_name.replace("strain=", "");
        if !organism_name.contains(&strain) && !strain.is_empty() {
            organism_name.push(' ');
            organism_name.push_str(&record.infraspecific_name);
        }
    }
    let desc = format!("{} {}", record.assembly_accession, organism_name);

    let mut record_outcome = RecordOutcome::new();

    for fna_type in &config.fna_types {
        let url = build_download_url(ftp_path, &asm_basename, fna_type);
        let bfname = format!("{}_{}", asm_basename, fna_type);
        let fname = format!("{}.fna", bfname);
        let gz_path = download_dir.join(format!("{}.gz", fname));
        let fna_path = download_dir.join(&fname);

        // Determine if we need to download
        let need_download = !fna_path.exists() || config.overwrite;

        let file_outcome = if need_download {
            // Download .gz
            match download_file(client, &url, &gz_path, config) {
                Ok(outcome) => outcome,
                Err(e) => {
                    eprintln!("  ERROR downloading {}: {}", url, e);
                    FileOutcome::Failed
                }
            }
        } else {
            // File already exists, skip download
            if config.verbose {
                eprintln!("  [SKIP] {} already exists", fna_path.display());
            }
            FileOutcome::Skipped
        };

        // If download succeeded or skipped, proceed with decompression and map generation
        match file_outcome {
            FileOutcome::Success { bytes } => {
                // Decompress .gz → .fna
                if let Err(e) = decompress_gzip(&gz_path, &fna_path, config) {
                    eprintln!("  ERROR decompressing {}: {}", gz_path.display(), e);
                    record_outcome.add_file_outcome(FileOutcome::Failed);
                    continue;
                }
                // Generate .map file
                if let Err(e) = generate_map_file(&fna_path, taxid, &desc, config) {
                    eprintln!("  ERROR generating map for {}: {}", fna_path.display(), e);
                    record_outcome.add_file_outcome(FileOutcome::Failed);
                    continue;
                }
                record_outcome.add_file_outcome(FileOutcome::Success { bytes });
            }
            FileOutcome::Skipped => {
                // File already exists, but we still need to ensure decompressed file exists
                // If gz file exists, decompress (maybe partial download)
                if gz_path.exists() {
                    if let Err(e) = decompress_gzip(&gz_path, &fna_path, config) {
                        eprintln!("  ERROR decompressing existing {}: {}", gz_path.display(), e);
                        record_outcome.add_file_outcome(FileOutcome::Failed);
                        continue;
                    }
                }
                // Generate map if missing
                let map_path_str = format!("{}.map", fna_path.display());
                let map_path = Path::new(&map_path_str);
                if !map_path.exists() || config.overwrite {
                    if let Err(e) = generate_map_file(&fna_path, taxid, &desc, config) {
                        eprintln!("  ERROR generating map for {}: {}", fna_path.display(), e);
                        record_outcome.add_file_outcome(FileOutcome::Failed);
                        continue;
                    }
                }
                record_outcome.add_file_outcome(FileOutcome::Skipped);
            }
            FileOutcome::Failed => {
                record_outcome.add_file_outcome(FileOutcome::Failed);
            }
        }
    }

    Ok(record_outcome)
}

/// Global statistics for the entire download session
#[derive(Debug, Default)]
struct SessionStats {
    total_files: usize,
    successful_files: usize,
    skipped_files: usize,
    failed_files: usize,
    total_bytes: u64,
    total_records: usize,
    successful_records: usize,
    failed_records: usize,
}

/// Main entry point: parse assembly records and download them concurrently
pub fn run_download(
    records: Vec<crate::assembly::AssemblyRecord>,
    config: &DownloadConfig,
) -> std::io::Result<(usize, usize)> {
    let output_path = Path::new(&config.output_dir);
    fs::create_dir_all(output_path)?;

    let total_records = records.len();
    eprintln!(
        "Starting download of {} genomes to {} ...",
        total_records, config.output_dir
    );
    eprintln!("Configuration: threads={}, fna_types={}, overwrite={}, verbose={}",
        config.threads,
        config.fna_types.join(","),
        config.overwrite,
        config.verbose
    );

    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(1800)) // 30 minutes timeout for large files
        .build()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;

    let n_threads = config.threads.max(1);

    let (tx, rx) = mpsc::sync_channel::<(usize, crate::assembly::AssemblyRecord)>(n_threads * 2);
    let rx_arc = Arc::new(Mutex::new(rx));

    // Shared statistics
    let stats = Arc::new(Mutex::new(SessionStats {
        total_records,
        ..SessionStats::default()
    }));

    let verbose = config.verbose;
    let total_copy = total_records;

    // Setup progress bars
    let multi_progress = MultiProgress::new();
    let overall_pb = multi_progress.add(ProgressBar::new(total_records as u64));
    overall_pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({percent}%) {msg}")
            .unwrap()
            .progress_chars("##-"),
    );
    overall_pb.set_message("Overall progress");

    // Determine if we should show detailed worker progress bars
    let show_worker_bars = !verbose;
    
    // Create a progress bar for each worker thread (only if not verbose)
    let worker_pbs: Arc<Mutex<Vec<ProgressBar>>> = Arc::new(Mutex::new(Vec::new()));
    if show_worker_bars {
        for i in 0..n_threads {
            let pb = multi_progress.add(ProgressBar::new_spinner());
            pb.set_style(
                ProgressStyle::default_spinner()
                    .template("{spinner} {msg}")
                    .unwrap()
            );
            pb.set_message(format!("Thread {}: idle", i + 1));
            pb.enable_steady_tick(std::time::Duration::from_millis(100));
            worker_pbs.lock().unwrap().push(pb);
        }
    }

    // Spawn worker threads using Arc<Mutex<Receiver>>
    let mut workers: Vec<thread::JoinHandle<()>> = Vec::with_capacity(n_threads);
    for ti in 0..n_threads {
        let rx_arc_clone = Arc::clone(&rx_arc);
        let client_clone = client.clone();
        let out_path_buf = output_path.to_path_buf();
        let cfg = config.clone();
        let stats_clone = Arc::clone(&stats);
        let worker_pbs_clone = Arc::clone(&worker_pbs);
        let overall_pb_clone = overall_pb.clone();
        let show_worker_bars_local = show_worker_bars;
        let handle = thread::spawn(move || {
            let thread_idx = ti;
            // Get worker progress bar if enabled
            let worker_pb_opt = if show_worker_bars_local {
                Some(worker_pbs_clone.lock().unwrap()[thread_idx].clone())
            } else {
                None
            };
            loop {
                let msg = {
                    let lock = rx_arc_clone.lock().unwrap();
                    lock.recv()
                };
                match msg {
                    Ok((idx, rec)) => {
                        // Update worker progress bar if enabled
                        if let Some(ref pb) = worker_pb_opt {
                            pb.set_message(format!("Thread {}: processing {} ({}/{})", 
                                thread_idx + 1, rec.assembly_accession, idx, total_copy));
                        }
                        
                        if verbose {
                            eprintln!(
                                "\n[{}/{}] Processing {}",
                                idx, total_copy, rec.assembly_accession
                            );
                        }
                        
                        let result = process_record(&client_clone, &rec, &out_path_buf, &cfg);
                        
                        // Update overall progress bar after processing
                        overall_pb_clone.inc(1);
                        
                        match result {
                            Ok(record_outcome) => {
                                // Update global statistics
                                let mut stats = stats_clone.lock().unwrap();
                                stats.total_files += record_outcome.successes + record_outcome.skips + record_outcome.failures;
                                stats.successful_files += record_outcome.successes;
                                stats.skipped_files += record_outcome.skips;
                                stats.failed_files += record_outcome.failures;
                                stats.total_bytes += record_outcome.total_bytes;
                                if record_outcome.failures == 0 && (record_outcome.successes > 0 || record_outcome.skips > 0) {
                                    stats.successful_records += 1;
                                } else if record_outcome.failures > 0 {
                                    stats.failed_records += 1;
                                }
                            }
                            Err(e) => {
                                eprintln!("\n  ERROR processing {}: {}", rec.assembly_accession, e);
                                let mut stats = stats_clone.lock().unwrap();
                                stats.failed_records += 1;
                                // Assume all files in this record failed
                                stats.total_files += cfg.fna_types.len();
                                stats.failed_files += cfg.fna_types.len();
                            }
                        }
                        
                        // Reset worker progress bar to idle if enabled
                        if let Some(ref pb) = worker_pb_opt {
                            pb.set_message(format!("Thread {}: idle", thread_idx + 1));
                        }
                    }
                    Err(_) => break, // Channel closed, no more work
                }
            }
            // Mark worker as finished if enabled
            if let Some(pb) = worker_pb_opt {
                pb.set_message(format!("Thread {}: finished", thread_idx + 1));
                pb.finish();
            }
        });
        workers.push(handle);
    }

    // Send all records to workers
    for (i, rec) in records.into_iter().enumerate() {
        if tx.send((i + 1, rec)).is_err() {
            break;
        }
    }
    drop(tx);

    // Wait for all workers to finish
    let mut thread_errors = 0usize;
    for worker in workers {
        if worker.join().is_err() {
            thread_errors += 1;
        }
    }

    // Finish overall progress bar
    overall_pb.finish_with_message("Download completed");

    // Clear worker progress bars if they were created
    if show_worker_bars {
        for pb in worker_pbs.lock().unwrap().iter() {
            pb.finish_and_clear();
        }
    }

    // Collect final statistics
    let final_stats = Arc::try_unwrap(stats).unwrap().into_inner().unwrap();

    // Print detailed summary
    eprintln!("\n\n=== Download Session Summary ===");
    eprintln!("Total assembly records processed: {}", final_stats.total_records);
    eprintln!("  Successfully processed records: {}", final_stats.successful_records);
    eprintln!("  Failed records:                {}", final_stats.failed_records);
    eprintln!("Total files (including all fna_types): {}", final_stats.total_files);
    eprintln!("  Successfully downloaded files: {}", final_stats.successful_files);
    eprintln!("  Skipped files (already exist): {}", final_stats.skipped_files);
    eprintln!("  Failed files:                  {}", final_stats.failed_files);
    eprintln!("Total data downloaded: {} bytes ({:.2} MB)",
        final_stats.total_bytes,
        final_stats.total_bytes as f64 / (1024.0 * 1024.0)
    );
    if final_stats.total_files > 0 {
        let avg_size = final_stats.total_bytes as f64 / final_stats.successful_files as f64;
        eprintln!("Average file size: {:.2} KB", avg_size / 1024.0);
    }
    if thread_errors > 0 {
        eprintln!("\nWarning: {} worker threads had errors.", thread_errors);
    }

    // Return approximate downloaded and skipped counts (for backward compatibility)
    // downloaded: number of records with at least one successful file download
    // skipped: number of records where all files were skipped
    // This is approximate but maintains compatibility with existing callers
    let downloaded = final_stats.successful_records;
    let skipped = final_stats.skipped_files.min(final_stats.total_records); // rough estimate
    Ok((downloaded, skipped))
}
