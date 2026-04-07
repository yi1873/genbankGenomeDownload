#![allow(dead_code, unused_variables)]
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Error, ErrorKind, Read, Write};
use std::path::Path;
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::time::Duration;
use regex::Regex;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use md5;

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

/// Download a file with retry, resume support, real-time progress bar, and optional MD5 verification.
fn download_file(
    client: &reqwest::blocking::Client,
    url: &str,
    dest_path: &Path,
    config: &DownloadConfig,
    expected_md5: Option<String>,
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

    let total_size = get_content_length(client, url).ok();
    let fname = dest_path.file_name().unwrap_or_default().to_string_lossy().to_string();

    // Create progress bar only in verbose mode
    let pb = if config.verbose {
        let bar = match total_size {
            Some(size) => ProgressBar::new(size),
            None => ProgressBar::new_spinner(),
        };
        bar.set_style(match total_size {
            Some(_) => ProgressStyle::default_bar()
                .template("{msg}\n{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) @ {bytes_per_sec}")
                .unwrap()
                .progress_chars("##-"),
            None => ProgressStyle::default_spinner()
                .template("{msg}\n{spinner:.green} [{elapsed_precise}] {bytes} @ {bytes_per_sec}")
                .unwrap(),
        });
        bar.set_message(format!("{}", fname));
        if total_size.is_none() {
            bar.enable_steady_tick(Duration::from_millis(100));
        }
        Some(bar)
    } else {
        None
    };

    // Use parallel download for large files when threads > 1
    const PARALLEL_THRESHOLD: u64 = 100 * 1024 * 1024;
    if config.threads > 1 {
        if let Some(size) = total_size {
            if size >= PARALLEL_THRESHOLD {
                let max_chunks = config.threads.min(16);
                let chunks = ((size / (50 * 1024 * 1024)) as usize + 1).min(max_chunks).max(2);
                return download_file_parallel(client, url, dest_path, config, chunks, expected_md5);
            }
        }
    }

    // Sequential download with retry and resume
    let mut last_err = String::new();

    for attempt in 1..=MAX_RETRIES {
        let resume_from = if dest_path.exists() {
            fs::metadata(dest_path).map(|m| m.len()).unwrap_or(0)
        } else {
            0
        };

        if config.verbose {
            if resume_from > 0 {
                eprintln!("  [RETRY {}/{}] Resuming from {} bytes", attempt, MAX_RETRIES, resume_from);
            } else {
                eprintln!("  [DOWNLOAD {}/{}] {}", attempt, MAX_RETRIES, url);
            }
        }

        let mut req = client.get(url);
        if resume_from > 0 {
            req = req.header("Range", format!("bytes={}-", resume_from));
        }

        match req.send() {
            Ok(mut resp) => {
                let status = resp.status();

                if resume_from > 0 && status == 200 {
                    if config.verbose {
                        eprintln!("  Server doesn't support Range header, restarting...");
                    }
                    let _ = fs::remove_file(dest_path);
                } else if status.as_u16() >= 400 {
                    last_err = format!("HTTP {}", status);
                    if attempt < MAX_RETRIES {
                        thread::sleep(Duration::from_secs(2u64.pow(attempt - 1)));
                    }
                    continue;
                }

                if let Some(ref p) = pb {
                    p.set_position(resume_from);
                }

                let bytes_downloaded = match &pb {
                    Some(p) => {
                        let mut reader = p.wrap_read(resp);
                        if resume_from > 0 && dest_path.exists() {
                            let mut file = BufWriter::new(OpenOptions::new().append(true).open(dest_path)?);
                            io_copy(&mut reader, &mut file).map_err(|e| Error::new(ErrorKind::Other, e.to_string()))?
                        } else {
                            let mut file = BufWriter::new(File::create(dest_path)?);
                            io_copy(&mut reader, &mut file).map_err(|e| Error::new(ErrorKind::Other, e.to_string()))?
                        }
                    }
                    None => {
                        if resume_from > 0 && dest_path.exists() {
                            let mut file = BufWriter::new(OpenOptions::new().append(true).open(dest_path)?);
                            std::io::copy(&mut resp, &mut file).map_err(|e| Error::new(ErrorKind::Other, e.to_string()))?
                        } else {
                            let mut file = BufWriter::new(File::create(dest_path)?);
                            std::io::copy(&mut resp, &mut file).map_err(|e| Error::new(ErrorKind::Other, e.to_string()))?
                        }
                    }
                };

                if let Some(ref p) = pb {
                    p.finish_and_clear();
                }
                
                // Verify MD5 if expected checksum is provided
                if let Some(ref expected) = expected_md5 {
                    if config.verbose {
                        eprintln!("  Verifying MD5 checksum...");
                    }
                    match calculate_md5(dest_path) {
                        Ok(actual_md5) => {
                            if &actual_md5 != expected {
                                if config.verbose {
                                    eprintln!("  MD5 mismatch: expected {}, got {}", expected, actual_md5);
                                }
                                let _ = fs::remove_file(dest_path);
                                last_err = format!("MD5 checksum mismatch: expected {}, got {}", expected, actual_md5);
                                if attempt < MAX_RETRIES {
                                    thread::sleep(Duration::from_secs(2u64.pow(attempt - 1)));
                                }
                                continue;
                            } else {
                                if config.verbose {
                                    eprintln!("  MD5 verification passed");
                                }
                            }
                        }
                        Err(e) => {
                            if config.verbose {
                                eprintln!("  Failed to calculate MD5: {}", e);
                            }
                            let _ = fs::remove_file(dest_path);
                            last_err = format!("Failed to calculate MD5: {}", e);
                            if attempt < MAX_RETRIES {
                                thread::sleep(Duration::from_secs(2u64.pow(attempt - 1)));
                            }
                            continue;
                        }
                    }
                }
                
                if config.verbose {
                    eprintln!("  OK ({} bytes)", bytes_downloaded);
                }
                return Ok(FileOutcome::Success { bytes: bytes_downloaded });
            }
            Err(e) => {
                last_err = e.to_string();
                if attempt < MAX_RETRIES {
                    thread::sleep(Duration::from_secs(2u64.pow(attempt - 1)));
                }
            }
        }
    }

    Err(Error::new(ErrorKind::Other, format!("Failed after {} attempts: {}", MAX_RETRIES, last_err)))
}

// Shorthand to avoid verbose type annotation
fn io_copy<R: Read, W: Write>(r: &mut R, w: &mut W) -> io::Result<u64> {
    std::io::copy(r, w)
}

/// Calculate MD5 checksum of a file
fn calculate_md5(file_path: &Path) -> std::io::Result<String> {
    let mut file = File::open(file_path)?;
    let mut context = md5::Context::new();
    let mut buffer = [0u8; 8192];
    
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        context.consume(&buffer[..n]);
    }
    
    let digest = context.compute();
    Ok(format!("{:x}", digest))
}

/// Parse md5checksums.txt content, returning mapping from filename to md5
fn parse_md5_checksums(content: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        // Format: "md5sum  filename" or "md5sum *filename"
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let md5sum = parts[0].to_string();
            let filename = parts[1..].join(" ");
            // Remove leading '*' if present (some formats use *filename)
            let filename = filename.trim_start_matches('*').to_string();
            // Remove leading './' or '.\' if present
            let filename = filename.trim_start_matches("./").to_string();
            // Also try to get basename in case filename contains path
            let basename = Path::new(&filename)
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or(&filename)
                .to_string();
            // Insert both full filename (without ./) and basename for flexible matching
            map.insert(filename.clone(), md5sum.clone());
            if basename != filename {
                map.insert(basename, md5sum);
            }
        }
    }
    map
}

/// Download and parse md5checksums.txt from a given FTP directory
fn fetch_md5_checksums(
    client: &reqwest::blocking::Client,
    ftp_dir_url: &str,
    config: &DownloadConfig,
) -> std::io::Result<std::collections::HashMap<String, String>> {
    let md5_url = format!("{}/md5checksums.txt", ftp_dir_url);
    if config.verbose {
        eprintln!("  Fetching MD5 checksums from: {}", md5_url);
    }
    
    let response = client.get(&md5_url).send()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
    
    if !response.status().is_success() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("Failed to download md5checksums.txt: HTTP {}", response.status())
        ));
    }
    
    let content = response.text()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
    
    let map = parse_md5_checksums(&content);
    if config.verbose {
        eprintln!("  Fetched {} MD5 checksums from {}", map.len(), md5_url);
        // Log only .fna.gz related entries to avoid clutter
        let fna_entries: Vec<(&String, &String)> = map.iter()
            .filter(|(filename, _)| filename.contains(".fna.gz"))
            .collect();
        if !fna_entries.is_empty() {
            let total_fna_entries = fna_entries.len();
            eprintln!("  Relevant .fna.gz entries ({}):", total_fna_entries);
            let mut count = 0;
            for (filename, md5) in &fna_entries {
                eprintln!("    {} -> {}", filename, md5);
                count += 1;
                if count >= 5 && total_fna_entries > 5 {
                    eprintln!("    ... and {} more", total_fna_entries - 5);
                    break;
                }
            }
        } else {
            eprintln!("  No .fna.gz entries found in MD5 checksums");
        }
    }
    Ok(map)
}

/// Get content length of a remote file via HEAD request
fn get_content_length(client: &reqwest::blocking::Client, url: &str) -> std::io::Result<u64> {
    let resp = client.head(url).send()
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
    if !resp.status().is_success() {
        return Err(std::io::Error::new(std::io::ErrorKind::Other, 
            format!("HEAD request failed with status: {}", resp.status())));
    }
    let len = resp.headers()
        .get(reqwest::header::CONTENT_LENGTH)
        .and_then(|v| v.to_str().ok())
        .and_then(|s| s.parse::<u64>().ok())
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::Other, 
            "Content-Length header missing or invalid"))?;
    Ok(len)
}

/// Download a specific byte range of a remote file with retry logic
fn download_range(
    client: &reqwest::blocking::Client,
    url: &str,
    range: (u64, u64),
    dest_path: &Path,
    config: &DownloadConfig,
) -> std::io::Result<u64> {
    let (start, end) = range;
    let range_header = format!("bytes={}-{}", start, end);
    let mut last_err = String::new();
    
    for attempt in 1..=MAX_RETRIES {
        if config.verbose {
            eprintln!("  [RANGE {}/{}] {} ({})", attempt, MAX_RETRIES, range_header, dest_path.display());
        }
        
        let mut resp = match client.get(url)
            .header("Range", &range_header)
            .send() {
                Ok(r) => r,
                Err(e) => {
                    last_err = e.to_string();
                    if attempt < MAX_RETRIES {
                        let wait = Duration::from_secs(2u64.pow(attempt - 1));
                        eprintln!("  Error: {}, retrying after {:?} ({}/{})", last_err, wait, attempt, MAX_RETRIES);
                        thread::sleep(wait);
                    }
                    continue;
                }
            };

        let status = resp.status();
        if !status.is_success() {
            last_err = format!("HTTP {}", status);
            if attempt < MAX_RETRIES {
                let wait = Duration::from_secs(2u64.pow(attempt - 1));
                eprintln!("  Error: {}, retrying after {:?} ({}/{})", last_err, wait, attempt, MAX_RETRIES);
                thread::sleep(wait);
            }
            continue;
        }

        // Check that server supports range requests
        if status != 206 {
            last_err = format!("Server does not support range requests (status {})", status);
            if attempt < MAX_RETRIES {
                let wait = Duration::from_secs(2u64.pow(attempt - 1));
                eprintln!("  Error: {}, retrying after {:?} ({}/{})", last_err, wait, attempt, MAX_RETRIES);
                thread::sleep(wait);
            }
            continue;
        }
        
        // Stream the response directly to file with buffering
        let mut file = std::io::BufWriter::new(File::create(dest_path)?);
        let bytes_downloaded = std::io::copy(&mut resp, &mut file)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        file.flush()?;
        
        if config.verbose {
            eprintln!("  OK ({} bytes)", bytes_downloaded);
        }
        return Ok(bytes_downloaded);
    }
    
    Err(std::io::Error::new(
        std::io::ErrorKind::Other,
        format!("Failed after {} attempts: {}", MAX_RETRIES, last_err),
    ))
}

/// Download a file in parallel using multiple threads (chunked download)
fn download_file_parallel(
    client: &reqwest::blocking::Client,
    url: &str,
    dest_path: &Path,
    config: &DownloadConfig,
    threads: usize,
    expected_md5: Option<String>,
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

    let total_size = match get_content_length(client, url) {
        Ok(size) => size,
        Err(e) => {
            eprintln!("  WARNING: Could not get content length for {}: {}", url, e);
            return download_file(client, url, dest_path, config, None);
        }
    };

    if total_size < 100 * 1024 * 1024 || threads < 2 {
        return download_file(client, url, dest_path, config, None);
    }

    let fname = dest_path.file_name().unwrap_or_default().to_string_lossy().to_string();

    // Create progress bar only in verbose mode
    let pb = if config.verbose {
        let bar = ProgressBar::new(total_size);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("{msg}\n{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta}) @ {bytes_per_sec}")
                .unwrap()
                .progress_chars("##-"),
        );
        bar.set_message(format!("{}", fname));
        Some(bar)
    } else {
        None
    };

    // Calculate chunks
    let chunk_size = total_size / threads as u64;
    let mut ranges = Vec::with_capacity(threads);
    for i in 0..threads {
        let start = i as u64 * chunk_size;
        let end = if i == threads - 1 { total_size - 1 } else { start + chunk_size - 1 };
        ranges.push((start, end));
    }

    // Create temp dir for chunks
    let temp_dir = dest_path.parent().unwrap().join(format!("{}.parts", fname));
    if temp_dir.exists() { fs::remove_dir_all(&temp_dir)?; }
    fs::create_dir_all(&temp_dir)?;

    // Download all chunks in parallel
    let (tx, rx) = mpsc::channel();
    for (i, range) in ranges.into_iter().enumerate() {
        let tx_c = tx.clone();
        let c = client.clone();
        let u = url.to_string();
        let p = temp_dir.join(format!("chunk_{}.part", i));
        let cfg = config.clone();
        thread::spawn(move || {
            tx_c.send(download_range(&c, &u, range, &p, &cfg)).unwrap();
        });
    }
    drop(tx);

    // Collect results and update progress in real-time
    let mut total_downloaded = 0u64;
    let mut has_error = false;
    for _ in 0..threads {
        match rx.recv().unwrap() {
            Ok(bytes) => {
                total_downloaded += bytes;
                if let Some(ref p) = pb { p.set_position(total_downloaded); }
            }
            Err(e) => { eprintln!("  ERROR downloading chunk: {}", e); has_error = true; }
        }
    }

    if has_error {
        let _ = fs::remove_dir_all(&temp_dir);
        return Err(Error::new(ErrorKind::Other, "One or more chunks failed to download"));
    }

    // Merge chunks into final file
    let mut final_file = File::create(dest_path)?;
    for i in 0..threads {
        let cp = temp_dir.join(format!("chunk_{}.part", i));
        let mut cf = File::open(&cp)?;
        std::io::copy(&mut cf, &mut final_file)?;
        fs::remove_file(&cp)?;
    }
    fs::remove_dir_all(&temp_dir)?;

    if let Some(ref p) = pb { p.finish_and_clear(); }
    
    // Verify MD5 if expected checksum is provided
    if let Some(ref expected) = expected_md5 {
        if config.verbose {
            eprintln!("  Verifying MD5 checksum...");
        }
        match calculate_md5(dest_path) {
            Ok(actual_md5) => {
                if &actual_md5 != expected {
                    let _ = fs::remove_file(dest_path);
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::Other,
                        format!("MD5 checksum mismatch: expected {}, got {}", expected, actual_md5)
                    ));
                } else {
                    if config.verbose {
                        eprintln!("  MD5 verification passed");
                    }
                }
            }
            Err(e) => {
                let _ = fs::remove_file(dest_path);
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Failed to calculate MD5: {}", e)
                ));
            }
        }
    }
    
    if config.verbose {
        eprintln!("  [PARALLEL] Completed: {} bytes", total_downloaded);
    }
    Ok(FileOutcome::Success { bytes: total_downloaded })
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
            // Count as failure for all fna_types
            let mut outcome = RecordOutcome::new();
            for _ in 0..config.fna_types.len() {
                outcome.add_file_outcome(FileOutcome::Failed);
            }
            return Ok(outcome);
        }
    };

    let subdir = match record.download_subdir() {
        Some(d) => d,
        None => {
            eprintln!("  Warning: cannot extract subdir from ftp_path '{}'", ftp_path);
            // Count as failure for all fna_types
            let mut outcome = RecordOutcome::new();
            for _ in 0..config.fna_types.len() {
                outcome.add_file_outcome(FileOutcome::Failed);
            }
            return Ok(outcome);
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

    // Try to fetch MD5 checksums for this FTP directory
    let md5_map = match fetch_md5_checksums(client, ftp_path, config) {
        Ok(map) => Some(map),
        Err(e) => {
            eprintln!("  Warning: could not fetch MD5 checksums for {}: {}", ftp_path, e);
            eprintln!("  Proceeding without MD5 verification (downloads may be corrupted)");
            None
        }
    };

    for fna_type in &config.fna_types {
        let url = build_download_url(ftp_path, &asm_basename, fna_type);
        let bfname = format!("{}_{}", asm_basename, fna_type);
        let fname = format!("{}.fna", bfname);
        let gz_path = download_dir.join(format!("{}.gz", fname));
        let fna_path = download_dir.join(&fname);

        // Determine expected MD5 checksum for this file
        let expected_md5 = md5_map.as_ref().and_then(|map| {
            let filename = format!("{}.fna.gz", bfname);
            map.get(&filename).cloned()
        });
        if config.verbose {
            if let Some(ref md5) = expected_md5 {
                eprintln!("  MD5 checksum for {}: {}", bfname, md5);
            } else if md5_map.is_some() {
                eprintln!("  WARNING: No MD5 checksum found for {}", bfname);
                // Debug: list only .fna.gz related keys in the map
                let map = md5_map.as_ref().unwrap();
                let fna_keys: Vec<&String> = map.keys()
                    .filter(|k| k.contains(".fna.gz"))
                    .collect();
                eprintln!("  Relevant .fna.gz keys in map ({} out of {} total):", fna_keys.len(), map.len());
                let mut sorted_keys: Vec<&String> = fna_keys;
                sorted_keys.sort();
                for key in sorted_keys.iter().take(10) {
                    eprintln!("    - {}", key);
                }
                if sorted_keys.len() > 10 {
                    eprintln!("    ... and {} more", sorted_keys.len() - 10);
                }
            } else {
                eprintln!("  WARNING: No MD5 checksums available for this directory");
            }
        }
        
        // Determine if we need to download
        let need_download = !fna_path.exists() || config.overwrite;
        
        let file_outcome = if need_download {
            // Check if .gz file already exists with correct MD5
            let mut skip_download = false;
            if gz_path.exists() {
                if let Some(ref expected) = expected_md5 {
                    match calculate_md5(&gz_path) {
                        Ok(actual_md5) => {
                            if &actual_md5 == expected {
                                if config.verbose {
                                    eprintln!("  [SKIP] {} already exists with correct MD5", gz_path.display());
                                }
                                skip_download = true;
                            } else {
                                if config.verbose {
                                    eprintln!("  [MD5 MISMATCH] {}: expected {}, got {}", 
                                        gz_path.display(), expected, actual_md5);
                                }
                                let _ = fs::remove_file(&gz_path);
                            }
                        }
                        Err(e) => {
                            if config.verbose {
                                eprintln!("  [WARN] Failed to calculate MD5 for {}: {}", gz_path.display(), e);
                            }
                        }
                    }
                }
            }
            
            if skip_download {
                // File exists with correct MD5, skip download
                FileOutcome::Skipped
            } else {
                // Download .gz with MD5 verification
                match download_file(client, &url, &gz_path, config, expected_md5) {
                    Ok(outcome) => outcome,
                    Err(e) => {
                        eprintln!("  ERROR downloading {}: {}", url, e);
                        FileOutcome::Failed
                    }
                }
            }
        } else {
            // .fna file already exists, skip download
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
        .tcp_nodelay(true) // Disable Nagle's algorithm for lower latency
        .pool_idle_timeout(std::time::Duration::from_secs(90)) // Keep connections alive longer
        .pool_max_idle_per_host(config.threads.max(4)) // Allow multiple idle connections per host for parallel downloads
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

    let final_stats = Arc::try_unwrap(stats).unwrap().into_inner().unwrap();

    eprintln!("\n=== Download Summary ===");
    eprintln!("  Records processed:   {}", final_stats.total_records);
    eprintln!("  Successful:          {}", final_stats.successful_records);
    eprintln!("  Failed:              {}", final_stats.failed_records);
    eprintln!("  Files downloaded:    {}", final_stats.successful_files);
    eprintln!("  Files skipped:       {}", final_stats.skipped_files);
    eprintln!("  Files failed:        {}", final_stats.failed_files);

    if final_stats.total_bytes > 0 {
        if final_stats.total_bytes >= 1024 * 1024 * 1024 {
            eprintln!(
                "  Total data:          {:.2} GB",
                final_stats.total_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
            );
        } else {
            eprintln!(
                "  Total data:          {:.2} MB",
                final_stats.total_bytes as f64 / (1024.0 * 1024.0)
            );
        }
    }

    Ok((final_stats.successful_records, final_stats.skipped_files))
}
