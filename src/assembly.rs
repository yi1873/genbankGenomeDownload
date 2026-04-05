/// Assembly record parsed from assembly_summary.txt
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AssemblyRecord {
    pub assembly_accession: String,
    pub bioproject: String,
    pub biosample: String,
    pub wgs_master: String,
    pub refseq_category: String,
    pub taxid: String,
    pub species_taxid: String,
    pub organism_name: String,
    pub infraspecific_name: String,
    pub isolate: String,
    pub version_status: String,
    pub assembly_level: String,
    pub release_type: String,
    pub genome_rep: String,
    pub seq_rel_date: String,
    pub asm_name: String,
    pub submitter: String,
    pub gbrs_paired_asm: String,
    pub paired_asm_comp: String,
    pub ftp_path: String,
    // Column indices per NCBI assembly_summary.txt spec
}

impl AssemblyRecord {
    /// Parse a single tab-separated line from assembly_summary.txt
    /// Returns None for comment lines or malformed rows
    pub fn parse_line(line: &str) -> Option<Self> {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return None;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 22 {
            eprintln!("Warning: skipping line with {} columns (expected >= 22)", fields.len());
            return None;
        }

        let ftp_path = fields[19].trim();
        if ftp_path.is_empty() {
            eprintln!("Warning: skipping line with empty ftp_path (assembly accession: {})", fields[0]);
            return None;
        }
        // Ensure ftp_path contains "/all/" which is required for NCBI directory structure
        if !ftp_path.contains("/all/") {
            eprintln!("Warning: ftp_path does not contain '/all/' pattern, may cause download issues: {}", ftp_path);
        }

        // Validate taxid is numeric (optional)
        let taxid = fields[5].trim();
        if !taxid.is_empty() && !taxid.chars().all(|c| c.is_ascii_digit()) {
            eprintln!("Warning: taxid '{}' is not numeric (assembly accession: {})", taxid, fields[0]);
        }

        Some(AssemblyRecord {
            assembly_accession: fields[0].to_string(),
            bioproject: fields[1].to_string(),
            biosample: fields[2].to_string(),
            wgs_master: fields[3].to_string(),
            refseq_category: fields[4].to_string(),
            taxid: taxid.to_string(),
            species_taxid: fields[6].to_string(),
            organism_name: fields[7].to_string(),
            infraspecific_name: fields[8].to_string(),
            isolate: fields[9].to_string(),
            version_status: fields[10].to_string(),
            assembly_level: fields[11].to_string().replace(' ', "_"),
            release_type: fields[12].to_string(),
            genome_rep: fields[13].to_string(),
            seq_rel_date: fields[14].to_string(),
            asm_name: fields[15].to_string(),
            submitter: fields[16].to_string(),
            gbrs_paired_asm: fields[17].to_string(),
            paired_asm_comp: fields[18].to_string(),
            ftp_path: ftp_path.to_string(),
        })
    }

    /// Extract the NCBI-style directory component from ftp_path.
    ///
    /// E.g. ftp_path = ".../all/GCA/022/175/585/GCA_022175585.2_ASM2217558v2"
    ///       => "GCA/022/175/585/GCA_022175585.2_ASM2217558v2"
    pub fn download_subdir(&self) -> Option<String> {
        if self.ftp_path.is_empty() {
            return None;
        }
        // Split on "/all/" to get the part after it
        self.ftp_path.split("/all/").nth(1).map(|s| s.to_string())
    }

    /// Get the base name of the assembly from ftp_path (last path component)
    pub fn asm_basename(&self) -> Option<String> {
        if self.ftp_path.is_empty() {
            return None;
        }
        self.ftp_path.split('/').last().map(|s| s.to_string())
    }
}

/// Parse all valid records from an assembly summary file
pub fn parse_assembly_file(path: &str) -> std::io::Result<Vec<AssemblyRecord>> {
    let content = std::fs::read_to_string(path)?;
    let mut records = Vec::new();
    let mut seen = std::collections::HashSet::new();

    for line in content.lines() {
        if let Some(rec) = AssemblyRecord::parse_line(line) {
            // Skip duplicates by assembly_accession
            if seen.insert(rec.assembly_accession.clone()) {
                records.push(rec);
            }
        }
    }
    Ok(records)
}
