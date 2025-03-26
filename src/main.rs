use chrono::prelude::*;
use clap::Parser;
use csv;
use flate2::read::MultiGzDecoder;
use indicatif::{ProgressBar, ProgressStyle};
use num_cpus;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use noodles::vcf::variant::record::info::field::Value;
use noodles::vcf::io::indexed_reader::IndexedReader;
use noodles::vcf::variant::record::AlternateBases;
use noodles::core::region::Region;
use noodles::csi;

/// Command-line arguments
#[derive(Parser)]
#[command(
    name = "pathogenic",
    about = "Identify pathogenic variants from a VCF using official ClinVar plus 1000 Genomes frequencies (GRCh37 or GRCh38)"
)]
struct Args {
    /// Genome build: must be "GRCh37" (hg19) or "GRCh38" (hg38)
    #[arg(short, long)]
    build: String,

    /// Path to the input VCF (can be uncompressed or gzipped)
    #[arg(short, long, value_name = "FILE")]
    input: PathBuf,
    
    /// Include variants of uncertain significance (VUS)
    #[arg(short = 'v', long, visible_alias = "vus")]
    include_vus: bool,
    
    /// Include benign variants
    #[arg(short = 'n', long, visible_alias = "benign")]
    include_benign: bool,
    
    /// Generate markdown report (enabled by default, use --markdown-report=false to disable)
    #[arg(long, visible_alias = "md-report", default_value_t = true, action = clap::ArgAction::Set)]
    markdown_report: bool,
}

/// Custom error type for downloads
#[derive(Debug)]
enum DownloadError {
    Io(std::io::Error),
    Reqwest(reqwest::Error),
}

impl fmt::Display for DownloadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DownloadError::Io(e) => write!(f, "IO error: {e}"),
            DownloadError::Reqwest(e) => write!(f, "Reqwest error: {e}"),
        }
    }
}
impl Error for DownloadError {}
impl From<std::io::Error> for DownloadError {
    fn from(e: std::io::Error) -> Self {
        DownloadError::Io(e)
    }
}
impl From<reqwest::Error> for DownloadError {
    fn from(e: reqwest::Error) -> Self {
        DownloadError::Reqwest(e)
    }
}

/// Download a remote file with a progress bar using parallel chunk downloads,
/// if the server supports HTTP range requests. Otherwise, fall back to a single-threaded download.
fn download_file(
    url: &str,
    out_path: &Path,
    log_file: &mut File,
) -> Result<(), DownloadError> {
    use std::sync::{
        atomic::{AtomicU64, Ordering},
        Arc,
    };
    use std::thread;

    println!("  -> Starting download from {url}");
    writeln!(log_file, "  -> Starting download from {url}")?;

    let client = reqwest::blocking::Client::new();
    let head_resp = client.head(url).send()?;
    let total_size = head_resp
        .headers()
        .get(reqwest::header::CONTENT_LENGTH)
        .and_then(|s| s.to_str().ok())
        .and_then(|s| s.parse::<u64>().ok())
        .unwrap_or(0);
    let accept_ranges = head_resp
        .headers()
        .get(reqwest::header::ACCEPT_RANGES)
        .and_then(|s| s.to_str().ok())
        .unwrap_or("");

    if total_size == 0 || accept_ranges != "bytes" {
        println!("  -> Server does not support parallel downloads; falling back to single-threaded download.");
        writeln!(log_file, "  -> Server does not support parallel downloads; falling back to single-threaded download.")?;
        let mut response = client.get(url).send()?;
        let pb = ProgressBar::new(total_size);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
                .unwrap()
                .progress_chars("=>-"),
        );
        let mut file = File::create(out_path)?;
        let mut downloaded: u64 = 0;
        let mut buffer = [0u8; 8192];
        loop {
            let bytes_read = response.read(&mut buffer)?;
            if bytes_read == 0 { break; }
            file.write_all(&buffer[..bytes_read])?;
            downloaded += bytes_read as u64;
            pb.set_position(downloaded);
        }
        pb.finish_with_message("Download complete");
        return Ok(());
    }

    let num_chunks = num_cpus::get();
    let chunk_size = total_size / num_chunks as u64;
    let mut ranges = Vec::new();
    for i in 0..num_chunks {
        let start = i as u64 * chunk_size;
        let end = if i == num_chunks - 1 { total_size - 1 } else { (i as u64 + 1) * chunk_size - 1 };
        ranges.push((start, end));
    }

    let progress = Arc::new(AtomicU64::new(0));
    let mut handles = Vec::new();

    for (i, (start, end)) in ranges.into_iter().enumerate() {
        let client = client.clone();
        let url = url.to_string();
        let progress = Arc::clone(&progress);
        let handle = thread::spawn(move || -> Result<(usize, Vec<u8>), DownloadError> {
            let range_header = format!("bytes={}-{}", start, end);
            let mut resp = client.get(&url)
                .header(reqwest::header::RANGE, range_header)
                .send()?;
            let mut buf = Vec::new();
            let mut local_buf = [0u8; 8192];
            loop {
                let n = resp.read(&mut local_buf)?;
                if n == 0 { break; }
                buf.extend_from_slice(&local_buf[..n]);
                progress.fetch_add(n as u64, Ordering::Relaxed);
            }
            Ok((i, buf))
        });
        handles.push(handle);
    }

    let pb = ProgressBar::new(total_size);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{total_bytes} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    while progress.load(Ordering::Relaxed) < total_size {
        pb.set_position(progress.load(Ordering::Relaxed));
        std::thread::sleep(std::time::Duration::from_millis(100));
    }
    pb.finish_with_message("Download complete");

    let mut chunks: Vec<(usize, Vec<u8>)> = Vec::with_capacity(handles.len());
    for handle in handles {
        let res = handle.join().map_err(|_| {
            DownloadError::Io(std::io::Error::new(std::io::ErrorKind::Other, "Thread join error"))
        })??;
        chunks.push(res);
    }
    chunks.sort_by_key(|(i, _)| *i);
    let mut file = File::create(out_path)?;
    for (_i, data) in chunks {
        file.write_all(&data)?;
    }
    Ok(())
}

/// Parse the semicolon-delimited INFO field into a HashMap
fn parse_info_field(info_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for item in info_str.split(';') {
        if item.is_empty() {
            continue;
        }
        let mut kv = item.splitn(2, '=');
        let key = kv.next().unwrap_or("").to_string();
        let val = kv.next().unwrap_or("").to_string();
        map.insert(key, val);
    }
    map
}

/// Convert ClinVar review status to a star rating
fn review_status_to_stars(status: Option<&str>) -> u8 {
    match status {
        Some(s) if s.contains("practice guideline") => 4,
        Some(s) if s.contains("reviewed by expert panel") => 3,
        Some(s) if s.contains("criteria provided, multiple submitters, no conflicts") => 2,
        Some(s) if s.contains("criteria provided, single submitter") => 1,
        _ => 0,
    }
}


/// Custom error to unify reading different VCFs (ClinVar or 1000G)
enum VcfReadError {
    Io(std::io::Error),
    Parse(std::num::ParseIntError),
}

impl From<std::io::Error> for VcfReadError {
    fn from(e: std::io::Error) -> Self {
        VcfReadError::Io(e)
    }
}

impl From<std::num::ParseIntError> for VcfReadError {
    fn from(e: std::num::ParseIntError) -> Self {
        VcfReadError::Parse(e)
    }
}

impl fmt::Display for VcfReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VcfReadError::Io(e) => write!(f, "IO error during VCF read: {}", e),
            VcfReadError::Parse(e) => write!(f, "Parse error during VCF read: {}", e),
        }
    }
}

/// A specialized type for records from the ClinVar VCF
#[derive(Debug, Clone)]
struct ClinVarRecord {
    chr: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clnsig: String,
    is_alt_pathogenic: bool,
    gene: Option<String>,
    allele_id: Option<i32>,
    clnrevstat: Option<String>,
    clndn: Option<String>,
    af_esp: Option<f64>,
    af_exac: Option<f64>,
    af_tgp: Option<f64>,
    clnsig_category: String, // "pathogenic", "benign", "vus", or "conflicting"
}

/// Container for ClinVar variants keyed by (chr, pos, ref, alt)
type ClinVarMap = HashMap<(String, u32, String, String), ClinVarRecord>;

/// Parse a single line from the ClinVar VCF
fn parse_clinvar_line(line: &str, include_vus: bool, include_benign: bool) -> Option<Vec<ClinVarRecord>> {
    if line.starts_with('#') || line.trim().is_empty() {
        return None;
    }

    let mut fields = line.split('\t');
    let chrom = fields.next()?.to_string();
    let pos_str = fields.next()?;
    let _ = fields.next(); // ID
    let ref_allele = fields.next()?;
    let alt_allele = fields.next()?;
    let _ = fields.next(); // QUAL
    let _ = fields.next(); // FILTER
    let info_str = fields.next()?;

    let pos_num = pos_str.parse::<u32>().ok()?;

    let chr_norm = if chrom.eq_ignore_ascii_case("chrM") || chrom.eq_ignore_ascii_case("chrMT") {
        "MT".to_string()
    } else if let Some(stripped) = chrom.strip_prefix("chr") {
        stripped.to_string()
    } else {
        chrom.clone()
    };

    let alt_list: Vec<&str> = alt_allele.split(',').collect();
    let info_map = parse_info_field(info_str);

    let clnsig_opt = info_map.get("CLNSIG");
    if clnsig_opt.is_none() {
        return None;
    }
    let clnsig_str = clnsig_opt.unwrap();

    // Filter variants based on clinical significance and command-line flags
    let contains_pathogenic = clnsig_str.contains("Pathogenic") || clnsig_str.contains("Likely_pathogenic");
    let contains_benign = clnsig_str.contains("Benign") || clnsig_str.contains("Likely_benign");
    let contains_vus = clnsig_str.contains("Uncertain_significance");
    let contains_conflicting = clnsig_str.contains("Conflicting_interpretations_of_pathogenicity");
    
    // Determine the variant classification category
    let clnsig_category = if contains_pathogenic && !contains_benign {
        "pathogenic"
    } else if contains_benign && !contains_pathogenic {
        "benign"
    } else if contains_vus {
        "vus"
    } else if contains_conflicting {
        "conflicting"
    } else {
        "other"
    };
    
    // Skip variants that don't match our inclusion criteria
    if !((contains_pathogenic) || 
         (include_vus && (contains_vus || contains_conflicting)) || 
         (include_benign && contains_benign))
    {
        return None;
    }

    // Determine if the variant is pathogenic for filtering in the match phase
    let is_alt_pathogenic = contains_pathogenic && !(contains_benign);

    let gene_opt = info_map
        .get("GENEINFO")
        .map(|g| g.split(':').next().unwrap_or(g).to_string());
    let allele_id_opt = info_map.get("ALLELEID").and_then(|a| a.parse::<i32>().ok());
    let clnrevstat = info_map.get("CLNREVSTAT").map(|s| s.to_string());
    let clndn = info_map.get("CLNDN").cloned();
    let af_esp = info_map.get("AF_ESP").and_then(|s| s.parse::<f64>().ok());
    let af_exac = info_map.get("AF_EXAC").and_then(|s| s.parse::<f64>().ok());
    let af_tgp = info_map.get("AF_TGP").and_then(|s| s.parse::<f64>().ok());

    let mut recs = Vec::with_capacity(alt_list.len());
    for alt_a in alt_list {
        recs.push(ClinVarRecord {
            chr: chr_norm.clone(),
            pos: pos_num,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_a.to_string(),
            clnsig: clnsig_str.to_string(),
            is_alt_pathogenic,
            gene: gene_opt.clone(),
            allele_id: allele_id_opt,
            clnrevstat: clnrevstat.clone(),
            clndn: clndn.clone(),
            af_esp,
            af_exac,
            af_tgp,
            clnsig_category: clnsig_category.to_string(),
        });
    }
    Some(recs)
}

/// Parse the ClinVar VCF in parallel
fn parse_clinvar_vcf_gz(
    path_gz: &Path,
    include_vus: bool,
    include_benign: bool,
    log_file: &mut File,
) -> Result<(ClinVarMap, String), Box<dyn Error>> {
    println!("[STEP] Parsing ClinVar VCF: {}", path_gz.display());
    writeln!(log_file, "[STEP] Parsing ClinVar VCF: {}", path_gz.display())?;

    // Extract file date from header if present
    let mut file_date = String::new();
    let file = File::open(&path_gz)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut lines = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;

        if line.starts_with("##fileDate=") {
            file_date = line.trim().replace("##fileDate=", "");
        } else if line.starts_with('#') {
            // Skip other header lines
            continue;
        } else {
            // Non-header line, add to batch for parallel processing
            lines.push(line);
        }
    }

    println!("  -> Found {} non-header lines.", lines.len());
    writeln!(log_file, "  -> Found {} non-header lines.", lines.len())?;

    println!("  -> Parsing in parallel with {} threads...", num_cpus::get());
    writeln!(
        log_file,
        "  -> Parsing in parallel with {} threads...",
        num_cpus::get()
    )?;

    let pb = ProgressBar::new(lines.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let chunk_maps: Vec<ClinVarMap> = lines
        .par_iter()
        .map(|line| {
            pb.inc(1);
            match parse_clinvar_line(line, include_vus, include_benign) {
                None => HashMap::new(),
                Some(records) => {
                    let mut local_map = HashMap::with_capacity(records.len());
                    for r in records {
                        let key = (r.chr.clone(), r.pos, r.ref_allele.clone(), r.alt_allele.clone());
                        local_map.insert(key, r);
                    }
                    local_map
                }
            }
        })
        .collect();

    pb.finish_with_message("Parsing complete.");

    println!("  -> Merging {} map chunks...", chunk_maps.len());
    writeln!(log_file, "  -> Merging {} map chunks...", chunk_maps.len())?;

    // Combine maps into a single map
    let mut final_map = HashMap::new();
    for cm in chunk_maps {
        final_map.extend(cm);
    }

    println!("  -> Final ClinVar map size: {}", final_map.len());
    writeln!(log_file, "  -> Final ClinVar map size: {}", final_map.len())?;
    Ok((final_map, file_date))
}

/// A specialized type for user input variants
#[derive(Debug)]
struct InputVariant {
    chr: String,
    pos: u32,
    ref_allele: String,
    alts: Vec<(String, bool)>, // (alt_allele, is_present_in_genotype)
    genotype: String,
}

/// Attempt to parse a single user line from VCF
fn parse_input_line(line: &str) -> Option<(String, InputVariant)> {
    if line.starts_with('#') || line.trim().is_empty() {
        return None;
    }

    let mut fields = line.split('\t');
    let chrom = fields.next()?.to_string();
    let pos_str = fields.next()?;
    let _ = fields.next(); // ID
    let ref_allele = fields.next()?;
    let alt_allele = fields.next()?;
    let _ = fields.next(); // QUAL
    let _ = fields.next(); // FILTER
    let _ = fields.next(); // INFO

    let mut rest_cols = Vec::new();
    for c in fields {
        rest_cols.push(c);
    }

    let pos_num = pos_str.parse::<u32>().ok()?;
    let chr_fixed = if chrom.eq_ignore_ascii_case("chrM") || chrom.eq_ignore_ascii_case("chrMT") {
        "MT".to_string()
    } else if let Some(stripped) = chrom.strip_prefix("chr") {
        stripped.to_string()
    } else {
        chrom.to_string()
    };

    let alt_list: Vec<String> = alt_allele.split(',').map(|s| s.to_string()).collect();
    let mut present_flags = HashSet::new();

    // Attempt to parse genotype from the next columns
    // Usually the first of rest_cols is the FORMAT, second is the sample data
    let genotype = if !rest_cols.is_empty() {
        let format_str = rest_cols[0];
        let format_items: Vec<&str> = format_str.split(':').collect();
        if let Some(gt_index) = format_items.iter().position(|&f| f == "GT") {
            if rest_cols.len() > 1 {
                let sample_str = rest_cols[1];
                let sample_items: Vec<&str> = sample_str.split(':').collect();
                if gt_index < sample_items.len() {
                    sample_items[gt_index].to_string()
                } else {
                    "1/1".to_string() 
                }
            } else {
                "1/1".to_string()
            }
        } else {
            "1/1".to_string()
        }
    } else {
        "1/1".to_string()
    };

    if !genotype.is_empty() {
        let split_gt: Vec<&str> = genotype.split(&['/', '|'][..]).collect();
        for g in split_gt {
            if let Ok(idx) = g.parse::<usize>() {
                if idx >= 1 {
                    present_flags.insert(idx);
                }
            }
        }
    }

    let mut alts = Vec::with_capacity(alt_list.len());
    for (i, alt_a) in alt_list.into_iter().enumerate() {
        let alt_idx = i + 1;
        let is_present = present_flags.contains(&alt_idx);
        alts.push((alt_a, is_present));
    }

    Some((line.to_string(), InputVariant {
        chr: chr_fixed,
        pos: pos_num,
        ref_allele: ref_allele.to_string(),
        alts,
        genotype,
    }))
}

/// Parse the user input VCF (can be gzipped or uncompressed)
fn parse_input_vcf(path: &Path, log_file: &mut File) -> Result<Vec<(String, InputVariant)>, Box<dyn Error>> {
    println!(
        "[STEP] Parsing user input VCF: {}",
        path.display()
    );
    writeln!(
        log_file,
        "[STEP] Parsing user input VCF: {}",
        path.display()
    )?;

    if !path.exists() {
        return Err(format!("User input VCF not found: {}", path.display()).into());
    }

    // Check if path is gzipped by extension
    let mut lines = Vec::new();
    if let Some(ext) = path.extension() {
        if ext == "gz" {
            let file = File::open(path)?;
            let decoder = MultiGzDecoder::new(file);
            let reader = BufReader::new(decoder);
            for l in reader.lines() {
                lines.push(l?);
            }
        } else {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            for l in reader.lines() {
                lines.push(l?);
            }
        }
    } else {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for l in reader.lines() {
            lines.push(l?);
        }
    }

    println!(
        "  -> Loaded user VCF ({} lines). Using {} CPU cores.",
        lines.len(),
        num_cpus::get()
    );
    writeln!(
        log_file,
        "  -> Loaded user VCF ({} lines). Using {} CPU cores.",
        lines.len(),
        num_cpus::get()
    )?;

    let total = lines.len() as u64;
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let chunk_variants: Vec<Vec<(String, InputVariant)>> = lines
        .par_iter()
        .map(|line| {
            pb.inc(1);
            match parse_input_line(line) {
                None => vec![],
                Some((l, iv)) => vec![(l, iv)],
            }
        })
        .collect();

    pb.finish_with_message("User input parse complete.");

    let final_list: Vec<(String, InputVariant)> = chunk_variants.into_iter().flatten().collect();
    println!("  -> Parsed {} variants from user input.", final_list.len());
    writeln!(
        log_file,
        "  -> Parsed {} variants from user input.",
        final_list.len()
    )?;
    Ok(final_list)
}

/// A specialized type for records from the ClinVar pathogenic summary TSV.
/// These fields are used for annotation in the final output.
#[derive(Debug, Clone)]
struct ClinVarTsvRecord {
    molecular_consequence: Option<String>,
    functional_consequence: Option<String>,
    mode_of_inheritance: Option<String>,
    preferred_values: Option<String>,
    citations: Option<String>,
    comments: Option<String>,
    family_data: Option<String>,
    record_status: Option<String>,
    description: Option<String>,
    date_last_evaluated: Option<String>,
}

/// Parse the ClinVar pathogenic summary TSV for deeper annotation
fn parse_clinvar_tsv(
    tsv_path: &Path,
    build: &str,
    log_file: &mut File,
) -> Result<HashMap<(String, u32, String, String), ClinVarTsvRecord>, Box<dyn Error>> {
    println!("[STEP] Parsing ClinVar summary TSV: {}", tsv_path.display());
    writeln!(
        log_file,
        "[STEP] Parsing ClinVar summary TSV: {}",
        tsv_path.display()
    )?;

    let file = File::open(tsv_path)?;
    let reader: Box<dyn BufRead> = if let Some(ext) = tsv_path.extension() {
        if ext == "xz" {
            Box::new(BufReader::new(xz2::read::XzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        }
    } else {
        Box::new(BufReader::new(file))
    };

    let mut line_count = 0;
    let mut map = HashMap::new();
    let mut header_indices: Option<HashMap<String, usize>> = None;

    // We handle GRCh37 or GRCh38 columns
    let want_37 = build == "GRCH37";

    for line_result in reader.lines() {
        let line = line_result?;
        line_count += 1;

        if line_count == 1 {
            // Header row
            let headers: Vec<&str> = line.split('\t').collect();
            let mut idx_map = HashMap::new();
            for (i, h) in headers.iter().enumerate() {
                idx_map.insert(h.to_string(), i);
            }
            header_indices = Some(idx_map);
            continue;
        }
        let headers_map = header_indices.as_ref().ok_or("No header row found")?;
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 28 {
            continue;
        }

        let c_chr = if want_37 {
            headers_map.get("GRCh37_Chromosome")
        } else {
            headers_map.get("GRCh38_Chromosome")
        };
        let c_pos = if want_37 {
            headers_map.get("GRCh37_PositionVCF")
        } else {
            headers_map.get("GRCh38_PositionVCF")
        };
        let c_ref = if want_37 {
            headers_map.get("GRCh37_ReferenceAlleleVCF")
        } else {
            headers_map.get("GRCh38_ReferenceAlleleVCF")
        };
        let c_alt = if want_37 {
            headers_map.get("GRCh37_AlternateAlleleVCF")
        } else {
            headers_map.get("GRCh38_AlternateAlleleVCF")
        };

        if c_chr.is_none() || c_pos.is_none() || c_ref.is_none() || c_alt.is_none() {
            continue;
        }
        let chr_val = cols[*c_chr.unwrap()].trim();
        let pos_val = cols[*c_pos.unwrap()].trim();
        let ref_val = cols[*c_ref.unwrap()].trim();
        let alt_val = cols[*c_alt.unwrap()].trim();

        let pos_parsed = pos_val.parse::<u32>().unwrap_or(0);
        if chr_val.is_empty() || pos_parsed == 0 || ref_val.is_empty() || alt_val.is_empty() {
            continue;
        }

        // Normalize chromosome
        let chr_norm = if chr_val.eq_ignore_ascii_case("chrM")
            || chr_val.eq_ignore_ascii_case("chrMT")
        {
            "MT".to_string()
        } else if let Some(stripped) = chr_val.strip_prefix("chr") {
            stripped.to_string()
        } else {
            chr_val.to_string()
        };

        let rec = ClinVarTsvRecord {
            molecular_consequence: headers_map
                .get("MolecularConsequence")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            functional_consequence: headers_map
                .get("FunctionalConsequence")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            mode_of_inheritance: headers_map
                .get("ModeOfInheritance")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            preferred_values: headers_map
                .get("PreferredValues")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            citations: headers_map
                .get("Citations")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            comments: headers_map
                .get("Comments")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            family_data: headers_map
                .get("FamilyData")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            record_status: headers_map
                .get("RecordStatus")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            description: headers_map
                .get("Description")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
            date_last_evaluated: headers_map
                .get("DateLastEvaluated")
                .and_then(|idx| cols.get(*idx).map(|s| s.to_string()))
                .filter(|s| !s.is_empty()),
        };

        let key = (chr_norm, pos_parsed, ref_val.to_string(), alt_val.to_string());
        map.insert(key, rec);
    }

    println!("  -> Parsed {} lines from summary TSV.", line_count - 1);
    writeln!(
        log_file,
        "  -> Parsed {} lines from summary TSV.",
        line_count - 1
    )?;
    Ok(map)
}

/// Our final output record combining user input, ClinVar VCF, ClinVar TSV, and 1000G frequencies.
#[derive(Debug)]
struct FinalRecord {
    chr: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clnsig: String,
    is_alt_pathogenic: bool,
    gene: Option<String>,
    allele_id: Option<i32>,
    genotype: String,
    review_stars: u8,
    af_esp: Option<f64>,
    af_exac: Option<f64>,
    af_tgp: Option<f64>,
    clndn: Option<String>,
    molecular_consequence: Option<String>,
    functional_consequence: Option<String>,
    mode_of_inheritance: Option<String>,
    preferred_values: Option<String>,
    citations: Option<String>,
    comments: Option<String>,
    family_data: Option<String>,
    record_status: Option<String>,
    description: Option<String>,
    date_last_evaluated: Option<String>,
    af_afr: Option<f64>,
    af_amr: Option<f64>,
    af_eas: Option<f64>,
    af_eur: Option<f64>,
    af_sas: Option<f64>,
    clnsig_category: String,
}

/// Get the numeric order for sorting chromosomes
fn get_chromosome_order(chr: &str) -> usize {
    if chr == "X" { 
        return 23; 
    } else if chr == "Y" { 
        return 24; 
    } else if chr == "MT" { 
        return 25; 
    } else if let Ok(num) = chr.parse::<usize>() {
        return num;
    } else {
        return 100; // Any other chromosome types will be sorted last
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut log_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open("pathogenic.log")?;

    println!("=== Pathogenic Variant Finder ===");
    writeln!(log_file, "=== Pathogenic Variant Finder ===")?;

    let now: DateTime<Utc> = Utc::now();
    println!("[LOG] Timestamp: {}", now.to_rfc3339());
    writeln!(log_file, "[LOG] Timestamp: {}", now.to_rfc3339())?;

    let args = Args::parse();
    let build = args.build.to_uppercase();
    let input_path = args.input.clone();
    let include_vus = args.include_vus;
    let include_benign = args.include_benign;

    println!("[LOG] Genome Build: {}", build);
    writeln!(log_file, "[LOG] Genome Build: {}", build)?;
    println!("[LOG] Input File: {}", input_path.display());
    writeln!(log_file, "[LOG] Input File: {}", input_path.display())?;
    println!("[LOG] Include VUS: {}", include_vus);
    writeln!(log_file, "[LOG] Include VUS: {}", include_vus)?;
    println!("[LOG] Include Benign: {}", include_benign);
    writeln!(log_file, "[LOG] Include Benign: {}", include_benign)?;

    println!("[STEP] Checking arguments...");
    writeln!(log_file, "[STEP] Checking arguments...")?;
    if build != "GRCH37" && build != "GRCH38" {
        return Err("Genome build must be 'GRCh37' or 'GRCh38'".into());
    }
    if !input_path.exists() {
        return Err(format!("Input VCF not found: {}", input_path.display()).into());
    }

    // Create reports directory if it doesn't exist
    let reports_dir = Path::new("reports");
    fs::create_dir_all(reports_dir)?;
    println!("[STEP] Ensuring reports directory exists: {}", reports_dir.display());
    writeln!(log_file, "[STEP] Ensuring reports directory exists: {}", reports_dir.display())?;

    // 1) CLINVAR VCF / TBI
    let (clinvar_url, clinvar_tbi_url) = match build.as_str() {
        "GRCH37" => (
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi",
        ),
        "GRCH38" => (
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
            "https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi",
        ),
        _ => unreachable!(),
    };

    let clinvar_dir = Path::new("clinvar_data");
    fs::create_dir_all(&clinvar_dir)?;
    let clinvar_vcf_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz", build));
    let clinvar_tbi_path = clinvar_dir.join(format!("clinvar_{}.vcf.gz.tbi", build));

    // 2) CLINVAR TSV for deeper annotation
    // We'll store it in the same directory
    let tsv_path = clinvar_dir.join("clinvar_summary_pathogenic.tsv.xz");
    let tsv_url = "https://raw.githubusercontent.com/SauersML/pathogenic/main/clinvar_summary_pathogenic.tsv.xz";

    // 3) 1000 GENOMES FREQUENCY FILE
    let onekg_url = match build.as_str() {
        "GRCH37" => "https://ftp.ensembl.org/pub/grch37/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz",
        "GRCH38" => "https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz",
        _ => unreachable!(),
    };
    let onekg_file_path = clinvar_dir.join(format!("1000GENOMES-phase_3_{}.vcf.gz", build));

    println!("[STEP] Ensuring necessary files exist locally...");
    writeln!(log_file, "[STEP] Ensuring necessary files exist locally...")?;

    // Download CLINVAR VCF if not found or corrupt
    if !clinvar_vcf_path.exists() {
        println!("  -> Missing ClinVar VCF for {build}, downloading...");
        writeln!(
            log_file,
            "  -> Missing ClinVar VCF for {build}, downloading..."
        )?;
        download_file(clinvar_url, &clinvar_vcf_path, &mut log_file)?;
    } else {
        println!("  -> Found local ClinVar VCF: {}", clinvar_vcf_path.display());
        writeln!(log_file, "  -> Found local ClinVar VCF: {}", clinvar_vcf_path.display())?;
        // Check if corrupt
        let testf = File::open(&clinvar_vcf_path)?;
        let mut test_decoder = MultiGzDecoder::new(testf);
        let mut buffer = [0u8; 1024];
        if let Err(e) = test_decoder.read(&mut buffer) {
            println!("  -> ClinVar VCF is corrupt ({e}), re-downloading...");
            writeln!(
                log_file,
                "  -> ClinVar VCF is corrupt ({e}), re-downloading..."
            )?;
            fs::remove_file(&clinvar_vcf_path)?;
            download_file(clinvar_url, &clinvar_vcf_path, &mut log_file)?;
        }
    }

    // Download CLINVAR TBI if missing
    if !clinvar_tbi_path.exists() {
        println!("  -> Missing ClinVar TBI, downloading...");
        writeln!(log_file, "  -> Missing ClinVar TBI, downloading...")?;
        download_file(clinvar_tbi_url, &clinvar_tbi_path, &mut log_file)?;
    } else {
        println!("  -> Found local ClinVar TBI: {}", clinvar_tbi_path.display());
        writeln!(log_file, "  -> Found local ClinVar TBI: {}", clinvar_tbi_path.display())?;
    }

    // Download the 1000G frequency file if missing
    if !onekg_file_path.exists() {
        println!("  -> Missing 1000 Genomes frequency file for {build}, downloading...");
        writeln!(
            log_file,
            "  -> Missing 1000 Genomes frequency file for {build}, downloading..."
        )?;
        download_file(onekg_url, &onekg_file_path, &mut log_file)?;
    } else {
        println!("  -> Found local 1000 Genomes: {}", onekg_file_path.display());
        writeln!(log_file, "  -> Found local 1000 Genomes: {}", onekg_file_path.display())?;
        let testf = File::open(&onekg_file_path)?;
        let mut test_decoder = MultiGzDecoder::new(testf);
        let mut buffer = [0u8; 1024];
        if let Err(e) = test_decoder.read(&mut buffer) {
            println!("  -> 1000 Genomes freq file is corrupt ({e}), re-downloading...");
            writeln!(
                log_file,
                "  -> 1000 Genomes freq file is corrupt ({e}), re-downloading..."
            )?;
            fs::remove_file(&onekg_file_path)?;
            download_file(onekg_url, &onekg_file_path, &mut log_file)?;
        }
    }

    // Download ClinVar summary TSV if missing
    if !tsv_path.exists() {
        println!("  -> Missing ClinVar summary TSV, downloading...");
        writeln!(log_file, "  -> Missing ClinVar summary TSV, downloading...")?;
        download_file(tsv_url, &tsv_path, &mut log_file)?;
    } else {
        println!("  -> Found local ClinVar summary TSV: {}", tsv_path.display());
        writeln!(
            log_file,
            "  -> Found local ClinVar summary TSV: {}",
            tsv_path.display()
        )?;
    }

    // Parse ClinVar VCF for pathogenic variants
    let (clinvar_map, _file_date) = parse_clinvar_vcf_gz(&clinvar_vcf_path, include_vus, include_benign, &mut log_file)?;
    
    // Parse the ClinVar TSV for deeper annotation
    let tsv_map = parse_clinvar_tsv(&tsv_path, &build, &mut log_file)?;

    // Parse the user VCF input
    let input_variants = parse_input_vcf(&input_path, &mut log_file)?;

    println!("[STEP] Matching user variants with ClinVar and 1000G...");
    writeln!(
        log_file,
        "[STEP] Matching user variants with ClinVar and 1000G..."
    )?;

    // Structure to hold temporary records during matching process
    struct TempRecord {
        chr: String,
        pos: u32,
        ref_allele: String,
        alt_allele: String,
        clnsig: String,
        is_alt_pathogenic: bool,
        gene: Option<String>,
        allele_id: Option<i32>,
        genotype: String,
        review_stars: u8,
        af_esp: Option<f64>,
        af_exac: Option<f64>,
        af_tgp: Option<f64>,
        clndn: Option<String>,
        clnsig_category: String,
    }

    // Match user variants with ClinVar
    println!("[STEP] Matching user variants to ClinVar...");
    writeln!(log_file, "[STEP] Matching user variants to ClinVar...")?;
    
    let pb = ProgressBar::new(input_variants.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    let temp_results: Vec<TempRecord> = input_variants
        .par_iter()
        .flat_map_iter(|(_, iv)| {
            let mut local_found = Vec::new();
            for (alt_a, is_present) in &iv.alts {
                if !is_present {
                    continue;
                }
                let key = (iv.chr.clone(), iv.pos, iv.ref_allele.clone(), alt_a.clone());
                if let Some(cv) = clinvar_map.get(&key) {
                    // Include all variants that match our criteria (controlled by flags when parsing)
                    // rather than checking is_alt_pathogenic here, as we did before
                    let genotype = iv.genotype.clone();
                    let review_stars = review_status_to_stars(cv.clnrevstat.as_deref());
                    local_found.push(TempRecord {
                        chr: cv.chr.clone(),
                        pos: cv.pos,
                        ref_allele: cv.ref_allele.clone(),
                        alt_allele: cv.alt_allele.clone(),
                        clnsig: cv.clnsig.clone(),
                        is_alt_pathogenic: cv.is_alt_pathogenic,
                        gene: cv.gene.clone(),
                        allele_id: cv.allele_id,
                        genotype,
                        review_stars,
                        af_esp: cv.af_esp,
                        af_exac: cv.af_exac,
                        af_tgp: cv.af_tgp,
                        clndn: cv.clndn.clone(),
                        clnsig_category: cv.clnsig_category.clone(),
                    });
                }
            }
            pb.inc(1);
            local_found
        })
        .collect();

    pb.finish_with_message("ClinVar matching complete.");

    println!("  -> Matched {} variants in ClinVar", temp_results.len());
    writeln!(
        log_file,
        "  -> Matched {} variants in ClinVar",
        temp_results.len()
    )?;

    // Collect keys of interest from temp_results for 1000 Genomes frequency lookup
    let keys_of_interest: HashSet<(String, u32, String, String)> = temp_results
        .iter()
        .map(|r| (r.chr.clone(), r.pos, r.ref_allele.clone(), r.alt_allele.clone()))
        .collect();

    // Stream 1000 Genomes VCF and extract frequencies for keys_of_interest
    println!("[STEP] Querying 1000 Genomes frequency VCF using CSI index: {}", onekg_file_path.display());
    writeln!(log_file, "[STEP] Querying 1000 Genomes frequency VCF using CSI index: {}", onekg_file_path.display())?;
    let onekg_index_path = clinvar_dir.join(format!("1000GENOMES-phase_3_{}_vcf.gz.csi", build));
    if !onekg_index_path.exists() {
        println!("  -> Missing 1000 Genomes CSI index for {}. Downloading...", build);
        writeln!(log_file, "  -> Missing 1000 Genomes CSI index for {}. Downloading...", build)?;
        download_file(
            if build == "GRCH37" {
                "https://ftp.ensembl.org/pub/grch37/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi"
            } else {
                "https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz.csi"
            },
            &onekg_index_path,
            &mut log_file
        )?;
    }

    // Improved error handling
    println!("  -> Opening 1000 Genomes file: {}", onekg_file_path.display());
    let file = match File::open(&onekg_file_path) {
        Ok(f) => f,
        Err(e) => {
            let error_msg = format!("Failed to open 1000 Genomes file: {}", e);
            println!("ERROR: {}", error_msg);
            writeln!(log_file, "ERROR: {}", error_msg)?;
            return Err(error_msg.into());
        }
    };

    println!("  -> Opening 1000 Genomes CSI index: {}", onekg_index_path.display());
    let index_file = match File::open(&onekg_index_path) {
        Ok(f) => f,
        Err(e) => {
            let error_msg = format!("Failed to open 1000 Genomes CSI index: {}", e);
            println!("ERROR: {}", error_msg);
            writeln!(log_file, "ERROR: {}", error_msg)?;
            return Err(error_msg.into());
        }
    };

    println!("  -> Reading CSI index...");
    let mut csi_reader = csi::io::Reader::new(index_file);
    let index = match csi_reader.read_index() {
        Ok(idx) => idx,
        Err(e) => {
            let error_msg = format!("Failed to read CSI index: {}", e);
            println!("ERROR: {}", error_msg);
            writeln!(log_file, "ERROR: {}", error_msg)?;
            return Err(error_msg.into());
        }
    };

    println!("  -> Creating indexed reader...");
    let mut reader = IndexedReader::new(file, index);

    println!("  -> Reading VCF header...");
    let header = match reader.read_header() {
        Ok(h) => h,
        Err(e) => {
            let error_msg = format!("Failed to read VCF header: {}", e);
            println!("ERROR: {}", error_msg);
            writeln!(log_file, "ERROR: {}", error_msg)?;
            return Err(error_msg.into());
        }
    };

    println!("  -> Setting up frequency mapping...");
    let mut onekg_freqs: HashMap<(String, u32, String, String), (Option<f64>, Option<f64>, Option<f64>, Option<f64>, Option<f64>)> = HashMap::with_capacity(keys_of_interest.len());
    let mut unique_regions: HashSet<(String, u32)> = HashSet::new();
    for (chr, pos, _ref, _alt) in keys_of_interest.iter() {
        unique_regions.insert((chr.clone(), *pos));
    }

    println!("  -> Will query {} unique genomic positions", unique_regions.len());
    writeln!(log_file, "  -> Will query {} unique genomic positions", unique_regions.len())?;

    // Store the size before moving unique_regions
    let unique_regions_count = unique_regions.len();

    // Progress bar for region querying
    let pb = ProgressBar::new(unique_regions_count as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-"),
    );

    // Add success/error counters
    let mut successful_queries = 0;
    let mut failed_queries = 0;

    for (chr, pos) in unique_regions {
        pb.inc(1);
        let region_str = format!("{}:{}-{}", chr, pos, pos);
        let region_obj = match region_str.parse::<Region>() {
            Ok(r) => r,
            Err(e) => {
                // Just log and continue to the next region
                writeln!(log_file, "  -> Failed to parse region {}: {}", region_str, e)?;
                failed_queries += 1;
                continue;
            }
        };
        
        match reader.query(&header, &region_obj) {
            Ok(mut query) => {
                successful_queries += 1;
                while let Some(record_result) = query.next() {
                    match record_result {
                        Ok(record) => {
                            // Process the record as before
                            let record_chr = record.reference_sequence_name().to_string();
                            let record_pos = match record.variant_start() {
                                Some(Ok(pos)) => pos,
                                _ => continue, // Skip if we can't get position
                            };
                            let record_ref = record.reference_bases().to_string();
                            
                            // Get alternative alleles
                            let record_alts: Vec<String> = record.alternate_bases()
                                .iter()
                                .filter_map(Result::ok)
                                .map(|alt| alt.to_string())
                                .collect();
                            
                            let info = record.info();
                            
                            for (alt_idx, alt) in record_alts.iter().enumerate() {
                                let key = (record_chr.clone(), record_pos.get() as u32, record_ref.clone(), alt.clone());
                                if keys_of_interest.contains(&key) {
                                    // Helper function to extract the frequency for a specific alternate allele from the INFO field
                                    let get_freq = |key: &str, alt_idx: usize| -> Option<f64> {
                                        match info.get(&header, key) {
                                            Some(Ok(Some(Value::Array(array)))) => {
                                                match array {
                                                    noodles::vcf::variant::record::info::field::value::Array::Float(values) => {
                                                        let mut iter = values.iter();
                                                        if let Some(Ok(Some(f))) = iter.nth(alt_idx) {
                                                            // Convert f32 to f64 without dereferencing
                                                            Some(f64::from(f))
                                                        } else {
                                                            None
                                                        }
                                                    },
                                                    _ => None, // Handle other Array variants (e.g., Integer, String)
                                                }
                                            },
                                            _ => None,
                                        }
                                    };
                                    let afr = get_freq("AFR", alt_idx);
                                    let amr = get_freq("AMR", alt_idx);
                                    let eas = get_freq("EAS", alt_idx);
                                    let eur = get_freq("EUR", alt_idx);
                                    let sas = get_freq("SAS", alt_idx);
                                    onekg_freqs.insert(key, (afr, amr, eas, eur, sas));
                                }
                            }
                        },
                        Err(e) => {
                            // Just log error and continue
                            writeln!(log_file, "  -> Error processing record for region {}: {}", region_str, e)?;
                        }
                    }
                }
            },
            Err(e) => {
                // Log the error and continue
                writeln!(log_file, "  -> Failed to query region {}: {}", region_str, e)?;
                failed_queries += 1;
            }
        }
    }

    pb.finish_with_message("1000 Genomes querying complete");

    println!("  -> Successfully queried {}/{} regions", successful_queries, unique_regions_count);
    println!("  -> Failed to query {} regions", failed_queries);
    writeln!(log_file, "  -> Successfully queried {}/{} regions", successful_queries, unique_regions_count)?;
    writeln!(log_file, "  -> Failed to query {} regions", failed_queries)?;
    println!("  -> Extracted frequencies for {} variants", onekg_freqs.len());
    writeln!(log_file, "  -> Extracted frequencies for {} variants", onekg_freqs.len())?;

    // Next, we combine with 1000G frequency data
    let mut final_records: Vec<FinalRecord> = Vec::new();
    for r in &temp_results {
        let key = (r.chr.clone(), r.pos, r.ref_allele.clone(), r.alt_allele.clone());
        let (af_afr, af_amr, af_eas, af_eur, af_sas) = match onekg_freqs.get(&key) {
            None => (None, None, None, None, None),
            Some((afr, amr, eas, eur, sas)) => (*afr, *amr, *eas, *eur, *sas),
        };
        let annotation = tsv_map.get(&key);
        let final_rec = FinalRecord {
            chr: r.chr.clone(),
            pos: r.pos,
            ref_allele: r.ref_allele.clone(),
            alt_allele: r.alt_allele.clone(),
            clnsig: r.clnsig.clone(),
            is_alt_pathogenic: r.is_alt_pathogenic,
            gene: r.gene.clone(),
            allele_id: r.allele_id,
            genotype: r.genotype.clone(),
            review_stars: r.review_stars,
            af_esp: r.af_esp,
            af_exac: r.af_exac,
            af_tgp: r.af_tgp,
            clndn: r.clndn.clone(),
            molecular_consequence: annotation.and_then(|ann| ann.molecular_consequence.clone()),
            functional_consequence: annotation.and_then(|ann| ann.functional_consequence.clone()),
            mode_of_inheritance: annotation.and_then(|ann| ann.mode_of_inheritance.clone()),
            preferred_values: annotation.and_then(|ann| ann.preferred_values.clone()),
            citations: annotation.and_then(|ann| ann.citations.clone()),
            comments: annotation.and_then(|ann| ann.comments.clone()),
            family_data: annotation.and_then(|ann| ann.family_data.clone()),
            record_status: annotation.and_then(|ann| ann.record_status.clone()),
            description: annotation.and_then(|ann| ann.description.clone()),
            date_last_evaluated: annotation.and_then(|ann| ann.date_last_evaluated.clone()),
            af_afr,
            af_amr,
            af_eas,
            af_eur,
            af_sas,
            clnsig_category: r.clnsig_category.clone(),
        };
        final_records.push(final_rec);
    }

    // Sort final records by chromosome, then by position
    final_records.sort_by(|a, b| {
        let a_order = get_chromosome_order(&a.chr);
        let b_order = get_chromosome_order(&b.chr);
        
        match a_order.cmp(&b_order) {
            std::cmp::Ordering::Equal => a.pos.cmp(&b.pos),
            other => other
        }
    });

    // Get input filename without extension for output naming
    let input_filename = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("input")
        .to_string()
        .replace(".vcf", "")  // Remove any additional .vcf if present in the stem
        .replace(".gz", "");  // Remove any additional .gz if present in the stem

    // Determine the analysis type for the filename
    let analysis_type = match (include_vus, include_benign) {
        (false, false) => "pathogenic_only",
        (true, false) => "pathogenic_vus",
        (false, true) => "pathogenic_benign",
        (true, true) => "pathogenic_vus_benign",
    };

    // Create the output filename using the input file name, analysis type, and timestamp
    let timestamp = now.format("%Y%m%d_%H%M%S");
    let output_filename = format!(
        "{}_{}_{}.csv",
        input_filename,
        analysis_type,
        timestamp
    );
    
    let output_path = reports_dir.join(output_filename);

    // Create a statistics file with the same base name
    let stats_filename = format!(
        "{}_{}_{}_stats.txt",
        input_filename,
        analysis_type,
        timestamp
    );
    let stats_path = reports_dir.join(stats_filename);
    
    println!("[STEP] Writing final CSV to: {}", output_path.display());
    writeln!(log_file, "[STEP] Writing final CSV to: {}", output_path.display())?;
    let mut wtr = csv::Writer::from_path(&output_path)?;

    wtr.write_record(&[
        "Chromosome",
        "Position",
        "Reference Allele",
        "Alternate Allele",
        "Clinical Significance",
        "Is Alt Pathogenic",
        "Significance Category",
        "Gene",
        "ClinVar Allele ID",
        "CLNDN",
        "Genotype",
        "Review Stars",
        "AF_ESP",
        "AF_EXAC",
        "AF_TGP",
        "Molecular Consequence",
        "Functional Consequence",
        "Mode of Inheritance",
        "Preferred Values",
        "Citations",
        "Comments",
        "Family Data",
        "Record Status",
        "Description",
        "Date Last Evaluated",
        "AF_AFR",
        "AF_AMR",
        "AF_EAS",
        "AF_EUR",
        "AF_SAS",
    ])?;

    for rec in &final_records {
        wtr.write_record(&[
            &rec.chr,
            &rec.pos.to_string(),
            &rec.ref_allele,
            &rec.alt_allele,
            &rec.clnsig,
            &rec.is_alt_pathogenic.to_string(),
            &rec.clnsig_category,
            rec.gene.as_deref().unwrap_or(""),
            &rec.allele_id.map(|id| id.to_string()).unwrap_or_default(),
            rec.clndn.as_deref().unwrap_or(""),
            &rec.genotype,
            &rec.review_stars.to_string(),
            &rec.af_esp.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_exac.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_tgp.map(|f| f.to_string()).unwrap_or_default(),
            rec.molecular_consequence.as_deref().unwrap_or(""),
            rec.functional_consequence.as_deref().unwrap_or(""),
            rec.mode_of_inheritance.as_deref().unwrap_or(""),
            rec.preferred_values.as_deref().unwrap_or(""),
            rec.citations.as_deref().unwrap_or(""),
            rec.comments.as_deref().unwrap_or(""),
            rec.family_data.as_deref().unwrap_or(""),
            rec.record_status.as_deref().unwrap_or(""),
            rec.description.as_deref().unwrap_or(""),
            rec.date_last_evaluated.as_deref().unwrap_or(""),
            &rec.af_afr.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_amr.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_eas.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_eur.map(|f| f.to_string()).unwrap_or_default(),
            &rec.af_sas.map(|f| f.to_string()).unwrap_or_default(),
        ])?;
    }

    wtr.flush()?;
    println!(
        "Done. Wrote {} variants to {}",
        final_records.len(),
        output_path.display()
    );
    writeln!(
        log_file,
        "Done. Wrote {} variants to {}",
        final_records.len(),
        output_path.display()
    )?;

    // Generate and write statistics file
    println!("[STEP] Writing statistics file to: {}", stats_path.display());
    writeln!(log_file, "[STEP] Writing statistics file to: {}", stats_path.display())?;
    
    // Collect statistics
    let mut unique_genes = HashSet::new();
    let mut category_counts = HashMap::new();
    let mut af_ranges = HashMap::new();
    
    for rec in &final_records {
        // Count unique genes
        if let Some(gene) = &rec.gene {
            unique_genes.insert(gene.clone());
        }
        
        // Count by classification category
        *category_counts.entry(rec.clnsig_category.clone()).or_insert(0) += 1;
        
        // Group by allele frequency range (for EAS, EUR, AFR, AMR, SAS)
        let add_af_range = |af: Option<f64>, population: &str, ranges: &mut HashMap<String, i32>| {
            if let Some(freq) = af {
                let range = if freq < 0.001 {
                    "< 0.1%"
                } else if freq < 0.01 {
                    "0.1% - 1%"
                } else if freq < 0.05 {
                    "1% - 5%"
                } else if freq < 0.10 {
                    "5% - 10%"
                } else {
                    "> 10%"
                };
                
                let key = format!("{}_{}", population, range);
                *ranges.entry(key).or_insert(0) += 1;
            }
        };
        
        add_af_range(rec.af_afr, "AFR", &mut af_ranges);
        add_af_range(rec.af_amr, "AMR", &mut af_ranges);
        add_af_range(rec.af_eas, "EAS", &mut af_ranges);
        add_af_range(rec.af_eur, "EUR", &mut af_ranges);
        add_af_range(rec.af_sas, "SAS", &mut af_ranges);
    }
    
    // Write statistics file
    let stats_path_clone = stats_path.clone();
    let mut stats_file = File::create(&stats_path)?;
    
    // Collect the command run for inclusion in the report and statistics
    let args_vec: Vec<String> = std::env::args().collect();
    
    // Create a cleaned command for display that uses "pathogenic" as the command name
    let mut clean_args = Vec::new();
    
    // Skip the binary path (first argument) and replace with "pathogenic"
    clean_args.push("pathogenic".to_string());
    
    // Add all command line arguments
    if args_vec.len() > 1 {
        for arg in &args_vec[1..] {
            clean_args.push(arg.clone());
        }
    }
    
    let command_run = clean_args.join(" ");
    
    writeln!(stats_file, "=== Pathogenic Variant Finder: Analysis Report ===")?;
    writeln!(stats_file, "Date/Time: {}", now.to_rfc3339())?;
    writeln!(stats_file, "\n=== Analysis Settings ===")?;
    writeln!(stats_file, "Input File: {}", input_path.display())?;
    writeln!(stats_file, "Genome Build: {}", build)?;
    writeln!(stats_file, "Include VUS: {}", include_vus)?;
    writeln!(stats_file, "Include Benign: {}", include_benign)?;
    writeln!(stats_file, "Command Used: {}", command_run)?;
    
    writeln!(stats_file, "\n=== Analysis Results ===")?;
    writeln!(stats_file, "Total Variants Processed: {}", input_variants.len())?;
    writeln!(stats_file, "Total Variants Reported: {}", final_records.len())?;
    writeln!(stats_file, "Unique Genes: {}", unique_genes.len())?;
    
    writeln!(stats_file, "\n=== Variant Classifications ===")?;
    for (category, count) in category_counts.iter() {
        writeln!(stats_file, "{}: {}", category, count)?;
    }
    
    writeln!(stats_file, "\n=== Allele Frequency Distribution ===")?;
    for (range, count) in af_ranges.iter() {
        writeln!(stats_file, "{}: {}", range, count)?;
    }
    
    writeln!(stats_file, "\n=== Population Coverage ===")?;
    let mut pop_coverage = HashMap::new();
    // Count variants with any frequency data for each population
    for rec in &final_records {
        if rec.af_afr.is_some() { *pop_coverage.entry("AFR").or_insert(0) += 1; }
        if rec.af_amr.is_some() { *pop_coverage.entry("AMR").or_insert(0) += 1; }
        if rec.af_eas.is_some() { *pop_coverage.entry("EAS").or_insert(0) += 1; }
        if rec.af_eur.is_some() { *pop_coverage.entry("EUR").or_insert(0) += 1; }
        if rec.af_sas.is_some() { *pop_coverage.entry("SAS").or_insert(0) += 1; }
    }
    
    for (pop, count) in pop_coverage.iter() {
        let percentage = if final_records.is_empty() { 
            0.0 
        } else { 
            (*count as f64 / final_records.len() as f64) * 100.0 
        };
        writeln!(stats_file, "{}: {} variants ({:.1}%)", pop, count, percentage)?;
    }
    
    println!(
        "Done. Wrote statistics to {}",
        stats_path_clone.display()
    );
    writeln!(
        log_file,
        "Done. Wrote statistics to {}",
        stats_path_clone.display()
    )?;

    // Generate markdown report if enabled
    if args.markdown_report {
        // Create the markdown filename
        let markdown_filename = format!(
            "{}_{}_{}_report.md",
            input_filename,
            analysis_type,
            timestamp
        );
        let markdown_path = reports_dir.join(markdown_filename);

        // Generate the markdown report
        println!("[STEP] Generating markdown report: {}", markdown_path.display());
        writeln!(log_file, "[STEP] Generating markdown report: {}", markdown_path.display())?;

        generate_markdown_report(
            &final_records,
            &markdown_path,
            &args,
            &category_counts,
            &unique_genes,
            &input_path,
            &build,
            input_variants.len(),
            &timestamp.to_string(),
            &command_run,
            &mut log_file,
        )?;

        println!(
            "Done. Wrote markdown report to {}",
            markdown_path.display()
        );
        writeln!(
            log_file,
            "Done. Wrote markdown report to {}",
            markdown_path.display()
        )?;
    }

    Ok(())
}

/// Generate markdown report of found variants with dynamic sections based on analysis settings
fn generate_markdown_report(
    final_records: &[FinalRecord],
    out_path: &Path,
    args: &Args,
    category_counts: &HashMap<String, usize>,
    unique_genes: &HashSet<String>,
    input_path: &Path,
    build: &str,
    total_variants: usize,
    _timestamp: &str,
    command_run: &str,
    log_file: &mut File,
) -> Result<(), Box<dyn Error>> {
    println!("[STEP] Generating markdown report: {}", out_path.display());
    writeln!(log_file, "[STEP] Generating markdown report: {}", out_path.display())?;
    
    let mut md_file = File::create(out_path)?;
    
    // Current date/time for report
    let now: DateTime<Utc> = Utc::now();
    let datetime = now.format("%B %d, %Y @ %H:%M:%S").to_string();
    
    // Group variants by classification and gene for better organization
    let mut gene_pathogenic: HashMap<String, Vec<&FinalRecord>> = HashMap::new();       // Strictly Pathogenic
    let mut gene_likely_pathogenic: HashMap<String, Vec<&FinalRecord>> = HashMap::new(); // Likely Pathogenic
    let mut gene_vus: HashMap<String, Vec<&FinalRecord>> = HashMap::new();              // Uncertain Significance
    let mut gene_conflicting: HashMap<String, Vec<&FinalRecord>> = HashMap::new();      // Conflicting Interpretations
    let mut gene_benign: HashMap<String, Vec<&FinalRecord>> = HashMap::new();           // Strictly Benign
    let mut gene_likely_benign: HashMap<String, Vec<&FinalRecord>> = HashMap::new();    // Likely Benign
    
    // Track which variant types are actually present in the results
    let mut has_pathogenic = false;
    let mut has_likely_pathogenic = false;
    let mut has_vus = false;
    let mut has_conflicting = false;
    let mut has_benign = false;
    let mut has_likely_benign = false;
    
    // Organize variants by type and gene
    for record in final_records {
        let gene_name = record.gene.clone().unwrap_or_else(|| "Unknown".to_string());
        let clnsig = &record.clnsig;
        
        if record.clnsig_category == "pathogenic" {
            // Split pathogenic variants based on their specific clnsig value
            if clnsig.contains("Likely_pathogenic") && !clnsig.contains("Pathogenic") {
                has_likely_pathogenic = true;
                gene_likely_pathogenic.entry(gene_name).or_default().push(record);
            } else {
                // Either strictly Pathogenic or Pathogenic/Likely_pathogenic
                has_pathogenic = true;
                gene_pathogenic.entry(gene_name).or_default().push(record);
            }
        } else if record.clnsig_category == "vus" && args.include_vus {
            has_vus = true;
            gene_vus.entry(gene_name).or_default().push(record);
        } else if record.clnsig_category == "benign" && args.include_benign {
            // Split benign variants based on their specific clnsig value
            if clnsig.contains("Likely_benign") && !clnsig.contains("Benign") {
                has_likely_benign = true;
                gene_likely_benign.entry(gene_name).or_default().push(record);
            } else {
                // Either strictly Benign or Benign/Likely_benign
                has_benign = true;
                gene_benign.entry(gene_name).or_default().push(record);
            }
        } else if record.clnsig_category == "conflicting" && args.include_vus {
            has_conflicting = true;
            gene_conflicting.entry(gene_name).or_default().push(record);
        }
    }
    
    // Create a list of sections for Table of Contents and report generation
    struct Section<'a> {
        id: String,
        title: String,
        present: bool,
        variants: HashMap<String, Vec<&'a FinalRecord>>,
    }
    
    let mut sections = vec![
        Section {
            id: "disclaimer".to_string(),
            title: "Disclaimer".to_string(),
            present: true, // Always include disclaimer
            variants: HashMap::new(),
        },
        Section {
            id: "summary".to_string(),
            title: "Summary".to_string(),
            present: true, // Always include summary
            variants: HashMap::new(),
        }
    ];
    
    // Add pathogenic sections if present
    if has_pathogenic {
        sections.push(Section {
            id: "pathogenic-variants".to_string(),
            title: "Pathogenic Variants".to_string(),
            present: true,
            variants: gene_pathogenic,
        });
    }
    
    if has_likely_pathogenic {
        sections.push(Section {
            id: "likely-pathogenic-variants".to_string(),
            title: "Likely Pathogenic Variants".to_string(),
            present: true,
            variants: gene_likely_pathogenic,
        });
    }
    
    // Dynamically add sections based on what's included and present
    if has_vus && args.include_vus {
        sections.push(Section {
            id: "uncertain-significance-variants".to_string(),
            title: "Variants of Uncertain Significance".to_string(),
            present: true,
            variants: gene_vus,
        });
    }
    
    if has_conflicting && args.include_vus {
        sections.push(Section {
            id: "conflicting-variants".to_string(),
            title: "Variants with Conflicting Interpretations".to_string(),
            present: true,
            variants: gene_conflicting,
        });
    }
    
    // Add benign sections if present and included
    if has_benign && args.include_benign {
        sections.push(Section {
            id: "benign-variants".to_string(),
            title: "Benign Variants".to_string(),
            present: true,
            variants: gene_benign,
        });
    }
    
    if has_likely_benign && args.include_benign {
        sections.push(Section {
            id: "likely-benign-variants".to_string(),
            title: "Likely Benign Variants".to_string(),
            present: true,
            variants: gene_likely_benign,
        });
    }
    
    // Always add understanding section at the end
    sections.push(Section {
        id: "understanding-this-report".to_string(),
        title: "Understanding This Report".to_string(),
        present: true,
        variants: HashMap::new(),
    });
    
    // Write the header
    writeln!(md_file, "# Genetic Variant Report")?;
    writeln!(md_file)?;
    writeln!(md_file, "**Date / Time Generated:** {}", datetime)?;
    
    // Write the variant types analyzed
    let mut variant_types = Vec::new();
    if has_pathogenic {
        variant_types.push("Pathogenic");
    }
    if has_likely_pathogenic {
        variant_types.push("Likely Pathogenic");
    }
    
    if args.include_vus {
        variant_types.push("Uncertain Significance");
        variant_types.push("Conflicting Interpretations");
    }
    
    if args.include_benign {
        if has_benign {
            variant_types.push("Benign");
        }
        if has_likely_benign {
            variant_types.push("Likely Benign");
        }
    }
    
    writeln!(md_file)?;
    writeln!(md_file, "**Variant Types Analyzed:** {}", variant_types.join(", "))?;
    writeln!(md_file)?;
    
    // Write analysis settings
    writeln!(md_file, "### Analysis Settings:")?;
    writeln!(md_file, "- **Input File:** `{}`", input_path.display())?;
    writeln!(md_file, "- **Genome Build:** `{}`", build)?;
    writeln!(md_file, "- **Include VUS:** `{}`", args.include_vus)?;
    writeln!(md_file, "- **Include Benign:** `{}`", args.include_benign)?;
    writeln!(md_file)?;
    
    // Write command run
    writeln!(md_file, "**Command Used:**")?;
    writeln!(md_file, "```bash")?;
    writeln!(md_file, "{}", command_run)?;
    writeln!(md_file, "```")?;
    
    // Generate dynamically correct table of contents
    writeln!(md_file)?;
    writeln!(md_file, "---")?;
    writeln!(md_file)?;
    writeln!(md_file, "## Table of Contents")?;
    writeln!(md_file)?;
    
    // Only include sections that are present in the results
    let present_sections: Vec<&Section> = sections.iter()
        .filter(|s| s.present)
        .collect();
    
    for (i, section) in present_sections.iter().enumerate() {
        writeln!(md_file, "{}. [{}](#{})", i + 1, section.title, section.id)?;
    }
    
    writeln!(md_file)?;
    writeln!(md_file, "---")?;
    
    // Generate each section
    for section in sections.iter().filter(|s| s.present) {
        if section.id == "disclaimer" {
            write_disclaimer_section(&mut md_file)?;
        } else if section.id == "summary" {
            write_summary_section(&mut md_file, category_counts, args, final_records.len(), total_variants, unique_genes.len())?;
        } else if section.id == "understanding-this-report" {
            write_understanding_section(&mut md_file)?;
        } else {
            // This is a variant section - use appropriate handler
            write_variant_section(&mut md_file, &section.variants, &section.id, &section.title)?;
        }
    }
    
    println!("  -> Markdown report generated successfully.");
    writeln!(log_file, "  -> Markdown report generated successfully.")?;
    
    Ok(())
}

/// Write the disclaimer section
fn write_disclaimer_section(md_file: &mut File) -> Result<(), Box<dyn Error>> {
    writeln!(md_file, "<a id=\"disclaimer\"></a>")?;
    writeln!(md_file, "# Disclaimer")?;
    writeln!(md_file)?;
    writeln!(md_file, "> **Important:** This report is provided for educational purposes only. Results generated or reported by the Pathogenic Variant Finder should not be used for medical diagnosis or to inform clinical decisions without proper validation by qualified healthcare professionals. The developers make no warranties, express or implied, regarding the accuracy, completeness, or reliability of the information provided by this tool. Users should be aware that variant classifications and pathogenicity assessments may change over time as new evidence emerges. Always consult with certified genetic counselors, clinical geneticists, or other appropriate healthcare providers for the interpretation of genetic variants in a medical context.")?;
    writeln!(md_file)?;
    writeln!(md_file, "---")?;
    
    Ok(())
}

/// Write the summary section
fn write_summary_section(
    md_file: &mut File,
    category_counts: &HashMap<String, usize>,
    args: &Args,
    variants_reported: usize,
    variants_processed: usize,
    unique_genes_count: usize,
) -> Result<(), Box<dyn Error>> {
    writeln!(md_file, "<a id=\"summary\"></a>")?;
    writeln!(md_file, "## Summary")?;
    writeln!(md_file)?;
    writeln!(md_file, "This report contains genetic variants found in your genome that match clinical databases. The variants are categorized by clinical significance:")?;
    writeln!(md_file)?;

    // Additional summary stats
    writeln!(md_file)?;
    writeln!(md_file, "**Total Number of Genetic Variants Processed:** {}", variants_processed)?;
    writeln!(md_file)?;
    writeln!(md_file, "**Number of Unique Genes with Identified Variants:** {}", unique_genes_count)?;
    writeln!(md_file)?;
    
    // Summary table
    writeln!(md_file, "| Category | Count |")?;
    writeln!(md_file, "|---------|-------|")?;
    
    let pathogenic_count = category_counts.get("pathogenic").unwrap_or(&0);
    let likely_pathogenic_count = category_counts.get("likely_pathogenic").unwrap_or(&0);
    let conflicting_count = category_counts.get("conflicting").unwrap_or(&0);
    let vus_count = category_counts.get("vus").unwrap_or(&0);
    let likely_benign_count = category_counts.get("likely_benign").unwrap_or(&0);
    let benign_count = category_counts.get("benign").unwrap_or(&0);
    
    writeln!(md_file, "| Pathogenic | {} |", pathogenic_count)?;
    writeln!(md_file, "| Likely Pathogenic | {} |", likely_pathogenic_count)?;
    
    if args.include_vus {
        writeln!(md_file, "| Conflicting | {} |", conflicting_count)?;
        writeln!(md_file, "| Uncertain Significance | {} |", vus_count)?;
    } else {
        writeln!(md_file, "| Conflicting | Not Analyzed |")?;
        writeln!(md_file, "| Uncertain Significance | Not Analyzed |")?;
    }
    
    if args.include_benign {
        writeln!(md_file, "| Likely Benign | {} |", likely_benign_count)?;
        writeln!(md_file, "| Benign | {} |", benign_count)?;
    } else {
        writeln!(md_file, "| Likely Benign | Not Analyzed |")?;
        writeln!(md_file, "| Benign | Not Analyzed |")?;
    }
    
    writeln!(md_file, "| **Total Variants Reported** | **{}** |", variants_reported)?;
    writeln!(md_file)?;
    writeln!(md_file, "---")?;
    
    Ok(())
}

/// Write the variant section (pathogenic, VUS, conflicting, or benign)
fn write_variant_section(
    md_file: &mut File,
    gene_variants: &HashMap<String, Vec<&FinalRecord>>,
    section_id: &str,
    section_title: &str,
) -> Result<(), Box<dyn Error>> {
    writeln!(md_file, "<a id=\"{}\"></a>", section_id)?;
    writeln!(md_file, "## {}", section_title)?;
    writeln!(md_file)?;
    
    if gene_variants.is_empty() {
        writeln!(md_file, "No {} were found in this analysis.", section_title.to_lowercase())?;
        writeln!(md_file)?;
        writeln!(md_file, "---")?;
        return Ok(());
    }
    
    // Summary table
    writeln!(md_file, "### Summary Table")?;
    writeln!(md_file)?;
    writeln!(md_file, "| Gene | Variant Location | DNA Change | Condition | Genotype, Zygosity | Clinical Significance | African Population Frequency | American Population Frequency | East Asian Population Frequency | European Population Frequency | South Asian Population Frequency |")?;
    writeln!(md_file, "|:-----|:-------|:----------|:----------|:---------|:---------|:---------------------|:---------------------|:---------------------|:---------------------|:---------------------|")?;
    
    // Track variant IDs for the details section
    let mut variant_id = 1;
    let mut summary_entries = Vec::new();
    
    // Collect all variants from all genes into a single list
    let mut all_variants: Vec<(&String, &FinalRecord)> = Vec::new();
    for (gene, variants) in gene_variants {
        for variant in variants {
            all_variants.push((gene, *variant));
        }
    }
    
    // Sort all variants by chromosome and position first
    all_variants.sort_by(|(_, a), (_, b)| {
        let a_order = get_chromosome_order(&a.chr);
        let b_order = get_chromosome_order(&b.chr);
        
        match a_order.cmp(&b_order) {
            std::cmp::Ordering::Equal => a.pos.cmp(&b.pos),
            other => other
        }
    });
    
    // Process the sorted variants
    for (gene, variant) in all_variants {
        let chr = &variant.chr;
        let pos = variant.pos;
        let ref_allele = &variant.ref_allele;
        let alt_allele = &variant.alt_allele;
        
        let dna_change = format!("{} -> {}", ref_allele, alt_allele);
        
        // Clean condition field by replacing pipe characters with commas and underscores with spaces
        let condition = match &variant.clndn {
            Some(c) => c.replace('|', ", ").replace('_', " "),
            None => "Not specified".to_string(),
        };
        
        // Determine zygosity from genotype
        let genotype = &variant.genotype;
        let zygosity = if genotype == "1/1" {
            "Homozygous"
        } else if genotype == "0/1" || genotype == "1/0" {
            "Heterozygous"
        } else {
            "Unknown"
        };
        
        let genotype_text = format!("{}, {}", genotype, zygosity);
        
        // Clean clinical significance - just use the classification type without extra data
        let clean_clnsig = variant.clnsig.replace('_', " ");
        
        // Format frequencies as percentages with 2 decimal places
        let format_freq = |f: Option<f64>| -> String {
            match f {
                Some(val) => format!("{:.2}%", val * 100.0),
                None => "N/A".to_string(),
            }
        };
        
        let afr_freq = format_freq(variant.af_afr);
        let amr_freq = format_freq(variant.af_amr);
        let eas_freq = format_freq(variant.af_eas);
        let eur_freq = format_freq(variant.af_eur);
        let sas_freq = format_freq(variant.af_sas);
        
        writeln!(
            md_file,
            "| [{}](#variant-{}) | Chr {}:{} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            gene, variant_id, chr, pos, dna_change, condition, genotype_text, 
            clean_clnsig, afr_freq, amr_freq, eas_freq, eur_freq, sas_freq
        )?;
        
        summary_entries.push((gene.clone(), variant, variant_id));
        variant_id += 1;
    }
    
    // Detailed descriptions
    writeln!(md_file)?;
    writeln!(md_file, "### Detailed Descriptions")?;
    writeln!(md_file)?;
    
    for (gene, variant, id) in summary_entries {
        write_variant_details(md_file, &gene, variant, id)?;
        writeln!(md_file)?;
        writeln!(md_file, "---")?;
        writeln!(md_file)?;
    }
    
    Ok(())
}

/// Write details for a single variant
fn write_variant_details(
    md_file: &mut File,
    gene: &str,
    variant: &FinalRecord,
    variant_id: usize,
) -> Result<(), Box<dyn Error>> {
    writeln!(md_file, "<a id=\"variant-{}\"></a>", variant_id)?;
    writeln!(md_file, "### {}. {} ", variant_id, gene)?;
    writeln!(md_file)?;
    writeln!(md_file, "<details open>")?;
    writeln!(md_file, "<summary><strong>Variant Details</strong></summary>")?;
    writeln!(md_file)?;
    
    // Clean condition field by replacing pipe characters with commas and underscores with spaces
    let condition = match &variant.clndn {
        Some(c) => c.replace('|', ", ").replace('_', " "),
        None => "Not specified".to_string(),
    };
    
    writeln!(md_file, "**Associated Condition:** {}", condition)?;
    writeln!(md_file)?;
    
    // Location
    writeln!(md_file, "**Location:** Chromosome {}, Position {}", variant.chr, variant.pos)?;
    writeln!(md_file)?;
    
    // DNA Change
    writeln!(
        md_file,
        "**DNA Change:** Your DNA has `{}` where the reference genome has `{}`",
        variant.alt_allele, variant.ref_allele
    )?;
    writeln!(md_file)?;
    
    // Genotype / Zygosity
    let zygosity = if variant.genotype == "1/1" {
        "Homozygous (both copies of the gene have this variant)"
    } else if variant.genotype == "0/1" || variant.genotype == "1/0" {
        "Heterozygous (one copy of the gene has this variant)"
    } else {
        ""
    };
    
    writeln!(
        md_file,
        "**Genotype / Zygosity:** {}, {}",
        variant.genotype, zygosity
    )?;
    writeln!(md_file)?;
    
    // Clinical significance - replace underscores with spaces
    let clean_clnsig = variant.clnsig.replace('_', " ");
    writeln!(md_file, "**Clinical Significance of Variant:** {}", clean_clnsig)?;
    writeln!(md_file)?;
    
    // Molecular effects (if available)
    if let Some(mol_effect) = &variant.molecular_consequence {
        if !mol_effect.is_empty() {
            writeln!(md_file, "**Molecular Effects:**")?;
            writeln!(md_file, "- Type: {}", mol_effect.replace('_', " "))?;
            writeln!(md_file)?;
        }
    }
    
    // Population frequencies
    let has_any_freq = variant.af_afr.is_some() || variant.af_amr.is_some() || 
                       variant.af_eas.is_some() || variant.af_eur.is_some() ||
                       variant.af_sas.is_some();
    
    if has_any_freq {
        writeln!(md_file, "**Population Frequencies:** (How common this variant is in different populations)")?;
        writeln!(md_file, "| Population | Frequency |")?;
        writeln!(md_file, "|:-----------|:----------|")?;
        
        // Format frequency helper
        let format_freq = |f: Option<f64>| -> String {
            match f {
                Some(val) => format!("{:.2}%", val * 100.0),
                None => "Not available".to_string(),
            }
        };
        
        if variant.af_afr.is_some() {
            writeln!(md_file, "| African | {} |", format_freq(variant.af_afr))?;
        }
        if variant.af_amr.is_some() {
            writeln!(md_file, "| American | {} |", format_freq(variant.af_amr))?;
        }
        if variant.af_eas.is_some() {
            writeln!(md_file, "| East Asian | {} |", format_freq(variant.af_eas))?;
        }
        if variant.af_eur.is_some() {
            writeln!(md_file, "| European | {} |", format_freq(variant.af_eur))?;
        }
        if variant.af_sas.is_some() {
            writeln!(md_file, "| South Asian | {} |", format_freq(variant.af_sas))?;
        }
        writeln!(md_file)?;
    }
    
    // Review status
    let stars = variant.review_stars;
    let review_status = if stars == 0 {
        " - No assertion criteria provided"
    } else if stars == 1 {
        " ★ - Criteria provided, single submitter"
    } else if stars == 2 {
        " ★★ - Criteria provided, multiple submitters, no conflicts"
    } else if stars == 3 {
        " ★★★ - Criteria provided, reviewed by expert panel"
    } else {
        " ★★★★ - Practice guideline"
    };
    
    writeln!(md_file, "> **Review Status:** {}", review_status)?;
    writeln!(md_file)?;
    
    // Additional fields if available
    if let Some(preferred) = &variant.preferred_values {
        if !preferred.is_empty() {
            writeln!(md_file, "> **Preferred Values:** {}", preferred.replace('_', " "))?;
            writeln!(md_file)?;
        }
    }
    
    if let Some(desc) = &variant.description {
        if !desc.is_empty() {
            writeln!(md_file, "> **Description:** {}", desc.replace('_', " "))?;
            writeln!(md_file)?;
        }
    }
    
    if let Some(citations) = &variant.citations {
        if !citations.is_empty() {
            writeln!(md_file, "> **Research Citations:** {}", citations.replace('_', " "))?;
            writeln!(md_file)?;
        }
    }
    
    if let Some(status) = &variant.record_status {
        if !status.is_empty() {
            writeln!(md_file, "> **Record Status:** {}", status.replace('_', " "))?;
            writeln!(md_file)?;
        }
    }
    
    if let Some(evaluated) = &variant.date_last_evaluated {
        if !evaluated.is_empty() {
            writeln!(md_file, "> **Date Last Evaluated:** {}", evaluated.replace('_', " "))?;
            writeln!(md_file)?;
        }
    }
    
    writeln!(md_file, "</details>")?;
    
    Ok(())
}

/// Write the understanding section of the markdown report
fn write_understanding_section(md_file: &mut File) -> Result<(), Box<dyn Error>> {
    writeln!(md_file, "<a id=\"understanding-this-report\"></a>")?;
    writeln!(md_file, "## Understanding This Report")?;
    writeln!(md_file)?;
    
    writeln!(md_file, "**Clinical Significance Categories:**")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Pathogenic variants** are genetic changes that have strong evidence for causing disease or health conditions.")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Likely Pathogenic variants** have good but not definitive evidence suggesting they cause disease.")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Conflicting Interpretations** are variants where different labs or researchers have come to different conclusions about their significance.")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Variants of Uncertain Significance (VUS)** are genetic changes where there is currently not enough evidence to determine if they are potentially harmful or harmless.")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Likely Benign variants** have good evidence suggesting they do not cause disease.")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Benign variants** are genetic changes that are known to be harmless based on strong evidence.")?;
    writeln!(md_file)?;
    
    writeln!(md_file, "**Important Note:** Having a pathogenic or likely pathogenic variant doesn't necessarily mean you will develop the condition, as other factors like environment, lifestyle, and additional genetic factors also play important roles.")?;
    writeln!(md_file)?;
    
    writeln!(md_file, "**Other Terms:**")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Genotype / Zygosity:**")?;
    writeln!(md_file)?;
    writeln!(md_file, "  - **1/1, Homozygous** means the variant is present on both copies of the gene.")?;
    writeln!(md_file)?;
    writeln!(md_file, "  - **0/1, Heterozygous** means the variant is present on only one copy of the gene.")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Population Frequencies** show how common the variant is in different populations. Rare variants (low percentage) may be more significant than common ones.")?;
    writeln!(md_file)?;
    writeln!(md_file, "- **Review Stars** indicate the level of review in ClinVar, with more stars representing more thorough evaluation.")?;
    
    Ok(())
}
