# Pathogenic Variant Finder

A high-performance Rust tool to identify pathogenic variants from VCF files using ClinVar and 1000 Genomes data.

## Overview

This tool analyzes VCF (Variant Call Format) files to identify potentially pathogenic genetic variants by cross-referencing them with the official ClinVar database and incorporating population frequency data from the 1000 Genomes Project. It supports both GRCh37 (hg19) and GRCh38 (hg38) genome builds.

## Features

- **Comprehensive Variant Detection**: Identifies variants classified as pathogenic, benign, or of uncertain significance in ClinVar
- **Population Frequency Integration**: Incorporates allele frequencies from 1000 Genomes (global and by population)
- **High Performance**:
  - Parallel processing using Rayon's work-stealing thread pool
  - Efficient file handling with memory-mapped files
  - Parallel chunk downloads for reference databases
- **Rich Annotation**: Provides comprehensive annotations including:
  - Clinical significance classifications
  - Gene information
  - Molecular and functional consequences
  - Mode of inheritance if available
  - Citations and additional evidence
- **Automatic Reference Data Management**: Downloads and maintains necessary reference databases
- **Advanced Reporting**: Generates detailed CSV reports and summary statistics files
- **Flexible Variant Selection**: Configurable inclusion of pathogenic, VUS, and benign variants

## Installation

There are multiple ways to install the Pathogenic Variant Finder:

### Option 1: Using the installer script (recommended)

```bash
# Clone the repository
git clone https://github.com/SauersML/pathogenic.git
cd pathogenic

# Run the installer script (you might need to make it executable first)
chmod +x install.sh
./install.sh

# Or run it directly with bash
bash install.sh
```

The installer will:
1. Build the release version of the tool
2. Give you options to create a symlink in a directory in your PATH
3. Allow you to run the tool simply as `pathogenic` from anywhere

### Option 2: Manual installation

```bash
# Clone the repository
git clone https://github.com/SauersML/pathogenic.git
cd pathogenic

# Build the project
cargo build --release

# The binary will be available at target/release/pathogenic_variant_finder
# You can run it directly:
./target/release/pathogenic_variant_finder -b GRCh38 -i your_variants.vcf

# Optional: Create a symlink for easier access
sudo ln -sf "$(pwd)/target/release/pathogenic_variant_finder" /usr/local/bin/pathogenic
```

### Development Container Users

If you're using the provided Dev Container, the tool will be automatically built and made available as `pathogenic` in your PATH when the container starts.

## Usage

Once installed, you can use the tool simply as:

```bash
# Basic usage (pathogenic variants only)
pathogenic -b GRCh38 -i your_variants.vcf

# Include Variants of Uncertain Significance (VUS)
pathogenic -b GRCh38 -i your_variants.vcf -v

# Include benign variants
pathogenic -b GRCh38 -i your_variants.vcf -n

# Include both VUS and benign variants
pathogenic -b GRCh38 -i your_variants.vcf -v -n

# Disable markdown report generation
pathogenic -b GRCh38 -i your_variants.vcf --markdown-report=false
```

### Command-line Arguments

- `--build`, `-b`: Genome build, must be either "GRCh37" (hg19) or "GRCh38" (hg38)
- `--input`, `-i`: Path to the input VCF file (can be uncompressed or gzipped)
- `--include-vus`, `-v`, `--vus`: Include variants of uncertain significance in the output
- `--include-benign`, `-n`, `--benign`: Include benign variants in the output
- `--markdown-report`, `--md-report`: Generate markdown report (enabled by default, use `--markdown-report=false` to disable)

## Output

The tool generates the following output files in the `reports/` directory:

### 1. CSV Report File

The CSV report contains detailed information about all identified variants with the following columns:

- Chromosome, Position, Reference Allele, Alternate Allele
- Clinical Significance
- Is Alt Pathogenic
- Significance Category (pathogenic, benign, vus, conflicting)
- Gene
- ClinVar Allele ID
- Clinical Disease Name (CLNDN)
- Genotype
- Review Stars (0-4 star rating system)
- Allele frequencies:
  - ESP, ExAC, and 1000 Genomes Project (TGP)
  - Population-specific frequencies (AFR, AMR, EAS, EUR, SAS)
- Detailed annotations:
  - Molecular and Functional Consequences
  - Mode of Inheritance
  - Preferred Values
  - Citations
  - Comments and Family Data
  - Record Status
  - Description
  - Date Last Evaluated
 
Depending on available data, many of these fields may not be available.

### 2. Statistics Text File

A companion statistics file summarizing the analysis settings and results:

- Analysis settings (input file, genome build, variant types included)
- Command used to run the analysis
- Total variants processed and reported
- Number of unique genes
- Counts of each variant classification type
- Distribution of allele frequencies across populations

### 3. Markdown Report File

A comprehensive, human-readable markdown report that includes:

- Analysis settings and command used
- Summary statistics with variant counts by category
- Table of contents linking to different sections
- Detailed variant information organized by clinical significance
- Variants grouped by gene with sortable tables
- Detailed information about each variant including:
  - Location, DNA change, and genotype
  - Clinical significance and disease associations
  - Population frequencies
  - Mode of inheritance and other annotations
- Understanding section with explanations of terms

The markdown report is generated by default but can be disabled using `--markdown-report=false`.

### Output File Naming

Output files follow this naming convention:
```
[input_filename]_[analysis_type]_[timestamp].csv
[input_filename]_[analysis_type]_[timestamp]_stats.txt
[input_filename]_[analysis_type]_[timestamp].md
```

Where `analysis_type` indicates which variant types were included in the report.

## How It Works

1. **Data Collection**: The tool automatically downloads necessary reference databases:
   - ClinVar VCF (GRCh37 or GRCh38)
   - ClinVar summary TSV for deeper annotation
   - 1000 Genomes frequency data

2. **Variant Processing**:
   - Parses the user's VCF input
   - Filters for variants present in the sample's genotype
   - Matches variants against ClinVar based on selected variant types
   - Integrates 1000 Genomes allele frequency data
   - Adds detailed annotations from ClinVar summary data

3. **Output Generation**:
   - Sorts variants by chromosome and position
   - Outputs comprehensive CSV with all annotations
   - Generates a statistics file summarizing the analysis
   - Creates a detailed markdown report with interactive sections (unless disabled)

## Logging

The tool maintains a log file (`pathogenic.log`) that captures all processing steps and can be useful for troubleshooting.

## Documentation

For more detailed information, see the documentation in the `docs/` directory:

- [Implementation Details](docs/implementation_details.md)
- [Noodles Integration](docs/noodles_integration.md)
- [Parallel Processing](docs/parallel_processing.md)
- [Reporting Features](docs/reporting_features.md)
- [1000 Genome Frequency Extraction](docs/1000genome_frequency_extraction.md)

## Usage
Do not use this for clinical purposes.
- I might have written a bug in the code.
- Your input file might be wrong.
- Your input file may be right, but different somehow than the files I tested on.
- I did not test this extensively. In fact, there are no unit tests. (Though I did run an integration test via other tooling.)
- You might have done something wrong. For example, you entered the wrong genome build.
- The data in the database might be wrong about pathogenicity of some variant.
- You don't know how to interpret the results.
- You don't know how to interpret mode of inheritance.
- You don't know how to assess evidence for or confidence in pathogenicity.
- Your input file is imputed or has spurious genotype calls.
- Many other reasons not listed here.

What can you use this for?
- Finding potentially pathogenic variants to follow up on via other means.

## Disclaimer
This software is provided for educational purposes only. Results generated by the Pathogenic Variant Finder should not be used for medical diagnosis or to inform clinical decisions without proper validation by qualified healthcare professionals. The developers make no warranties, express or implied, regarding the accuracy, completeness, or reliability of the information provided by this tool. Users should be aware that variant classifications and pathogenicity assessments may change over time as new evidence emerges. Always consult with certified genetic counselors, clinical geneticists, or other appropriate healthcare providers for the interpretation of genetic variants in a medical context.
