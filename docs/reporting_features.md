# Reporting Features in Pathogenic Variant Finder

This document details the output reporting features of the Pathogenic Variant Finder tool, including the naming conventions, file formats, and statistical information provided in the output files.

## Output Files Overview

The Pathogenic Variant Finder generates three primary output files for each analysis run:

1. **CSV Report File**: Contains detailed information about all variants detected that match the filtering criteria.
2. **Statistics Text File**: Provides a summary of the analysis settings, results, and distribution statistics.
3. **Markdown Report File**: Offers a comprehensive, human-readable report with organized variant information and explanations.

All files are automatically generated in the `reports/` directory upon completion of an analysis. The markdown report can be disabled using the `--markdown-report=false` command-line option.

## Naming Conventions

The output files follow a specific naming convention to help with organization and identification:

```
[input_filename]_[analysis_type]_[timestamp].(csv|txt|md)
```

Where:
- `input_filename`: The base name of the input VCF file (without extension)
- `analysis_type`: Indicates the types of variants included in the analysis
  - `pathogenic`: Only pathogenic variants
  - `pathogenic_vus`: Pathogenic and variants of uncertain significance
  - `pathogenic_benign`: Pathogenic and benign variants
  - `pathogenic_vus_benign`: All variant types
- `timestamp`: Date and time of the analysis in format `YYYYMMDD_HHMMSS`
- File extensions: `.csv` for the variant report, `_stats.txt` for the statistics file, and `.md` for the markdown report

## CSV Report Format

The CSV report contains the following columns:

| Column | Description |
|--------|-------------|
| Chromosome | Chromosome identifier |
| Position | Genomic position of the variant |
| Reference Allele | The reference allele |
| Alternate Allele | The alternate allele |
| Clinical Significance | ClinVar classification of clinical significance |
| Is Alt Pathogenic | Boolean indicating if the alternate allele is pathogenic |
| Significance Category | Simplified category (pathogenic, benign, vus, conflicting) |
| Gene | Gene symbol |
| ClinVar Allele ID | ClinVar identifier for the allele |
| CLNDN | Clinical disease name associated with the variant |
| Genotype | Genotype of the sample (e.g., 0/1, 1/1) |
| Review Stars | ClinVar review status represented as stars (0-4) |
| AF_ESP | Allele frequency from Exome Sequencing Project |
| AF_EXAC | Allele frequency from Exome Aggregation Consortium |
| AF_TGP | Allele frequency from 1000 Genomes Project |
| AF_AFR | African population allele frequency (1000 Genomes) |
| AF_AMR | American population allele frequency (1000 Genomes) |
| AF_EAS | East Asian population allele frequency (1000 Genomes) |
| AF_EUR | European population allele frequency (1000 Genomes) |
| AF_SAS | South Asian population allele frequency (1000 Genomes) |
| Molecular Consequence | Effect on DNA/RNA |
| Functional Consequence | Effect on protein function |
| Mode of Inheritance | Pattern of inheritance |
| Preferred Values | Preferred variants |
| Citations | Academic citations |
| Comments | Additional notes |
| Family Data | Information about family genetics |
| Record Status | Status of the record |
| Description | Detailed description of the variant |
| Date Last Evaluated | Date of last evaluation in ClinVar |

Variants in the CSV file are sorted by chromosome and position for ease of analysis.

## Statistics File Format

The statistics file provides a comprehensive summary of the analysis in a readable text format. It includes:

### Analysis Settings Section
- Input file path
- Genome build used (e.g., GRCH38)
- Which variant types were included (VUS, benign)
- Command used to run the analysis
- Date and time of analysis

### Analysis Results Section
- Total number of variants processed
- Total number of variants reported
- Number of unique genes covered

### Variant Classifications Section
- Count of each variant classification:
  - Pathogenic
  - Benign
  - VUS (Variants of Uncertain Significance)

### Allele Frequency Distribution Section
- Distribution of variants across different allele frequency ranges for each population:
  - < 0.1% (rare variants)
  - 0.1% - 1% (uncommon variants)
  - 1% - 5% (low-frequency variants)
  - 5% - 10% (common variants)
  - > 10% (high-frequency variants)
- Population codes:
  - AFR: African
  - AMR: American
  - EAS: East Asian
  - EUR: European
  - SAS: South Asian

## Markdown Report Format

The markdown report presents a comprehensive, human-readable summary of the analysis with the following sections:

### Header and Metadata
- Report title and generation timestamp
- Analysis settings: input file, genome build, variant types analyzed
- Command used to generate the report

### Table of Contents
- Dynamic navigation links to all sections present in the report

### Disclaimer
- Clear statement about the limitations of the tool and appropriate usage

### Summary Section
- Overview of total variants processed and reported
- Count of unique genes with identified variants
- Summary table with counts for each variant classification type

### Variant Sections (Dynamically Generated)
- Separate sections for each variant classification present in the results:
  - Pathogenic Variants
  - Likely Pathogenic Variants
  - Variants of Uncertain Significance (if `--include-vus` is specified)
  - Variants with Conflicting Interpretations (if `--include-vus` is specified)
  - Benign Variants (if `--include-benign` is specified)
  - Likely Benign Variants (if `--include-benign` is specified)

### Variant Organization
- Variants are grouped by gene and sorted by chromosome and position
- Each variant section includes:
  - Summary table with key information
  - Detailed information for each variant

### Understanding Section
- Explanations of key terms and classifications
- Information about data sources and interpretation guidelines
- Guidance on how to use the report

The markdown report is designed to be both readable as plain text and rendered properly with markdown viewers for enhanced visualization.

## Example Output

### Example Statistics File
```
=== Pathogenic Variant Finder: Analysis Report ===
Date/Time: 2025-03-25T23:40:18.591768597+00:00

=== Analysis Settings ===
Input File: /workspaces/pathogenic/data/TrevorCampbell-SQ63A788-30x-WGS-Sequencing_com-02-22-25.snp-indel.genome.vcf.gz
Genome Build: GRCH38
Include VUS: true
Include Benign: true
Command Used: pathogenic -b GRCh38 -i data/TrevorCampbell-SQ63A788-30x-WGS-Sequencing_com-02-22-25.snp-indel.genome.vcf.gz -v -n

=== Analysis Results ===
Total Variants Processed: 13027212
Total Variants Reported: 36853
Unique Genes: 4882

=== Variant Classifications ===
benign: 36647
pathogenic: 11
vus: 195

=== Allele Frequency Distribution ===
EAS_0.1% - 1%: 1339
AMR_1% - 5%: 1165
AFR_5% - 10%: 2162
...
```

## Usage Recommendations

1. **Archiving Reports**: All reports have timestamps to prevent overwriting, making it safe to run multiple analyses on the same file with different settings.

2. **Comparing Analyses**: The naming convention makes it easy to compare different analyses of the same input file.

3. **Filtering Results**: The CSV format allows for easy filtering and further analysis in spreadsheet software or programming languages.

4. **Statistical Analysis**: The statistics file provides a quick overview without having to parse the entire CSV file.

## Implementation Details

The reporting functionality is implemented in the `main.rs` file. Key components include:

1. **Filename Generation**: Input filename extraction, analysis type determination, and timestamp generation.

2. **CSV Writing**: Sorted variant data is written to a CSV file with appropriate headers.

3. **Statistics Collection**: Counts of various metrics are collected during processing and summarized.

4. **Statistics File Generation**: A formatted text file is created with sections for different statistical categories.

## Future Enhancements

Potential improvements to reporting features:

1. **HTML Report Generation**: Interactive HTML reports with visualization of variant distributions.

2. **Advanced Filtering Options**: Ability to specify custom filters for report generation.

3. **Integration with External Resources**: Links to external databases for each variant.

4. **Comparison Reports**: Automated comparison between different analysis runs.

5. **Exportable Visualizations**: Generation of charts and graphs for key statistics. 