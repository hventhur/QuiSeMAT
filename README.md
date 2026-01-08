# QuiSeMAT - Quick Sequence Match Among Transcriptomes
## Pipeline Implementation: All Steps 1-4 Complete ✅

---

## 🎯 Project Overview

**QuiSeMAT** is a bioinformatics pipeline written in Python that facilitates the search of already-annotated genes from RNA sequencing projects (transcriptomes) by comparing them against a database of known genes with homology-based search.

### Current Implementation Status
- ✅ **Step 1**: File Validation 
- ✅ **Step 2**: Homology Search using BLAST
- ✅ **Step 3**: Isoform Filtering & Best Match Selection
- ✅ **Step 4**: Final Results Report Generation

---

## 📋 Quick Start

### Prerequisites
- Python 3.9+
- NCBI BLAST (installed via Homebrew on macOS)
- Virtual environment (included in project)

### Installation

1. **Clone/navigate to the project**:
```bash
cd /Users/herbertventhur/Desktop/dev/gennotate
```

2. **Activate virtual environment** (already set up):
```bash
source venv/bin/activate
```

3. **Verify dependencies**:
```bash
pip install -r requirements.txt
```

### Basic Usage

```bash
# Run complete pipeline (all 4 steps)
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4

# With verbose output (shows details and top matches)
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4 --verbose

# Run specific steps only
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 2  # Steps 1-2 only
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 1  # Step 1 only (validation)
```

# Step 1 only (validation)
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 1

# Custom BLAST method
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --blast-type blastn

# Custom output directory
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --output-dir my_results
```

---

## 🔧 Command-Line Interface

```
usage: QuiSeMAT [-h] -t TRANSCRIPTOME -d DATABASE [-v] [--step {1,2,3,4}]
                [--blast-type {blastx,blastn,blastp,auto}] [--output-dir OUTPUT_DIR] [--version]

Quick Sequence Match Among Transcriptomes (QuiSeMAT) - 
A bioinformatics pipeline for searching annotated genes in transcriptome assemblies.

Required Arguments:
  -t, --transcriptome TRANSCRIPTOME
                        Path to transcriptome assembly FASTA file (file one)
  -d, --database DATABASE
                        Path to database of annotated genes FASTA file (file two)

Optional Arguments:
  -v, --verbose         Enable verbose output
  --step {1,2,3,4}      Steps to execute: 1 (validation), 2 (+ BLAST search), 
                        3 (+ isoform filtering), 4 (complete pipeline).
                        Default: 4
  --blast-type {blastx,blastn,blastp,auto}
                        BLAST method to use. 'auto' detects automatically.
                        Default: auto
  --output-dir OUTPUT_DIR
                        Directory for output files. Default: blast_results
  -h, --help            Show help message and exit
  --version             Show program version and exit
```

---

## 📊 Step 1: File Validation

### What It Does
- Validates FASTA format of both input files
- Checks file existence, readability, and non-empty content
- Counts sequences and analyzes sequence lengths
- Provides comprehensive file statistics

### Output
```
Transcriptome Assembly (Gmel2025_v2.fasta):
  • Total sequences: 36,343
  • Total bases: 54,279,974
  • Sequence length range: 200 - 48,084 bp

Database (GmelOBPs.txt):
  • Total sequences: 25
  • Total bases: 4,147
  • Sequence length range: 121 - 322 bp
```

### Example
```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 1
```

---

## 🔍 Step 2: Homology Search using BLAST

### What It Does
1. **Sequence Type Detection**: Automatically identifies nucleotide vs. protein sequences
2. **BLAST Method Selection**: Chooses appropriate tool (BLASTx, BLASTn, or BLASTp)
3. **Database Creation**: Builds indexed BLAST database from annotated genes
4. **Homology Search**: Runs BLAST search with statistical thresholds
5. **Results Parsing**: Extracts alignment information and statistics
6. **Results Export**: Saves in TSV and JSON formats

### Output Files

**In `blast_results/` directory:**

1. **blast_results.tsv** (11 KB)
   - Raw BLAST output in tabular format (12 columns)
   - Columns: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
   - Ready for further processing

2. **blast_results.json** (43 KB)
   - Structured data for programmatic access
   - Contains: total_hits, unique_queries, unique_subjects, matches array
   - Full alignment details for each match

3. **BLAST Database Files** (250 KB)
   - BLAST index files for protein/nucleotide sequences
   - Can be reused for additional searches

### Test Results
```
Running BLASTx on test data:
  Total BLAST hits found: 150
  Unique query sequences matched: 40
  Unique database sequences matched: 25 (100%)

Top Match:
  NonamEVm006891t1 → GmelOBP7aL
  Identity: 98.75% | Length: 319 bp | E-value: 1.72e-174
```

### Example
```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 2
```

---

## 🧬 Step 3: Filter Isoforms and Select Best Matches

### What It Does
1. **Load Sequences**: Retrieves all transcriptome sequences and their lengths
2. **Group Matches**: Organizes BLAST hits by database gene
3. **Apply Identity Threshold**: Filters matches with ≥98% sequence identity
4. **Select Best Isoform**: Picks longest transcript among high-identity matches for each database gene

### Key Algorithm
- For each database gene, finds all matching transcriptome sequences (potential isoforms)
- Filters to only those with ≥98% sequence identity
- Among filtered matches, selects the longest to represent the most complete gene

### Test Results
```
Total database genes: 25
Passed filtering (≥98% identity): 24
Filtered out (no ≥98% match): 1

Statistics:
  • 23 genes with 100% identity
  • 1 gene with 99.33% identity
  • 1 gene with 98.94% identity
  • 24 genes exceeded threshold
  • GmelOBP17L filtered out (best: 97.40% identity)
```

### Example
```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 3
```

---

## 📊 Step 4: Generate Final Results Report

### What It Does
1. **Create TSV Report**: Generates comprehensive alignment information
2. **Create Summary CSV**: Produces simplified quick-reference file
3. **Generate Statistics**: Calculates and reports quality metrics
4. **Organize Results**: Ensures each database gene has exactly one match

### Output Files

**In `blast_results/` directory:**

1. **final_annotated_genes.tsv** (Main Results - 3 KB)
   - Tab-separated values with complete alignment details
   - 13 columns: Database_Gene, Matched_Transcript, Sequence_Identity_%, Alignment_Length_bp, Mismatches, Gap_Opens, Query_Start, Query_End, Subject_Start, Subject_End, E_value, Bit_Score, Transcript_Length_bp
   - One match per database gene, sorted alphabetically
   - Ready for downstream analysis

2. **final_annotated_genes_summary.csv** (Quick Reference - 1 KB)
   - Simplified CSV format (5 columns)
   - Columns: Database_Gene, Matched_Transcript, Identity_%, E-value, Transcript_Length_bp
   - Easy to open in spreadsheet applications

3. **results_statistics.txt** (Quality Report - 2 KB)
   - Comprehensive statistics including:
     - Match rate (96% for test data)
     - Sequence identity statistics (min/max/mean/median)
     - E-value distribution analysis
     - Bit score statistics
     - Transcript length analysis

### Test Results Summary
```
Database genes: 25
Matched genes: 24 (96.0% match rate)

Sequence Identity:
  Mean: 99.82%
  Min: 98.75%
  Max: 100.00%
  Median: 100.00%

E-value (Statistical Significance):
  Best: 1.72e-174
  Worst: 1.93e-86
  (All < 1e-80: Extremely significant)

Transcript Length:
  Mean: 814 bp
  Min: 446 bp
  Max: 1719 bp
```

### Example
```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4
```

---

## 📁 Project Structure

```
gennotate/
├── main.py                              # Main pipeline script
├── venv/                                # Virtual environment
├── requirements.txt                     # Python dependencies
├── INSTRUCTIONS.md                      # Original project requirements
├── STEP_TWO_IMPLEMENTATION.md          # Step 2 detailed documentation
├── STEP_THREE_FOUR_IMPLEMENTATION.md   # Steps 3-4 detailed documentation
├── README.md                            # This file
│
├── Gmel2025_v2.fasta                   # Test data: Transcriptome assembly
├── GmelOBPs.txt                        # Test data: Annotated genes database
├── GxH_TPM.tsv                         # Gene expression data (not used yet)
│
└── blast_results/                      # Output directory
    ├── blast_results.tsv                           # Step 2: BLAST results (tabular)
    ├── blast_results.json                         # Step 2: BLAST results (JSON)
    ├── annotated_genes_db.*                       # Step 2: BLAST database files
    ├── final_annotated_genes.tsv                  # Step 4: Final results (detailed)
    ├── final_annotated_genes_summary.csv          # Step 4: Final results (summary)
    └── results_statistics.txt                     # Step 4: Quality statistics
```

---

## 🛠️ Implementation Details

### Technologies Used
- **Python 3.9**: Core language
- **Biopython 1.83**: FASTA parsing and sequence handling
- **NCBI BLAST 2.17**: Homology search engine
- **JSON**: Structured data export

### Key Features
1. **Robust Error Handling**: Graceful handling of missing files, invalid formats, and BLAST errors
2. **Automatic Detection**: Detects sequence types and selects appropriate BLAST method
3. **Multi-threading**: Uses 4 parallel threads for faster BLAST execution
4. **Configurable Thresholds**: E-value threshold (default: 1e-5) can be modified
5. **Comprehensive Logging**: Detailed output with optional verbose mode

### BLAST Parameters
- **E-value threshold**: 1e-5 (statistically significant alignments)
- **Output format**: Tabular with all columns
- **Threads**: 4 (configurable)
- **Database type**: Protein or Nucleotide

---

## 📝 Usage Examples

### Example 1: Complete Pipeline Run
```bash
cd /Users/herbertventhur/Desktop/dev/gennotate
source venv/bin/activate
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt
```

### Example 2: With Verbose Output
```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt -v
```

Shows:
- BLAST commands being executed
- Top 10 matches with identity percentages
- Detailed sequence type detection
- Database creation commands

### Example 3: Validation Only
```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 1
```

### Example 4: Custom BLAST Method
```bash
# Force BLASTn (nucleotide vs nucleotide)
python main.py -t assembly.fasta -d genes.fasta --blast-type blastn
```

### Example 5: Custom Output Directory
```bash
python main.py -t assembly.fasta -d genes.fasta --output-dir results_v2
```

---

## 🔧 Troubleshooting

### BLAST Not Found
```
Error: BLAST tools are not installed or not in PATH.
```

**Solution**:
```bash
brew install blast
```

### Invalid FASTA File
```
Error: [filename] is not valid FASTA format.
```

**Solution**: Ensure file is valid FASTA with:
- `>header` lines starting sequences
- Sequence data on following lines
- No empty sequences

### Virtual Environment Issues
```bash
# Recreate virtual environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## 🚀 Next Steps (Towards Step 3)

### Step 3: Isoform Filtering & Header Matching
Will implement:
1. **Group transcripts by database gene**: Organize multiple matches per gene
2. **Filter by sequence identity**: Keep sequences ≥98% identical
3. **Select longest isoform**: Among filtered sequences, choose the longest
4. **Output CSV mapping**: Query → Subject matches

### Step 4: Final Output
Will produce:
1. **Final CSV file**: With filtered, unique matches
2. **Quality metrics**: Coverage, identity, alignment length
3. **Summary report**: Statistics about matched vs unmatched genes

---

## 📚 References

### Input File Format
- **Transcriptome FASTA**: Standard FASTA format with sequence headers and nucleotide sequences
- **Database FASTA**: Can contain nucleotide or protein sequences with annotated headers

### BLAST Documentation
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs)
- [BLAST command-line guide](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

---

## 📄 License & Usage
This is a bioinformatics research tool. Use according to your institution's guidelines and cite appropriately.

---

## 👤 Development Info
- **Project**: QuiSeMAT (Quick Sequence Match Among Transcriptomes)
- **Version**: 0.2.0 (Steps 1 & 2 Complete)
- **Python**: 3.9+
- **Last Updated**: January 6, 2026

---

## 💡 Tips & Best Practices

1. **Always use virtual environment** to avoid dependency conflicts
2. **Validate input files** before running BLAST (use `--step 1`)
3. **Use verbose mode** (`-v`) for first runs to see what's happening
4. **Check output files** in `blast_results/` directory after completion
5. **Save JSON results** for Step 3 processing
6. **Keep BLAST database** if running multiple searches on same data

---

**Questions? Check STEP_TWO_IMPLEMENTATION.md for detailed Step 2 documentation.**
