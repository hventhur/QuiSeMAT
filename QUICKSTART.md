# QuiSeMAT Quick Start Guide

## 🚀 Complete Analysis Workflow

### Prerequisites
- Python 3.9 or higher
- NCBI BLAST installed (via `brew install blast` on macOS)

### Prerequisites Setup

```bash
# Activate virtual environment
cd /Users/herbertventhur/Desktop/dev/gennotate
source venv/bin/activate
```

---

## 📋 Step-by-Step Complete Analysis

Follow these steps in order for a complete analysis:

### **Step 1: Main Pipeline - Gene Matching & Filtering**
```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4
```
**What it does:**
- Validates input files (transcriptome assembly + annotated gene database)
- Performs BLAST homology search
- Filters isoforms (keeps longest sequence per gene at ~98% identity)
- Generates matched gene results

**Output files:**
- `blast_results/final_annotated_genes.tsv` - Complete alignment details
- `blast_results/final_annotated_genes_summary.csv` - Quick reference (5 columns)
- `blast_results/results_statistics.txt` - Quality metrics

---

### **Step 2 (Optional): Analyze Unmatched Genes**
```bash
python analyze_unmatched_genes.py
```
**When to use:** To investigate genes from the database that didn't match because of high % of sequence identity (98% as default)

**What it does:**
- Identifies genes that failed to match or didn't meet filtering criteria
- Provides BLAST statistics for unmatched genes
- Helps assess data quality

**Output file:**
- `blast_results/unmatched_genes_analysis.txt` - Detailed analysis report

---

### **Step 3 (Optional): Extract TPM Expression Values**
```bash
python extract_tpm_values.py
```
**When to use:** If you have TPM expression data in `GxH_TPM.tsv`

**What it does:**
- Matches transcripts with their TPM expression values
- Creates enriched output with expression data

**Output file:**
- `blast_results/annotated_genes_tpm_values.tsv` - Genes with TPM values

---

## ✅ Quick Run: Full Analysis in One Command

```bash
# Run main pipeline
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4 && \

# Analyze unmatched genes
python analyze_unmatched_genes.py && \

# Extract TPM values (if GxH_TPM.tsv exists)
python extract_tpm_values.py
```

---

## 📊 Expected Results

For the test data (Gmel2025_v2.fasta + GmelOBPs.txt):

| Metric | Value |
|--------|-------|
| Database genes | 25 |
| Matched genes | 24 (96%) |
| Mean identity | 99.82% |
| Perfect matches (100%) | 17 (70.8%) |
| Weakest match | 98.75% |
| Strongest E-value | 1.72e-174 |

---

## 📖 What Each Output File Contains

### final_annotated_genes.tsv
Complete alignment information for each matched gene:
```
Database_Gene    Matched_Transcript    Sequence_Identity_%    E_value    ...
GmelGOBP2L       NonamEVm010587t1      100.000               4.15e-123  ...
GmelGOBP56aL     NonamEVm011615t1      100.000               6.27e-95   ...
GmelGOBP72L      NonamEVm011325t1      100.000               5.59e-101  ...
```

### final_annotated_genes_summary.csv
Quick reference (5 columns):
```csv
Database_Gene,Matched_Transcript,Identity_%,E-value,Transcript_Length_bp
GmelGOBP2L,NonamEVm010587t1,100.00,4.15e-123,570
GmelGOBP56aL,NonamEVm011615t1,100.00,6.27e-95,1161
```

### results_statistics.txt
Quality metrics:
```
Matched database genes: 24/25 (96.0%)
Mean sequence identity: 99.82%
Mean bit score: 314.3
All matches: E-value < 1e-80 (extremely significant)
```

---

## 🎛️ Advanced Command Options (main.py)

```bash
# Help information
python main.py --help

# Verbose output (shows detailed logs)
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4 --verbose

# Run specific pipeline steps only
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 1  # Validation only
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 2  # + BLAST
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 3  # + Filtering
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4  # Complete pipeline

# Custom output directory
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --output-dir my_results --step 4
```

---

## � Extracting TPM Values (Optional)

After running the main pipeline, you can extract TPM (Transcripts Per Million) expression values for your matched genes.

### When to Use This

Use `extract_tpm_values.py` when you have:
- Completed the main pipeline (`python main.py --step 4`)
- A TPM data file (`GxH_TPM.tsv`) with expression values for your transcripts

### How to Run

```bash
# After running the main pipeline
python extract_tpm_values.py
```

### What It Does

- Reads matched transcripts from `blast_results/final_annotated_genes_summary.csv`
- Extracts TPM values from `GxH_TPM.tsv` for each matched transcript
- Creates `blast_results/annotated_genes_tpm_values.tsv` with combined data

### Output Format

Creates a TSV file with the following columns:
```
Database_Gene    Matched_Transcript    GR1    GR2    GR3    HR1    HR2    HR3
GmelGOBP2L       NonamEVm010587t1      0.5    1.2    0.8    2.1    1.9    2.3
GmelGOBP56aL     NonamEVm011615t1      3.4    4.1    3.8    5.2    4.9    5.1
```

### Requirements

- `GxH_TPM.tsv` file must be in the working directory
- Must have run the main pipeline first (creates `final_annotated_genes_summary.csv`)

---

## �🐛 Troubleshooting

**BLAST not found?**
```bash
# Check if installed
which blastx

# Install on macOS
brew install blast
```

**File not found?**
- Check full path to FASTA files
- Ensure files are readable: `ls -l file.fasta`

**Different results?**
- Results may vary slightly with different BLAST versions
- The algorithm is deterministic for same inputs

---

## 📚 Learn More

See detailed documentation:
- `README.md` - Complete project overview
- `IMPLEMENTATION_SUMMARY.md` - Technical implementation details
- `STEP_THREE_FOUR_IMPLEMENTATION.md` - Steps 3-4 algorithms
- `INSTRUCTIONS.md` - Original project requirements

---

## ✅ You're All Set!

The QuiSeMAT pipeline is ready to use. Run your first analysis:

```bash
python main.py -t Gmel2025_v2.fasta -d GmelOBPs.txt --step 4 --verbose
```

Questions? Check the detailed documentation files!
