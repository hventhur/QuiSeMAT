"""
Extract TPM values for matched transcripts from annotation results.
Creates a new TSV file with transcript IDs and their corresponding TPM values.
"""

import csv
import sys

def extract_tpm_values():
    """Read CSV annotations and extract TPM values from TSV data."""
    
    # Read GxH_TPM.tsv into a dictionary
    tpm_data = {}
    tpm_columns = []  # To store column names dynamically
    try:
        with open('GxH_TPM.tsv', 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            # Get column names, excluding the first column (transcript ID)
            if reader.fieldnames:
                # First column is transcript ID, rest are condition/sample columns
                tpm_columns = [col for col in reader.fieldnames if col]  # Filter empty column names
                # The first column (if empty name) is the transcript ID column
                transcript_id_col = reader.fieldnames[0]
                
            for row in reader:
                # Get transcript ID from the first column (may have empty name)
                transcript_id = row[transcript_id_col] if transcript_id_col in row else list(row.values())[0]
                
                # Store all TPM values dynamically
                tpm_data[transcript_id] = {}
                for col in tpm_columns:
                    if col and col in row:
                        tpm_data[transcript_id][col] = row[col]
                        
    except Exception as e:
        print(f"Error reading GxH_TPM.tsv: {e}", file=sys.stderr)
        return False
    
    # Read CSV and extract matched transcripts
    output_data = []
    try:
        with open('blast_results/final_annotated_genes_summary.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row_num, row in enumerate(reader, start=2):  # Start from row 2
                transcript = row['Matched_Transcript']
                database_gene = row['Database_Gene']
                
                if transcript in tpm_data:
                    # Build output row with database gene, transcript, and all TPM values
                    output_row = {
                        'Database_Gene': database_gene,
                        'Matched_Transcript': transcript
                    }
                    # Add all TPM columns dynamically
                    output_row.update(tpm_data[transcript])
                    output_data.append(output_row)
                else:
                    print(f"Warning: Transcript '{transcript}' not found in TPM data", file=sys.stderr)
    except Exception as e:
        print(f"Error reading final_annotated_genes_summary.csv: {e}", file=sys.stderr)
        return False
    
    # Write output TSV
    try:
        with open('blast_results/annotated_genes_tpm_values.tsv', 'w', newline='') as f:
            # Build fieldnames: Database_Gene, Matched_Transcript, then all TPM columns
            fieldnames = ['Database_Gene', 'Matched_Transcript'] + tpm_columns
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(output_data)
        
        print(f"✓ Successfully created annotated_genes_tpm_values.tsv with {len(output_data)} entries")
        return True
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        return False

if __name__ == '__main__':
    success = extract_tpm_values()
    sys.exit(0 if success else 1)
