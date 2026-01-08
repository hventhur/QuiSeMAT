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
    try:
        with open('GxH_TPM.tsv', 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                transcript_id = row['']  # First column (transcript ID)
                tpm_data[transcript_id] = {
                    'GR1': row['GR1'],
                    'GR2': row['GR2'],
                    'GR3': row['GR3'],
                    'HR1': row['HR1'],
                    'HR2': row['HR2'],
                    'HR3': row['HR3']
                }
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
                    output_data.append({
                        'Database_Gene': database_gene,
                        'Matched_Transcript': transcript,
                        'GR1': tpm_data[transcript]['GR1'],
                        'GR2': tpm_data[transcript]['GR2'],
                        'GR3': tpm_data[transcript]['GR3'],
                        'HR1': tpm_data[transcript]['HR1'],
                        'HR2': tpm_data[transcript]['HR2'],
                        'HR3': tpm_data[transcript]['HR3']
                    })
                else:
                    print(f"Warning: Transcript '{transcript}' not found in TPM data", file=sys.stderr)
    except Exception as e:
        print(f"Error reading final_annotated_genes_summary.csv: {e}", file=sys.stderr)
        return False
    
    # Write output TSV
    try:
        with open('blast_results/annotated_genes_tpm_values.tsv', 'w', newline='') as f:
            fieldnames = ['Database_Gene', 'Matched_Transcript', 'GR1', 'GR2', 'GR3', 'HR1', 'HR2', 'HR3']
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
