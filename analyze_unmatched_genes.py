"""
Analyze unmatched database genes and their BLAST statistics.
Provides detailed information about genes that didn't meet filtering criteria.
Allows user to optionally include unmatched genes in final summary.
"""

import csv
import sys
import os
from collections import defaultdict

def analyze_unmatched_genes():
    """Analyze genes that didn't pass the filtering criteria."""
    
    # Read matched genes from the CSV
    matched_genes = set()
    try:
        with open('blast_results/final_annotated_genes_summary.csv', 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                matched_genes.add(row['Database_Gene'])
    except Exception as e:
        print(f"Error reading matched genes: {e}", file=sys.stderr)
        return False
    
    # Read database genes
    database_genes = set()
    try:
        with open('GmelOBPs.txt', 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    gene_id = line[1:].split()[0]  # Get ID without '>'
                    database_genes.add(gene_id)
    except Exception as e:
        print(f"Error reading database: {e}", file=sys.stderr)
        return False
    
    # Find unmatched genes
    unmatched_genes = database_genes - matched_genes
    
    if not unmatched_genes:
        print("All database genes were matched!")
        return True
    
    # Read BLAST results and find best matches for unmatched genes
    unmatched_stats = {}
    try:
        with open('blast_results/blast_results.tsv', 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 3:
                    continue
                sseqid = row[1]  # Database gene (subject)
                
                if sseqid in unmatched_genes:
                    if sseqid not in unmatched_stats:
                        unmatched_stats[sseqid] = {
                            'best_match': row[0],  # Query transcript
                            'identity': float(row[2]),
                            'evalue': float(row[10]),
                            'bitscore': float(row[11]),
                            'all_matches': []
                        }
                    
                    # Store all matches for this gene
                    unmatched_stats[sseqid]['all_matches'].append({
                        'query': row[0],
                        'identity': float(row[2]),
                        'evalue': float(row[10]),
                        'bitscore': float(row[11])
                    })
                    
                    # Update best match if this one has higher identity
                    if float(row[2]) > unmatched_stats[sseqid]['identity']:
                        unmatched_stats[sseqid]['best_match'] = row[0]
                        unmatched_stats[sseqid]['identity'] = float(row[2])
                        unmatched_stats[sseqid]['evalue'] = float(row[10])
                        unmatched_stats[sseqid]['bitscore'] = float(row[11])
    except Exception as e:
        print(f"Error reading BLAST results: {e}", file=sys.stderr)
        return False
    
    # Print report
    print("\n" + "="*70)
    print("UNMATCHED DATABASE GENES ANALYSIS")
    print("="*70)
    print(f"\nTotal unmatched genes: {len(unmatched_genes)}")
    print(f"Unmatched genes: {', '.join(sorted(unmatched_genes))}")
    
    if unmatched_stats:
        print("\n" + "-"*70)
        print("SEQUENCE IDENTITY STATISTICS FOR UNMATCHED GENES")
        print("-"*70)
        
        identities = [stats['identity'] for stats in unmatched_stats.values()]
        
        print(f"\nMean identity: {sum(identities)/len(identities):.2f}%")
        print(f"Min identity: {min(identities):.2f}%")
        print(f"Max identity: {max(identities):.2f}%")
        print(f"Median identity: {sorted(identities)[len(identities)//2]:.2f}%")
        
        print("\n" + "-"*70)
        print("DETAILED RESULTS FOR EACH UNMATCHED GENE")
        print("-"*70)
        
        for gene in sorted(unmatched_stats.keys()):
            stats = unmatched_stats[gene]
            print(f"\n{gene}:")
            print(f"  Best Match Transcript: {stats['best_match']}")
            print(f"  Identity: {stats['identity']:.2f}%")
            print(f"  E-value: {stats['evalue']:.2e}")
            print(f"  Bit Score: {stats['bitscore']:.1f}")
            print(f"  Total BLAST hits: {len(stats['all_matches'])}")
            
            # Show why it didn't match (below 98% threshold)
            if stats['identity'] < 98.0:
                print(f"  ⚠ Reason not matched: Identity {stats['identity']:.2f}% < 98% threshold")
    else:
        print("\nNo BLAST hits found for unmatched genes")
    
    print("\n" + "="*70)
    
    # Ask user if they want to include unmatched genes in final summary
    if unmatched_stats:
        print("INCLUDE UNMATCHED GENES IN FINAL SUMMARY?")
        print("="*70)
        
        for gene in sorted(unmatched_stats.keys()):
            stats = unmatched_stats[gene]
            identity_pct = stats['identity']
            best_transcript = stats['best_match']
            
            print(f"\n⚠️  Unmatched Gene: {gene}")
            print(f"   Your unmatched database gene(s) have {identity_pct:.2f}% of sequence identity")
            print(f"   (Best BLAST Hit) with {best_transcript} assembly header.")
            print(f"   Would you like to include it into your final_annotated_genes_summary?")
        
        print("\n" + "-"*70)
        response = input("\nEnter 'yes' or 'y' to include unmatched genes, any other key to skip: ").strip().lower()
        
        if response in ['yes', 'y']:
            success = add_unmatched_to_summary(unmatched_stats, matched_genes)
            return success
    
    print("\n" + "="*70 + "\n")
    return True


def add_unmatched_to_summary(unmatched_stats, matched_genes):
    """Add unmatched genes to the final summary CSV file."""
    
    csv_file = 'blast_results/final_annotated_genes_summary.csv'
    
    if not os.path.exists(csv_file):
        print(f"\n✗ Error: {csv_file} not found", file=sys.stderr)
        return False
    
    # Read existing data
    try:
        existing_rows = []
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            existing_rows = list(reader)
    except Exception as e:
        print(f"\n✗ Error reading CSV: {e}", file=sys.stderr)
        return False
    
    # Prepare new rows for unmatched genes
    # We need to read BLAST results to get full alignment details
    new_rows = []
    try:
        with open('blast_results/blast_results.tsv', 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 12:
                    continue
                sseqid = row[1]  # Database gene
                
                if sseqid in unmatched_stats:
                    # Check if we haven't already added this gene
                    if sseqid not in [r['Database_Gene'] for r in new_rows]:
                        stats = unmatched_stats[sseqid]
                        new_row = {
                            'Database_Gene': sseqid,
                            'Matched_Transcript': stats['best_match'],
                            'Identity_%': f"{stats['identity']:.2f}",
                            'E-value': f"{stats['evalue']:.2e}",
                            'Transcript_Length_bp': 'N/A'  # Would need to look up in transcriptome
                        }
                        new_rows.append(new_row)
    except Exception as e:
        print(f"\n✗ Error reading BLAST results: {e}", file=sys.stderr)
        return False
    
    if not new_rows:
        print("\n✓ No new unmatched genes to add")
        return True
    
    # Combine existing and new rows
    all_rows = existing_rows + new_rows
    
    # Sort by Database_Gene for consistency
    all_rows.sort(key=lambda x: x['Database_Gene'])
    
    # Write updated CSV
    try:
        with open(csv_file, 'w', newline='') as f:
            fieldnames = ['Database_Gene', 'Matched_Transcript', 'Identity_%', 'E-value', 'Transcript_Length_bp']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        
        print(f"\n✓ Successfully updated {csv_file}")
        print(f"  • Added {len(new_rows)} unmatched gene(s)")
        print(f"  • Total genes now: {len(all_rows)}")
        
        # Also update the results statistics file
        update_statistics_file(len(all_rows))
        
        return True
    except Exception as e:
        print(f"\n✗ Error writing CSV: {e}", file=sys.stderr)
        return False


def update_statistics_file(total_matched):
    """Update results_statistics.txt with new match count."""
    
    stats_file = 'blast_results/results_statistics.txt'
    
    try:
        with open(stats_file, 'r') as f:
            content = f.read()
        
        # Update match rate and count
        old_count = "Matched database genes: 24"
        new_count = f"Matched database genes: {total_matched}"
        content = content.replace(old_count, new_count)
        
        # Recalculate match rate (out of 25 total genes)
        match_rate = (total_matched / 25) * 100
        old_rate = "Match rate: 96.0%"
        new_rate = f"Match rate: {match_rate:.1f}%"
        content = content.replace(old_rate, new_rate)
        
        with open(stats_file, 'w') as f:
            f.write(content)
        
        print(f"  • Updated {stats_file} with new statistics")
        
    except Exception as e:
        print(f"  ⚠ Warning: Could not update statistics file: {e}", file=sys.stderr)

if __name__ == '__main__':
    success = analyze_unmatched_genes()
    sys.exit(0 if success else 1)
