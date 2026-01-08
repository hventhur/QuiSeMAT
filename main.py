#!/usr/bin/env python3
"""
Quick Sequence Match Among Transcriptomes (QuiSeMAT)
A bioinformatics pipeline for searching annotated genes in transcriptome assemblies.

This script facilitates the search of already annotated genes from RNA sequencing 
projects (transcriptomes) by comparing them against a database of known genes.
"""

import argparse
import sys
import os
import subprocess
import tempfile
import json
from pathlib import Path
from typing import Tuple, Dict, List, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class FileValidator:
    """Handles validation of input FASTA files."""
    
    @staticmethod
    def validate_fasta_file(filepath: str, file_type: str) -> Tuple[bool, str]:
        """
        Validate that a file exists and is in valid FASTA format.
        
        Args:
            filepath: Path to the file to validate
            file_type: Description of the file type (e.g., 'transcriptome assembly')
            
        Returns:
            Tuple of (is_valid: bool, message: str)
        """
        # Check if file exists
        if not os.path.exists(filepath):
            return False, f"Error: {file_type} file '{filepath}' does not exist."
        
        # Check if file is readable
        if not os.access(filepath, os.R_OK):
            return False, f"Error: {file_type} file '{filepath}' is not readable."
        
        # Check if file is empty
        if os.path.getsize(filepath) == 0:
            return False, f"Error: {file_type} file '{filepath}' is empty."
        
        # Try to parse as FASTA
        try:
            record_count = 0
            for record in SeqIO.parse(filepath, "fasta"):
                if len(record.seq) == 0:
                    return False, f"Error: Record '{record.id}' in {file_type} file has empty sequence."
                record_count += 1
            
            if record_count == 0:
                return False, f"Error: {file_type} file '{filepath}' contains no valid FASTA records."
            
            return True, f"✓ {file_type} file '{filepath}' is valid ({record_count} sequences found)"
        
        except Exception as e:
            return False, f"Error: Failed to parse {file_type} file '{filepath}': {str(e)}"
    
    @staticmethod
    def get_fasta_stats(filepath: str) -> Dict[str, any]:
        """
        Get statistics about a FASTA file.
        
        Args:
            filepath: Path to the FASTA file
            
        Returns:
            Dictionary containing sequence count, min/max lengths, and total bases
        """
        stats = {
            'sequence_count': 0,
            'total_bases': 0,
            'min_length': float('inf'),
            'max_length': 0,
            'headers': []
        }
        
        for record in SeqIO.parse(filepath, "fasta"):
            stats['sequence_count'] += 1
            seq_len = len(record.seq)
            stats['total_bases'] += seq_len
            stats['min_length'] = min(stats['min_length'], seq_len)
            stats['max_length'] = max(stats['max_length'], seq_len)
            stats['headers'].append(record.id)
        
        # Handle case where file is empty
        if stats['sequence_count'] == 0:
            stats['min_length'] = 0
        
        return stats


class BLASTManager:
    """Handles BLAST database creation and homology searches."""
    
    def __init__(self, verbose: bool = False):
        """
        Initialize BLASTManager.
        
        Args:
            verbose: Enable verbose output
        """
        self.verbose = verbose
        self._check_blast_installation()
    
    def _check_blast_installation(self) -> bool:
        """
        Check if BLAST tools are installed.
        
        Returns:
            True if BLAST is installed, raises exception otherwise
        """
        try:
            subprocess.run(
                ['which', 'makeblastdb'],
                capture_output=True,
                check=True
            )
            subprocess.run(
                ['which', 'blastx'],
                capture_output=True,
                check=True
            )
            if self.verbose:
                print("✓ BLAST tools are available")
            return True
        except subprocess.CalledProcessError:
            raise RuntimeError(
                "Error: BLAST tools are not installed or not in PATH.\n"
                "Please install NCBI BLAST:\n"
                "  On macOS: brew install blast\n"
                "  On Linux: sudo apt-get install ncbi-blast+\n"
                "Or download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
            )
    
    def detect_sequence_type(self, filepath: str) -> str:
        """
        Detect if sequences are protein or nucleotide.
        
        Args:
            filepath: Path to FASTA file
            
        Returns:
            'protein' or 'nucleotide'
        """
        nucleotide_chars = set('ATGCNatgcn')
        protein_specific = set('EFIPQZ')
        
        for record in SeqIO.parse(filepath, "fasta"):
            seq_str = str(record.seq).upper()
            non_standard = set(seq_str) - nucleotide_chars
            
            # If we find protein-specific characters, it's protein
            if any(char in seq_str for char in protein_specific):
                return 'protein'
            
            # If sequences are predominantly not nucleotides, assume protein
            if len(non_standard) / len(seq_str) > 0.1:
                return 'protein'
            
            # Check first 3 sequences, then decide
            break
        
        return 'nucleotide'
    
    def create_blast_database(self, fasta_file: str, db_name: str, db_type: str) -> Tuple[bool, str]:
        """
        Create BLAST database from FASTA file.
        
        Args:
            fasta_file: Path to input FASTA file
            db_name: Name/path for the database (without extension)
            db_type: 'nucl' for nucleotide or 'prot' for protein
            
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            cmd = [
                'makeblastdb',
                '-in', fasta_file,
                '-dbtype', db_type,
                '-out', db_name,
                '-title', os.path.basename(db_name)
            ]
            
            if self.verbose:
                print(f"  Creating BLAST database: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )
            
            if result.returncode == 0:
                return True, f"✓ BLAST database created: {db_name}"
            else:
                return False, f"Error creating BLAST database: {result.stderr}"
        
        except Exception as e:
            return False, f"Error: {str(e)}"
    
    def run_blast_search(self, query_file: str, db_path: str, blast_type: str, 
                        output_file: str, evalue: float = 1e-5) -> Tuple[bool, str]:
        """
        Run BLAST search.
        
        Args:
            query_file: Path to query sequences
            db_path: Path to BLAST database
            blast_type: 'blastx' or 'blastn'
            output_file: Path to output file
            evalue: E-value threshold
            
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            cmd = [
                blast_type,
                '-query', query_file,
                '-db', db_path,
                '-out', output_file,
                '-outfmt', '6',  # Tabular output
                '-evalue', str(evalue),
                '-num_threads', '4'
            ]
            
            if self.verbose:
                print(f"  Running BLAST: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )
            
            if result.returncode == 0:
                return True, f"✓ BLAST search completed: {output_file}"
            else:
                return False, f"Error in BLAST search: {result.stderr}"
        
        except Exception as e:
            return False, f"Error: {str(e)}"
    
    def parse_blast_output(self, blast_output_file: str) -> Dict:
        """
        Parse BLAST tabular output.
        
        Args:
            blast_output_file: Path to BLAST output file
            
        Returns:
            Dictionary with parsed results
        """
        results = {
            'total_hits': 0,
            'unique_queries': set(),
            'unique_subjects': set(),
            'matches': []
        }
        
        try:
            with open(blast_output_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) >= 12:
                        match = {
                            'qseqid': fields[0],      # Query ID
                            'sseqid': fields[1],      # Subject ID
                            'pident': float(fields[2]),    # Percentage identity
                            'length': int(fields[3]),      # Alignment length
                            'mismatch': int(fields[4]),    # Mismatches
                            'gapopen': int(fields[5]),     # Gap opens
                            'qstart': int(fields[6]),      # Query start
                            'qend': int(fields[7]),        # Query end
                            'sstart': int(fields[8]),      # Subject start
                            'send': int(fields[9]),        # Subject end
                            'evalue': float(fields[10]),   # E-value
                            'bitscore': float(fields[11])  # Bit score
                        }
                        results['matches'].append(match)
                        results['unique_queries'].add(fields[0])
                        results['unique_subjects'].add(fields[1])
            
            results['total_hits'] = len(results['matches'])
            results['unique_queries'] = len(results['unique_queries'])
            results['unique_subjects'] = len(results['unique_subjects'])
            
            return results
        
        except Exception as e:
            raise RuntimeError(f"Error parsing BLAST output: {str(e)}")


class QuiSeMAT:
    """Main class for QuiSeMAT pipeline."""
    
    def __init__(self, transcriptome_file: str, database_file: str, verbose: bool = False,
                 blast_type: str = 'auto', output_dir: str = 'blast_results'):
        """
        Initialize QuiSeMAT with input files.
        
        Args:
            transcriptome_file: Path to transcriptome assembly FASTA file
            database_file: Path to database of annotated genes FASTA file
            verbose: Enable verbose output
            blast_type: Type of BLAST to use ('blastx', 'blastn', or 'auto')
            output_dir: Directory for output files
        """
        self.transcriptome_file = transcriptome_file
        self.database_file = database_file
        self.verbose = verbose
        self.blast_type = blast_type
        self.output_dir = output_dir
        self.validator = FileValidator()
        self.blast_manager = BLASTManager(verbose=verbose)
        self.transcriptome_stats = None
        self.database_stats = None
        self.blast_results = None
        self.transcriptome_sequences = {}  # Store transcriptome sequences for filtering
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
    
    def step_one_validate_files(self) -> bool:
        """
        STEP ONE: Validate input files.
        
        This step validates that:
        1. The transcriptome assembly file (FASTA format) exists and is readable
        2. The database file with annotated genes (FASTA format) exists and is readable
        3. Both files contain valid FASTA records with non-empty sequences
        
        Returns:
            True if both files are valid, False otherwise
        """
        print("\n" + "="*70)
        print("STEP ONE: File Validation")
        print("="*70)
        
        # Validate transcriptome file
        print("\nValidating transcriptome assembly file...")
        is_valid_trans, msg_trans = self.validator.validate_fasta_file(
            self.transcriptome_file, 
            "Transcriptome assembly"
        )
        print(msg_trans)
        
        if not is_valid_trans:
            return False
        
        # Validate database file
        print("\nValidating database file (annotated genes)...")
        is_valid_db, msg_db = self.validator.validate_fasta_file(
            self.database_file,
            "Database"
        )
        print(msg_db)
        
        if not is_valid_db:
            return False
        
        # Get statistics for both files
        print("\n" + "-"*70)
        print("File Statistics")
        print("-"*70)
        
        self.transcriptome_stats = self.validator.get_fasta_stats(self.transcriptome_file)
        self.database_stats = self.validator.get_fasta_stats(self.database_file)
        
        # Print transcriptome statistics
        print(f"\nTranscriptome Assembly ({self.transcriptome_file}):")
        print(f"  • Total sequences: {self.transcriptome_stats['sequence_count']}")
        print(f"  • Total bases: {self.transcriptome_stats['total_bases']:,}")
        print(f"  • Sequence length range: {self.transcriptome_stats['min_length']:,} - {self.transcriptome_stats['max_length']:,} bp")
        
        # Print database statistics
        print(f"\nDatabase ({self.database_file}):")
        print(f"  • Total sequences: {self.database_stats['sequence_count']}")
        print(f"  • Total bases: {self.database_stats['total_bases']:,}")
        print(f"  • Sequence length range: {self.database_stats['min_length']:,} - {self.database_stats['max_length']:,} bp")
        
        if self.verbose:
            print("\n" + "-"*70)
            print("Database Headers (first 10):")
            print("-"*70)
            for i, header in enumerate(self.database_stats['headers'][:10], 1):
                print(f"  {i}. {header}")
            
            if len(self.database_stats['headers']) > 10:
                print(f"  ... and {len(self.database_stats['headers']) - 10} more")
        
        print("\n" + "="*70)
        print("✓ STEP ONE COMPLETED SUCCESSFULLY")
        print("="*70)
        return True
    
    def step_two_homology_search(self) -> bool:
        """
        STEP TWO: Perform homology search using BLAST.
        
        This step:
        1. Detects sequence types (nucleotide or protein)
        2. Creates BLAST database from annotated genes
        3. Performs BLASTx (protein-translated) or BLASTn search
        4. Parses and stores results
        
        Returns:
            True if search is successful, False otherwise
        """
        print("\n" + "="*70)
        print("STEP TWO: Homology Search using BLAST")
        print("="*70)
        
        # Step 2.1: Detect sequence types
        print("\nDetecting sequence types...")
        trans_type = self.blast_manager.detect_sequence_type(self.transcriptome_file)
        db_type = self.blast_manager.detect_sequence_type(self.database_file)
        
        print(f"  • Transcriptome sequences: {trans_type}")
        print(f"  • Database sequences: {db_type}")
        
        # Determine BLAST type to use
        if self.blast_type == 'auto':
            if trans_type == 'nucleotide' and db_type == 'protein':
                blast_to_use = 'blastx'
                db_blast_type = 'prot'
            elif trans_type == 'nucleotide' and db_type == 'nucleotide':
                blast_to_use = 'blastn'
                db_blast_type = 'nucl'
            elif trans_type == 'protein' and db_type == 'protein':
                blast_to_use = 'blastp'
                db_blast_type = 'prot'
            else:
                # Default to BLASTx if mixed
                blast_to_use = 'blastx'
                db_blast_type = 'prot' if db_type == 'protein' else 'nucl'
        else:
            blast_to_use = self.blast_type
            db_blast_type = 'prot' if db_type == 'protein' else 'nucl'
        
        print(f"  • Selected BLAST method: {blast_to_use}")
        
        # Step 2.2: Create BLAST database
        print("\nCreating BLAST database from annotated genes...")
        db_path = os.path.join(self.output_dir, 'annotated_genes_db')
        success, msg = self.blast_manager.create_blast_database(
            self.database_file,
            db_path,
            db_blast_type
        )
        print(msg)
        
        if not success:
            print("Error: Failed to create BLAST database")
            return False
        
        # Step 2.3: Run BLAST search
        print(f"\nRunning {blast_to_use} search...")
        blast_output = os.path.join(self.output_dir, 'blast_results.tsv')
        
        success, msg = self.blast_manager.run_blast_search(
            self.transcriptome_file,
            db_path,
            blast_to_use,
            blast_output,
            evalue=1e-5
        )
        print(msg)
        
        if not success:
            print("Error: BLAST search failed")
            return False
        
        # Step 2.4: Parse BLAST results
        print("\nParsing BLAST results...")
        try:
            self.blast_results = self.blast_manager.parse_blast_output(blast_output)
            
            print(f"  • Total hits found: {self.blast_results['total_hits']}")
            print(f"  • Unique query sequences matched: {self.blast_results['unique_queries']}")
            print(f"  • Unique database sequences matched: {self.blast_results['unique_subjects']}")
            
            # Show top matches
            if self.verbose and self.blast_results['matches']:
                print("\n" + "-"*70)
                print("Top 10 BLAST Hits (by bit score):")
                print("-"*70)
                
                sorted_matches = sorted(
                    self.blast_results['matches'],
                    key=lambda x: x['bitscore'],
                    reverse=True
                )[:10]
                
                for i, match in enumerate(sorted_matches, 1):
                    print(f"\n  {i}. Query: {match['qseqid']} → Subject: {match['sseqid']}")
                    print(f"     Identity: {match['pident']:.2f}% | Length: {match['length']} bp | E-value: {match['evalue']:.2e}")
            
        except Exception as e:
            print(f"Error parsing BLAST results: {str(e)}")
            return False
        
        # Save BLAST results as JSON for Step Three
        results_json = os.path.join(self.output_dir, 'blast_results.json')
        try:
            json_data = {
                'total_hits': self.blast_results['total_hits'],
                'unique_queries': self.blast_results['unique_queries'],
                'unique_subjects': self.blast_results['unique_subjects'],
                'matches': self.blast_results['matches'],
                'blast_method': blast_to_use,
                'database_file': self.database_file,
                'transcriptome_file': self.transcriptome_file
            }
            
            with open(results_json, 'w') as f:
                json.dump(json_data, f, indent=2)
            
            print(f"\n✓ BLAST results saved to: {results_json}")
        
        except Exception as e:
            print(f"Warning: Could not save JSON results: {str(e)}")
        
        print("\n" + "="*70)
        print("✓ STEP TWO COMPLETED SUCCESSFULLY")
        print("="*70)
        return True
    
    def step_three_filter_isoforms(self) -> bool:
        """
        STEP THREE: Filter isoforms and select best matches.
        
        This step:
        1. Groups BLAST results by database gene (subject)
        2. For each database gene, identifies all transcriptome matches (potential isoforms)
        3. Filters matches: keeps only those ≥ 98% sequence identity
        4. Among filtered matches, selects the longest transcriptome sequence
        5. Stores filtered results for Step Four
        
        Returns:
            True if filtering is successful, False otherwise
        """
        print("\n" + "="*70)
        print("STEP THREE: Filter Isoforms and Select Best Matches")
        print("="*70)
        
        if not self.blast_results or not self.blast_results['matches']:
            print("Error: No BLAST results to filter. Run Step Two first.")
            return False
        
        # Step 3.1: Load transcriptome sequences to get their lengths
        print("\nLoading transcriptome sequence lengths...")
        try:
            for record in SeqIO.parse(self.transcriptome_file, "fasta"):
                self.transcriptome_sequences[record.id] = len(record.seq)
            print(f"  • Loaded {len(self.transcriptome_sequences)} transcriptome sequences")
        except Exception as e:
            print(f"Error loading transcriptome sequences: {str(e)}")
            return False
        
        # Step 3.2: Group matches by database gene (sseqid)
        print("\nGrouping matches by database gene...")
        db_gene_groups = {}
        
        for match in self.blast_results['matches']:
            db_gene = match['sseqid']
            if db_gene not in db_gene_groups:
                db_gene_groups[db_gene] = []
            db_gene_groups[db_gene].append(match)
        
        print(f"  • Found {len(db_gene_groups)} unique database genes with matches")
        
        # Step 3.3: Filter and select best match per database gene
        print("\nFiltering matches (≥98% sequence identity) and selecting longest transcripts...")
        
        IDENTITY_THRESHOLD = 98.0  # 98% as per instructions
        filtered_results = {}
        filtering_stats = {
            'total_groups': len(db_gene_groups),
            'passed_threshold': 0,
            'filtered_out': 0
        }
        
        for db_gene, matches in db_gene_groups.items():
            # Filter by identity threshold
            high_identity_matches = [m for m in matches if m['pident'] >= IDENTITY_THRESHOLD]
            
            if high_identity_matches:
                # Sort by query sequence length (descending) to get longest
                sorted_by_length = sorted(
                    high_identity_matches,
                    key=lambda m: self.transcriptome_sequences.get(m['qseqid'], 0),
                    reverse=True
                )
                
                # Select the longest match
                best_match = sorted_by_length[0]
                filtered_results[db_gene] = best_match
                filtering_stats['passed_threshold'] += 1
                
                if self.verbose:
                    transcript_len = self.transcriptome_sequences.get(best_match['qseqid'], 0)
                    print(f"\n  {db_gene}:")
                    print(f"    Matches found: {len(matches)}")
                    print(f"    High identity (≥98%): {len(high_identity_matches)}")
                    print(f"    Selected: {best_match['qseqid']} (length: {transcript_len} bp)")
                    print(f"    Identity: {best_match['pident']:.2f}% | E-value: {best_match['evalue']:.2e}")
            else:
                filtering_stats['filtered_out'] += 1
                if self.verbose:
                    print(f"\n  {db_gene}: No matches with ≥98% identity (best: {max(m['pident'] for m in matches):.2f}%)")
        
        # Step 3.4: Print filtering summary
        print("\n" + "-"*70)
        print("Filtering Summary")
        print("-"*70)
        print(f"Total database genes: {filtering_stats['total_groups']}")
        print(f"Passed filtering (≥98% identity): {filtering_stats['passed_threshold']}")
        print(f"Filtered out (no ≥98% identity match): {filtering_stats['filtered_out']}")
        
        # Store filtered results for Step Four
        self.filtered_results = filtered_results
        
        print("\n" + "="*70)
        print("✓ STEP THREE COMPLETED SUCCESSFULLY")
        print("="*70)
        return True
    
    def step_four_generate_final_report(self) -> bool:
        """
        STEP FOUR: Generate final results report.
        
        This step:
        1. Takes filtered results from Step Three
        2. Generates a final TSV file with matched headers
        3. Ensures each database gene is represented exactly once (as per instructions)
        4. Includes alignment statistics and sequence information
        
        Returns:
            True if report generation is successful, False otherwise
        """
        print("\n" + "="*70)
        print("STEP FOUR: Generate Final Results Report")
        print("="*70)
        
        if not hasattr(self, 'filtered_results'):
            print("Error: No filtered results available. Run Step Three first.")
            return False
        
        if not self.database_stats:
            print("Error: Database stats not available. Run Step One first.")
            return False
        
        # Step 4.1: Check completeness
        print("\nGenerating final report...")
        
        num_db_genes = self.database_stats['sequence_count']
        num_matched = len(self.filtered_results)
        
        print(f"  • Database genes: {num_db_genes}")
        print(f"  • Matched genes: {num_matched}")
        
        if num_matched < num_db_genes:
            unmatched = set(self.database_stats['headers']) - set(self.filtered_results.keys())
            print(f"  • Unmatched genes: {num_db_genes - num_matched}")
            if self.verbose:
                print("  • Unmatched database genes:")
                for gene in sorted(unmatched):
                    print(f"    - {gene}")
        
        # Step 4.2: Generate final TSV file
        print("\nWriting final results TSV file...")
        
        output_file = os.path.join(self.output_dir, 'final_annotated_genes.tsv')
        
        try:
            with open(output_file, 'w') as f:
                # Write header
                header = [
                    'Database_Gene',
                    'Matched_Transcript',
                    'Sequence_Identity_%',
                    'Alignment_Length_bp',
                    'Mismatches',
                    'Gap_Opens',
                    'Query_Start',
                    'Query_End',
                    'Subject_Start',
                    'Subject_End',
                    'E_value',
                    'Bit_Score',
                    'Transcript_Length_bp'
                ]
                f.write('\t'.join(header) + '\n')
                
                # Write data rows, sorted by database gene name
                for db_gene in sorted(self.filtered_results.keys()):
                    match = self.filtered_results[db_gene]
                    transcript_len = self.transcriptome_sequences.get(match['qseqid'], 0)
                    
                    row = [
                        db_gene,
                        match['qseqid'],
                        f"{match['pident']:.3f}",
                        str(match['length']),
                        str(match['mismatch']),
                        str(match['gapopen']),
                        str(match['qstart']),
                        str(match['qend']),
                        str(match['sstart']),
                        str(match['send']),
                        f"{match['evalue']:.2e}",
                        f"{match['bitscore']:.1f}",
                        str(transcript_len)
                    ]
                    f.write('\t'.join(row) + '\n')
            
            print(f"✓ Final results written to: {output_file}")
            
        except Exception as e:
            print(f"Error writing results file: {str(e)}")
            return False
        
        # Step 4.3: Generate summary CSV file (optional, for easier viewing)
        print("\nWriting summary CSV file...")
        
        summary_file = os.path.join(self.output_dir, 'final_annotated_genes_summary.csv')
        
        try:
            with open(summary_file, 'w') as f:
                # Write header
                header = 'Database_Gene,Matched_Transcript,Identity_%,E-value,Transcript_Length_bp\n'
                f.write(header)
                
                # Write data rows
                for db_gene in sorted(self.filtered_results.keys()):
                    match = self.filtered_results[db_gene]
                    transcript_len = self.transcriptome_sequences.get(match['qseqid'], 0)
                    
                    row = f"{db_gene},{match['qseqid']},{match['pident']:.2f},{match['evalue']:.2e},{transcript_len}\n"
                    f.write(row)
            
            print(f"✓ Summary CSV written to: {summary_file}")
            
        except Exception as e:
            print(f"Warning: Could not write summary CSV: {str(e)}")
        
        # Step 4.4: Generate detailed statistics report
        print("\nGenerating detailed statistics report...")
        
        stats_file = os.path.join(self.output_dir, 'results_statistics.txt')
        
        try:
            identities = [m['pident'] for m in self.filtered_results.values()]
            evalues = [m['evalue'] for m in self.filtered_results.values()]
            bitscores = [m['bitscore'] for m in self.filtered_results.values()]
            transcript_lengths = [
                self.transcriptome_sequences.get(m['qseqid'], 0) 
                for m in self.filtered_results.values()
            ]
            
            with open(stats_file, 'w') as f:
                f.write("="*70 + "\n")
                f.write("QuiSeMAT Pipeline Results Statistics\n")
                f.write("="*70 + "\n\n")
                
                f.write("INPUT FILES\n")
                f.write("-"*70 + "\n")
                f.write(f"Transcriptome Assembly: {self.transcriptome_file}\n")
                f.write(f"Database (Annotated Genes): {self.database_file}\n\n")
                
                f.write("FILE STATISTICS\n")
                f.write("-"*70 + "\n")
                f.write(f"Transcriptome sequences: {self.transcriptome_stats['sequence_count']}\n")
                f.write(f"Database sequences: {self.database_stats['sequence_count']}\n\n")
                
                f.write("MATCHING RESULTS\n")
                f.write("-"*70 + "\n")
                f.write(f"Matched database genes: {len(self.filtered_results)}\n")
                f.write(f"Unmatched database genes: {num_db_genes - num_matched}\n")
                f.write(f"Match rate: {(num_matched/num_db_genes)*100:.1f}%\n\n")
                
                f.write("SEQUENCE IDENTITY STATISTICS\n")
                f.write("-"*70 + "\n")
                f.write(f"Mean identity: {sum(identities)/len(identities):.2f}%\n")
                f.write(f"Min identity: {min(identities):.2f}%\n")
                f.write(f"Max identity: {max(identities):.2f}%\n")
                f.write(f"Median identity: {sorted(identities)[len(identities)//2]:.2f}%\n\n")
                
                f.write("E-VALUE STATISTICS\n")
                f.write("-"*70 + "\n")
                f.write(f"Mean e-value: {sum(evalues)/len(evalues):.2e}\n")
                f.write(f"Min e-value: {min(evalues):.2e}\n")
                f.write(f"Max e-value: {max(evalues):.2e}\n\n")
                
                f.write("BIT SCORE STATISTICS\n")
                f.write("-"*70 + "\n")
                f.write(f"Mean bit score: {sum(bitscores)/len(bitscores):.1f}\n")
                f.write(f"Min bit score: {min(bitscores):.1f}\n")
                f.write(f"Max bit score: {max(bitscores):.1f}\n\n")
                
                f.write("TRANSCRIPT LENGTH STATISTICS\n")
                f.write("-"*70 + "\n")
                f.write(f"Mean length: {sum(transcript_lengths)/len(transcript_lengths):.0f} bp\n")
                f.write(f"Min length: {min(transcript_lengths)} bp\n")
                f.write(f"Max length: {max(transcript_lengths)} bp\n\n")
            
            print(f"✓ Statistics report written to: {stats_file}")
            
        except Exception as e:
            print(f"Warning: Could not write statistics report: {str(e)}")
        
        # Step 4.5: Print final summary
        print("\n" + "-"*70)
        print("FINAL RESULTS SUMMARY")
        print("-"*70)
        print(f"✓ Total matched genes: {len(self.filtered_results)}/{num_db_genes}")
        print(f"✓ Output files generated:")
        print(f"  1. {output_file}")
        print(f"  2. {summary_file}")
        print(f"  3. {stats_file}")
        
        print("\n" + "="*70)
        print("✓ STEP FOUR COMPLETED SUCCESSFULLY")
        print("="*70)
        return True


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser.
    
    Returns:
        Configured ArgumentParser object
    """
    parser = argparse.ArgumentParser(
        prog='QuiSeMAT',
        description='''Quick Sequence Match Among Transcriptomes (QuiSeMAT) - 
A bioinformatics pipeline for searching annotated genes in transcriptome assemblies.''',
        epilog='''Examples:
  python main.py -t assembly.fasta -d genes.fasta
  python main.py --transcriptome assembly.fasta --database genes.fasta --verbose''',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '-t', '--transcriptome',
        required=True,
        type=str,
        help='Path to transcriptome assembly FASTA file (file one)'
    )
    
    parser.add_argument(
        '-d', '--database',
        required=True,
        type=str,
        help='Path to database of annotated genes FASTA file (file two)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    parser.add_argument(
        '--step',
        type=int,
        choices=[1, 2, 3, 4],
        default=4,
        help='Steps to execute: 1 (validation only), 2 (+ BLAST search), 3 (+ isoform filtering), or 4 (full pipeline). Default: 4'
    )
    
    parser.add_argument(
        '--blast-type',
        type=str,
        choices=['blastx', 'blastn', 'blastp', 'auto'],
        default='auto',
        help='BLAST method to use. auto=automatically detect. Default: auto'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='blast_results',
        help='Directory for output files. Default: blast_results'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s 0.4.0 (Steps One, Two, Three & Four Implementation)'
    )
    
    return parser


def main():
    """Main entry point for QuiSeMAT."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    try:
        # Initialize QuiSeMAT
        quismat = QuiSeMAT(
            transcriptome_file=args.transcriptome,
            database_file=args.database,
            verbose=args.verbose,
            blast_type=args.blast_type,
            output_dir=args.output_dir
        )
        
        # Execute Step One: Validate Files
        if not quismat.step_one_validate_files():
            sys.exit(1)
        
        # Execute Step Two if requested
        if args.step >= 2:
            if not quismat.step_two_homology_search():
                sys.exit(1)
        
        # Execute Step Three if requested
        if args.step >= 3:
            if not quismat.step_three_filter_isoforms():
                sys.exit(1)
        
        # Execute Step Four if requested
        if args.step >= 4:
            if not quismat.step_four_generate_final_report():
                sys.exit(1)
        
        sys.exit(0)
    
    except KeyboardInterrupt:
        print("\n\nExecution interrupted by user.")
        sys.exit(130)
    
    except Exception as e:
        print(f"\nFatal error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
