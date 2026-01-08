#!/usr/bin/env python3
"""
Generate TPM Expression Heatmap
================================
Generates professional heatmaps from TPM (Transcripts Per Million) expression data
with log-scale transformation for better visualization of wide value ranges.

Usage:
    python generate_tpm_heatmap.py [OPTIONS]

Options:
    --input FILE            Input TSV file (default: blast_results/annotated_genes_tpm_values.tsv)
    --output FILE           Output image file (default: blast_results/tpm_heatmap.png)
    --format FORMAT         Export format: png, svg, or both (default: png)
    --cmap COLORMAP         Color palette (default: viridis)
    --dpi DPI               Image resolution in DPI (default: 300, PNG only)
    --figsize WIDTHxHEIGHT  Figure size (default: 12x14)
    --annot                 Show values in cells
    --cluster               Enable hierarchical clustering
    --log-scale             Apply log2 transformation (default: True)
    --verbose               Show detailed output
    --help                  Show this help message

Supported colormaps:
    viridis, plasma, inferno, magma (perceptually uniform)
    coolwarm, RdYlGn, RdBu, PiYG (diverging)
    YlOrRd, YlGnBu, Blues, Greens (sequential)

Examples:
    # Default: PNG format with clustering
    python generate_tpm_heatmap.py

    # Export in both PNG and SVG formats
    python generate_tpm_heatmap.py --format both

    # SVG format only (for publications)
    python generate_tpm_heatmap.py --format svg --verbose

    # Custom colormap with both formats
    python generate_tpm_heatmap.py --cmap coolwarm --figsize 14x16 --format both
"""

import argparse
import sys
import os
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def validate_file(file_path):
    """Validate that input file exists and is readable."""
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    if not path.is_file():
        raise ValueError(f"Path is not a file: {file_path}")
    return path


def load_tpm_data(input_file, verbose=False):
    """Load and prepare TPM data for heatmap visualization."""
    if verbose:
        print(f"📖 Loading TPM data from: {input_file}")
    
    df = pd.read_csv(input_file, sep='\t')
    
    if verbose:
        print(f"   ✓ Loaded {len(df)} genes")
        print(f"   ✓ Columns: {list(df.columns)}")
    
    # Set Database_Gene as index for row labels
    df_heatmap = df.set_index('Database_Gene')
    
    # Remove the Matched_Transcript column
    if 'Matched_Transcript' in df_heatmap.columns:
        df_heatmap = df_heatmap.drop('Matched_Transcript', axis=1)
    
    if verbose:
        print(f"   ✓ TPM samples: {list(df_heatmap.columns)}")
        print(f"   ✓ TPM range: {df_heatmap.values.min():.2f} - {df_heatmap.values.max():.2f}")
    
    return df_heatmap


def apply_log_transformation(df_heatmap, verbose=False):
    """Apply log2 transformation with pseudocount to handle zero values."""
    if verbose:
        print("🔄 Applying log2 transformation (adding pseudocount of 1)...")
    
    df_log = np.log2(df_heatmap + 1)
    
    if verbose:
        print(f"   ✓ Transformed range: {df_log.values.min():.2f} - {df_log.values.max():.2f}")
    
    return df_log


def save_heatmap(fig, output_file, formats='png', dpi=300, verbose=False):
    """Save heatmap figure in specified format(s)."""
    # Normalize format input
    if isinstance(formats, str):
        formats = [f.lower().strip() for f in formats.split(',')]
    
    valid_formats = {'png', 'svg', 'both'}
    if 'both' in formats:
        formats = ['png', 'svg']
    
    # Validate formats
    for fmt in formats:
        if fmt not in valid_formats and fmt != 'both':
            raise ValueError(f"Unsupported format: {fmt}. Choose from: png, svg, both")
    
    # Remove duplicates and sort for consistent output
    formats = sorted(set(formats))
    
    output_path = Path(output_file)
    output_stem = output_path.stem
    output_dir = output_path.parent
    
    saved_files = []
    
    for fmt in formats:
        if fmt == 'png':
            output_file_fmt = output_dir / f"{output_stem}.png"
            if verbose:
                print(f"💾 Saving PNG to: {output_file_fmt}")
            fig.savefig(str(output_file_fmt), dpi=dpi, bbox_inches='tight', format='png')
            saved_files.append(output_file_fmt)
            
            if verbose:
                file_size = output_file_fmt.stat().st_size / (1024 * 1024)
                print(f"   ✓ PNG saved successfully ({file_size:.2f} MB)")
        
        elif fmt == 'svg':
            output_file_fmt = output_dir / f"{output_stem}.svg"
            if verbose:
                print(f"💾 Saving SVG to: {output_file_fmt}")
            fig.savefig(str(output_file_fmt), bbox_inches='tight', format='svg')
            saved_files.append(output_file_fmt)
            
            if verbose:
                file_size = output_file_fmt.stat().st_size / (1024)
                print(f"   ✓ SVG saved successfully ({file_size:.2f} KB)")
    
    return saved_files


def create_heatmap(df_data, output_file, cmap='viridis', figsize=(12, 14), 
                   dpi=300, annot=False, cluster=True, formats='png', verbose=False):
    """Create and save heatmap visualization in specified format(s)."""
    if verbose:
        print("🎨 Creating heatmap...")
        print(f"   ✓ Colormap: {cmap}")
        print(f"   ✓ Figure size: {figsize[0]}x{figsize[1]} inches")
        print(f"   ✓ DPI: {dpi}")
        print(f"   ✓ Clustering: {cluster}")
        print(f"   ✓ Show values: {annot}")
        print(f"   ✓ Export formats: {formats}")
    
    # Set style
    sns.set_style("whitegrid")
    
    if cluster:
        # Use clustermap for hierarchical clustering
        if verbose:
            print("   ✓ Performing hierarchical clustering...")
        
        g = sns.clustermap(
            df_data,
            cmap=cmap,
            figsize=figsize,
            cbar_kws={'label': 'log2(TPM + 1)'},
            linewidths=0.5,
            linecolor='white',
            method='ward',  # Linkage method for clustering
            metric='euclidean',  # Distance metric
            annot=annot,
            fmt='.2f' if annot else '',
            xticklabels=True,
            yticklabels=True,
            dendrogram_ratio=(.1, .2)  # Dendrogram size
        )
        
        # Improve labels and title
        g.ax_heatmap.set_xlabel('Sample', fontsize=12, fontweight='bold')
        g.ax_heatmap.set_ylabel('Database Gene', fontsize=12, fontweight='bold')
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(),
            rotation=45,
            ha='right',
            fontsize=10
        )
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_yticklabels(),
            rotation=0,
            fontsize=9
        )
        
        plt.suptitle(
            'Gene Expression Heatmap (log2-transformed TPM with Hierarchical Clustering)',
            fontsize=14,
            fontweight='bold',
            y=0.98
        )
        
        fig = g.fig
    else:
        # Use standard heatmap without clustering
        fig, ax = plt.subplots(figsize=figsize)
        
        sns.heatmap(
            df_data,
            cmap=cmap,
            cbar_kws={'label': 'log2(TPM + 1)', 'shrink': 0.8},
            linewidths=0.5,
            linecolor='white',
            annot=annot,
            fmt='.2f' if annot else '',
            ax=ax,
            square=False
        )
        
        ax.set_xlabel('Sample', fontsize=12, fontweight='bold')
        ax.set_ylabel('Database Gene', fontsize=12, fontweight='bold')
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            ha='right',
            fontsize=10
        )
        ax.set_yticklabels(
            ax.get_yticklabels(),
            rotation=0,
            fontsize=9
        )
        
        plt.title(
            'Gene Expression Heatmap (log2-transformed TPM)',
            fontsize=14,
            fontweight='bold',
            pad=20
        )
    
    # Adjust layout
    plt.tight_layout()
    
    # Save in specified format(s)
    saved_files = save_heatmap(fig, output_file, formats=formats, dpi=dpi, verbose=verbose)
    
    plt.close()
    
    return saved_files


def print_summary(df_heatmap, df_log, verbose=False):
    """Print summary statistics of the TPM data."""
    if not verbose:
        return
    
    print("\n📊 TPM Data Summary:")
    print(f"   Genes: {len(df_heatmap)}")
    print(f"   Samples: {len(df_heatmap.columns)}")
    print(f"\n   Original TPM values:")
    print(f"     Min: {df_heatmap.values.min():.2f}")
    print(f"     Max: {df_heatmap.values.max():.2f}")
    print(f"     Mean: {df_heatmap.values.mean():.2f}")
    print(f"     Median: {np.median(df_heatmap.values):.2f}")
    print(f"\n   Log2-transformed values:")
    print(f"     Min: {df_log.values.min():.2f}")
    print(f"     Max: {df_log.values.max():.2f}")
    print(f"     Mean: {df_log.values.mean():.2f}")
    print(f"     Median: {np.median(df_log.values):.2f}")


def parse_figsize(figsize_str):
    """Parse figsize string in format 'WIDTHxHEIGHT'."""
    try:
        width, height = map(float, figsize_str.lower().split('x'))
        return (width, height)
    except (ValueError, AttributeError):
        raise ValueError(f"Invalid figsize format: {figsize_str}. Use format: WIDTHxHEIGHT (e.g., 12x14)")


def main():
    parser = argparse.ArgumentParser(
        description='Generate professional TPM expression heatmaps with log-scale transformation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_tpm_heatmap.py
  python generate_tpm_heatmap.py --cmap coolwarm --figsize 14x16 --verbose
  python generate_tpm_heatmap.py --annot --no-cluster --dpi 150
        """
    )
    
    parser.add_argument(
        '--input',
        type=str,
        default='blast_results/annotated_genes_tpm_values.tsv',
        help='Input TSV file with TPM values (default: blast_results/annotated_genes_tpm_values.tsv)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='blast_results/tpm_heatmap.png',
        help='Output image file path (default: blast_results/tpm_heatmap.png)'
    )
    
    parser.add_argument(
        '--format',
        type=str,
        default='png',
        help='Export format: png, svg, or both (default: png)'
    )
    
    parser.add_argument(
        '--cmap',
        type=str,
        default='viridis',
        help='Color palette (default: viridis). Options: viridis, plasma, inferno, magma, coolwarm, RdYlGn, RdBu, PiYG, YlOrRd, YlGnBu, Blues, Greens'
    )
    
    parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='Image resolution in DPI (default: 300)'
    )
    
    parser.add_argument(
        '--figsize',
        type=str,
        default='12x14',
        help='Figure size as WIDTHxHEIGHT (default: 12x14)'
    )
    
    parser.add_argument(
        '--annot',
        action='store_true',
        help='Show TPM values in heatmap cells'
    )
    
    parser.add_argument(
        '--no-cluster',
        action='store_true',
        help='Disable hierarchical clustering'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show detailed output messages'
    )
    
    args = parser.parse_args()
    
    try:
        # Validate input file
        input_file = validate_file(args.input)
        
        # Parse figsize
        figsize = parse_figsize(args.figsize)
        
        # Create output directory if needed
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if args.verbose:
            print("=" * 60)
            print("🔬 TPM Expression Heatmap Generator")
            print("=" * 60)
        
        # Load data
        df_heatmap = load_tpm_data(str(input_file), verbose=args.verbose)
        
        # Apply transformation
        df_log = apply_log_transformation(df_heatmap, verbose=args.verbose)
        
        # Print summary
        print_summary(df_heatmap, df_log, verbose=args.verbose)
        
        # Create heatmap
        saved_files = create_heatmap(
            df_log,
            str(output_path),
            cmap=args.cmap,
            figsize=figsize,
            dpi=args.dpi,
            annot=args.annot,
            cluster=not args.no_cluster,
            formats=args.format,
            verbose=args.verbose
        )
        
        if args.verbose:
            print(f"\n📁 Generated files:")
            for f in saved_files:
                print(f"   ✓ {f}")
            print("\n✅ Heatmap generation completed successfully!")
            print("=" * 60)
        else:
            print(f"✅ Heatmap saved successfully!")
            for f in saved_files:
                print(f"   {f}")
    
    except FileNotFoundError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
