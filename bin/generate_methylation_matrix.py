#!/usr/bin/env python3

"""
Generate a methylation matrix from bedGraph files.

This script merges multiple bedGraph files based on genomic regions from a BED file
and creates a consolidated methylation rate matrix.
"""

import os
import sys
import pandas as pd
import click
from rich.console import Console
from rich.panel import Panel
from tqdm import tqdm

console = Console()

@click.command()
@click.option('--bedgraph_list', required=True, help='CSV file with columns: sample,bedgraph_file_path')
@click.option('--bed', required=True, help='BED file with genomic regions')
@click.option('--output', required=True, help='Output file name')
def generate_methylation_matrix(bedgraph_list, bed, output):
    """
    Generate a methylation rate matrix from multiple bedGraph files.
    
    Args:
        bedgraph_list: Path to CSV file containing sample names and bedGraph file paths.
        bed: Path to BED file containing genomic regions of interest.
        output: Path to output file for the methylation matrix.
    """
    try:
        # Check if input files exist
        for file_path in [bedgraph_list, bed]:
            if not os.path.exists(file_path):
                console.print(Panel(f"[bold red]Error: File not found: {file_path}[/bold red]"))
                sys.exit(1)
        
        # Read the BED file
        console.print("[bold blue]Reading BED file...[/bold blue]")
        try:
            bed_df = pd.read_csv(bed, sep='\t', header=None, usecols=[0, 1, 2])
            bed_df.columns = ['chr', 'start', 'end']
        except Exception as e:
            console.print(Panel(f"[bold red]Error reading BED file: {str(e)}[/bold red]"))
            sys.exit(1)
        
        # Read the bedgraph list file
        console.print("[bold blue]Reading bedGraph list file...[/bold blue]")
        try:
            bedgraph_list_df = pd.read_csv(bedgraph_list)
            if 'sample' not in bedgraph_list_df.columns or 'bedgraph_file_path' not in bedgraph_list_df.columns:
                console.print(Panel("[bold red]Error: bedgraph_list must have 'sample' and 'bedgraph_file_path' columns[/bold red]"))
                sys.exit(1)
        except Exception as e:
            console.print(Panel(f"[bold red]Error reading bedGraph list file: {str(e)}[/bold red]"))
            sys.exit(1)
        
        # Initialize merged_df with the bed_df
        merged_df = bed_df.copy()
        
        # Process each bedGraph file
        console.print("[bold blue]Processing bedGraph files...[/bold blue]")
        for _, row in tqdm(bedgraph_list_df.iterrows(), total=len(bedgraph_list_df), desc="Processing samples"):
            sample = row['sample']
            bedgraph_path = row['bedgraph_file_path']
            
            if not os.path.exists(bedgraph_path):
                console.print(f"[yellow]Warning: bedGraph file not found: {bedgraph_path}. Skipping sample {sample}.[/yellow]")
                continue
            
            try:
                # Read bedGraph file
                bg_df = pd.read_csv(bedgraph_path, sep='\t', header=None, usecols=[0, 1, 2, 3])
                bg_df.columns = ['chr', 'start', 'end', 'meth_rate']
                
                # Merge with bed_df
                # For bedGraph, we need exact matches on chr, start, and end
                sample_df = pd.merge(
                    bed_df, 
                    bg_df, 
                    on=['chr', 'start', 'end'],
                    how='left'
                )
                
                # Rename the meth_rate column to the sample name
                sample_df.rename(columns={'meth_rate': sample}, inplace=True)
                
                # Merge with the main dataframe
                if len(merged_df.columns) == 3:  # Only has chr, start, end
                    merged_df = sample_df
                else:
                    merged_df = pd.merge(
                        merged_df,
                        sample_df[['chr', 'start', 'end', sample]],
                        on=['chr', 'start', 'end'],
                        how='outer'
                    )
            except Exception as e:
                console.print(f"[yellow]Error processing {bedgraph_path} for sample {sample}: {str(e)}. Skipping this sample.[/yellow]")
                continue
        
        # Write the merged matrix to file
        console.print("[bold blue]Writing output file...[/bold blue]")
        try:
            merged_df.to_csv(output, sep='\t', index=False, na_rep='NA')
            console.print(f"[bold green]Methylation matrix successfully written to {output}[/bold green]")
        except Exception as e:
            console.print(Panel(f"[bold red]Error writing output file: {str(e)}[/bold red]"))
            sys.exit(1)
            
    except Exception as e:
        console.print(Panel(f"[bold red]Unexpected error: {str(e)}[/bold red]"))
        sys.exit(1)

if __name__ == "__main__":
    generate_methylation_matrix()
