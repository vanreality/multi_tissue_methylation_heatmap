#!/usr/bin/env python3
"""
Calculate average methylation rates for DMR regions from a bedGraph file.

This script processes a bedGraph file containing methylation rates for CpG sites
along with their corresponding DMR regions, and calculates the average
methylation rate for each unique DMR region.
"""

import os
import sys
import click
import pandas as pd
import numpy as np
from rich.console import Console

console = Console()


@click.command()
@click.option('--bedgraph', required=True, help='Input bedGraph file with methylation rates and DMR regions')
@click.option('--output', required=True, help='Output prefix for the result file')
def main(bedgraph, output):
    """
    Calculate average methylation rates for DMR regions.

    Args:
        bedgraph: Path to the input bedGraph file containing methylation rates and DMR regions.
        output: Prefix for the output file.
    """
    try:
        console.print(f"[bold blue]Processing bedGraph file: {bedgraph}[/bold blue]")
        
        # Read bedGraph file
        # Columns: chr, start, end, meth_rate, chr_dmr, start_dmr, end_dmr
        column_names = ['chr', 'start', 'end', 'meth_rate', 'chr_dmr', 'start_dmr', 'end_dmr']
        
        # Check if file exists
        if not os.path.exists(bedgraph):
            console.print(f"[bold red]Error: Input file {bedgraph} does not exist.[/bold red]")
            sys.exit(1)
        
        # Read the file with proper column names
        with console.status("[bold green]Reading bedGraph file..."):
            try:
                df = pd.read_csv(bedgraph, sep='\t', header=None, names=column_names, usecols=[0, 1, 2, 3, 4, 5, 6])
            except Exception as e:
                console.print(f"[bold red]Error reading bedGraph file: {str(e)}[/bold red]")
                sys.exit(1)
        
        # Validate numeric columns
        numeric_cols = ['start', 'end', 'meth_rate', 'start_dmr', 'end_dmr']
        for col in numeric_cols:
            if not pd.api.types.is_numeric_dtype(df[col]):
                console.print(f"[bold yellow]Warning: Column '{col}' contains non-numeric values. Attempting to convert...[/bold yellow]")
                try:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                    df = df.dropna(subset=[col])
                except:
                    console.print(f"[bold red]Error: Failed to convert column '{col}' to numeric type.[/bold red]")
                    sys.exit(1)
        
        # Calculate average methylation rate for each DMR
        console.print("[bold blue]Calculating average methylation rates for DMRs...[/bold blue]")
        
        # Group by DMR regions and calculate mean methylation rate
        result_df = df.groupby(['chr_dmr', 'start_dmr', 'end_dmr']).agg({
            'meth_rate': 'mean'
        }).reset_index()
        
        # Write results to output file
        result_df.to_csv(output, sep='\t', header=False, index=False)
            
    except Exception as e:
        console.print(f"[bold red]Unexpected error: {str(e)}[/bold red]")
        sys.exit(1)


if __name__ == "__main__":
    main()
