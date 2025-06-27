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
@click.option('--mode', type=click.Choice(['rate', 'count']), default='rate', help="Calculation mode: 'rate' for average of rates, 'count' for sum of counts.")
def main(bedgraph, output, mode):
    """
    Calculate average methylation rates for DMR regions.

    Args:
        bedgraph: Path to the input bedGraph file containing methylation rates and DMR regions.
        output: Prefix for the output file.
        mode: Calculation mode. 'rate' calculates the average of methylation rates.
              'count' calculates the rate from methylated and unmethylated counts.
    """
    try:
        console.print(f"[bold blue]Processing bedGraph file: {bedgraph} (mode: {mode})[/bold blue]")
        
        # Define columns and numeric columns based on mode
        if mode == 'rate':
            column_names = ['chr', 'start', 'end', 'meth_rate', 'chr_dmr', 'start_dmr', 'end_dmr']
            use_cols = [0, 1, 2, 3, 4, 5, 6]
            numeric_cols = ['start', 'end', 'meth_rate', 'start_dmr', 'end_dmr']
        else:  # mode == 'count'
            column_names = ['chr', 'start', 'end', 'meth_count', 'unmeth_count', 'chr_dmr', 'start_dmr', 'end_dmr']
            use_cols = [0, 1, 2, 4, 5, 6, 7, 8]
            numeric_cols = ['start', 'end', 'meth_count', 'unmeth_count', 'start_dmr', 'end_dmr']

        # Check if file exists
        if not os.path.exists(bedgraph):
            console.print(f"[bold red]Error: Input file {bedgraph} does not exist.[/bold red]")
            sys.exit(1)
        
        # Read the file with proper column names
        with console.status("[bold green]Reading bedGraph file..."):
            try:
                df = pd.read_csv(bedgraph, sep='\t', header=None, names=column_names, usecols=use_cols)
            except Exception as e:
                console.print(f"[bold red]Error reading bedGraph file: {str(e)}[/bold red]")
                sys.exit(1)
        
        # Validate numeric columns
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
        
        if mode == 'rate':
            # Group by DMR regions and calculate mean methylation rate
            result_df = df.groupby(['chr_dmr', 'start_dmr', 'end_dmr']).agg({
                'meth_rate': 'mean'
            }).reset_index()
        else: # mode == 'count'
            # Group by DMR and sum counts
            summed_counts = df.groupby(['chr_dmr', 'start_dmr', 'end_dmr']).agg({
                'meth_count': 'sum',
                'unmeth_count': 'sum'
            }).reset_index()

            # Calculate methylation rate
            total_counts = summed_counts['meth_count'] + summed_counts['unmeth_count']
            # Avoid division by zero
            summed_counts['meth_rate'] = np.where(total_counts > 0, summed_counts['meth_count'] / total_counts, 0)
            
            # Select final columns
            result_df = summed_counts[['chr_dmr', 'start_dmr', 'end_dmr', 'meth_rate']]

        # Write results to output file
        result_df.to_csv(output, sep='\t', header=False, index=False)
            
    except Exception as e:
        console.print(f"[bold red]Unexpected error: {str(e)}[/bold red]")
        sys.exit(1)


if __name__ == "__main__":
    main()
