#!/usr/bin/env python3
"""
Generate methylation matrix from multiple bedgraph files.

This script reads a CSV file that contains sample names and paths to bedgraph files,
combines the methylation data into a single matrix, and outputs it as a TSV file.
"""

import os
import sys
import pandas as pd
import click
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from tqdm import tqdm

console = Console()

@click.command()
@click.option('--bedgraph_list', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True))
@click.option('--output', type=click.Path(file_okay=True, dir_okay=False, writable=True))
@click.option('--chunk-size', default=100000, help='Number of rows to process at once for large files.')
@click.option('--verbose', is_flag=True, help='Display detailed processing information.')
def generate_matrix(bedgraph_list, output, chunk_size, verbose):
    """
    Generate a methylation matrix from multiple bedgraph files.
    
    BEDGRAPH_LIST: CSV file with columns 'sample' and 'bedgraph_file_path'
    
    OUTPUT: Path to save the resulting methylation matrix in TSV format
    """
    try:
        # Read input CSV file
        console.print(Panel(f"[bold green]Reading input CSV: {bedgraph_list}[/]"))
        samples_df = pd.read_csv(bedgraph_list)
        
        # Validate input columns
        required_columns = ['sample', 'bedgraph_file_path']
        missing_columns = [col for col in required_columns if col not in samples_df.columns]
        if missing_columns:
            console.print(f"[bold red]Error: Missing required columns in input CSV: {', '.join(missing_columns)}[/]")
            sys.exit(1)
        
        # Validate bedgraph file paths
        invalid_paths = []
        for idx, row in samples_df.iterrows():
            if not os.path.isfile(row['bedgraph_file_path']):
                invalid_paths.append(f"Sample {row['sample']}: {row['bedgraph_file_path']}")
        
        if invalid_paths:
            console.print("[bold red]Error: The following bedgraph files do not exist:[/]")
            for path in invalid_paths:
                console.print(f"  - {path}")
            sys.exit(1)
        
        # Process bedgraph files to build the matrix
        console.print(Panel("[bold blue]Processing bedgraph files...[/]"))
        
        # First pass: collect all unique genomic positions
        all_positions = set()
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn()
        ) as progress:
            task = progress.add_task("[cyan]Collecting genomic positions...", total=len(samples_df))
            
            for idx, row in samples_df.iterrows():
                sample = row['sample']
                bedgraph_path = row['bedgraph_file_path']
                
                if verbose:
                    console.print(f"Processing {sample}: {bedgraph_path}")
                
                # Read bedgraph in chunks to handle large files
                for chunk in pd.read_csv(bedgraph_path, sep='\t', header=None, 
                                        names=['chr', 'start', 'end', 'meth_rate'],
                                        chunksize=chunk_size):
                    # Create position tuples and add to set
                    positions = set(zip(chunk['chr'], chunk['start'], chunk['end']))
                    all_positions.update(positions)
                
                progress.update(task, advance=1)
        
        console.print(f"[green]Found {len(all_positions):,} unique genomic positions.[/]")
        
        # Convert positions set to DataFrame for the matrix index
        positions_list = list(all_positions)
        positions_df = pd.DataFrame(positions_list, columns=['chr', 'start', 'end'])
        
        # Sort by chromosome and position
        positions_df = positions_df.sort_values(['chr', 'start', 'end'])
        
        # Create the matrix with the positions as index
        matrix = positions_df.copy()
        
        # Second pass: fill in methylation rates for each sample
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn()
        ) as progress:
            task = progress.add_task("[cyan]Building methylation matrix...", total=len(samples_df))
            
            for idx, row in samples_df.iterrows():
                sample = row['sample']
                bedgraph_path = row['bedgraph_file_path']
                
                if verbose:
                    console.print(f"Adding methylation data for {sample}")
                
                # Read bedgraph file
                bedgraph_df = pd.read_csv(bedgraph_path, sep='\t', header=None, 
                                        names=['chr', 'start', 'end', 'meth_rate'])
                
                # Create a lookup dictionary for faster access
                meth_dict = {(chrom, start, end): rate for chrom, start, end, rate 
                            in zip(bedgraph_df['chr'], bedgraph_df['start'], 
                                bedgraph_df['end'], bedgraph_df['meth_rate'])}
                
                # Add methylation rates to the matrix
                matrix[sample] = matrix.apply(
                    lambda row: meth_dict.get((row['chr'], row['start'], row['end']), float('nan')), 
                    axis=1
                )
                
                progress.update(task, advance=1)
        
        # Save the matrix to TSV file
        console.print(Panel(f"[bold green]Saving methylation matrix to: {output}[/]"))
        matrix.to_csv(output, sep='\t', index=False, na_rep='NA')
        console.print(f"[bold green]✓[/] Successfully generated methylation matrix with dimensions: {matrix.shape}")
        
    except Exception as e:
        console.print(f"[bold red]Error: {str(e)}[/]")
        if verbose:
            console.print_exception()
        sys.exit(1)

if __name__ == "__main__":
    generate_matrix()
