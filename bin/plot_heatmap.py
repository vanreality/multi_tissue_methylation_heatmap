import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import re
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.colors import LinearSegmentedColormap
import click
from rich.console import Console
from rich.panel import Panel
from rich import print as rprint
from tqdm import tqdm

# Initialize rich console
console = Console()

def preprocess_data(input_file, target_label, target_label_na_threshold, labels_to_drop):
    """Preprocess methylation data for heatmap visualization.
    
    Args:
        input_file (str): Path to the input TSV file with methylation data.
        target_label (str): Label to use for filtering NA threshold.
        target_label_na_threshold (float): Maximum acceptable NA percentage for the target label.
        labels_to_drop (list): List of label patterns to exclude.
        
    Returns:
        pandas.DataFrame: Processed dataframe ready for visualization.
        
    Raises:
        FileNotFoundError: If the input file doesn't exist.
        ValueError: If the processed data is empty after filtering.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    console.print(f"[bold blue]Reading data from {input_file}[/bold blue]")
    
    # Data preprocess
    df = pd.read_csv(input_file, sep='\t')
    
    # Row filter, only keep autosomes
    with console.status("[bold green]Filtering autosomes...[/bold green]"):
        autosome_list = [f'chr{i}' for i in range(1, 23)]
        df = df[df['chr'].isin(autosome_list)].reset_index(drop=True)
    
    console.print(f"[green]Filtered to {len(df)} autosomal regions[/green]")
    
    with console.status("[bold green]Creating region IDs...[/bold green]"):
        region_id = (
            df[['chr', 'start', 'end']]
            .astype(str)
            .agg('-'.join, axis=1)
        )
        
        data_only = df.drop(columns=['chr', 'start', 'end'])
        data_only.index = region_id
        
        tdf = data_only.T.reset_index().rename(columns={'index': 'sample'})
        
        tissue_number = tdf['sample'].str.split('_\d', n=1, expand=True)
        tdf['label'] = tissue_number[0]
        
        tdf = (
            tdf.sort_values('label')
               .reset_index(drop=True)
        )
        
        tdf = tdf.drop(columns=['sample'])
    
    # Remove columns based on NA threshold in target label
    console.print(f"[bold green]Filtering columns with >{target_label_na_threshold*100}% NA in {target_label}...[/bold green]")
    initial_cols = tdf.shape[1] - 1  # Subtract 'label' column
    tdf = tdf.loc[:, (tdf[tdf['label'] == target_label].isna().mean() <= target_label_na_threshold) | (tdf.columns == 'label')]
    filtered_cols = tdf.shape[1] - 1  # Subtract 'label' column
    console.print(f"[green]Removed {initial_cols - filtered_cols} columns with high NA percentage[/green]")
    
    # Remove labels based on pattern matching
    if labels_to_drop:
        initial_rows = len(tdf)
        for pattern in labels_to_drop:
            tdf = tdf[~tdf['label'].str.contains(pattern)]
        console.print(f"[green]Removed {initial_rows - len(tdf)} rows matching patterns: {', '.join(labels_to_drop)}[/green]")
    
    if len(tdf) == 0 or tdf.shape[1] <= 1:
        raise ValueError("No data left after filtering. Please check your filter criteria.")
    
    return tdf

def plot_heatmap(df, title, output_file=None):
    """Generate a heatmap visualization of methylation data.
    
    Args:
        df (pandas.DataFrame): Processed dataframe with methylation data.
        title (str): Title for the heatmap.
        output_file (str, optional): Path to save the heatmap image. If None, displays the plot.
        
    Returns:
        None
    """
    with console.status("[bold green]Generating heatmap...[/bold green]"):
        data = df.drop(columns='label')
        row_labels = df['label'].values
        
        # Perform hierarchical clustering on columns
        col_linkage = linkage(data.T.fillna(0), method='average')
        col_order = leaves_list(col_linkage)
        data = data.iloc[:, col_order]
        
        # Create mask for NA values
        mask = data.isna()
        
        # Identify group boundaries for horizontal lines
        group_boundaries = np.where(row_labels[:-1] != row_labels[1:])[0] + 1
        
        # Create plot
        plt.figure(figsize=(24, 16))
        cmap = LinearSegmentedColormap.from_list('blue_yellow', ['#0000FF', '#FFD700'])
        ax = sns.heatmap(data, cmap=cmap, cbar=True, mask=mask, vmin=0, vmax=1)
        ax.imshow(mask, aspect='auto', cmap=LinearSegmentedColormap.from_list('gray', ['gray', 'gray']), alpha=0.5)
        
        # Add horizontal lines between tissue groups
        for pos in group_boundaries:
            plt.axhline(pos, color='black', linewidth=2)
        
        plt.title(title, fontsize=24)
        plt.xlabel('')
        plt.ylabel('')
        plt.xticks([])
        plt.yticks([])
        
        # Add tissue labels
        group_names = pd.Series(row_labels).groupby(row_labels).apply(lambda x: x.index[0])
        for group, idx in group_names.items():
            plt.text(-5, idx + len(df[df['label'] == group]) / 2 - 0.5, group, va='center', ha='right', fontsize=18)
        
        plt.tight_layout()
        
        # Save or display the plot
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            console.print(f"[bold green]Heatmap saved to {output_file}[/bold green]")
        else:
            plt.show()

def filter_columns(df, threshold=0.8, target_label='placenta'):
    """Filter columns based on methylation difference between target label and other tissues.
    
    Args:
        df (pandas.DataFrame): Processed dataframe with methylation data.
        threshold (float): Minimum absolute difference in methylation between target and other tissues.
        target_label (str): The tissue label to compare against others.
        
    Returns:
        pandas.DataFrame: Filtered dataframe containing only columns with significant differences.
    """
    console.print(f"[bold blue]Filtering columns with difference > {threshold}...[/bold blue]")
    initial_cols = df.shape[1] - 1  # Subtract 'label' column
    
    group_means = df.groupby(df['label'] == target_label).mean(numeric_only=True)
    diffs = (group_means.loc[True] - group_means.loc[False]).abs()
    selected_cols = diffs[diffs > threshold].index.tolist()
    
    result = df[['label'] + selected_cols]
    console.print(f"[green]Selected {len(selected_cols)} columns with difference > {threshold}[/green]")
    
    return result

def col_to_bed(df, output):
    """Convert column names in 'chr-start-end' format to BED file.
    
    Args:
        df (pandas.DataFrame): Dataframe with column names in 'chr-start-end' format.
        output (str): Path to save the BED file.
        
    Returns:
        None
    """
    col_list = list(df.drop(columns=['label']).columns)
    
    # col_list format 'chr-start-end' to bed file
    bed = []
    for col in tqdm(col_list, desc="Converting to BED format"):
        chrom, start, end = col.split('-')
        bed.append([chrom, start, end])
    
    pd.DataFrame(bed).to_csv(output, sep='\t', header=None, index=False)
    console.print(f"[bold green]BED file saved to {output}[/bold green]")

@click.command()
@click.option('--input', required=True, help='Path to input methylation matrix TSV file')
@click.option('--target-label', default='placenta', help='Target tissue label for filtering')
@click.option('--na-threshold', default=0.1, type=float, help='Maximum NA percentage allowed in target label')
@click.option('--drop-labels', multiple=True, default=['fetal'], help='Patterns of labels to drop (can be used multiple times)')
@click.option('--title', help='Heatmap title (defaults to auto-generated title based on input)')
@click.option('--output-plot', help='Path to save heatmap image (displays plot if not provided)')
@click.option('--filter-threshold', type=float, help='Filter columns by methylation difference threshold')
@click.option('--output-bed', help='Path to save filtered regions as BED file')
def main(input, target_label, na_threshold, drop_labels, title, output_plot, filter_threshold, output_bed):
    """Generate methylation heatmaps from multi-tissue methylation data.
    
    This tool processes methylation data, filters it based on various criteria,
    and generates visualizations. It can also output filtered regions in BED format.
    """
    try:
        rprint(Panel.fit(
            "[bold blue]Multi-tissue Methylation Heatmap Generator[/bold blue]",
            border_style="blue"
        ))
        
        # Process the data
        tdf = preprocess_data(
            input_file=input,
            target_label=target_label,
            target_label_na_threshold=na_threshold,
            labels_to_drop=drop_labels
        )
        
        console.print(f"[bold green]Processed data: {len(tdf)} rows, {tdf.shape[1] - 1} columns[/bold green]")
        
        # Apply filtering if threshold is provided
        if filter_threshold is not None:
            filtered_df = filter_columns(tdf, threshold=filter_threshold, target_label=target_label)
            plot_df = filtered_df
            
            # Generate BED file if output path is provided
            if output_bed:
                col_to_bed(filtered_df, output_bed)
        else:
            plot_df = tdf
        
        # Generate plot title if not provided
        if not title:
            if filter_threshold is not None:
                title = f'Heatmap of Filtered Multi-tissue Methylation Rate, n={plot_df.shape[1]-1}, min diff = {filter_threshold}'
            else:
                title = f'Heatmap of Multi-tissue Methylation Rate, n={plot_df.shape[1]-1}'
        
        # Generate the heatmap
        plot_heatmap(plot_df, title=title, output_file=output_plot)
        
        rprint(Panel.fit(
            "[bold green]Processing complete![/bold green]",
            border_style="green"
        ))
        
    except Exception as e:
        console.print(f"[bold red]Error: {str(e)}[/bold red]")
        raise click.Abort()

if __name__ == "__main__":
    main()

