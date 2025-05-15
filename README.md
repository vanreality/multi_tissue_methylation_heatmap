# MTDMR-Heatmap

A Nextflow pipeline for analyzing multi-tissue DNA methylation data and generating heatmaps of Differentially Methylated Regions (DMRs).

## Description

MTDMR-Heatmap processes DNA methylation data across multiple tissue samples to identify and visualize methylation patterns in DMRs. The pipeline automates the extraction of methylation sites, calculation of methylation levels at DMRs, generation of methylation matrices, and creation of publication-ready heatmaps.

## Features

- Extract methylation sites from bedgraph files using bedtools
- Calculate methylation levels for specific regions defined in BED files
- Generate consolidated methylation matrices across multiple samples
- Create publication-quality heatmaps for visualizing multi-tissue methylation patterns
- Containerized workflow using Singularity for reproducibility
- Optimized for high-performance computing environments with SLURM

## Installation

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (v24.10.3 or later)
- [Singularity](https://sylabs.io/singularity/) (v3.0 or later)
- SLURM job scheduler (for HPC execution)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/vanreality/multi_tissue_methylation_heatmap.git
   cd MTDMR-Heatmap
   ```

2. The pipeline uses a Singularity container with required tools. The container definition is located in the `images` directory.
   ```bash
   singularity pull --arch amd64 library://syfan/mtdmr/common_tools:latest
   ``` 


## Usage

The basic command to run the pipeline:

```bash
nextflow run plot.nf \
  --input_csv samples.csv \
  --input_bed regions.bed \
  --output_dir results \
  -profile singularity
```

### Parameters

- `--input_csv`: CSV file with sample information (required)
- `--input_bed`: BED file containing regions of interest (required)
- `--output_dir`: Directory for output files (required)

## Examples

### Sample Input CSV Format

Create a CSV file with sample information:

```csv
sample,bedgraph_file_path
liver,/path/to/liver_methylation.bedgraph
brain,/path/to/brain_methylation.bedgraph
heart,/path/to/heart_methylation.bedgraph
```

### Complete Example

```bash
# Prepare sample data
cat > samples.csv << EOF
sample,bedgraph_file_path
liver,/path/to/liver_methylation.bedgraph
brain,/path/to/brain_methylation.bedgraph
heart,/path/to/heart_methylation.bedgraph
EOF

# Run the pipeline
nextflow run plot.nf \
  --input_csv samples.csv \
  --input_bed regions.bed \
  --output_dir results \
  -profile singularity
```

## Input/Output

### Input Files

1. **Sample CSV file**:
   - Headers: `sample,bedgraph_file_path`
   - Each row represents a tissue sample with its methylation bedgraph file path

2. **BED file**:
   - Standard BED format containing regions of interest (e.g., DMRs)
   - At minimum, contains chromosome, start, and end positions

3. **Bedgraph files**:
   - Four-column format: chromosome, start, end, methylation value
   - Methylation values typically range from 0-1 or 0-100%

### Output Files

The pipeline generates the following outputs in the specified output directory:

1. **Extracted methylation sites**:
   - Methylation values for each position in the regions of interest

2. **DMR methylation calculations**:
   - Aggregated methylation values for each region in each sample

3. **Methylation matrix**:
   - Combined matrix of all samples with methylation values for each region

4. **Heatmap visualizations**:
   - PNG/PDF heatmaps showing methylation patterns across tissues and regions

## Requirements

- Nextflow (v24.10.3 or later)
- Singularity (v3.0 or later)
- SLURM job scheduler
- Python 3.6+ with libraries:
  - numpy
  - pandas
  - matplotlib
  - seaborn
- bedtools (included in the Singularity container)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

Developed by vanreality.

For questions or issues, please open an issue on the GitHub repository.
