


# RNA MAP

[![Docker Linux Build](https://github.com/YesselmanLab/rna_map/actions/workflows/docker_linux_build.yml/badge.svg)](https://github.com/YesselmanLab/rna_map/actions/workflows/docker_linux_build.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![linting: flake8](https://img.shields.io/badge/linting-flake8-greenyellow)](https://github.com/PyCQA/flake8)
[![PYPI package](https://badge.fury.io/py/rna-map.png)](http://badge.fury.io/py/rna-map)

## Table of Contents

- [Overview](#overview)
- [Software Requirements](#software-requirements)
- [Installation](#installation)
  - [Using Conda (Recommended)](#using-conda-recommended)
  - [Using Docker](#using-docker)
  - [From Source](#from-source)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
  - [Basic Usage](#basic-usage)
  - [Advanced Usage](#advanced-usage)
  - [Working with Large Datasets](#working-with-large-datasets)
- [Command Line Reference](#command-line-reference)
- [Troubleshooting](#troubleshooting)
- [Development](#development)
- [Citation](#citation)

## Overview

**RNA MAP** is an open-source tool for rapid analysis of RNA mutational profiling (MaP) experiments. This tool was inspired by the DREEM algorithm developed by the Rouskin Lab and provides a comprehensive platform for analyzing DMS-reactivity of RNA molecules.

### What is RNA MaP?

RNA mutational profiling (MaP) is a powerful technique that uses dimethyl sulfate (DMS) to probe RNA structure by introducing mutations during reverse transcription. The resulting sequencing data reveals the structural state of RNA molecules at single-nucleotide resolution.

### Key Features

- **Fast and accurate analysis** of DMS-MaPseq experiments
- **Support for both single-end and paired-end** sequencing data
- **Quality control integration** with FastQC and Trim Galore
- **Flexible alignment options** using Bowtie2
- **Interactive visualizations** with Plotly
- **Batch processing** capabilities for large datasets
- **Docker support** for reproducible environments

### Input Requirements

- **FASTA file**: Reference RNA sequence(s) of interest
- **FASTQ file(s)**: Raw sequencing data from DMS-MaPseq experiment
- **Optional**: Dot-bracket structure file for enhanced visualization

## Software Requirements

### Core Dependencies
- **Python**: 3.8 or greater
- **Bowtie2**: 2.2.9+ (for sequence alignment)
- **FastQC**: 0.11.9+ (for quality control)
- **Trim Galore**: 0.6.6+ (for adapter trimming)
- **Cutadapt**: 1.18+ (for quality trimming)

### Optional Dependencies
- **Conda/Mamba**: For environment management
- **Docker**: For containerized deployment

> **💡 Recommendation**: If you're trying the software for the first time, we highly recommend using the Docker image for a hassle-free experience.

## Installation

### Using Conda (Recommended)

The easiest way to install RNA MAP is using the provided `environment.yml` file, which automatically installs all dependencies:

```bash
# Clone the repository
git clone https://github.com/YesselmanLab/rna_map
cd rna_map

# Create environment from environment.yml
conda env create -f environment.yml

# Activate the environment
conda activate rna-map

# Install the package in development mode
pip install -e .
```

Alternatively, you can install from PyPI:

```bash
# Create a new conda environment
conda create -n rna-map python=3.11

# Activate the environment
conda activate rna-map

# Install bioinformatics tools
conda install -c bioconda bowtie2 fastqc cutadapt trim-galore

# Install RNA MAP
pip install rna-map
```

### Using Docker

Docker provides the most reliable installation method with all dependencies pre-configured:

```bash
# Clone the repository
git clone https://github.com/YesselmanLab/rna_map
cd rna_map

# Build the Docker image
# For Linux and Intel Mac:
docker build -t rna-map -f docker/Dockerfile .

# For Apple Silicon Mac (ARM64):
docker build -t rna-map --platform linux/amd64 -f docker/Dockerfile .

# Run RNA MAP with Docker
rna-map -fa <fasta_file> -fq1 <fastq_file> --docker
```

### From Source

For development or if you need the latest features:

```bash
# Clone the repository
git clone https://github.com/YesselmanLab/rna_map
cd rna_map

# Install in development mode
pip install -e .

# Install development dependencies
pip install -e ".[dev]"
```

## Quick Start

After installation, the `rna-map` command-line tool will be available in your environment.

### Basic Analysis

```bash
# Single-end sequencing
rna-map -fa reference.fasta -fq1 reads.fastq

# Paired-end sequencing
rna-map -fa reference.fasta -fq1 reads_R1.fastq -fq2 reads_R2.fastq

# With secondary structure information
rna-map -fa reference.fasta -fq1 reads.fastq --dot-bracket structures.csv
```

## Usage Examples

### Basic Usage

#### Single-End Analysis
```bash
rna-map -fa my_rna.fasta -fq1 my_reads.fastq
```

#### Paired-End Analysis
```bash
rna-map -fa my_rna.fasta -fq1 reads_R1.fastq -fq2 reads_R2.fastq
```

#### With Secondary Structure
```bash
rna-map -fa my_rna.fasta -fq1 my_reads.fastq --dot-bracket structures.csv
```

### Advanced Usage

#### Using Parameter Files
```bash
# Create a custom parameter file
rna-map -fa my_rna.fasta -fq1 my_reads.fastq -pf custom_params.yml

# Use preset parameters for barcoded libraries
rna-map -fa my_rna.fasta -fq1 my_reads.fastq -pp barcoded-library
```

#### Quality Control Options
```bash
# Skip quality control steps for faster processing
rna-map -fa my_rna.fasta -fq1 my_reads.fastq --skip-fastqc --skip-trim-galore

# Custom quality cutoff
rna-map -fa my_rna.fasta -fq1 my_reads.fastq --tg-q-cutoff 20
```

#### Bit Vector Filtering
```bash
# Apply strict filtering criteria
rna-map -fa my_rna.fasta -fq1 my_reads.fastq \
  --map-score-cutoff 20 \
  --qscore-cutoff 20 \
  --mutation-count-cutoff 2 \
  --percent-length-cutoff 0.8
```

### Working with Large Datasets

For large datasets with thousands of reference sequences:

```bash
# Generate summary only (no individual plots)
rna-map -fa large_dataset.fasta -fq1 reads.fastq --summary-output-only

# Skip bit vector generation for very large datasets
rna-map -fa large_dataset.fasta -fq1 reads.fastq --skip-bit-vector
```

### Docker Usage

```bash
# Run with Docker (Linux/Intel Mac)
rna-map -fa my_rna.fasta -fq1 my_reads.fastq --docker

# Run with Docker (Apple Silicon Mac)
rna-map -fa my_rna.fasta -fq1 my_reads.fastq --docker --docker-platform linux/amd64
```



## Command Line Reference

For a complete list of available options, run:

```bash
rna-map --help
```

### Main Arguments
- `-fa, --fasta PATH`: Reference sequences in FASTA format (required)
- `-fq1, --fastq1 PATH`: Single-end reads or first pair of paired-end reads (required)
- `-fq2, --fastq2 PATH`: Second pair of paired-end reads (optional)
- `--dot-bracket PATH`: CSV file with secondary structure information
- `-pf, --param-file PATH`: Custom parameter file (YAML format)
- `-pp, --param-preset TEXT`: Use preset parameters (e.g., 'barcoded-library')

### Quality Control Options
- `--skip-fastqc`: Skip FastQC quality control
- `--skip-trim-galore`: Skip Trim Galore adapter trimming
- `--tg-q-cutoff INTEGER`: Quality cutoff for Trim Galore (default: 20)

### Alignment Options
- `--bt2-alignment-args TEXT`: Custom Bowtie2 arguments (comma-separated)
- `--save-unaligned PATH`: Save unaligned reads to specified path

### Bit Vector Options
- `--skip-bit-vector`: Skip bit vector generation
- `--summary-output-only`: Generate summary only (no individual plots)
- `--map-score-cutoff INTEGER`: Minimum mapping score threshold
- `--qscore-cutoff INTEGER`: Quality score threshold for nucleotides
- `--mutation-count-cutoff INTEGER`: Maximum mutations per read
- `--percent-length-cutoff FLOAT`: Minimum read length percentage
- `--min-mut-distance INTEGER`: Minimum distance between mutations

### Docker Options
- `--docker`: Run in Docker container
- `--docker-image TEXT`: Specify Docker image
- `--docker-platform TEXT`: Specify platform (e.g., linux/amd64)

### Other Options
- `--overwrite`: Overwrite existing output directory
- `--debug`: Enable debug mode
- `--help`: Show help message

## Troubleshooting

### Common Issues

#### Installation Problems

**Problem**: `rna-map` command not found after installation
```bash
# Solution: Ensure the conda environment is activated
conda activate rna-map
which rna-map
```

**Problem**: Missing dependencies (bowtie2, fastqc, etc.)
```bash
# Solution: Install using conda
conda install -c bioconda bowtie2 fastqc cutadapt trim-galore
```

#### Runtime Errors

**Problem**: Low alignment rates
- Check that your FASTQ files are properly formatted
- Verify that the reference sequence matches your experimental setup
- Try adjusting Bowtie2 parameters: `--bt2-alignment-args "very-sensitive-local"`

**Problem**: Memory issues with large datasets
```bash
# Solution: Use summary-only mode for large datasets
rna-map -fa large_dataset.fasta -fq1 reads.fastq --summary-output-only
```

**Problem**: Docker permission errors
```bash
# Solution: Add user to docker group (Linux) or use sudo
sudo usermod -aG docker $USER
# Then log out and back in
```

#### Quality Control Issues

**Problem**: Poor quality scores
- Check your sequencing quality
- Adjust quality cutoff: `--tg-q-cutoff 15` (lower threshold)
- Skip quality control if data is known to be good: `--skip-fastqc --skip-trim-galore`

**Problem**: No mutations detected
- Verify DMS treatment was performed correctly
- Check that reference sequence is correct
- Lower mutation count cutoff: `--mutation-count-cutoff 1`

### Getting Help

1. **Check the logs**: RNA MAP provides detailed logging information
2. **Use debug mode**: Add `--debug` flag for verbose output
3. **Test with sample data**: Use the provided test files in `test/resources/`
4. **GitHub Issues**: Report bugs and ask questions on the [GitHub repository](https://github.com/YesselmanLab/rna_map/issues)

### Example Output

A successful run will show:
```
rna_map.CLI - INFO - Starting RNA MAP analysis...
rna_map.RUN - INFO - Found 1 valid reference sequences
rna_map.MAPPING - INFO - Bowtie2 2.4.5 detected!
rna_map.MAPPING - INFO - FastQC v0.11.9 detected!
rna_map.MAPPING - INFO - Trim Galore 0.6.6 detected!
rna_map.EXTERNAL_CMD - INFO - Running FastQC...
rna_map.EXTERNAL_CMD - INFO - FastQC completed successfully
rna_map.EXTERNAL_CMD - INFO - Running Trim Galore...
rna_map.EXTERNAL_CMD - INFO - Trim Galore completed successfully
rna_map.EXTERNAL_CMD - INFO - Running Bowtie2 alignment...
rna_map.EXTERNAL_CMD - INFO - Bowtie2 alignment completed successfully
rna_map.BIT_VECTOR - INFO - Starting bit vector generation...
rna_map.BIT_VECTOR - INFO - Analysis completed successfully!
```

## Development

### Setting up Development Environment

```bash
# Clone the repository
git clone https://github.com/YesselmanLab/rna_map
cd rna_map

# Create development environment
conda env create -f environment.yml
conda activate rna-map

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest

# Run linting
black rna_map/
flake8 rna_map/
```


### Code Style

This project uses:
- **Black** for code formatting
- **Flake8** for linting
- **MyPy** for type checking
- **Pytest** for testing

## Citation

If you use RNA MAP in your research, please cite:

```
Yesselman, J.D., et al. (2022). RNA mutational profiling (MaP) with DREEM reveals 
structural and functional insights into RNA structure. Nucleic Acids Research, 
50(12), e70. https://doi.org/10.1093/nar/gkac435
```
