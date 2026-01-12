# ProkaRegia Snakemake Pipeline

An automated long-read assembly and binning pipeline for prokaryotic genomes in contaminated datasets, implemented with Snakemake for easier dependency management and reproducibility.

## Installation

```bash
# Create and activate environment
conda create -n prokaregia -c conda-forge -c bioconda -c schmigle prokaregia
conda activate prokaregia

# Run the pipeline
prokaregia -i reads.fastq -s ont -t 16
```

### Quick Start (Bioconda Installation)

If you installed via bioconda, use the `prokaregia` command:
```bash
# ONT data
prokaregia -i reads.fastq -s ont -t 16

# PacBio data
prokaregia -i reads.fastq -s pacbio -t 16

# Custom output directory
prokaregia -i reads.fastq -s ont -t 16 -o MyResults
```
