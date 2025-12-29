# ProkaRegia Snakemake - Quick Start

## Installation (First Time Only)

```bash
# Install Miniconda if you don't have it
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create and activate snakemake environment
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

## Running ProkaRegia

### Method 1: Direct command (simplest!)
```bash
cd /path/to/prokaregia_snakemake
conda activate snakemake

# ONT data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont

# PacBio data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=pacbio

# Custom output directory
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont output_dir=MyResults
```

### Method 2: Using wrapper script
```bash
# ONT data
./run_prokaregia.sh -i /path/to/reads.fastq -s ont -t 16

# PacBio data
./run_prokaregia.sh -i /path/to/reads.fastq -s pacbio -t 16
```

### Method 3: Using config file
```bash
# Edit config.yaml first, then:
snakemake --use-conda --cores 16
```

## What Happens

1. First run will take longer as conda environments are created (one-time setup)
2. Subsequent runs use cached environments and are much faster
3. If interrupted, just run the same command again - it resumes automatically
4. Final genomes will be in `ProkaRegia/polished_bins/`

## Check Results

```bash
# View final quality report
less ProkaRegia/polished_checkm/quality_report.tsv

# View quality graphs
xdg-open ProkaRegia/checkm2_graphs/Completeness.png
xdg-open ProkaRegia/checkm2_graphs/Contamination.png

# List final genome bins
ls -lh ProkaRegia/polished_bins/
```

## Troubleshooting

### Pipeline fails
```bash
# Check the error message, then resume with:
snakemake --use-conda --cores 16
```

### Out of disk space
```bash
# Remove intermediate files (keeps final results)
rm -rf ProkaRegia/flye ProkaRegia/*_refined ProkaRegia/racon_temp ProkaRegia/ntlink_temp
```

### Start completely fresh
```bash
rm -rf ProkaRegia/
snakemake --use-conda --cores 16
```

## Expected Runtime

- Small dataset (1-5 GB): 4-8 hours
- Medium dataset (5-20 GB): 8-24 hours
- Large dataset (20+ GB): 24-72 hours

Runtime depends on data size, complexity, and available compute resources.

## Key Advantages Over Singularity

✅ No huge container to build/download
✅ Easier to debug individual steps
✅ Automatic parallelization
✅ Resume from failures
✅ Lighter on disk space
✅ More flexible and customizable

## Need More Help?

- See [CHEATSHEET.md](CHEATSHEET.md) for one-line command reference
- See [USAGE.md](USAGE.md) for detailed usage patterns
- See [README.md](README.md) for full documentation
- Check Snakemake docs: https://snakemake.readthedocs.io/
