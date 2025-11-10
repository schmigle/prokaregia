# ProkaRegia Snakemake Pipeline

An automated long-read assembly and binning pipeline for prokaryotic genomes in contaminated datasets, implemented with Snakemake for easier dependency management and reproducibility.

## Overview

This Snakemake implementation replaces the Singularity/Apptainer-based workflow with a more flexible conda-based approach. The pipeline performs:

1. **Quality Control** - Filters long reads with Filtlong
2. **Assembly** - De novo assembly with Flye
3. **Classification** - Taxonomic filtering with Tiara to remove eukaryotic contamination
4. **Binning** - Multiple binning algorithms (Rosella, MetaCoAG, MetaBAT, MaxBin, SemiBin2)
5. **Refinement** - Cross-bin refinement with DAS Tool and Rosella
6. **Scaffolding** - Gap filling with ntLink
7. **Polishing** - Two-pass consensus polishing with Racon
8. **Quality Assessment** - Comprehensive QC with CheckM2 and visualization

## Installation

### Option 1: Install via Bioconda (Recommended)

```bash
# Create and activate environment
conda create -n prokaregia -c conda-forge -c bioconda prokaregia
conda activate prokaregia

# Run the pipeline
prokaregia -i reads.fastq -s ont -t 16
```

### Option 2: Manual Installation from Source

**Prerequisites:**
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/) (recommended)

**Setup:**

1. **Clone this repository**:
   ```bash
   git clone https://github.com/yourusername/prokaregia_snakemake.git
   cd prokaregia_snakemake
   ```

2. **Create the prokaregia conda environment**:
   ```bash
   conda env create -f environment.yml
   conda activate prokaregia
   ```

   All tool-specific dependencies will be automatically installed by Snakemake in isolated conda environments when you first run the pipeline.

### Database Setup

The CheckM2 database (~3.5GB) will be automatically downloaded to `~/databases/CheckM2_database/` the first time you run the pipeline. This is a one-time download that will be reused for all future runs.

## Usage

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

### Quick Start (Manual Installation)

If you cloned from source, use snakemake directly:
```bash
# ONT data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont

# PacBio data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=pacbio

# Custom output directory (use absolute path)
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont output_dir=/path/to/MyResults
```

**Alternative: Use config file**
1. Edit [config.yaml](config.yaml):
   ```yaml
   input_fastq: "path/to/your/reads.fastq"
   output_dir: "ProkaRegia"
   seq_tech: "ont"  # or "pacbio"
   threads: 16
   ```

2. Run the pipeline:
   ```bash
   snakemake --use-conda --cores 16
   ```

**See [CHEATSHEET.md](CHEATSHEET.md) for quick reference of all commands!**

### Advanced Usage

#### Using Mamba for faster dependency resolution
```bash
snakemake --use-conda --conda-frontend mamba --cores 16
```

#### Run specific steps only
```bash
# Run only up to assembly
snakemake --use-conda --cores 16 ProkaRegia/flye/assembly.fasta

# Run only binning steps
snakemake --use-conda --cores 16 ProkaRegia/rosella_output/bins
```

#### Generate workflow visualization
```bash
snakemake --dag | dot -Tpng > workflow.png
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

#### Run on a cluster (SLURM example)
```bash
snakemake --use-conda --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {threads}" \
    --cluster-config cluster.yaml --jobs 10
```

## Directory Structure

```
prokaregia_snakemake/
├── Snakefile                    # Main workflow definition
├── config.yaml                  # Configuration file
├── README.md                    # This file
├── envs/                        # Conda environment definitions
│   ├── filtlong.yml
│   ├── flye.yml
│   ├── tiara.yml
│   ├── minimap_racon.yml
│   ├── samtools.yml
│   ├── rosella.yml
│   ├── metacoag.yml
│   ├── metabat.yml
│   ├── maxbin.yml
│   ├── semibin.yml
│   ├── gfastats.yml
│   ├── das_tool.yml
│   ├── ntLink.yml
│   └── checkm2.yml
└── scripts/                     # Helper scripts
    ├── Fasta_to_Contig2Bin.sh  # Converts bins to TSV format
    └── prokaregia_graph.py      # Generates quality visualization plots
```

## Output Structure

By default, outputs are created in the current working directory:

```
ProkaRegia/  (or your specified output directory)
├── QC.fastq                           # Filtered reads
├── flye/                              # Flye assembly outputs
├── initial_assembly_flye.fasta        # Initial assembly before filtering
├── tiara_classified.txt               # Taxonomic classification results
├── prokarya_contigs.fasta            # Filtered prokaryotic contigs
├── aln.sorted.bam                    # Read alignment to contigs
├── rosella_output/                    # Rosella binning results
├── metacoag/                          # MetaCoAG binning results
├── metabat/bins/                      # MetaBAT bins
├── maxbin/bins/                       # MaxBin bins
├── semibin/bins/                      # SemiBin bins
├── *_refined/bins/                    # Individual refined bins
├── *_checkm/                          # CheckM2 reports for each binner
├── DAStool_bins/                      # DAS Tool integrated bins
├── DAStool_refined/bins/              # Final refined bins (pre-polish)
├── polished_bins/                     # Final polished genomes
├── polished_checkm/quality_report.tsv # Quality metrics for final bins
├── checkm2_combined/                  # Combined quality reports
└── checkm2_graphs/                    # Quality visualization plots
    ├── Completeness.png
    └── Contamination.png
```

## Key Features

### Compared to Singularity Implementation

**Advantages:**
- ✅ No need to build large container images
- ✅ Easier to modify and debug individual steps
- ✅ Better dependency management per tool
- ✅ Parallel execution of independent steps
- ✅ Automatic resume on failure
- ✅ Better resource management
- ✅ Works on any system with conda

**Considerations:**
- First run takes time to set up conda environments
- Environments are cached for subsequent runs
- May need more disk space for multiple conda environments

## Pipeline Steps

### 1. Quality Control
- **Tool**: Filtlong
- **Purpose**: Filters reads to minimum 1000bp, keeps top 90% by quality

### 2. Initial Assembly
- **Tool**: Flye
- **Purpose**: Metagenomic assembly optimized for ONT or PacBio

### 3. Taxonomic Classification
- **Tool**: Tiara
- **Purpose**: Identifies and removes eukaryotic contamination

### 4. Read Mapping
- **Tools**: minimap2, samtools
- **Purpose**: Maps reads to prokaryotic contigs for coverage calculation

### 5. Binning (Multiple Algorithms)
- **Rosella**: Long-read-specific binning with coverage
- **MetaCoAG**: Graph-based binning using assembly graphs
- **MetaBAT**: Coverage-based binning
- **MaxBin**: Expectation-maximization binning
- **SemiBin2**: Deep learning-based binning

### 6. Individual Refinement
- **Tool**: Rosella refine
- **Purpose**: Refines bins from each algorithm independently

### 7. Cross-Bin Integration
- **Tool**: DAS Tool
- **Purpose**: Integrates predictions from all binners, selecting best bins

### 8. Final Refinement
- **Tool**: Rosella refine
- **Purpose**: Final bin refinement on integrated results

### 9. Scaffolding
- **Tool**: ntLink
- **Purpose**: Scaffolds contigs and fills gaps

### 10. Polishing
- **Tool**: Racon (2 rounds)
- **Purpose**: Consensus polishing to improve accuracy

### 11. Quality Assessment
- **Tool**: CheckM2
- **Purpose**: Assesses completeness and contamination of all bins

### 12. Visualization
- **Tool**: Custom Python script
- **Purpose**: Generates comparative quality plots

## Troubleshooting

### Conda environment creation fails
```bash
# Try using mamba instead
conda install -c conda-forge mamba
snakemake --use-conda --conda-frontend mamba
```

### CheckM2 database not found
```bash
# CheckM2 database should auto-download to ~/databases/CheckM2_database/
# If it fails, manually download:
checkm2 database --download
```

### Out of memory errors
- Reduce the number of parallel jobs: `--cores 8`
- Increase memory allocation if using cluster submission
- Consider running binners individually

### Pipeline stops unexpectedly
```bash
# Snakemake will automatically resume from last checkpoint
snakemake --use-conda --cores 16
```

## Dependencies

All dependencies are managed via conda. The pipeline uses:
- Filtlong (read filtering)
- Flye (assembly)
- Tiara (classification)
- seqtk (sequence extraction)
- minimap2 (mapping)
- samtools (BAM processing)
- Rosella (binning & refinement)
- MetaCoAG (binning)
- MetaBAT (binning)
- MaxBin (binning)
- SemiBin2 (binning)
- gfastats (graph processing)
- DAS Tool (bin integration)
- ntLink (scaffolding)
- Racon (polishing)
- CheckM2 (quality assessment)
- Python 3 with matplotlib, pandas, etc. (visualization)

## Citation

If you use this pipeline, please cite the original tools:
- **Flye**: Kolmogorov et al., Nature Biotechnology 2019
- **Tiara**: Karlicki et al., Bioinformatics 2022
- **Rosella**: Roach et al., bioRxiv 2023
- **MetaCoAG**: Mallawaarachchi et al., Genome Biology 2022
- **MetaBAT**: Kang et al., PeerJ 2019
- **MaxBin**: Wu et al., Microbiome 2016
- **SemiBin**: Pan et al., Nature Communications 2022
- **DAS Tool**: Sieber et al., Nature Microbiology 2018
- **CheckM2**: Chklovski et al., Nature Methods 2023

## License

This pipeline implementation is provided as-is. Please refer to individual tool licenses.

## Support

For issues specific to this Snakemake implementation, please check:
- Snakemake documentation: https://snakemake.readthedocs.io/
- Original ProkaRegia repository: https://github.com/schmigle/prokaregia

## Differences from Original Script

This Snakemake version differs from the original bash script in a few ways:

1. **Dependency Management**: Uses individual conda environments per tool instead of a monolithic Singularity container
2. **Parallelization**: Snakemake automatically parallelizes independent steps
3. **Resumability**: Can resume from any failed step without rerunning completed steps
4. **Modularity**: Each step is a separate rule, making it easier to modify or debug
5. **Sequence Extraction**: Uses seqtk (via its own conda environment) for robust extraction of prokaryotic contigs

## Quick Reference Commands

```bash
# Run complete pipeline
snakemake --use-conda --cores 16

# Dry run (check what will be executed)
snakemake --use-conda -n

# Force re-run specific rule
snakemake --use-conda --cores 16 --forcerun tiara_classify

# Clean up all outputs
snakemake --delete-all-output

# Generate workflow report
snakemake --report report.html

# Unlock directory (if previous run was interrupted)
snakemake --unlock
```
