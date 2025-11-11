# ProkaRegia Snakemake Pipeline

An automated long-read assembly and binning pipeline for prokaryotic genomes in contaminated datasets.

## Installation

Install with bioconda, then run:

```bash
mamba create -n prokaregia -c conda-forge -c bioconda prokaregia
mamba activate prokaregia

prokaregia -i [READS FILE] -s [SEQUENCING TECH] [OPTIONS]
```
## Usage

### Quick Start (Bioconda Installation)

If you installed via bioconda, use the `prokaregia` command:
```bash
# ONT data
prokaregia -i reads.fastq -s ont

# PacBio data
prokaregia -i reads.fastq -s pacbio

# Custom output directory
prokaregia -i reads.fastq -s ont -t -o Prokaregia_Out
```

#### Run specific steps only
```bash
# Run only up to assembly
snakemake --use-conda --cores 16 ProkaRegia/flye/assembly.fasta

# Run only binning steps
snakemake --use-conda --cores 16 ProkaRegia/rosella_output/bins

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

## License

This pipeline implementation is provided as-is. Please refer to individual tool licenses.
