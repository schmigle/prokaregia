# ProkaRegia Snakemake - Usage Guide

## Basic Setup

**Activate snakemake environment**:
```bash
conda activate snakemake
cd /path/to/prokaregia_snakemake
```

## Running the Pipeline

### Option 1: Run the entire pipeline (command line config)
```bash
# ONT data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont

# PacBio data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=pacbio

# With custom output directory
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont output_dir=MyOutput
```

**Or use config file**: Edit `config.yaml` with your settings, then run:
```bash
snakemake --use-conda --cores 16
```

This will run everything from start to finish, automatically managing dependencies.

### Option 2: Run one step at a time
```bash
# Add your config parameters, then specify the target file

# Run just the quality control step
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/QC.fastq

# Then run assembly
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/flye/assembly.fasta

# Then run classification
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/tiara_classified.txt

# And so on...
```

### Option 3: Run from a specific step forward
```bash
# Start from binning and run everything after it
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/rosella_output/bins

# Or start from polishing
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/polished_bins
```

Snakemake will automatically figure out what needs to be run and skip any steps that are already complete.

## Useful Commands

### Check what will run (dry run)
```bash
snakemake --use-conda -n --config input_fastq=reads.fastq seq_tech=ont
```

### Resume after interruption
If the pipeline stops or fails, just run the same command again:
```bash
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont
```
It will automatically resume from where it left off.

### Force re-run a specific step
```bash
# Re-run just the assembly step
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont --forcerun flye_assembly
```

### Clean up and start fresh
```bash
# Remove all output files
rm -rf ProkaRegia/

# Or use snakemake's built-in cleanup (add your config params)
snakemake --delete-all-output --config input_fastq=reads.fastq seq_tech=ont
```

## Key Output Files

The main results you'll want:
- **Final genomes**: `ProkaRegia/polished_bins/*.fasta`
- **Quality metrics**: `ProkaRegia/polished_checkm/quality_report.tsv`
- **Quality graphs**: `ProkaRegia/checkm2_graphs/*.png`

## Workflow Steps (in order)

1. filtlong → Quality control
2. flye_assembly → Initial assembly
3. tiara_classify → Remove eukaryotic contamination
4. minimap2_mapping → Map reads to contigs
5. rosella_binning → First binning method
6. metacoag_binning → Second binning method
7. metabat_binning → Third binning method
8. maxbin_binning → Fourth binning method
9. semibin_binning → Fifth binning method
10. refine_bins → Refine each binner's output
11. checkm_individual → Quality check individual binners
12. create_bins_tsv → Prepare for integration
13. dastool_aggregate → Integrate all bins
14. final_refinement → Final refinement
15. ntlink_scaffold → Scaffolding
16. racon_polish_round1 → First polish
17. racon_polish_round2 → Second polish
18. checkm_polished → Final quality check
19. generate_graphs → Create visualization

You can run up to any of these steps by specifying the output file you want!
