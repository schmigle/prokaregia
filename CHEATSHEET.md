# ProkaRegia Snakemake - Command Cheat Sheet

## One-Line Commands

### Run complete pipeline
```bash
# ONT data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont

# PacBio data
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=pacbio

# Custom output directory
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont output_dir=MyResults
```

### Dry run (see what will execute)
```bash
snakemake --use-conda -n --config input_fastq=reads.fastq seq_tech=ont
```

### Run specific step only
```bash
# Just quality control
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/QC.fastq

# Just assembly
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/flye/assembly.fasta

# Just binning
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/DAStool_bins

# Just polishing
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont ProkaRegia/polished_bins
```

### Resume after failure
```bash
# Just run the same command again - it resumes automatically!
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont
```

### Force re-run a rule
```bash
snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont --forcerun rule_name
```

### Check results
```bash
# View quality report
less ProkaRegia/polished_checkm/quality_report.tsv

# List final bins
ls -lh ProkaRegia/polished_bins/

# View graphs
xdg-open ProkaRegia/checkm2_graphs/Completeness.png
```

### Clean up
```bash
# Remove all outputs
rm -rf ProkaRegia/

# Remove intermediate files only (keep final results)
rm -rf ProkaRegia/flye ProkaRegia/*_refined ProkaRegia/racon_temp ProkaRegia/ntlink_temp
```

## Parameter Reference

| Parameter | Required | Values | Default | Description |
|-----------|----------|--------|---------|-------------|
| `input_fastq` | ✅ Yes | file path | - | Path to input FASTQ file |
| `seq_tech` | ✅ Yes | `ont` or `pacbio` | - | Sequencing technology |
| `output_dir` | ❌ No | directory name | `ProkaRegia` | Output directory name |
| `threads` | ❌ No | integer | `--cores` value | Threads per job (uses `--cores` by default) |

## Quick Examples

### Minimal command
```bash
snakemake --use-conda --cores 16 --config input_fastq=data.fastq seq_tech=ont
```

### All parameters
```bash
snakemake --use-conda --cores 32 --config input_fastq=/path/to/data.fastq seq_tech=pacbio output_dir=MyRun threads=32
```

### With absolute path
```bash
snakemake --use-conda --cores 16 --config input_fastq=/home/user/data/reads.fastq seq_tech=ont output_dir=Run_2024
```

### Faster with mamba
```bash
snakemake --use-conda --conda-frontend mamba --cores 16 --config input_fastq=reads.fastq seq_tech=ont
```

## Key Output Files

| File/Directory | Description |
|----------------|-------------|
| `ProkaRegia/polished_bins/*.fasta` | Final polished genome bins (YOUR RESULTS!) |
| `ProkaRegia/polished_checkm/quality_report.tsv` | Quality metrics (completeness/contamination) |
| `ProkaRegia/checkm2_graphs/*.png` | Quality visualization plots |
| `ProkaRegia/DAStool_refined/bins/*.fna` | Unpolished bins (if you want them) |

## Workflow Steps (Target Files)

Run up to any step by specifying its output file:

```bash
ProkaRegia/QC.fastq                          # 1. Quality control
ProkaRegia/flye/assembly.fasta               # 2. Assembly
ProkaRegia/tiara_classified.txt              # 3. Classification
ProkaRegia/prokarya_contigs.fasta            # 4. Filtered contigs
ProkaRegia/aln.sorted.bam                    # 5. Read mapping
ProkaRegia/rosella_output/bins               # 6. Binning
ProkaRegia/DAStool_bins                      # 7. Bin integration
ProkaRegia/DAStool_refined/bins              # 8. Refinement
ProkaRegia/polished_bins                     # 9. Polishing
ProkaRegia/polished_checkm/quality_report.tsv # 10. Final QC
ProkaRegia/checkm2_graphs/Completeness.png   # 11. Graphs
```

## Tips

- **First run is slow**: Conda environments are created (one-time)
- **Subsequent runs are fast**: Environments are cached
- **Auto-resume**: Just rerun the command if interrupted
- **Parallel execution**: Snakemake runs independent steps in parallel
- **Check progress**: Use `--dryrun` or `-n` to preview
