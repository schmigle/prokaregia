# ProkaRegia

An automated long-read assembly and binning pipeline for prokaryotic genomes in contaminated datasets.

## Installation

```bash
# Create and activate environment
conda create -n prokaregia -c conda-forge -c bioconda -c schmigle prokaregia
conda activate prokaregia

# Run the pipeline
prokaregia -i reads.fastq -s ont -t 16
```

### Quick Start

Use the `prokaregia` command:
```bash
# ONT data
prokaregia -i reads.fastq -s ont -t 16

# PacBio data
prokaregia -i reads.fastq -s pacbio -t 16

# Custom output directory
prokaregia -i reads.fastq -s ont -t 16 -o MyResults
```

## Citations
In addition to citing ProkaRegia, please cite the following tools:

- **CheckM2**: Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods*. 2023;20(8):1203-1212. doi:10.1038/s41592-023-01940-w

- **CoverM**: Woodcroft BJ. CoverM: Read coverage calculator for metagenomics. GitHub. https://github.com/wwood/CoverM

- **DAS Tool**: Sieber CMK, Probst AJ, Sharrar A, Thomas BC, Hess M, Tringe SG, Banfield JF. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology*. 2018;3(7):836-843. doi:10.1038/s41564-018-0171-1

- **Filtlong**: Wick RR. Filtlong: Quality filtering tool for long reads. GitHub. https://github.com/rrwick/Filtlong

- **Flye**: Kolmogorov M, Yuan J, Lin Y, Pevzner PA. Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*. 2019;37(5):540-546. doi:10.1038/s41587-019-0072-8

- **gfatools**: Li H. gfatools: Tools for manipulating sequence graphs in the GFA and rGFA formats. GitHub. https://github.com/lh3/gfatools

- **MaxBin 2.0**: Wu YW, Simmons BA, Singer SW. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics*. 2016;32(4):605-607. doi:10.1093/bioinformatics/btv638

- **MetaBAT 2**: Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ*. 2019;7:e7359. doi:10.7717/peerj.7359

- **MetaCoAG**: Mallawaarachchi V, Lin Y. MetaCoAG: Binning metagenomic contigs via composition, coverage and assembly graphs. *Genome Biology*. 2022;23:162. doi:10.1186/s13059-022-02730-z

- **ntLink**: Coombe L, Li JX, Lo T, Wong J, Nikolic V, Warren RL, Birol I. ntLink: A toolkit for de novo genome assembly scaffolding and mapping using long reads. *Current Protocols*. 2023;3(4):e733. doi:10.1002/cpz1.733

- **Racon**: Vaser R, Sović I, Nagarajan N, Šikić M. Fast and accurate de novo genome assembly from long uncorrected reads. *Genome Research*. 2017;27(5):737-746. doi:10.1101/gr.214270.116

- **Rosella**: Roach MJ, Pierce-Ward NT, Suchecki R, Mallawaarachchi V, Papudeshi B, Handley SA, Brown CT, Watson-Haigh NS, Edwards RA. Rosella: A high-quality, reference-guided, metagenomic binning tool. *bioRxiv*. 2023. doi:10.1101/2023.05.30.542897

- **SemiBin**: Pan S, Zhu C, Zhao XM, Coelho LP. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. *Nature Communications*. 2022;13:2326. doi:10.1038/s41467-022-29843-y

- **Tiara**: Karlicki M, Antonowicz S, Karnkowska A. Tiara: deep learning-based classification system for eukaryotic sequences. *Bioinformatics*. 2022;38(2):344-350. doi:10.1093/bioinformatics/btab672
