 # Prokaregia
 Prokaregia is an end-to-end long-read assembly and binning suite aimed at prokaryotic genomes in datasets with heavy eukaryotic contamination. It is designed to require minimal troubleshooting--you can generally plug in raw reads at one end, and you will get clean, polished bins out the other.

 ## Usage
 Prokaregia is provided as a singularity container. It can be downloaded and then run with the following commands:
 ```bash
singularity pull library://schmiggle/prokaregia/prokaregia
singularity run -v [WORKING DIRECTORY FULL PATH] -w [WORKING DIRECTORY] prokaregia.sif [OPTIONS]
```
Or run all at once as follows:
```bash
singularity run -v [WORKING DIRECTORY FULL PATH] -w [WORKING DIRECTORY] library://schmiggle/prokaregia/prokaregia [OPTIONS]
```
In  addition to 

## Options

-i, --input         Input file. Must be a long-read file in FASTQ format.
-o, --output        Output directory (default: ProkaRegia).
-r, --retain        Retain all intermediate files.
-t, --threads       Number of threads (default: all).
-s, --seq_tech      Technology used to generate long read data (accepts: ont, pacbio).
-h, --help          Display this help page

