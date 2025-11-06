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
You can download a relatively small test file here: https://zenodo.org/records/17540419
### Output Files and Directories
maxbin, metabat, metacoag, rosella, semibin: single binning algorithm outputs, refined with Rosella.

DASTool_refined: DASTool is a refinement tool that extracts the best bins from across different algorithms; these refined algorithms are further refined with Rosella.

polished_bins: DASTool bins further polished with NTLink and two rounds of Racon.

checkm2_graphs: quality graphs displaying completeness and contamination as assessed by Checkm2 for the above directories. On small datasets, frequently some or all binning algorithms do not produce any output with high enough quality to appear in these graphs.

NOTE: Checkm2 often poorly assesses highly reduced symbionts; in my testing, it reports complete Wolbachia, Liberibacter, and Carsonella genomes as less than 40% complete. If you believe these are present in your sample, I recommend checking individual bins with BLAST or a 16S search to see if it is present, then perhaps comparing the bin to a known reference genome to estimate completeness. Automatic implementation of this functionality will be available in the future.

## Options
```bash
-i, --input         Input file. Must be a long-read file in FASTQ format.
-o, --output        Output directory (default: ProkaRegia).
-r, --retain        Retain all intermediate files.
-t, --threads       Number of threads (default: all).
-s, --seq_tech      Technology used to generate long read data (accepts: ont, pacbio).
-h, --help          Display this help page
``` 
