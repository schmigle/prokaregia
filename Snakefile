"""
ProkaRegia Snakemake Pipeline
Automated long-read assembly and binning pipeline for prokaryotic genomes

Usage:
    snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=ont
    snakemake --use-conda --cores 16 --config input_fastq=reads.fastq seq_tech=pacbio output_dir=MyOutput
"""

import os
import sys
import glob
from shutil import which

# Get parameters from command line config or config file
# Command line takes precedence over config file
if os.path.exists("config.yaml"):
    configfile: "config.yaml"

# Get parameters with defaults
INPUT_FASTQ = config.get("input_fastq", None)
SEQ_TECH = config.get("seq_tech", None)
THREADS = config.get("threads", workflow.cores if workflow.cores else 8)

# Output directory: default to current working directory
output_dir_config = config.get("output_dir", "ProkaRegia")

if os.path.isabs(output_dir_config):
    # Absolute path, use as-is
    OUTPUT_DIR = output_dir_config
elif "/" in output_dir_config or "\\" in output_dir_config:
    # Relative path with slashes - not supported
    sys.exit(f"ERROR: output_dir '{output_dir_config}' contains path separators but is not an absolute path. "
             f"Please use either a simple directory name (e.g., 'ProkaRegia') or an absolute path (e.g., '/full/path/to/output').")
else:
    # Just a name, create in current working directory
    OUTPUT_DIR = os.path.abspath(output_dir_config)

# Script directory: relative to Snakefile location
SCRIPT_DIR = os.path.join(workflow.basedir, "scripts")

# Database location: use CheckM2 default location (~/databases/)
# This is where CheckM2 expects the database by default
CHECKM2_DB_DIR = os.path.expanduser("~/databases")
CHECKM2_DB_FILE = os.path.join(CHECKM2_DB_DIR, "CheckM2_database", "uniref100.KO.1.dmnd")

# Validate required parameters
if not INPUT_FASTQ:
    sys.exit("ERROR: input_fastq is required. Use: --config input_fastq=reads.fastq")
if not SEQ_TECH:
    sys.exit("ERROR: seq_tech is required. Use: --config seq_tech=ont (or pacbio)")
if SEQ_TECH not in ["ont", "pacbio"]:
    sys.exit(f"ERROR: seq_tech must be 'ont' or 'pacbio', got '{SEQ_TECH}'")

# Define output files
rule all:
    input:admittedly
        expand("{outdir}/polished_checkm/quality_report.tsv", outdir=OUTPUT_DIR),
        expand("{outdir}/checkm2_graphs/Completeness.png", outdir=OUTPUT_DIR),
        expand("{outdir}/checkm2_graphs/Contamination.png", outdir=OUTPUT_DIR)

# Step 0: Download CheckM2 database (one-time setup)
rule download_checkm2_db:
    output:
        db_file = CHECKM2_DB_FILE
    conda:
        "envs/checkm2.yml"
    shell:
        """
        export PYTHONNOUSERSITE=1
        unset PYTHONPATH

        echo "Downloading CheckM2 database to default location: {CHECKM2_DB_FILE}"
        echo "This is a large download and may take some time..."

        # Download database if it doesn't exist (uses CheckM2 default location)
        if [ ! -f {output.db_file} ]; then
            checkm2 database --download

            # Verify download succeeded
            if [ ! -f {output.db_file} ] || [ ! -s {output.db_file} ]; then
                echo "ERROR: CheckM2 database download failed!"
                echo "The database file was not created at: {output.db_file}"
                rm -rf {CHECKM2_DB_DIR}/CheckM2_database 2>/dev/null || true
                exit 1
            fi
        fi

        echo "CheckM2 database ready at: {CHECKM2_DB_FILE}/"
        """

# Step 1: Quality control with filtlong
rule filtlong:
    input:
        fastq = INPUT_FASTQ
    output:
        qc = "{outdir}/QC.fastq"
    conda:
        "envs/filtlong.yml"
    threads: 1
    shell:
        "filtlong --min_length 1000 --keep_percent 90 {input.fastq} > {output.qc}"

# Step 2: Initial assembly with Flye
rule flye_assembly:
    input:
        qc = "{outdir}/QC.fastq"
    output:
        assembly = "{outdir}/flye/assembly.fasta",
        graph = "{outdir}/flye/assembly_graph.gfa",
        info = "{outdir}/flye/assembly_info.txt"
    params:
        mode = "--nano-raw" if SEQ_TECH == "ont" else "--pacbio-raw"
    conda:
        "envs/flye.yml"
    threads: THREADS
    shell:
        "flye {params.mode} {input.qc} --meta -o $(dirname {output.assembly}) -t {threads}"

# Step 3: Copy initial assembly
rule copy_initial_assembly:
    input:
        "{outdir}/flye/assembly.fasta"
    output:
        "{outdir}/initial_assembly_flye.fasta"
    shell:
        "cp {input} {output}"

# Step 4: Taxonomic classification with Tiara
rule tiara_classify:
    input:
        assembly = "{outdir}/flye/assembly.fasta"
    output:
        classification = "{outdir}/tiara_classified.txt"
    conda:
        "envs/tiara.yml"
    threads: THREADS
    shell:
        "tiara -i {input.assembly} -o {output.classification} -t {threads} --tf bac arc pro"

# Step 5: Concatenate prokaryotic contigs (Tiara outputs them with --tf flag)
rule extract_prokarya:
    input:
        classification = "{outdir}/tiara_classified.txt"
    output:
        prokarya = "{outdir}/prokarya_contigs.fasta"
    shell:
        """
        # Tiara --tf flag creates separate fasta files for each classification
        # Concatenate all prokaryotic assemblies (bacteria, archaea, prokarya)
        cat {wildcards.outdir}/*assembly.fasta > {output.prokarya}
        """

# Step 6: Map reads to prokaryotic contigs
rule minimap2_mapping:
    input:
        contigs = "{outdir}/prokarya_contigs.fasta",
        reads = "{outdir}/QC.fastq"
    output:
        sam = "{outdir}/aln.sam"
    params:
        mode = "-ax map-ont" if SEQ_TECH == "ont" else "-ax map-pb"
    conda:
        "envs/minimap_racon.yml"
    threads: THREADS
    shell:
        "minimap2 {params.mode} -t {threads} {input.contigs} {input.reads} -o {output.sam}"

# Step 7: SAM to BAM conversion
rule samtools_view:
    input:
        sam = "{outdir}/aln.sam"
    output:
        bam = "{outdir}/aln.bam"
    conda:
        "envs/samtools.yml"
    threads: THREADS
    shell:
        "samtools view -@ {threads} -hbS {input.sam} -o {output.bam}"

# Step 8: Sort BAM file
rule samtools_sort:
    input:
        bam = "{outdir}/aln.bam"
    output:
        sorted_bam = "{outdir}/aln.sorted.bam"
    conda:
        "envs/samtools.yml"
    threads: THREADS
    shell:
        "samtools sort -@ {threads} {input.bam} -o {output.sorted_bam}"

# Step 9: Index sorted BAM
rule samtools_index:
    input:
        sorted_bam = "{outdir}/aln.sorted.bam"
    output:
        index = "{outdir}/aln.sorted.bam.bai"
    conda:
        "envs/samtools.yml"
    threads: THREADS
    shell:
        "samtools index -@ {threads} {input.sorted_bam}"

# Step 10: Create edits file for GFA
rule create_gfa_edits:
    input:
        prokarya = "{outdir}/prokarya_contigs.fasta",
        assembly = "{outdir}/flye/assembly.fasta"
    output:
        edits = "{outdir}/edits.sak"
    shell:
        """
        # Create array of prokaryotic contig names
        prokarya_names=$(grep "^>" {input.prokarya} | sed 's/^>//')

        # Find contigs in assembly not in prokarya and mark for exclusion
        grep "^>" {input.assembly} | sed 's/^>//' | while read name; do
            if ! echo "$prokarya_names" | grep -q "^${{name}}$"; then
                echo -e "EXCLUDE\\t${{name}}"
            fi
        done > {output.edits}
        """

# Step 11: Filter GFA graph
rule filter_gfa:
    input:
        gfa = "{outdir}/flye/assembly_graph.gfa",
        edits = "{outdir}/edits.sak"
    output:
        filtered_gfa = "{outdir}/prokarya_graph.gfa"
    conda:
        "envs/gfastats.yml"
    shell:
        "gfastats {input.gfa} -k {input.edits} -o gfa > {output.filtered_gfa}"

# Step 12: Filter assembly info
rule filter_assembly_info:
    input:
        info = "{outdir}/flye/assembly_info.txt",
        prokarya = "{outdir}/prokarya_contigs.fasta"
    output:
        filtered_info = "{outdir}/prokarya_info.txt"
    shell:
        """
        # Extract contig names into array
        prokarya_names=$(grep "^>" {input.prokarya} | sed 's/^>//')

        # Keep header
        head -n 1 {input.info} > {output.filtered_info}

        # Filter lines - keep only those matching prokaryotic contigs
        tail -n +2 {input.info} | while IFS= read -r line; do
            name=$(echo "$line" | awk '{{print $1}}')
            if echo "$prokarya_names" | grep -q "^${{name}}$"; then
                echo "$line"
            fi
        done >> {output.filtered_info}
        """

# Step 13: Rosella binning
rule rosella_binning:
    input:
        contigs = "{outdir}/prokarya_contigs.fasta",
        reads = "{outdir}/QC.fastq",
        bam = "{outdir}/aln.sorted.bam",
        bai = "{outdir}/aln.sorted.bam.bai",
        filtered_gfa = "{outdir}/prokarya_graph.gfa",
        filtered_info = "{outdir}/prokarya_info.txt"
    output:
        coverage = "{outdir}/rosella_output/coverage.tsv",
        bins_dir = directory("{outdir}/rosella_output/bins")
    params:
        mapper = "minimap2-ont" if SEQ_TECH == "ont" else "minimap2-pb"
    conda:
        "envs/rosella.yml"
    threads: THREADS
    shell:
        """
        rosella recover --longread-mapper {params.mapper} \
            -r {input.contigs} -t {threads} -o {wildcards.outdir}/rosella_output \
            -b {input.bam}
        mkdir -p {output.bins_dir}
        mv {wildcards.outdir}/rosella_output/rosella*fna {output.bins_dir}/ 2>/dev/null || true
        """

# Step 14: Prepare coverage for MetaCoAG
rule prepare_coverm:
    input:
        coverage = "{outdir}/rosella_output/coverage.tsv"
    output:
        coverm = "{outdir}/coverm.tsv"
    shell:
        "sed '1d' {input.coverage} > {output.coverm}"

# Step 15: MetaCoAG binning
rule metacoag_binning:
    input:
        gfa = "{outdir}/prokarya_graph.gfa",
        contigs = "{outdir}/prokarya_contigs.fasta",
        info = "{outdir}/prokarya_info.txt",
        coverage = "{outdir}/coverm.tsv"
    output:
        bins_dir = directory("{outdir}/metacoag")
    conda:
        "envs/metacoag.yml"
    threads: THREADS
    shell:
        """
        # Remove stale intermediate files from previous runs
        rm -f {wildcards.outdir}/*.hmmout* {wildcards.outdir}/*.frag*

        mkdir -p {output.bins_dir}
        metacoag --assembler flye --graph {input.gfa} --contigs {input.contigs} \
            --paths {input.info} --abundance {input.coverage} \
            --output {output.bins_dir} --nthreads {threads}
        """

# Step 16: MetaBAT binning (also generates depth file for MaxBin)
rule metabat_binning:
    input:
        contigs = "{outdir}/prokarya_contigs.fasta",
        bam = "{outdir}/aln.sorted.bam",
        bai = "{outdir}/aln.sorted.bam.bai"
    output:
        bins_dir = directory("{outdir}/metabat/bins"),
        depth = "{outdir}/prokarya_contigs.fasta.depth.txt"
    conda:
        "envs/metabat.yml"
    threads: THREADS
    shell:
        """
        cd {wildcards.outdir}
        # Remove old depth file to force regeneration
        rm -f prokarya_contigs.fasta.depth.txt
        runMetaBat.sh -s 100000 -t {threads} prokarya_contigs.fasta aln.sorted.bam
        mkdir -p metabat/bins
        mv metabat*bin*.f* metabat/bins/ 2>/dev/null || true
        """

# Step 17: MaxBin binning (uses depth file from MetaBAT)
rule maxbin_binning:
    input:
        contigs = "{outdir}/prokarya_contigs.fasta",
        depth = "{outdir}/prokarya_contigs.fasta.depth.txt"
    output:
        bins_dir = directory("{outdir}/maxbin/bins")
    conda:
        "envs/maxbin.yml"
    threads: THREADS
    shell:
        """
        run_MaxBin.pl -abund {input.depth} \
            -contig {input.contigs} -out {wildcards.outdir}/maxbin

        mkdir -p {output.bins_dir}
        mv {wildcards.outdir}/maxbin*fasta {output.bins_dir}/ 2>/dev/null || true
        """

# Step 18: SemiBin binning
rule semibin_binning:
    input:
        contigs = "{outdir}/prokarya_contigs.fasta",
        bam = "{outdir}/aln.sorted.bam",
        bai = "{outdir}/aln.sorted.bam.bai"
    output:
        bins_dir = directory("{outdir}/semibin/bins")
    conda:
        "envs/semibin.yml"
    threads: THREADS
    shell:
        """
        SemiBin2 single_easy_bin --environment global -t {threads} \
            --sequencing-type long_read -i {input.contigs} \
            -b {input.bam} -o {wildcards.outdir}/semibin

        mkdir -p {output.bins_dir}
        if [ -d {wildcards.outdir}/semibin/output_bins ]; then
            mv {wildcards.outdir}/semibin/output_bins {output.bins_dir}
            gzip -d {output.bins_dir}/*fa.gz 2>/dev/null || true
        fi
        """

# Step 19a: Refine Rosella bins
rule refine_rosella_bins:
    input:
        bins_dir = "{outdir}/rosella_output/bins",
        contigs = "{outdir}/prokarya_contigs.fasta",
        coverage = "{outdir}/rosella_output/coverage.tsv"
    output:
        refined_dir = directory("{outdir}/rosella_refined/bins")
    conda:
        "envs/rosella.yml"
    threads: THREADS
    shell:
        """
        for f in {input.bins_dir}/*.fa {input.bins_dir}/*.fasta; do
            [ -f "$f" ] && mv "$f" "${{f%.fa*}}.fna"
        done 2>/dev/null || true
        if [ -n "$(ls -A {input.bins_dir}/*.fna 2>/dev/null)" ]; then
            rosella refine -r {input.contigs} -o {wildcards.outdir}/rosella_refined \
                -C {input.coverage} -t {threads} --genome-fasta-directory {input.bins_dir}
            mkdir -p {output.refined_dir}
            mv {wildcards.outdir}/rosella_refined/*_bins/*fna {output.refined_dir}/ 2>/dev/null || true
        else
            mkdir -p {output.refined_dir}
        fi
        """

# Step 19b: Refine MetaCoAG bins
rule refine_metacoag_bins:
    input:
        bins_dir = "{outdir}/metacoag",
        contigs = "{outdir}/prokarya_contigs.fasta",
        coverage = "{outdir}/rosella_output/coverage.tsv"
    output:
        refined_dir = directory("{outdir}/metacoag_refined/bins")
    conda:
        "envs/rosella.yml"
    threads: THREADS
    shell:
        """
        for f in {input.bins_dir}/*.fa {input.bins_dir}/*.fasta; do
            [ -f "$f" ] && mv "$f" "${{f%.fa*}}.fna"
        done 2>/dev/null || true
        if [ -n "$(ls -A {input.bins_dir}/*.fna 2>/dev/null)" ]; then
            rosella refine -r {input.contigs} -o {wildcards.outdir}/metacoag_refined \
                -C {input.coverage} -t {threads} --genome-fasta-directory {input.bins_dir}
            mkdir -p {output.refined_dir}
            mv {wildcards.outdir}/metacoag_refined/*_bins/*fna {output.refined_dir}/ 2>/dev/null || true
        else
            mkdir -p {output.refined_dir}
        fi
        """

# Step 19c: Refine MetaBAT bins
rule refine_metabat_bins:
    input:
        bins_dir = "{outdir}/metabat/bins",
        contigs = "{outdir}/prokarya_contigs.fasta",
        coverage = "{outdir}/rosella_output/coverage.tsv"
    output:
        refined_dir = directory("{outdir}/metabat_refined/bins")
    conda:
        "envs/rosella.yml"
    threads: THREADS
    shell:
        """
        for f in {input.bins_dir}/*.fa {input.bins_dir}/*.fasta; do
            [ -f "$f" ] && mv "$f" "${{f%.fa*}}.fna"
        done 2>/dev/null || true
        if [ -n "$(ls -A {input.bins_dir}/*.fna 2>/dev/null)" ]; then
            rosella refine -r {input.contigs} -o {wildcards.outdir}/metabat_refined \
                -C {input.coverage} -t {threads} --genome-fasta-directory {input.bins_dir}
            mkdir -p {output.refined_dir}
            mv {wildcards.outdir}/metabat_refined/*_bins/*fna {output.refined_dir}/ 2>/dev/null || true
        else
            mkdir -p {output.refined_dir}
        fi
        """

# Step 19d: Refine MaxBin bins
rule refine_maxbin_bins:
    input:
        bins_dir = "{outdir}/maxbin/bins",
        contigs = "{outdir}/prokarya_contigs.fasta",
        coverage = "{outdir}/rosella_output/coverage.tsv"
    output:
        refined_dir = directory("{outdir}/maxbin_refined/bins")
    conda:
        "envs/rosella.yml"
    threads: THREADS
    shell:
        """
        for f in {input.bins_dir}/*.fa {input.bins_dir}/*.fasta; do
            [ -f "$f" ] && mv "$f" "${{f%.fa*}}.fna"
        done 2>/dev/null || true
        if [ -n "$(ls -A {input.bins_dir}/*.fna 2>/dev/null)" ]; then
            rosella refine -r {input.contigs} -o {wildcards.outdir}/maxbin_refined \
                -C {input.coverage} -t {threads} --genome-fasta-directory {input.bins_dir}
            mkdir -p {output.refined_dir}
            mv {wildcards.outdir}/maxbin_refined/*_bins/*fna {output.refined_dir}/ 2>/dev/null || true
        else
            mkdir -p {output.refined_dir}
        fi
        """

# Step 19e: Refine SemiBin bins
rule refine_semibin_bins:
    input:
        bins_dir = "{outdir}/semibin/bins",
        contigs = "{outdir}/prokarya_contigs.fasta",
        coverage = "{outdir}/rosella_output/coverage.tsv"
    output:
        refined_dir = directory("{outdir}/semibin_refined/bins")
    conda:
        "envs/rosella.yml"
    threads: THREADS
    shell:
        """
        for f in {input.bins_dir}/*.fa {input.bins_dir}/*.fasta; do
            [ -f "$f" ] && mv "$f" "${{f%.fa*}}.fna"
        done 2>/dev/null || true
        if [ -n "$(ls -A {input.bins_dir}/*.fna 2>/dev/null)" ]; then
            rosella refine -r {input.contigs} -o {wildcards.outdir}/semibin_refined \
                -C {input.coverage} -t {threads} --genome-fasta-directory {input.bins_dir}
            mkdir -p {output.refined_dir}
            mv {wildcards.outdir}/semibin_refined/*_bins/*fna {output.refined_dir}/ 2>/dev/null || true
        else
            mkdir -p {output.refined_dir}
        fi
        """

# Step 20: CheckM2 quality check for individual binners
rule checkm_individual:
    input:
        bins_dir = "{outdir}/{binner}_refined/bins",
        db_file = str(CHECKM2_DB_FILE)
    output:
        checkm_dir = directory("{outdir}/{binner}_checkm")
    conda:
        "envs/checkm2.yml"
    threads: THREADS
    shell:
        """
        export PYTHONNOUSERSITE=1
        unset PYTHONPATH

        if [ -n "$(ls -A {input.bins_dir}/*.fna 2>/dev/null)" ]; then
            checkm2 predict --threads {threads} --input {input.bins_dir} \
                --output-directory {output.checkm_dir} -x fna
        else
            mkdir -p {output.checkm_dir}
        fi
        """

# Step 21: Create bins TSV files for DAS Tool
rule create_bins_tsv:
    input:
        bins_dir = "{outdir}/{binner}_refined/bins"
    output:
        tsv = "{outdir}/{binner}_refined_bins.tsv"
    params:
        script_dir = SCRIPT_DIR
    shell:
        """
        if [ -n "$(ls -A {input.bins_dir}/*.fna 2>/dev/null)" ]; then
            bash {params.script_dir}/Fasta_to_Contig2Bin.sh -e fna -i {input.bins_dir} > {output.tsv}
        else
            touch {output.tsv}
        fi
        """

# Step 22: Aggregate bins with DAS Tool
rule dastool_aggregate:
    input:
        rosella_tsv = "{outdir}/rosella_refined_bins.tsv",
        metacoag_tsv = "{outdir}/metacoag_refined_bins.tsv",
        metabat_tsv = "{outdir}/metabat_refined_bins.tsv",
        maxbin_tsv = "{outdir}/maxbin_refined_bins.tsv",
        semibin_tsv = "{outdir}/semibin_refined_bins.tsv",
        contigs = "{outdir}/prokarya_contigs.fasta"
    output:
        bins_dir = directory("{outdir}/DAStool_bins")
    params:
        tsv_list = lambda wildcards, input: ",".join([
            f for f in [input.rosella_tsv, input.metacoag_tsv, input.metabat_tsv,
                       input.maxbin_tsv, input.semibin_tsv]
            if os.path.getsize(f) > 0
        ])
    conda:
        "envs/das_tool.yml"
    threads: THREADS
    shell:
        """
        DAS_Tool -i {params.tsv_list} -c {input.contigs} \
            -o {wildcards.outdir}/d --write_bins --write_unbinned

        mkdir -p {output.bins_dir}
        mv {wildcards.outdir}/d_DASTool_bins/* {output.bins_dir}/ 2>/dev/null || true
        rm {output.bins_dir}/*.k32.w100* 2>/dev/null || true

        # Rename to .fna
        for f in {output.bins_dir}/*.fa {output.bins_dir}/*.fasta; do
            [ -f "$f" ] && mv "$f" "${{f%.fa*}}.fna"
        done 2>/dev/null || true
        """

# Step 23: Final refinement with Rosella
rule final_refinement:
    input:
        bins_dir = "{outdir}/DAStool_bins",
        contigs = "{outdir}/prokarya_contigs.fasta",
        coverage = "{outdir}/rosella_output/coverage.tsv"
    output:
        refined_dir = directory("{outdir}/DAStool_refined/bins")
    conda:
        "envs/rosella.yml"
    threads: THREADS
    shell:
        """
        rosella refine -r {input.contigs} -o {wildcards.outdir}/DAStool_refined \
            -C {input.coverage} -t {threads} --genome-fasta-directory {input.bins_dir}

        mkdir -p {output.refined_dir}
        mv {wildcards.outdir}/DAStool_refined/*_bins/* {output.refined_dir}/ 2>/dev/null || true
        """

# Step 24: Scaffold with ntLink
rule ntlink_scaffold:
    input:
        bin_file = "{outdir}/DAStool_refined/bins/{bin}.fna",
        reads = "{outdir}/QC.fastq"
    output:
        scaffolded = "{outdir}/ntlink_temp/{bin}_scaffolded.fa"
    conda:
        "envs/ntLink.yml"
    threads: THREADS
    shell:
        """
        cd {wildcards.outdir}/ntlink_temp
        ntLink scaffold gap_fill target=../../{input.bin_file} reads=../../{input.reads} t={threads}

        # Find the output file
        output_file=$(ls *ntLink.scaffolds.gap_fill.fa 2>/dev/null | head -1)
        if [ -n "$output_file" ] && [ -s "$output_file" ]; then
            cp "$output_file" ../../{output.scaffolded}
        else
            cp ../../{input.bin_file} ../../{output.scaffolded}
        fi
        """

# Step 25: Polish with Racon (round 1)
rule racon_polish_round1:
    input:
        assembly = "{outdir}/ntlink_temp/{bin}_scaffolded.fa",
        reads = "{outdir}/QC.fastq"
    output:
        polished = "{outdir}/racon_temp/{bin}_racon1.fasta"
    params:
        mode = "-ax map-ont" if SEQ_TECH == "ont" else "-ax map-pb"
    conda:
        "envs/minimap_racon.yml"
    threads: THREADS
    shell:
        """
        minimap2 {params.mode} -t {threads} {input.assembly} {input.reads} \
            -o {wildcards.outdir}/racon_temp/{wildcards.bin}_aln1.sam

        racon -t {threads} -u {input.reads} \
            {wildcards.outdir}/racon_temp/{wildcards.bin}_aln1.sam \
            {input.assembly} > {output.polished}

        rm {wildcards.outdir}/racon_temp/{wildcards.bin}_aln1.sam
        """

# Step 26: Polish with Racon (round 2)
rule racon_polish_round2:
    input:
        assembly = "{outdir}/racon_temp/{bin}_racon1.fasta",
        reads = "{outdir}/QC.fastq"
    output:
        polished = "{outdir}/polished_bins/{bin}_polished.fasta"
    params:
        mode = "-ax map-ont" if SEQ_TECH == "ont" else "-ax map-pb"
    conda:
        "envs/minimap_racon.yml"
    threads: THREADS
    shell:
        """
        minimap2 {params.mode} -t {threads} {input.assembly} {input.reads} \
            -o {wildcards.outdir}/racon_temp/{wildcards.bin}_aln2.sam

        racon -t {threads} -u {input.reads} \
            {wildcards.outdir}/racon_temp/{wildcards.bin}_aln2.sam \
            {input.assembly} > {output.polished}

        rm {wildcards.outdir}/racon_temp/{wildcards.bin}_aln2.sam
        """

# Step 27a: Scaffold all bins with ntLink
rule scaffold_all_bins:
    input:
        refined_dir = "{outdir}/DAStool_refined/bins",
        reads = "{outdir}/QC.fastq"
    output:
        done = "{outdir}/scaffolded_bins/.done"
    params:
        outdir = "{outdir}"
    conda:
        "envs/ntLink.yml"
    threads: THREADS
    shell:
        """
        # Get absolute paths at the start (before any cd commands)
        READS_ABS=$(realpath {input.reads})
        REFINED_DIR_ABS=$(realpath {input.refined_dir})
        OUTDIR_ABS=$(realpath {params.outdir})

        mkdir -p "$OUTDIR_ABS"/scaffolded_bins
        mkdir -p "$OUTDIR_ABS"/ntlink_temp

        SCAFFOLDED_DIR_ABS="$OUTDIR_ABS"/scaffolded_bins
        DONE_MARKER="$OUTDIR_ABS"/scaffolded_bins/.done

        # Process each bin with ntLink
        for bin_file in "$REFINED_DIR_ABS"/*.fna; do
            if [ ! -f "$bin_file" ]; then
                continue
            fi

            bin_name=$(basename "$bin_file" .fna)
            echo "Scaffolding $bin_name with ntLink..."

            # Get absolute path to bin file
            bin_abs=$(realpath "$bin_file")

            # Create temp directory for this bin
            mkdir -p "$OUTDIR_ABS"/ntlink_temp/${{bin_name}}
            cd "$OUTDIR_ABS"/ntlink_temp/${{bin_name}}

            # Run ntLink with absolute paths
            ntLink scaffold gap_fill target="$bin_abs" reads="$READS_ABS" t={threads}

            # Check if ntLink succeeded and copy output
            if [ -f *.ntLink.scaffolds.gap_fill.fa ]; then
                cp *.ntLink.scaffolds.gap_fill.fa "$SCAFFOLDED_DIR_ABS/${{bin_name}}_scaffolded.fna"
            else
                # If ntLink failed, copy original bin
                echo "ntLink failed for $bin_name, using original bin"
                cp "$bin_abs" "$SCAFFOLDED_DIR_ABS/${{bin_name}}_scaffolded.fna"
            fi
        done

        touch "$DONE_MARKER"
        """

# Step 27b: Polish all scaffolded bins with Racon
rule polish_all_bins:
    input:
        done = "{outdir}/scaffolded_bins/.done",
        reads = "{outdir}/QC.fastq"
    output:
        done = "{outdir}/polished_bins/.done"
    params:
        outdir = "{outdir}",
        mode = "-ax map-ont" if SEQ_TECH == "ont" else "-ax map-pb"
    conda:
        "envs/minimap_racon.yml"
    threads: THREADS
    shell:
        """
        mkdir -p {params.outdir}/polished_bins
        mkdir -p {params.outdir}/racon_temp

        # Process each scaffolded bin
        for scaffolded_file in {params.outdir}/scaffolded_bins/*.fna; do
            if [ ! -f "$scaffolded_file" ]; then
                continue
            fi

            bin_name=$(basename "$scaffolded_file" _scaffolded.fna)
            echo "Polishing $bin_name with Racon..."

            # Racon round 1
            minimap2 {params.mode} -t {threads} \
                "$scaffolded_file" {input.reads} \
                -o {params.outdir}/racon_temp/${{bin_name}}_aln1.sam

            racon -t {threads} -u {input.reads} \
                {params.outdir}/racon_temp/${{bin_name}}_aln1.sam \
                "$scaffolded_file" > {params.outdir}/racon_temp/${{bin_name}}_racon1.fasta

            rm {params.outdir}/racon_temp/${{bin_name}}_aln1.sam

            # Racon round 2
            minimap2 {params.mode} -t {threads} \
                {params.outdir}/racon_temp/${{bin_name}}_racon1.fasta {input.reads} \
                -o {params.outdir}/racon_temp/${{bin_name}}_aln2.sam

            racon -t {threads} -u {input.reads} \
                {params.outdir}/racon_temp/${{bin_name}}_aln2.sam \
                {params.outdir}/racon_temp/${{bin_name}}_racon1.fasta \
                > {params.outdir}/polished_bins/${{bin_name}}_polished.fasta

            rm {params.outdir}/racon_temp/${{bin_name}}_aln2.sam
        done

        touch {output.done}
        """

# Step 28: CheckM2 on unpolished DAS Tool bins
rule checkm_dastool:
    input:
        bins_dir = "{outdir}/DAStool_refined/bins",
        db_file = str(CHECKM2_DB_FILE)
    output:
        report = "{outdir}/dastool_checkm/quality_report.tsv"
    conda:
        "envs/checkm2.yml"
    threads: THREADS
    shell:
        """
        export PYTHONNOUSERSITE=1
        unset PYTHONPATH
        checkm2 predict --threads {threads} --input {input.bins_dir} \
            --output-directory {wildcards.outdir}/dastool_checkm -x fna
        """

# Step 29: CheckM2 on polished bins
rule checkm_polished:
    input:
        done = "{outdir}/polished_bins/.done",
        db_file = str(CHECKM2_DB_FILE)
    output:
        report = "{outdir}/polished_checkm/quality_report.tsv"
    params:
        bins_dir = "{outdir}/polished_bins"
    conda:
        "envs/checkm2.yml"
    threads: THREADS
    shell:
        """
        export PYTHONNOUSERSITE=1
        unset PYTHONPATH
        checkm2 predict --threads {threads} --input {params.bins_dir} \
            --output-directory {wildcards.outdir}/polished_checkm -x fasta
        """

# Step 30: Combine CheckM2 reports
rule combine_checkm_reports:
    input:
        rosella = "{outdir}/rosella_checkm",
        metacoag = "{outdir}/metacoag_checkm",
        metabat = "{outdir}/metabat_checkm",
        maxbin = "{outdir}/maxbin_checkm",
        semibin = "{outdir}/semibin_checkm",
        dastool = "{outdir}/dastool_checkm/quality_report.tsv",
        polished = "{outdir}/polished_checkm/quality_report.tsv"
    output:
        combined_dir = directory("{outdir}/checkm2_combined")
    shell:
        """
        mkdir -p {output.combined_dir}

        for checkm_dir in {input.rosella} {input.metacoag} {input.metabat} \
                          {input.maxbin} {input.semibin}; do
            if [ -f "$checkm_dir/quality_report.tsv" ] && \
               [ "$(wc -l < "$checkm_dir/quality_report.tsv")" -ge 2 ]; then
                binner_name=$(basename "$checkm_dir" | sed 's/_checkm$//')
                cp "$checkm_dir/quality_report.tsv" "{output.combined_dir}/${{binner_name}}.tsv"
            fi
        done

        # Copy DAS Tool and polished reports
        if [ -f {input.dastool} ] && [ "$(wc -l < {input.dastool})" -ge 2 ]; then
            cp {input.dastool} {output.combined_dir}/dastool.tsv
        fi

        if [ -f {input.polished} ] && [ "$(wc -l < {input.polished})" -ge 2 ]; then
            cp {input.polished} {output.combined_dir}/polished.tsv
        fi
        """

# Step 30: Generate quality graphs
rule generate_graphs:
    input:
        combined_dir = "{outdir}/checkm2_combined"
    output:
        completeness = "{outdir}/checkm2_graphs/Completeness.png",
        contamination = "{outdir}/checkm2_graphs/Contamination.png"
    params:
        script_dir = SCRIPT_DIR
    conda:
        "envs/plotting.yml"
    shell:
        """
        cd {wildcards.outdir}
        python {params.script_dir}/prokaregia_graph.py
        """
