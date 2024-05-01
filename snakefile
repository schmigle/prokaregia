from os import mkdir, chdir, cpu_count
from os.path import exists, abspath
import argparse

def check_positive_int(value):
    try:
        int_value = int(value)
        if int_value <= 0:
            raise argparse.ArgumentTypeError("Threads must be a positive integer")
        return str(int_value)
    except ValueError:
        raise argparse.ArgumentTypeError("Threads must be a positive integer")


parser = argparse.ArgumentParser()

parser.add_argument("--seq_tech", choices=["ont", "pb"], type=str, help="Sequencing technology: 'ont' or 'pb'")
parser.add_argument("--input", type=str, help="Input file path")
parser.add_argument("--output", default="prokaregia", type=str, help="Output directory path (default: prokaregia)")
parser.add_argument("--threads", default=cpu_count(), type=check_positive_int, help="Number of threads (default: number of available cores)")

args = parser.parse_args()

# Set config values based on command-line arguments
config["seq_tech"] = args.seq_tech
config["input"] = abspath(args.input)
config["output"] = args.output
config["threads"] = args.threads

mkdir(args.output); chdir(args.output)

rule filtlong:
    input:
        "input"
    conda:
        "envs/filtlong.yml"
    shell:
        "filtlong --min_length 1000 --keep_percent 90 {input} > QC.fastq"   

rule tiara_fq:
    threads:
        "threads"
    shell:
        "singularity run tiara.sif tiara -i QC.fastq --tf pro -o tiara.txt -t {threads} --fq"

rule flye:
    params:
        "technology=lambda wildcards: "pacbio" if seq_tech == "pb" else "nano""
    threads:
        "threads"
    conda:
        "envs/flye.yml"
    shell:
        "flye --{params.technology}-raw prokarya_QC.fastq --meta -o flye -t {threads}"

rule tiara_fa:
    threads:
        "threads"
    shell:
        "singularity run tiara.sif tiara -i flye/assembly.fasta --tf pro -o tiara.txt -t {threads} --fa"

rule ntLink:
    threads:
        "threads"
    conda:
        "envs/ntLink.yml"
    shell:
        ntLink scaffold gap_fill target=prokarya_assembly.fasta reads=prokarya_QC.fastq
        mv *ntLink.scaffolds.fa ntlink_out.fasta && rm *ntLink*

rule minimap:
    params:
        technology=lambda wildcards: "pb" if seq_tech == "pb" else "nano"
    input:
        reads=lambda wildcards: ["ntlink_out.fasta"] if not exists("racon_intermediate.fasta")
                                else ["racon_intermediate.fasta"] if exists("racon_intermediate.fasta") and not exists("polished.fasta")
                                else ["polished.fasta"]
    threads:
        "threads"
    conda:
        "envs/minimap.yml"
    shell:
        "minimap2 -ax map-{params.technology} -t {threads} --sam-hit-only {input.reads} prokarya_QC.fastq > aln.sam"

rule racon:
    input:
        contigs=lambda wildcards: ["ntlink_out.fasta"] if not exists("racon_intermediate.fasta")
                                else ["racon_intermediate.fasta"]
    output:
        contigs=lambda wildcards: ["racon_intermediate.fasta"] if not exists("racon_intermediate.fasta")
                                else ["polished.fasta"]
    threads:
        "threads"
    conda:
        "envs/racon.yml"
    shell:
        "racon -t {threads} prokarya_QC.fastq aln.sam {input.contigs} > {output.contigs}"

rule LRBinner:
    threads:
        "threads"
    shell:
        "singularity run LRBinner.sif python /LRBinner contigs -c polished.fasta --reads-path prokarya_QC.fastq -t {threads} -o lrbinner -sep -k 3"