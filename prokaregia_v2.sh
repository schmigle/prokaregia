#!/usr/bin/env bash
#POSIX

input=''
outfolder='ProkaRegia'
seq_tech=''
threads=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print 2*$4}')

is_positive_integer() {
    re='^[0-9]+$'
    if ! [[ $1 =~ $re ]] ; then
       return 1
    fi
    if [ "$1" -le 0 ]; then
        return 1
    fi
    return 0
}

die() {
    printf '%s\n' "$1" >&2
    exit 1
}

usage() {
    echo "Usage: prokaregia.sh [-i infile] [-o output_file] [-t thread_count] [-s seq_tech]"
    echo "Options:"
    echo "  -i, --input         Input file"
    echo "  -o, --output        Output directory (default: ProkaRegia)"
    echo "  -p, --purge         Remove unnecessary intermediate files after run (optional)" 
    echo "  -t, --threads       Number of threads (default: all)"
    echo "  -s, --seq_tech      Technology used to generate long read data (accepts: ont, pacbio)"  
    echo "  -h, --help          Display this help page"
    exit 1
}

if [[ $# -eq 1 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            input="$(readlink -f "$2")"
            shift 2
            ;;
        -o|--output)
            output=$2
            shift 2
            ;;
        -p|--purge)
            purge_option=true
            shift
            ;;
        -s|--seq_tech)
            if [[ $2 == "ont" || $2 == "pacbio" ]]; then
                seq_tech=$2
                shift 2
            else
                echo "Error: Sequencing technology must be either 'ont' or 'pacbio'."
                exit 1
            fi
            ;;
        -t|--threads)
            if is_positive_integer "$2"; then
                threads=$2
                shift 2
            else
                echo "Error: Thread count must be a positive integer."
                exit 1
            fi
            ;;

        -h|-/?|--help)
            usage
            ;;
        *)
            break
            ;;
    esac
done

if [ -z "$input" ] || [ -z "$seq_tech" ]; then
    echo "Error: Missing required options."
    usage
    exit 1
fi

if [ ! -f "$input" ]; then
    echo "Error: Input file not found."
    exit 1
fi

if [ -d "$output" ]; then
    rm -r "$output"
fi

mkdir "$output"; cd "$output"

filtlong --min_length 1000 --keep_percent 90 "$input" > QC.fastq

tiara -i QC.fastq -o tiara_classified.txt -t "$threads" --tf pro --fq

if [ $seq_tech = 'ont' ]; then
    flye --nano-raw prokarya_QC.fastq --meta -o flye -t "$threads"
else
    flye --pacbio-raw prokarya_QC.fastq --meta -o flye -t "$threads"
fi

tiara -i flye/assembly.fasta -o tiara_classified.txt -t "$threads" --tf pro --fa
ntLink scaffold gap_fill target=prokarya_assembly.fasta reads=prokarya_QC.fastq
mv *ntLink.scaffolds.fa ntlink_out.fasta
 
if [ $seq_tech = 'ont' ]; then
    minimap2 -ax map-ont -t "$threads" --sam-hit-only ntlink_out.fasta prokarya_QC.fastq > aln.sam
else
    minimap2 -ax map-pb -t "$threads" --sam-hit-only ntlink_out.fasta prokarya_QC.fastq > aln.sam
fi
    
racon -t "$threads" prokarya_QC.fastq aln.sam ntlink_out.fasta > racon_intermediate.fasta

if [ $seq_tech = 'ont' ]; then
    minimap2 -ax map-ont -t "$threads" --sam-hit-only racon_intermediate.fasta prokarya_QC.fastq > aln.sam
else
    minimap2 -ax map-pb -t "$threads" --sam-hit-only racon_intermediate.fasta prokarya_QC.fastq > aln.sam
fi

racon -t "$threads" prokarya_QC.fastq aln.sam racon_intermediate.fasta > polished.fasta

if [ $seq_tech = 'ont' ]; then
    minimap2 -ax map-ont -t "$threads" racon_intermediate.fasta prokarya_QC.fastq > aln.sam
else
    minimap2 -ax map-pb -t "$threads" racon_intermediate.fasta prokarya_QC.fastq > aln.sam
fi

samtools view -Sb -@ "$threads" aln.sam > aln.bam && rm aln.sam
samtools sort aln.bam -@ "$threads" > aln.sorted.bam && rm aln.bam
 
runMetaBat.sh -t $threads polished.fasta aln.sorted.bam

rm -r QC.fastq aln.sorted.bam racon_intermediate.fasta 
