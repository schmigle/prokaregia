#!/usr/bin/env bash

# ProkaRegia wrapper script for conda installation
# Converts CLI flags to Snakemake config parameters

set -e

# Resolve symlinks to find the actual script location
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
    DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SCRIPT_DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"

SNAKEFILE="${SCRIPT_DIR}/Snakefile"

# Default values
INPUT_FASTQ=""
SEQ_TECH=""
OUTPUT_DIR="ProkaRegia"
THREADS=$(nproc)
USE_MAMBA=false

# Parse command line arguments
show_help() {
    cat << EOF
ProkaRegia - Automated long-read assembly and binning pipeline

Usage: prokaregia -i INPUT -s SEQ_TECH [OPTIONS]

Required arguments:
  -i, --input FILE       Input FASTQ file (ONT or PacBio reads)
  -s, --seq-tech TYPE    Sequencing technology: 'ont' or 'pacbio'

Optional arguments:
  -o, --output DIR       Output directory (default: ProkaRegia in current directory)
  -t, --threads NUM      Number of threads (default: all available cores)
  -m, --mamba            Use mamba instead of conda for faster dependency resolution
  -h, --help             Show this help message

Examples:
  prokaregia -i reads.fastq -s ont -t 16
  prokaregia -i reads.fastq -s pacbio -o MyResults -t 32 -m

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: -i/--input requires a value"
                show_help
                exit 1
            fi
            INPUT_FASTQ="$2"
            shift 2
            ;;
        -s|--seq-tech)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: -s/--seq-tech requires a value"
                show_help
                exit 1
            fi
            SEQ_TECH="$2"
            shift 2
            ;;
        -o|--output)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: -o/--output requires a value"
                show_help
                exit 1
            fi
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: -t/--threads requires a value"
                show_help
                exit 1
            fi
            THREADS="$2"
            shift 2
            ;;
        -m|--mamba)
            USE_MAMBA=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_FASTQ" ]]; then
    echo "Error: Input FASTQ file is required (-i/--input)"
    show_help
    exit 1
fi

if [[ -z "$SEQ_TECH" ]]; then
    echo "Error: Sequencing technology is required (-s/--seq-tech)"
    show_help
    exit 1
fi

if [[ ! -f "$INPUT_FASTQ" ]]; then
    echo "Error: Input file not found: $INPUT_FASTQ"
    exit 1
fi

if [[ "$SEQ_TECH" != "ont" ]] && [[ "$SEQ_TECH" != "pacbio" ]]; then
    echo "Error: seq-tech must be 'ont' or 'pacbio', got: $SEQ_TECH"
    exit 1
fi

# Build snakemake command
SNAKEMAKE_CMD="snakemake --snakefile $SNAKEFILE --use-conda --cores $THREADS"
SNAKEMAKE_CMD="$SNAKEMAKE_CMD --config input_fastq=$INPUT_FASTQ seq_tech=$SEQ_TECH output_dir=$OUTPUT_DIR"

if [[ "$USE_MAMBA" == true ]]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --conda-frontend mamba"
fi

# Display configuration
echo "=========================================="
echo "ProkaRegia Pipeline Configuration"
echo "=========================================="
echo "Input FASTQ:     $INPUT_FASTQ"
echo "Sequencing tech: $SEQ_TECH"
echo "Output dir:      $OUTPUT_DIR"
echo "Threads:         $THREADS"
echo "Use mamba:       $USE_MAMBA"
echo "=========================================="
echo ""
echo "Running: $SNAKEMAKE_CMD"
echo ""

# Run snakemake
exec $SNAKEMAKE_CMD
