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
STEP=""

# Optional step-specific inputs
ASSEMBLY=""
CONTIGS=""
BAM=""
COVERAGE=""
BINS_DIR=""
GFA=""
ASSEMBLY_INFO=""

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
  --step STEP            Run specific pipeline step (default: full pipeline)
  -h, --help             Show this help message

Step-specific input arguments (used with --step):
  --assembly FILE        Assembly FASTA file (for contamination step)
  --contigs FILE         Contigs FASTA file (for binning, refinement, dastool steps)
  --bam FILE             Sorted BAM alignment file (for binning step)
  --gfa FILE             Assembly graph GFA file (for binning step)
  --assembly-info FILE   Assembly info file from Flye (for binning step)
  --coverage FILE        Coverage TSV file (for refinement, dastool steps)
  --bins-dir DIR         Directory containing bins (for refinement, polish steps)

Pipeline steps (use with --step):
  assembly               QC and assembly only (filtlong + flye)
  contamination          Remove contamination only (tiara)
  binning                Run all binners only
  refinement             Refine all bins with CheckM2 only
  dastool                DAS Tool aggregation only
  polish                 Scaffold and polish bins only (ntLink + racon)
  final                  Full pipeline (all steps)

Note: Snakemake automatically skips completed steps, so you can resume or
rerun specific steps as needed.

Examples:
  # Full pipeline
  prokaregia -i reads.fastq -s ont -t 16
  prokaregia -i reads.fastq -s pacbio -o MyResults

  # Run specific steps with user-provided inputs
  prokaregia -i reads.fastq -s ont --step binning --contigs my_contigs.fasta --bam my_alignment.bam --gfa my_graph.gfa
  prokaregia -i reads.fastq -s pacbio --step contamination --assembly my_assembly.fasta
  prokaregia -i reads.fastq -s ont --step polish --bins-dir my_bins/

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
        --step)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --step requires a value"
                show_help
                exit 1
            fi
            STEP="$2"
            shift 2
            ;;
        --assembly)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --assembly requires a value"
                show_help
                exit 1
            fi
            ASSEMBLY="$2"
            shift 2
            ;;
        --contigs)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --contigs requires a value"
                show_help
                exit 1
            fi
            CONTIGS="$2"
            shift 2
            ;;
        --bam)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --bam requires a value"
                show_help
                exit 1
            fi
            BAM="$2"
            shift 2
            ;;
        --coverage)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --coverage requires a value"
                show_help
                exit 1
            fi
            COVERAGE="$2"
            shift 2
            ;;
        --bins-dir)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --bins-dir requires a value"
                show_help
                exit 1
            fi
            BINS_DIR="$2"
            shift 2
            ;;
        --gfa)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --gfa requires a value"
                show_help
                exit 1
            fi
            GFA="$2"
            shift 2
            ;;
        --assembly-info)
            if [[ -z "$2" ]] || [[ "$2" == -* ]]; then
                echo "Error: --assembly-info requires a value"
                show_help
                exit 1
            fi
            ASSEMBLY_INFO="$2"
            shift 2
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

# Validate step-specific inputs and prepare input files
if [[ -n "$STEP" ]]; then
    case "$STEP" in
        contamination)
            if [[ -n "$ASSEMBLY" ]]; then
                if [[ ! -f "$ASSEMBLY" ]]; then
                    echo "Error: Assembly file not found: $ASSEMBLY"
                    exit 1
                fi
                # Copy to expected location
                mkdir -p "$OUTPUT_DIR/flye"
                cp "$ASSEMBLY" "$OUTPUT_DIR/flye/assembly.fasta"
                echo "Using provided assembly: $ASSEMBLY"
            fi
            ;;
        binning)
            if [[ -n "$CONTIGS" ]]; then
                if [[ ! -f "$CONTIGS" ]]; then
                    echo "Error: Contigs file not found: $CONTIGS"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR"
                cp "$CONTIGS" "$OUTPUT_DIR/prokarya_contigs.fasta"
                echo "Using provided contigs: $CONTIGS"
            fi
            if [[ -n "$BAM" ]]; then
                if [[ ! -f "$BAM" ]]; then
                    echo "Error: BAM file not found: $BAM"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR"
                cp "$BAM" "$OUTPUT_DIR/aln.sorted.bam"
                # Also check for BAI file
                if [[ -f "${BAM}.bai" ]]; then
                    cp "${BAM}.bai" "$OUTPUT_DIR/aln.sorted.bam.bai"
                else
                    echo "Warning: BAM index file ${BAM}.bai not found. Will be generated."
                fi
                echo "Using provided BAM: $BAM"
            fi
            if [[ -n "$GFA" ]]; then
                if [[ ! -f "$GFA" ]]; then
                    echo "Error: GFA file not found: $GFA"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR"
                cp "$GFA" "$OUTPUT_DIR/prokarya_graph.gfa"
                echo "Using provided GFA: $GFA"
            fi
            if [[ -n "$ASSEMBLY_INFO" ]]; then
                if [[ ! -f "$ASSEMBLY_INFO" ]]; then
                    echo "Error: Assembly info file not found: $ASSEMBLY_INFO"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR"
                cp "$ASSEMBLY_INFO" "$OUTPUT_DIR/prokarya_info.txt"
                echo "Using provided assembly info: $ASSEMBLY_INFO"
            fi
            ;;
        refinement)
            if [[ -n "$CONTIGS" ]]; then
                if [[ ! -f "$CONTIGS" ]]; then
                    echo "Error: Contigs file not found: $CONTIGS"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR"
                cp "$CONTIGS" "$OUTPUT_DIR/prokarya_contigs.fasta"
                echo "Using provided contigs: $CONTIGS"
            fi
            if [[ -n "$COVERAGE" ]]; then
                if [[ ! -f "$COVERAGE" ]]; then
                    echo "Error: Coverage file not found: $COVERAGE"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR/rosella_output"
                cp "$COVERAGE" "$OUTPUT_DIR/rosella_output/coverage.tsv"
                echo "Using provided coverage: $COVERAGE"
            fi
            if [[ -n "$BINS_DIR" ]]; then
                if [[ ! -d "$BINS_DIR" ]]; then
                    echo "Error: Bins directory not found: $BINS_DIR"
                    exit 1
                fi
                # Copy bins to expected binner directories
                for binner in rosella metabat maxbin semibin metacoag; do
                    mkdir -p "$OUTPUT_DIR/${binner}/bins"
                    cp "$BINS_DIR"/* "$OUTPUT_DIR/${binner}/bins/" 2>/dev/null || true
                done
                echo "Using provided bins from: $BINS_DIR"
            fi
            ;;
        dastool)
            if [[ -n "$CONTIGS" ]]; then
                if [[ ! -f "$CONTIGS" ]]; then
                    echo "Error: Contigs file not found: $CONTIGS"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR"
                cp "$CONTIGS" "$OUTPUT_DIR/prokarya_contigs.fasta"
                echo "Using provided contigs: $CONTIGS"
            fi
            if [[ -n "$COVERAGE" ]]; then
                if [[ ! -f "$COVERAGE" ]]; then
                    echo "Error: Coverage file not found: $COVERAGE"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR/rosella_output"
                cp "$COVERAGE" "$OUTPUT_DIR/rosella_output/coverage.tsv"
                echo "Using provided coverage: $COVERAGE"
            fi
            ;;
        polish)
            if [[ -n "$BINS_DIR" ]]; then
                if [[ ! -d "$BINS_DIR" ]]; then
                    echo "Error: Bins directory not found: $BINS_DIR"
                    exit 1
                fi
                mkdir -p "$OUTPUT_DIR/DAStool_refined/bins"
                cp "$BINS_DIR"/* "$OUTPUT_DIR/DAStool_refined/bins/" 2>/dev/null || true
                echo "Using provided bins from: $BINS_DIR"
            fi
            ;;
    esac
fi

# Map step names to Snakemake rule targets
TARGET_RULES=""
if [[ -n "$STEP" ]]; then
    case "$STEP" in
        assembly)
            # QC and assembly
            TARGET_RULES="copy_initial_assembly"
            ;;
        contamination)
            # Tiara classification and extraction
            TARGET_RULES="extract_prokarya"
            ;;
        binning)
            # Run all binners
            TARGET_RULES="rosella_binning metabat_binning maxbin_binning semibin_binning metacoag_binning"
            ;;
        refinement)
            # Refine all bins with CheckM2
            TARGET_RULES="refine_rosella_bins refine_metabat_bins refine_maxbin_bins refine_semibin_bins refine_metacoag_bins"
            ;;
        dastool)
            # DAS Tool aggregation and final refinement
            TARGET_RULES="final_refinement"
            ;;
        polish)
            # Scaffold and polish all bins
            TARGET_RULES="polish_all_bins"
            ;;
        final)
            # Full pipeline (use 'all' rule)
            TARGET_RULES=""
            ;;
        *)
            echo "Error: Unknown step '$STEP'"
            echo "Valid steps: assembly, contamination, binning, refinement, dastool, polish, final"
            exit 1
            ;;
    esac
fi

# Build snakemake command
SNAKEMAKE_CMD="snakemake --snakefile $SNAKEFILE --use-conda --cores $THREADS"
SNAKEMAKE_CMD="$SNAKEMAKE_CMD --config input_fastq=$INPUT_FASTQ seq_tech=$SEQ_TECH output_dir=$OUTPUT_DIR"

# Add target rules if specified, otherwise run 'all'
if [[ -n "$TARGET_RULES" ]]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD $TARGET_RULES"
fi

# Display configuration
echo "=========================================="
echo "ProkaRegia Pipeline Configuration"
echo "=========================================="
echo "Input FASTQ:     $INPUT_FASTQ"
echo "Sequencing tech: $SEQ_TECH"
echo "Output dir:      $OUTPUT_DIR"
echo "Threads:         $THREADS"
if [[ -n "$STEP" ]]; then
    echo "Pipeline step:   $STEP"
    if [[ -n "$ASSEMBLY" ]]; then
        echo "  - Assembly:    $ASSEMBLY"
    fi
    if [[ -n "$CONTIGS" ]]; then
        echo "  - Contigs:     $CONTIGS"
    fi
    if [[ -n "$BAM" ]]; then
        echo "  - BAM:         $BAM"
    fi
    if [[ -n "$GFA" ]]; then
        echo "  - GFA:         $GFA"
    fi
    if [[ -n "$ASSEMBLY_INFO" ]]; then
        echo "  - Assembly info: $ASSEMBLY_INFO"
    fi
    if [[ -n "$COVERAGE" ]]; then
        echo "  - Coverage:    $COVERAGE"
    fi
    if [[ -n "$BINS_DIR" ]]; then
        echo "  - Bins dir:    $BINS_DIR"
    fi
fi
echo "=========================================="
echo ""
echo "Running: $SNAKEMAKE_CMD"
echo ""

# Run snakemake
exec $SNAKEMAKE_CMD
