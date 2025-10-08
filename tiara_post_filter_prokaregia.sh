#!/usr/bin/env bash
#POSIX

source /opt/miniconda/etc/profile.d/conda.sh

input=''
output='ProkaRegia'
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

conda run -n filtlong filtlong --min_length 1000 --keep_percent 90 "$input" > QC.fastq

#initial assembly and mapping
if [ $seq_tech = "ont" ]; then
    conda run -n flye flye --nano-raw QC.fastq --meta -o flye -t "$threads"
    cp flye/assembly.fasta initial_assembly_flye.fasta
    conda run -n tiara tiara -i flye/assembly.fasta -o tiara_classified.txt -t "$threads" --tf bac arc pro
    cat *assembly.fasta > prokarya_contigs.fasta
    conda run -n minimap_racon minimap2 -ax map-ont -t "$threads" prokarya_contigs.fasta QC.fastq -o aln.sam
else
    conda run -n flye flye --pacbio-raw QC.fastq --meta -o flye -t "$threads"
    cp flye/assembly.fasta initial_assembly_flye.fasta
    conda run -n tiara tiara -i flye/assembly.fasta -o tiara_classified.txt -t "$threads" --tf bac arc pro
    cat *assembly.fasta > prokarya_contigs.fasta 
    conda run -n minimap_racon minimap2 -ax map-pb -t "$threads" prokarya_contigs.fasta QC.fastq -o aln.sam
fi

#samtools processing and coverage
conda run -n samtools samtools view -@ "$threads" -hbS aln.sam -o aln.bam
conda run -n samtools samtools sort -@ "$threads" aln.bam -o aln.sorted.bam && rm aln.bam
conda run -n samtools samtools index -@ "$threads" aln.sorted.bam

# conda run -n coverm coverm contig -b aln.sorted.bam -t "$threads" -o coverm.tsv
# sed -i '1d' coverm.tsv 

#gfa file edits for metacoag
contig_names=($(grep "^>" prokarya_contigs.fasta | sed 's/^>//'))
for name in $(grep "^>" flye/assembly.fasta | sed 's/^>//'); do
    if ! [[ " ${contig_names[@]} " =~ " ${name} " ]]; then
        echo -e "EXCLUDE\t$name" >> edits.sak
    fi
done

conda run -n gfastats gfastats flye/assembly_graph.gfa -k edits.sak -o gfa > prokarya_graph.gfa

#info file edits for metacoag
head -n 1 flye/assembly_info.txt > prokarya_info.txt
tail -n +2 flye/assembly_info.txt | while IFS= read -r line; do
    name=$(echo "$line" | awk '{print $1}')
    if [[ " ${contig_names[@]} " =~ " ${name} " ]]; then
        echo "$line" >> prokarya_info.txt
    fi
done


#binning
if [ $seq_tech = "ont" ]; then
    conda run -n rosella rosella recover --longread-mapper minimap2-ont -r prokarya_contigs.fasta -t "$threads" -o rosella_output --longreads QC.fastq && \
        mkdir rosella_output/bins && mv rosella_output/rosella*fna rosella_output/bins
else
    conda run -n rosella rosella recover --longread-mapper minimap2-ont -r prokarya_contigs.fasta -t "$threads" -o rosella_output --longreads QC.fastq && \
        mkdir rosella_output/bins && mv rosella_output/rosella*fna rosella_output/bins
fi

sed '1d' rosella_output/coverage.tsv > coverm.tsv

mkdir metacoag; conda run -n metacoag metacoag --assembler flye --graph prokarya_graph.gfa --contigs prokarya_contigs.fasta --paths prokarya_info.txt --abundance coverm.tsv --output metacoag --nthreads "$threads"
conda run -n metabat runMetaBat.sh -s 100000 -t "$threads" prokarya_contigs.fasta aln.sorted.bam && mkdir metabat && mkdir metabat/bins && mv metabat*bin*.f* metabat/bins
conda run -n maxbin run_MaxBin.pl -abund prokarya_contigs.fasta.depth.txt -contig prokarya_contigs.fasta -out maxbin && mkdir maxbin && mkdir maxbin/bins && mv maxbin*fasta maxbin/bins
conda run -n semibin SemiBin2 single_easy_bin --environment global -t "$threads" --sequencing-type long_read -i prokarya_contigs.fasta -b aln.sorted.bam -o semibin && mv semibin/output_bins semibin/bins && gzip -d semibin/bins/*fa.gz
    
#bin tsv preparation 
for i in */bins; do
    if [ "$(ls -A $i)" ]; then
        for j in $i/*.fa $i/*.fasta; do mv $j "${j%.fa*}.fna"; done
        conda run -n rosella rosella refine -r prokarya_contigs.fasta -o "${i%/bins}_refined" -C rosella_output/coverage.tsv -t "$threads" --genome-fasta-directory $i && \
            cd "${i%/bins}_refined" && mkdir bins && mv *_bins/*fna bins && cd ..
    else
        continue
    fi
done

for i in *_refined/bins; do
    echo $i
    if [ ! -z "$(ls -A $i)" ]; then
        bash /Fasta_to_Contig2Bin.sh -e fna -i $i > "${i%/*}_bins.tsv"
        conda run -n checkm checkm2 predict --threads "$threads" --input $i --output-directory "${i%_refined*}_checkm" -x fna
        # directory="${i%/*}"; extension=$(ls $i/* | head -n1 | sed 's/.*\.//') 
        # bash /Fasta_to_Contig2Bin.sh -e $extension -i $i > "${directory}_bins.tsv"
    fi
done

bins=$(find . -maxdepth 1 -type f -name '*_bins.tsv' -printf '%P\n'); formatted_bins=$(printf "%s," $bins | sed 's/,$//')

#cross-bin refinement
conda run -n das_tool DAS_Tool -i "$formatted_bins" -c prokarya_contigs.fasta -o d --write_bins --write_unbinned
mv *DAS*bins DAStool_bins
rm DAStool_bins/*.k32.w100*

for i in DAStool_bins/*.fa DAStool_bins/*.fasta; do mv $i "${i%.fa*}.fna"; done

#final refinement
conda run -n rosella rosella refine -r prokarya_contigs.fasta -o DAStool_refined -C rosella_output/coverage.tsv -t "$threads" --genome-fasta-directory DAStool_bins
mkdir DAStool_refined/bins; mv DAStool_refined/*_bins/* DAStool_refined/bins

#final scaffolding and polishing
mkdir polished_bins

for i in DAStool_refined/bins/*.fna; do
    filename="${i##*/}"; filename="${filename%.*}"
    conda run -n ntLink ntLink scaffold gap_fill target="$i" reads=QC.fastq t="$threads"
    cp DAStool_refined/*ntLink.scaffolds.gap_fill.fa "ntlink_${filename}.fa" && infile="ntlink_${filename}.fa" 
    wc -c "ntlink_${filename}.fa" || infile=$i

    if [ $seq_tech = "ont" ]; then
        conda run -n minimap_racon minimap2 -ax map-ont -t "$threads" "$infile" QC.fastq -o aln.sam
    else
        conda run -n minimap_racon minimap2 -ax map-pb -t "$threads" "$infile" QC.fastq -o aln.sam
    fi

    conda run -n minimap_racon racon -t "$threads" -u QC.fastq aln.sam "$infile" > racon_intermediate.fasta
    
    if [ $seq_tech = "ont" ]; then
        conda run -n minimap_racon minimap2 -ax map-ont -t "$threads" racon_intermediate.fasta QC.fastq -o aln.sam
    else
        conda run -n minimap_racon minimap2 -ax map-pb -t "$threads" racon_intermediate.fasta QC.fastq -o aln.sam
    fi

    conda run -n minimap_racon racon -t "$threads" -u QC.fastq aln.sam racon_intermediate.fasta > "${filename}_polished.fasta" && mv "${filename}_polished.fasta" polished_bins

    rm DAStool_refined/*.k32.w100*
done

#checkm2 quality check and graphs: first without polish, then with
conda run -n checkm checkm2 predict --threads "$threads" --input DAStool_refined/bins --output-directory dastool_checkm -x fa
conda run -n checkm checkm2 predict --threads "$threads" --input polished_bins --output_directory polished_checkm -x fasta
mkdir checkm2_combined
for i in *_checkm; do
    if [ "$(wc -l < "$i/quality_report.tsv")" -ge 2 ]; then mv $i/quality_report.tsv "checkm2_combined/${i%_checkm}.tsv"; fi
done 
    # mv dastool_checkm/quality_report.tsv checkm2_combined/dastool.tsv; mv polished_checkm/quality_report.tsv checkm2_combined/polished.tsv
conda run -n checkm python /prokaregia_graph.py
