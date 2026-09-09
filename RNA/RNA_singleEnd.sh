#!/usr/bin/env bash
set -e

# ----------------------
# This script is intended to be used for the processing of RNA sequencing data.
# 
# Each step of this pipeline will skip if the 1st sample's 1st replicate's associated output
# file exists. If the previous run of the script failed in that step, please delete at least
# the 1st sample's 1st replicate associated output file. 
# ----------------------

# ----------------------
# This script expects a file structure like the following:
# /usr
# |-- sample_1_1
# |   |-- sample_1_1.fq.gz
# |-- sample_1_2
# |   |-- sample_1_2.fq.gz
# |-- sample_2_1
# |   |-- sample_2_1.fq.gz
# |-- sample_2_2
# |   |-- sample_2_2.fq.gz
# ----------------------

# ----------------------
# To run this script, enter the following into a Linux terminal while in your usr directory:
# bash RNA_singleEnd.sh \
#   -s {number_of_samples} {sample_name_1} {sample_name_2} {etc} \
#   -S {job_name} \
#   -i {STAR_index_directory} \
#   -a {genomic_annotation_file}
#
# For example:
# bash RNA_singleEnd.sh \
#   -s 2 MYC TP53 \
#   -S oncogene \
#   -i /star_index/human_v3 \
#   -a /anno/gencode.v38.annotation.gtf
# ----------------------

usage() {
    echo "Help: This script takes in -s <number_of_samples> <sample_name_1> <sample_name_2> <etc>,"
    echo "-S <job_name>, -i <STAR_index_directory>, -a <genomic_annotation_file>, and -h for help."
}

echo ""

echo "SETTINGS:"

while getopts "s:S:i:a:h" opt; do
    case $opt in
        s) 
            num_samples="$OPTARG"
            echo "Number of samples: $num_samples"

            for ((i=1;i<=num_samples; i++)); do
                sample_name="${!OPTIND}"

                if [[ -z "$sample_name" || "$sample_name" == -* ]]; then
                    echo "Error: Expected $num_samples sample names after -s, but only found $((i-1))."
                    usage
                    exit 1
                fi

                declare -g "sample_${i}=${sample_name}"
                declare -g "sample_${i}_1=${sample_name}_1"
                declare -g "sample_${i}_2=${sample_name}_2"

                OPTIND=$((OPTIND + 1))
            done
        S)
            job_name="$OPTARG"
            echo "Job name: $job_name"
            ;;
        i)
            index="$OPTARG"
            echo "STAR index directory path: $index."
            ;;
        a)
            annotation_path="$OPTARG"
            echo "Genomic annotation file path: $annotation_path"
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Error: Invalid option -$OPTARG."
            usage
            exit 1
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument."
            usage
            exit 1
            ;;
    esac
done

echo ""

echo -e "Beginning single-end RNA sequencing upstream analysis.\n"

# runs FastQC HTML generation for users to check sequencing quality
if compgen -G "${sample_1_1}/${sample_1_1}_fastqc.html" > /dev/null; then
    echo -e "Skipping FastQC HTML generation.\n"
else
    echo -e "Starting FastQC HTML generation.\n"

    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        if compgen -G "${!sample1}"/"${!sample1}".fq.gz > /dev/null; then
            mv "${!sample1}"/"${!sample1}".fastq.gz "${!sample1}"/"${!sample1}".fq.gz
            mv "${!sample2}"/"${!sample2}".fastq.gz "${!sample2}"/"${!sample2}".fq.gz
        fi

        fastqc -t 2 \
            ${sample1}.fq.gz \
            ${sample2}.fq.gz &

    done

    wait

    echo ""

    echo "FastQC HTMLs now available. Please take a moment to review."

    read -p "Press 'Enter' to continue. "
    
    echo ""
fi

if compgen -G "${sample_1_1}/trimmed_${sample_1_1}.fq.gz" > /dev/null; then
    echo -e "Skipping fastp trimming.\n"
else
    echo -e "Now trimming with fastp.\n"

    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        fastp \
            -i ${sample1}.fq.gz \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --cut_mean_quality 20 \
            -q 36 \
            -u 80 \
            --length_required 30 \
            -h ${sample1}.fastp.html \
            -j ${sample1}.fastp.json \
            -o trimmed_${sample1}.fq.gz &

        fastp \
            -i ${sample2}.fq.gz \
            --detect_adapter_for_pe \
            --cut_front \
            --cut_tail \
            --cut_mean_quality 20 \
            -q 36 \
            -u 80 \
            --length_required 30 \
            -h ${sample2}.fastp.html \
            -j ${sample2}.fastp.json \
            -o trimmed_${sample2}.fq.gz &

    done
    
    wait

    echo -e "Trimming now complete.\n"  
fi

# if quality before trimming was subpar (or had overrepresented sequences or high adapter content),
# users may recheck FastQC HTMLs
read -p "Would you like to re-run FastQC HTML generation? <yes/no>: "
if [[ $REPLY = "yes" ]]; then
    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        fastqc -t 2 \
            "${!sample1}"/trimmed_"${!sample1}".fq.gz \
            "${!sample2}"/trimmed_"${!sample2}".fq.gz &
            
    done

    wait

    echo "Re-run FastQC HTMLs now available. Please take a moment to review."
    read -p "Press 'Enter' to continue. "
fi

echo ""

# aligning reads with reference genomes using STAR
if [[ ! -d "${index}" ]]; then
    echo "Warning: STAR index directory matching "${index}" not found. Check the path."
    exit 1
else
    echo -e "Index files matching "${index}" found.\n"
fi

if compgen -G "${sample_1_1}/${sample_1_1}_Aligned.sortedByCoord.out.bam" > /dev/null; then
    echo -e "Skipping STAR alignment.\n"
else
    echo -e "Beginning STAR alignment.\n"

    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        STAR \
            --runThreadN 18 \
            --genomeDir ${index} \
            --readFilesIn trimmed_${sample1}.fq.gz \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD \
            --outFilterMultimapNmax 1 \
            --outFilterMismatchNmax 3 \
            --alignIntronMax 1000000 \
            --quantMode GeneCounts \
            --outBAMsortingThreadN 18 \
            --outFileNamePrefix ${sample1}_ &
        STAR \
            --runThreadN 18 \
            --genomeDir ${index} \
            --readFilesIn trimmed_${sample2}.fq.gz \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD \
            --outFilterMultimapNmax 1 \
            --outFilterMismatchNmax 3 \
            --alignIntronMax 1000000 \
            --quantMode GeneCounts \
            --outBAMsortingThreadN 18 \
    
    done

    wait

    echo ""
fi

echo -e "Alignment now complete. Moving on to feature count.\n"

# counting number of sequenced reaeds mapping to gene features
if [[ -f "${anno}" ]]; then

    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        computing_samples_files+="${!sample1}_Aligned.sortedByCoord.out.bam "
        computing_samples_files+="${!sample2}_Aligned.sortedByCoord.out.bam "
    done

    featureCounts \
        -T 16 \
        -p \
        --countReadPairs \
        -a ${annotation_path} \
        -o ${job_name}_counts.txt \
        ${computing_samples_files}
else
    echo "The annotation file path is not correct. Exiting now."
    exit 1
fi

echo ""

echo "Feature count now complete. Output in '${job_name}_counts.txt'."
echo -e "Upstream analysis now complete.\n"

exit 0