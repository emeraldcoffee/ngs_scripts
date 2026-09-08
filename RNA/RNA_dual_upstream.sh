#!/bin/bash
set -e

# arguments
usage() {
    echo "Help: This script takes in -s <sample_1_name>, -S <sample_2_name>,"
    echo "-i <STAR_index_directory>, -a <annotation.gtf>, and -h for help."
}

while getopts "s:S:i:a:h" opt; do
    case $opt in
        s) 
            sample_1="$OPTARG"
            echo "Sample 1: $sample_1"
            sample_1_1=${sample_1}_1
            sample_1_2=${sample_1}_2
            ;;
        S)
            sample_2="$OPTARG"
            echo "Sample 2: $sample_2"
            sample_2_1=${sample_2}_1
            sample_2_2=${sample_2}_2
            ;;
        i)
            index="$OPTARG"
            echo "STAR index directory path: $index."
            ;;
        a)
            annotation_path="$OPTARG"
            echo "Annotation .gtf file path: $annotation_path"
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

# RNA-seq analysis
echo "Beginning RNA sequencing upstream analysis."
echo ""

# FastQC quality control

if compgen -G "${sample_2}/*_R_fastqc.html" > /dev/null; then
    echo -e "Skipping FastQC HTML generation.\n"
else
    echo "Starting FastQC HTML generation."

    fastqc \
        -t 2 ${sample_1_1}_F.fq.gz ${sample_1_1}_R.fq.gz &
    fastqc \
        -t 2 ${sample_1_2}_F.fq.gz ${sample_1_2}_R.fq.gz &
    fastqc \
        -t 2 ${sample_2_1}_F.fq.gz ${sample_2_1}_R.fq.gz &
    fastqc \
        -t 2 ${sample_2_2}_F.fq.gz ${sample_2_2}_R.fq.gz &

    wait

    echo ""

    echo "FastQC HTMLs now available. Please take a moment to review."

    read -p "Press 'Enter' to continue. "
    echo ""
fi

echo "Now trimming using fastp."

# fastp trimming
fastp \
    -i ${sample_1_1}_F.fq.gz \
    -I ${sample_1_1}_R.fq.gz \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    -q 36 \
    -u 80 \
    --length_required 30 \
    -h ${sample_1_1}.fastp.html \
    -j ${sample_1_1}.fastp.json \
    -o ${sample_1_1}_trimmed_F.fq.gz \
    -O ${sample_1_1}_trimmed_R.fq.gz &
fastp \
    -i ${sample_1_2}_F.fq.gz \
    -I ${sample_1_2}_R.fq.gz \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    -q 36 \
    -u 80 \
    --length_required 30 \
    -h ${sample_1_2}.fastp.html \
    -j ${sample_1_2}.fastp.json \
    -o ${sample_1_2}_trimmed_F.fq.gz \
    -O ${sample_1_2}_trimmed_R.fq.gz &
fastp \
    -i ${sample_2_1}_F.fq.gz \
    -I ${sample_2_1}_R.fq.gz \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    -q 36 \
    -u 80 \
    --length_required 30 \
    -h ${sample_2_1}.fastp.html \
    -j ${sample_2_1}.fastp.json \
    -o ${sample_2_1}_trimmed_F.fq.gz \
    -O ${sample_2_1}_trimmed_R.fq.gz &
fastp \
    -i ${sample_2_2}_F.fq.gz \
    -I ${sample_2_2}_R.fq.gz \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    -q 36 \
    -u 80 \
    --length_required 30 \
    -h ${sample_2_2}.fastp.html \
    -j ${sample_2_2}.fastp.json \
    -o ${sample_2_2}_trimmed_F.fq.gz \
    -O ${sample_2_2}_trimmed_R.fq.gz &

wait

echo "Trimming now complete."
echo ""

read -p "Would you like to re-run FastQC HTML generation? <yes/no>: "
if [[ $REPLY = "yes" ]]; then

    # FastQC
    fastqc \
        -t 2 ${sample_1_1}_trimmed_F.fq.gz ${sample_1_1}_trimmed_R.fq.gz &
    fastqc \
        -t 2 ${sample_1_2}_trimmed_F.fq.gz ${sample_1_2}_trimmed_R.fq.gz &
    fastqc \
        -t 2 ${sample_2_1}_trimmed_F.fq.gz ${sample_2_1}_trimmed_R.fq.gz &
    fastqc \
        -t 2 ${sample_2_2}_trimmed_F.fq.gz ${sample_2_2}_trimmed_R.fq.gz &

    wait

    echo "Re-run FastQC HTMLs now available. Please take a moment to review."

    read -p "Press 'Enter' to continue. "
fi

echo ""
echo "Beginning alignment. "

# Alignment
STAR \
    --runThreadN 18 \
    --genomeDir ${index} \
    --readFilesIn ${sample_1_1}_trimmed_F.fq.gz ${sample_1_1}_trimmed_R.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 3 \
    --alignIntronMax 1000000 \
    --quantMode GeneCounts \
    --outBAMsortingThreadN 18 \
    --outFileNamePrefix ${sample_1_1}_ &
STAR \
    --runThreadN 18 \
    --genomeDir ${index} \
    --readFilesIn ${sample_1_2}_trimmed_F.fq.gz ${sample_1_2}_trimmed_R.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 3 \
    --alignIntronMax 1000000 \
    --quantMode GeneCounts \
    --outBAMsortingThreadN 18 \
    --outFileNamePrefix ${sample_1_2}_ &
STAR \
    --runThreadN 18 \
    --genomeDir ${index} \
    --readFilesIn ${sample_2_1}_trimmed_F.fq.gz ${sample_2_1}_trimmed_R.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 3 \
    --alignIntronMax 1000000 \
    --quantMode GeneCounts \
    --outBAMsortingThreadN 18 \
    --outFileNamePrefix ${sample_2_1}_ &
STAR \
    --runThreadN 18 \
    --genomeDir ${index} \
    --readFilesIn ${sample_2_2}_trimmed_F.fq.gz ${sample_2_2}_trimmed_R.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 3 \
    --alignIntronMax 1000000 \
    --quantMode GeneCounts \
    --outBAMsortingThreadN 18 \
    --outFileNamePrefix ${sample_2_2}_ &

wait

echo "Alignment now complete. Moving on to feature count."
echo ""

# featureCounts
featureCounts \
    -T 16 \
    -p \
    --countReadPairs \
    -a ${annotation_path} \
    -o counts.txt \
    ${sample_1_1}_Aligned.sortedByCoord.out.bam \
    ${sample_1_2}_Aligned.sortedByCoord.out.bam \
    ${sample_2_1}_Aligned.sortedByCoord.out.bam \
    ${sample_2_2}_Aligned.sortedByCoord.out.bam

echo "Feature count now complete. Output in 'counts.txt'."
