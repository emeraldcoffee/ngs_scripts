#!/bin/bash
set -e

# arguments
usage() {
    echo "Help: This script takes in -s <sample_1_name>, -S <sample_2_name>,"
    echo "-i <Bowtie2_index>, -n <normalization method: CPM, RPKM, RPGC, or spike-in>,"
    echo "-a <annotation path>, and -h for help."
}

while getopts "s:S:t:T:u:U:i:n:a:h" opt; do
    case $opt in
        s) 
            sample_1="$OPTARG"
            echo "Sample 1 name: $sample_1"
            sample_1_2=${sample_1}_2
            ;;
        S)
            sample_2="$OPTARG"
            echo "Sample 2 name: $sample_2"
            sample_2_1=${sample_2}_1
            sample_2_2=${sample_2}_2
            ;;
        t)
            sample_3="$OPTARG"
            echo "Sample 3 name: $sample_3"
            sample_3_2=${sample_3}_2
            ;;
        T)
            sample_4="$OPTARG"
            echo "Sample 4 name: $sample_4"
            sample_4_1=${sample_4}_1
            sample_4_2=${sample_4}_2
            ;;
        u)
            sample_5="$OPTARG"
            echo "Sample 5 name: $sample_5"
            sample_5_1=${sample_5}_1
            sample_5_2=${sample_5}_2
            ;;
        U)
            sample_6="$OPTARG"
            echo "Sample 6 name: $sample_6"
            sample_6_1=${sample_6}_1
            sample_6_2=${sample_6}_2
            ;;
        i)
            index="$OPTARG"
            echo "Bowtie2 index: $index. Ensure you have a soft-link to this index in this directory."
            ;;
        n)
            norm="$OPTARG"
            echo "Normalization method: $norm"
            ;;
        a)
            anno="$OPTARG"
            echo "Annotation path: $anno"
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

echo -e "Beginning dual end sequencing upstream analysis.\n"

# FastQC HTML generation to check sequencing quality
if compgen -G "${sample_6_2}/*_2_fastqc.html" > /dev/null; then
    echo -e "Skipping FastQC HTML generation.\n"
else
    echo "Starting FastQC HTML generation."

    fastqc -t 2 \
        ${sample_1_2}/*_1.fq.gz \
        ${sample_1_2}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_2_1}/*_1.fq.gz \
        ${sample_2_1}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_2_2}/*_1.fq.gz \
        ${sample_2_2}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_3_2}/*_1.fq.gz \
        ${sample_3_2}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_4_1}/*_1.fq.gz \
        ${sample_4_1}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_4_2}/*_1.fq.gz \
        ${sample_4_2}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_5_1}/*_1.fq.gz \
        ${sample_5_1}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_5_2}/*_1.fq.gz \
        ${sample_5_2}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_6_1}/*_1.fq.gz \
        ${sample_6_1}/*_2.fq.gz &
    fastqc -t 2 \
        ${sample_6_2}/*_1.fq.gz \
        ${sample_6_2}/*_2.fq.gz &

    wait

    echo ""

    echo "FastQC HTMLs now available. Please take a moment to review."

    read -p "Press 'Enter' to continue. "
    echo ""
fi

echo -e "Now trimming with cutadapt.\n"

# trimming adapters, poly-G tails
if [[ -f "${sample_2_2}/trimmed_${sample_2_2}_2.fastq.gz" ]]; then
    echo -e "Skipping cutadapt trimming.\n"
else
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_1_2}/trimmed_${sample_1_2}_1.fastq.gz \
        -p ${sample_1_2}/trimmed_${sample_1_2}_2.fastq.gz \
        ${sample_1_2}/*_1.fq.gz \
        ${sample_1_2}/*_2.fq.gz \
        -j 12 &
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_2_1}/trimmed_${sample_2_1}_1.fastq.gz \
        -p ${sample_2_1}/trimmed_${sample_2_1}_2.fastq.gz \
        ${sample_2_1}/*_1.fq.gz \
        ${sample_2_1}/*_2.fq.gz \
        -j 12 &
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_2_2}/trimmed_${sample_2_2}_1.fastq.gz \
        -p ${sample_2_2}/trimmed_${sample_2_2}_2.fastq.gz \
        ${sample_2_2}/*_1.fq.gz \
        ${sample_2_2}/*_2.fq.gz \
        -j 12 &
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_3_2}/trimmed_${sample_3_2}_1.fastq.gz \
        -p ${sample_3_2}/trimmed_${sample_3_2}_2.fastq.gz \
        ${sample_3_2}/*_1.fq.gz \
        ${sample_3_2}/*_2.fq.gz \
        -j 12 &
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_4_1}/trimmed_${sample_4_1}_1.fastq.gz \
        -p ${sample_4_1}/trimmed_${sample_4_1}_2.fastq.gz \
        ${sample_4_1}/*_1.fq.gz \
        ${sample_4_1}/*_2.fq.gz \
        -j 12 &
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_4_2}/trimmed_${sample_4_2}_1.fastq.gz \
        -p ${sample_4_2}/trimmed_${sample_4_2}_2.fastq.gz \
        ${sample_4_2}/*_1.fq.gz \
        ${sample_4_2}/*_2.fq.gz \
        -j 12 &
   cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_5_1}/trimmed_${sample_5_1}_1.fastq.gz \
        -p ${sample_5_1}/trimmed_${sample_5_1}_2.fastq.gz \
        ${sample_5_1}/*_1.fq.gz \
        ${sample_5_1}/*_2.fq.gz \
        -j 12 &
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_5_2}/trimmed_${sample_5_2}_1.fastq.gz \
        -p ${sample_5_2}/trimmed_${sample_5_2}_2.fastq.gz \
        ${sample_5_2}/*_1.fq.gz \
        ${sample_5_2}/*_2.fq.gz \
        -j 12 &
   cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_6_1}/trimmed_${sample_6_1}_1.fastq.gz \
        -p ${sample_6_1}/trimmed_${sample_6_1}_2.fastq.gz \
        ${sample_6_1}/*_1.fq.gz \
        ${sample_6_1}/*_2.fq.gz \
        -j 12 &
    cutadapt \
        --trim-n \
        -m 20 \
        -q 20 \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -a "G{20}" \
        -A "G{20}" \
        -o ${sample_6_2}/trimmed_${sample_6_2}_1.fastq.gz \
        -p ${sample_6_2}/trimmed_${sample_6_2}_2.fastq.gz \
        ${sample_6_2}/*_1.fq.gz \
        ${sample_6_2}/*_2.fq.gz \
        -j 12 &
    
    wait

    echo ""
fi

echo -e "Trimming now complete.\n"

# if quality before trimming was subpar (or had overrepresented sequences or high adapter content), recheck FastQC HTMLs
read -p "Would you like to re-run FastQC HTML generation? <yes/no>: "
if [[ ${REPLY} = "yes" ]]; then
    fastqc -t 2 \
        ${sample_1_2}/trimmed_${sample_1_2}_1.fastq.gz \
        ${sample_1_2}/trimmed_${sample_1_2}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_2_1}/trimmed_${sample_2_1}_1.fastq.gz \
        ${sample_2_1}/trimmed_${sample_2_1}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_2_2}/trimmed_${sample_2_2}_1.fastq.gz \
        ${sample_2_2}/trimmed_${sample_2_2}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_3_2}/trimmed_${sample_3_2}_1.fastq.gz \
        ${sample_3_2}/trimmed_${sample_3_2}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_4_1}/trimmed_${sample_4_1}_1.fastq.gz \
        ${sample_4_1}/trimmed_${sample_4_1}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_4_2}/trimmed_${sample_4_2}_1.fastq.gz \
        ${sample_4_2}/trimmed_${sample_4_2}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_5_1}/trimmed_${sample_5_1}_1.fastq.gz \
        ${sample_5_1}/trimmed_${sample_5_1}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_5_2}/trimmed_${sample_5_2}_1.fastq.gz \
        ${sample_5_2}/trimmed_${sample_5_2}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_6_1}/trimmed_${sample_6_1}_1.fastq.gz \
        ${sample_6_1}/trimmed_${sample_6_1}_2.fastq.gz &
    fastqc -t 2 \
        ${sample_6_2}/trimmed_${sample_6_2}_1.fastq.gz \
        ${sample_6_2}/trimmed_${sample_6_2}_2.fastq.gz &

    wait

    echo ""

    echo "Re-run FastQC HTMLs now available. Please take a moment to review."
    read -p "Press 'Enter' to continue. "
fi

echo ""

echo -e "Beginning alignment.\n"

# aligning with reference genomes using Bowtie2
if [[ ! -f "${index}.1.bt2" ]]; then
    echo "Warning: index files matching "${index}" not found. Check the path."
    exit 1
else
    echo -e "Index files matching "${index}" found. Continuing to Bowtie2.\n"
fi

if [[ -f "${sample_1_2}/mapped_${sample_1_2}.bam" && -f "${sample_6_2}/mapped_${sample_6_2}.bam" ]]; then
    echo -e "Skipping Bowtie2 alignment.\n"
else
    echo -e "Beginning Bowtie2 alignment.\n"
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_1_2}/trimmed_${sample_1_2}_1.fastq.gz \
        -2 ${sample_1_2}/trimmed_${sample_1_2}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_1_2}/mapped_${sample_1_2}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_2_1}/trimmed_${sample_2_1}_1.fastq.gz \
        -2 ${sample_2_1}/trimmed_${sample_2_1}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_2_1}/mapped_${sample_2_1}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_2_2}/trimmed_${sample_2_2}_1.fastq.gz \
        -2 ${sample_2_2}/trimmed_${sample_2_2}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_2_2}/mapped_${sample_2_2}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_3_2}/trimmed_${sample_3_2}_1.fastq.gz \
        -2 ${sample_3_2}/trimmed_${sample_3_2}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_3_2}/mapped_${sample_3_2}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_4_1}/trimmed_${sample_4_1}_1.fastq.gz \
        -2 ${sample_4_1}/trimmed_${sample_4_1}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_4_1}/mapped_${sample_4_1}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_4_2}/trimmed_${sample_4_2}_1.fastq.gz \
        -2 ${sample_4_2}/trimmed_${sample_4_2}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_4_2}/mapped_${sample_4_2}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_5_1}/trimmed_${sample_5_1}_1.fastq.gz \
        -2 ${sample_5_1}/trimmed_${sample_5_1}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_5_1}/mapped_${sample_5_1}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_5_2}/trimmed_${sample_5_2}_1.fastq.gz \
        -2 ${sample_5_2}/trimmed_${sample_5_2}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_5_2}/mapped_${sample_5_2}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_6_1}/trimmed_${sample_6_1}_1.fastq.gz \
        -2 ${sample_6_1}/trimmed_${sample_6_1}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_6_1}/mapped_${sample_6_1}.bam &
    bowtie2 \
        -p 12 \
        -x ${index} \
        -1 ${sample_6_2}/trimmed_${sample_6_2}_1.fastq.gz \
        -2 ${sample_6_2}/trimmed_${sample_6_2}_2.fastq.gz \
        | \
        samtools view -bS -> ${sample_6_2}/mapped_${sample_6_2}.bam &

    wait

    echo ""
fi

echo -e "Alignment complete. Continuing to sorting.\n"

# sorting & indexing files
if [[ -f "${sample_6_2}/sorted_mapped_${sample_6_2}.bam" && -f "${sample_1_2}/sorted_mapped_${sample_1_2}.bam" ]]; then
    echo -e "Skipping sorting.\n"
else
    if [[ -f "${sample_1_2}/mapped_${sample_1_2}.bam" && -f "${sample_6_2}/mapped_${sample_6_2}.bam" ]]; then
        echo -e "Now sorting .bam files.\n"
        samtools sort \
            -@ 50 \
            ${sample_1_2}/mapped_${sample_1_2}.bam \
            -o ${sample_1_2}/sorted_mapped_${sample_1_2}.bam &
        samtools sort \
            -@ 50 \
            ${sample_2_1}/mapped_${sample_2_1}.bam \
            -o ${sample_2_1}/sorted_mapped_${sample_2_1}.bam &
        samtools sort \
            -@ 50 \
            ${sample_2_2}/mapped_${sample_2_2}.bam \
            -o ${sample_2_2}/sorted_mapped_${sample_2_2}.bam &
        samtools sort \
            -@ 50 \
            ${sample_3_2}/mapped_${sample_3_2}.bam \
            -o ${sample_3_2}/sorted_mapped_${sample_3_2}.bam &
        samtools sort \
            -@ 50 \
            ${sample_4_1}/mapped_${sample_4_1}.bam \
            -o ${sample_4_1}/sorted_mapped_${sample_4_1}.bam &
        samtools sort \
            -@ 50 \
            ${sample_4_2}/mapped_${sample_4_2}.bam \
            -o ${sample_4_2}/sorted_mapped_${sample_4_2}.bam &
        samtools sort \
            -@ 50 \
            ${sample_5_1}/mapped_${sample_5_1}.bam \
            -o ${sample_5_1}/sorted_mapped_${sample_5_1}.bam &
        samtools sort \
            -@ 50 \
            ${sample_5_2}/mapped_${sample_5_2}.bam \
            -o ${sample_5_2}/sorted_mapped_${sample_5_2}.bam &
        samtools sort \
            -@ 50 \
            ${sample_6_1}/mapped_${sample_6_1}.bam \
            -o ${sample_6_1}/sorted_mapped_${sample_6_1}.bam &
        samtools sort \
            -@ 50 \
            ${sample_6_2}/mapped_${sample_6_2}.bam \
            -o ${sample_6_2}/sorted_mapped_${sample_6_2}.bam &

        wait

        echo -e "Now creating indices for .bam files.\n"
        samtools index \
            -b ${sample_1_2}/sorted_mapped_${sample_1_2}.bam &
        samtools index \
            -b ${sample_2_1}/sorted_mapped_${sample_2_1}.bam &
        samtools index \
            -b ${sample_2_2}/sorted_mapped_${sample_2_2}.bam &
        samtools index \
            -b ${sample_3_2}/sorted_mapped_${sample_3_2}.bam &
        samtools index \
            -b ${sample_4_1}/sorted_mapped_${sample_4_1}.bam &
        samtools index \
            -b ${sample_4_2}/sorted_mapped_${sample_4_2}.bam &
        samtools index \
            -b ${sample_5_1}/sorted_mapped_${sample_5_1}.bam &
        samtools index \
            -b ${sample_5_2}/sorted_mapped_${sample_5_2}.bam &
        samtools index \
            -b ${sample_6_1}/sorted_mapped_${sample_6_1}.bam &
        samtools index \
            -b ${sample_6_2}/sorted_mapped_${sample_6_2}.bam &

        wait

        echo -e "Sorting complete. Moving forward to normalization.\n"
    else
        echo "Error: Bowtie2 output does not seem to exist. Exiting now."
        exit 1
    fi
fi

# normalization of read numbers
if [[ -f "${sample_1_2}/sorted_mapped_${sample_1_2}.bw" && -f "${sample_6_2}/sorted_mapped_${sample_6_2}.bw" ]]; then
    echo -e "Skipping normalization.\n"
else
    if [[ "${norm}" == "RPKM" || "${norm}" == "rpkm" ]]; then
        echo -e "Normalizing using RPKM.\n"

        eval "$(conda shell.bash hook)"
        conda activate deeptools_env

        bamCoverage \
            -b ${sample_1_2}/sorted_mapped_${sample_1_2}.bam \
            -o ${sample_1_2}/sorted_mapped_${sample_1_2}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_2_1}/sorted_mapped_${sample_2_1}.bam \
            -o ${sample_2_1}/sorted_mapped_${sample_2_1}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_2_2}/sorted_mapped_${sample_2_2}.bam \
            -o ${sample_2_2}/sorted_mapped_${sample_2_2}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_3_2}/sorted_mapped_${sample_3_2}.bam \
            -o ${sample_3_2}/sorted_mapped_${sample_3_2}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_4_1}/sorted_mapped_${sample_4_1}.bam \
            -o ${sample_4_1}/sorted_mapped_${sample_4_1}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_4_2}/sorted_mapped_${sample_4_2}.bam \
            -o ${sample_4_2}/sorted_mapped_${sample_4_2}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_5_1}/sorted_mapped_${sample_5_1}.bam \
            -o ${sample_5_1}/sorted_mapped_${sample_5_1}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_5_2}/sorted_mapped_${sample_5_2}.bam \
            -o ${sample_5_2}/sorted_mapped_${sample_5_2}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_6_1}/sorted_mapped_${sample_6_1}.bam \
            -o ${sample_6_1}/sorted_mapped_${sample_6_1}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &
        bamCoverage \
            -b ${sample_6_2}/sorted_mapped_${sample_6_2}.bam \
            -o ${sample_6_2}/sorted_mapped_${sample_6_2}.bw \
            --normalizeUsing RPKM \
            --binSize 20 \
            -p 50 \
            --minMappingQuality 10 &

        wait

    # elif [[ "${norm}" == "CPM" || "${norm}" == "cpm" ]]; then
    #     # echo -e "Normalizing using CPM.\n"

        # eval "$(conda shell.bash hook)"
        # conda activate deeptools_env

        # bamCoverage \
        #     -b ${sample_1_1}/sorted_mapped_${sample_1_1}.bam \
        #     -o ${sample_1_1}/sorted_mapped_${sample_1_1}.bw \
        #     --normalizeUsing CPM \
        #     --binSize 20 \
        #     -p 12 \
        #     --minMappingQuality 10 \
        #     --smoothLength 150 &
        # bamCoverage \
        #     -b ${sample_1_2}/sorted_mapped_${sample_1_2}.bam \
        #     -o ${sample_1_2}/sorted_mapped_${sample_1_2}.bw \
        #     --normalizeUsing CPM \
        #     --binSize 20 \
        #     -p 12 \
        #     --minMappingQuality 10 \
        #     --smoothLength 150 &
        # bamCoverage \
        #     -b ${sample_2_1}/sorted_mapped_${sample_2_1}.bam \
        #     -o ${sample_2_1}/sorted_mapped_${sample_2_1}.bw \
        #     --normalizeUsing CPM \
        #     --binSize 20 \
        #     -p 12 \
        #     --minMappingQuality 10 \
        #     --smoothLength 150 &
        # bamCoverage \
        #     -b ${sample_2_2}/sorted_mapped_${sample_2_2}.bam \
        #     -o ${sample_2_2}/sorted_mapped_${sample_2_2}.bw \
        #     --normalizeUsing CPM \
        #     --binSize 20 \
        #     -p 12 \
        #     --minMappingQuality 10 \
        #     --smoothLength 150 &

        # wait

    elif [[ "${norm}" == "RPGC" || "${norm}" == "rpgc" ]]; then
            echo -e "Normalizing using RPGC.\n"

            eval "$(conda shell.bash hook)"
            conda activate deeptools_env

            bamCoverage \
                -b ${sample_1_2}/sorted_mapped_${sample_1_2}.bam \
                -o ${sample_1_2}/sorted_mapped_${sample_1_2}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_2_1}/sorted_mapped_${sample_2_1}.bam \
                -o ${sample_2_1}/sorted_mapped_${sample_2_1}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_2_2}/sorted_mapped_${sample_2_2}.bam \
                -o ${sample_2_2}/sorted_mapped_${sample_2_2}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_3_2}/sorted_mapped_${sample_3_2}.bam \
                -o ${sample_3_2}/sorted_mapped_${sample_3_2}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_4_1}/sorted_mapped_${sample_4_1}.bam \
                -o ${sample_4_1}/sorted_mapped_${sample_4_1}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_4_2}/sorted_mapped_${sample_4_2}.bam \
                -o ${sample_4_2}/sorted_mapped_${sample_4_2}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_5_1}/sorted_mapped_${sample_5_1}.bam \
                -o ${sample_5_1}/sorted_mapped_${sample_5_1}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_5_2}/sorted_mapped_${sample_5_2}.bam \
                -o ${sample_5_2}/sorted_mapped_${sample_5_2}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_6_1}/sorted_mapped_${sample_6_1}.bam \
                -o ${sample_6_1}/sorted_mapped_${sample_6_1}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &
            bamCoverage \
                -b ${sample_6_2}/sorted_mapped_${sample_6_2}.bam \
                -o ${sample_6_2}/sorted_mapped_${sample_6_2}.bw \
                --normalizeUsing RPGC \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 150 &

            wait

    elif [[ "${norm}" == "spike-in" ]]; then
        echo -e "Normalizing using spike-in control.\n"

        # if [[ ! -f spike_index ]]; then
        #     bowtie2-build \
        #         spike.fa spike_index
        # fi

        # bowtie2 \
        #     -p 24 \
        #     -x spike_index \
        #     -1 ${sample_1_1}/trimmed_${sample_1_1}_1.fastq.gz \
        #     -2 ${sample_1_1}/trimmed_${sample_1_1}_2.fastq.gz \
        #     -S ${sample_1_1}_spike.sam 2> ${sample_1_1}_spike_stats.txt &
        # bowtie2 \
        #     -p 24 \
        #     -x spike_index \
        #     -1 ${sample_1_2}/trimmed_${sample_1_2}_1.fastq.gz \
        #     -2 ${sample_1_2}/trimmed_${sample_1_2}_2.fastq.gz \
        #     -S ${sample_1_2}_spike.sam 2> ${sample_1_2}_spike_stats.txt &
        # bowtie2 \
        #     -p 24 \
        #     -x spike_index \
        #     -1 ${sample_2_1}/trimmed_${sample_2_1}_1.fastq.gz \
        #     -2 ${sample_2_1}/trimmed_${sample_2_1}_2.fastq.gz \
        #     -S ${sample_2_1}_spike.sam 2> ${sample_2_1}_spike_stats.txt &
        # bowtie2 \
        #     -p 24 \
        #     -x spike_index \
        #     -1 ${sample_2_2}/trimmed_${sample_2_2}_1.fastq.gz \
        #     -2 ${sample_2_2}/trimmed_${sample_2_2}_2.fastq.gz \
        #     -S ${sample_2_2}_spike.sam 2> ${sample_2_2}_spike_stats.txt &
        
        # wait

        # declare -A spike_counts
        # for f in ${sample_1_1} ${sample_1_2} ${sample_2_1} ${sample_2_2}; do
        #     concordant_1=$(grep "concordantly exactly 1 time" ${f}_spike_stats.txt | awk '{print $1}')
        #     concordant_multi=$(grep "concordantly >1 times" ${f}_spike_stats.txt | head -1 | awk '{print $1}')
        #     total_align=$((concordant_1 + concordant_multi))
        #     spike_counts[${f}]=${total_align}
        #     echo "${f} spike-in aligned reads: ${total_align}"
        # done
        
        # min_count=""
        # for f in ${sample_1_1} ${sample_1_2} ${sample_2_1} ${sample_2_2}; do
        #     if [[ -z "$min_count" || ${spike_counts[${f}]} -lt $min_count ]]; then
        #         min_count=${spike_counts[${f}]}
        #     fi
        # done

        # echo "min count is: ${min_count}"

        # declare -A sf
        # for f in ${sample_1_1} ${sample_1_2} ${sample_2_1} ${sample_2_2}; do
        #     sf[${f}]=$(awk -v c=${min_count} -v n=${spike_counts[${f}]} 'BEGIN{printf "%.6f", c/n}')
        #     echo "${f} scale factor: ${sf[${f}]}"
        # done

        # eval "$(conda shell.bash hook)"
        # conda activate deeptools_env

        # bamCoverage \
        #     -b ${sample_1_1}/sorted_mapped_${sample_1_1}.bam \
        #     -o ${sample_1_1}/sorted_mapped_${sample_1_1}.bw \
        #     --scaleFactor ${sf[${sample_1_1}]} \
        #     --binSize 20 \
        #     -p 24 \
        #     --minMappingQuality 10 &
        # bamCoverage \
        #     -b ${sample_1_2}/sorted_mapped_${sample_1_2}.bam \
        #     -o ${sample_1_2}/sorted_mapped_${sample_1_2}.bw \
        #     --scaleFactor ${sf[${sample_1_2}]} \
        #     --binSize 20 \
        #     -p 24 \
        #     --minMappingQuality 10 &
        # bamCoverage \
        #     -b ${sample_2_1}/sorted_mapped_${sample_2_1}.bam \
        #     -o ${sample_2_1}/sorted_mapped_${sample_2_1}.bw \
        #     --scaleFactor ${sf[${sample_2_1}]} \
        #     --binSize 20 \
        #     -p 24 \
        #     --minMappingQuality 10 &
        # bamCoverage \
        #     -b ${sample_2_2}/sorted_mapped_${sample_2_2}.bam \
        #     -o ${sample_2_2}/sorted_mapped_${sample_2_2}.bw \
        #     --scaleFactor ${sf[${sample_2_2}]} \
        #     --binSize 20 \
        #     -p 24 \
        #     --minMappingQuality 10 &
        
        # wait
    fi
fi

echo ""

echo -e "Normalization complete. Now merging replicates.\n"

# merging replicates of sample_1 and of sample_2
if [[ -f "${sample_2}_merged.bw" && -f "${sample_6}_merged.bw" ]]; then
    echo -e "Skipping replicate merging.\n"
else
    eval "$(conda shell.bash hook)"
    conda activate deeptools_env

    mv ${sample_1_2}/sorted_mapped_${sample_1_2}.bw ${sample_1}_merged.bw
    mv ${sample_3_2}/sorted_mapped_${sample_3_2}.bw ${sample_3}_merged.bw

    bigwigCompare \
        -b1 ${sample_2_1}/sorted_mapped_${sample_2_1}.bw \
        -b2 ${sample_2_2}/sorted_mapped_${sample_2_2}.bw \
        --operation mean \
        --binSize 50 \
        --outFileFormat bigwig \
        -p 24 \
        --outFileName ${sample_2}_merged.bw &
    bigwigCompare \
        -b1 ${sample_4_1}/sorted_mapped_${sample_4_1}.bw \
        -b2 ${sample_4_2}/sorted_mapped_${sample_4_2}.bw \
        --operation mean \
        --binSize 50 \
        --outFileFormat bigwig \
        -p 24 \
        --outFileName ${sample_4}_merged.bw &
    bigwigCompare \
        -b1 ${sample_5_1}/sorted_mapped_${sample_5_1}.bw \
        -b2 ${sample_5_2}/sorted_mapped_${sample_5_2}.bw \
        --operation mean \
        --binSize 50 \
        --outFileFormat bigwig \
        -p 24 \
        --outFileName ${sample_5}_merged.bw &
    bigwigCompare \
        -b1 ${sample_6_1}/sorted_mapped_${sample_6_1}.bw \
        -b2 ${sample_6_2}/sorted_mapped_${sample_6_2}.bw \
        --operation mean \
        --binSize 50 \
        --outFileFormat bigwig \
        -p 24 \
        --outFileName ${sample_6}_merged.bw &
        
    wait
fi

echo ""

echo -e "Replicate merging complete. Now creating matrix for plotting profiles and heatmaps.\n"

# computing a matrix for profile & heatmap
if [[ -f "${anno}" ]]; then

    eval "$(conda shell.bash hook)"
    conda activate deeptools_env

    computeMatrix scale-regions \
        -S ${sample_1}_merged.bw ${sample_2}_merged.bw ${sample_3}_merged.bw ${sample_4}_merged.bw ${sample_5}_merged.bw ${sample_6}_merged.bw \
        -R ${anno} \
        -o SCP1_log2_matrix \
        -m 5000 \
        --startLabel TSS \
        --endLabel TES \
        -b 2000 \
        -a 2000 \
        --binSize 50 \
        --samplesLabel ${sample_1} ${sample_2} ${sample_3} ${sample_4} ${sample_5} ${sample_6} \
        -p 20 \
        # --missingDataAsZero
else
    echo "The annotation file path is not correct. Exiting now."
    exit 1
fi

echo -e "Matrix creation complete. Now moving forward to plotting profiles and heatmaps.\n"

eval "$(conda shell.bash hook)"
conda activate deeptools_env

# creating a profile
plotProfile \
    -m SCP1_matrix \
    -o SCP1_profile.pdf \
    --averageType mean \
    --yAxisLabel Log2Input \
    --plotTitle SCP1_profile \
    --legendLocation best \
    --perGroup \
    --startLabel TSS \
    --endLabel TES 

eval "$(conda shell.bash hook)"
conda activate deeptools_env

# creating a heatmap
plotHeatmap \
    -m SCP1_matrix \
    -o SCP1_heatmap.pdf \
    --sortRegions descend \
    --linesAtTickMarks \
    --sortUsing mean \
    --averageTypeSummaryPlot mean \
    --plotType lines \
    --xAxisLabel "distance_from TSS (bp)" \
    --heatmapWidth 6 \
    --heatmapHeight 16 \
    --startLabel TSS \
    --endLabel TES \
    --refPointLabel TSS \
    --legendLocation best \
    --colorMap RdBu_r \
    --missingDataColor 1 \
    --yAxisLabel Average_Signal \
    --samplesLabel ${sample_1} ${sample_2} ${sample_3} ${sample_4} ${sample_5} ${sample_6}

echo ""

echo -e "Profile & heatmap plotting complete. Upstream analysis now complete."
exit 0
