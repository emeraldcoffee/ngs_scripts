#!/usr/bin/env bash
set -e

# ----------------------
# This script is intended to be used for the processing of CUT&Run, ChIP, or other protein-DNA
# interaction sequencing. It is intended to be used for experiments where there are no replicates
# for each sample.
# 
# Each step of this pipeline will skip if the 1st sample's 1st replicate's associated output
# file exists. If the previous run of the script failed in that step, please delete at least
# the 1st sample's 1st replicate associated output file. 
# 
# This script expects the user to have a conda environment available named "deeptools_env" with
# the deepTools suite installed.
# 
# Note that if the user is using spike-in normalization, either spike.fa or spike_index should
# be available in the usr directory.
# ----------------------

# ----------------------
# This script expects a file structure like the following:
# /usr
# |-- sample_1
# |   |-- sample_1_R1.fq.gz
# |   |-- sample_1_R2.fq.gz
# |-- sample_2
# |   |-- sample_2_R1.fq.gz
# |   |-- sample_2_R2.fq.gz
# ----------------------

# ----------------------
# To run this script, enter the following into a Linux terminal while in your usr directory:
# bash DNA_pairedEnd_noRep.sh \
#   -s {number_of_samples} {sample_name_1} {sample_name_2} {etc} \
#   -S {job_name} \
#   -i {Bowtie2_index_files} \
#   -n {normalization_method} \
#   -a {genomic_annotation_file}
#
# For example:
# bash DNA_pairedEnd_noRep.sh \
#   -s 3 HSF_1 E2F4 CEBPa \
#   -S H3K56 \
#   -i /bt2_index/human/HG38 \
#   -n RPKM \
#   -a /anno/gencode.v38.annotation.gtf
# ----------------------

usage() {
    echo "Help: This script takes in -s <number_of_samples> <sample_name_1> <sample_name_2> <etc>,"
    echo "-S <job_name>, -i <Bowtie2_index_files>, -n <normalization method: CPM, RPKM, or spike-in>,"
    echo "-a <genomic_annotation_file>, and -h for help."
}

echo ""

echo "SETTINGS:"

while getopts "s:S:i:n:a:h" opt; do
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

                OPTIND=$((OPTIND + 1))
            done
            ;;
        S)
            job_name="$OPTARG"
            echo "Job name: $job_name"
            ;;
        i)
            index="$OPTARG"
            echo "Bowtie2 index path: $index."
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

echo -e "Beginning paired-end DNA sequencing upstream analysis.\n"

# runs FastQC HTML generation for users to check sequencing quality
if compgen -G "${sample_1}/${sample_1}_R1_fastqc.html" > /dev/null; then
    echo -e "Skipping FastQC HTML generation.\n"
else
    echo -e "Starting FastQC HTML generation.\n"

    for ((i=1;i<=num_samples; i++)); do 
        sample="sample_${i}"

        # if sequencing data is in .fastq.gz format, changes the file extension to fq.gz
        if compgen -G "${!sample}"/"${!sample}"_R1.fq.gz > /dev/null; then
            mv "${!sample}"/"${!sample}"_R1.fastq.gz "${!sample}"/"${!sample}"_R1.fq.gz
            mv "${!sample}"/"${!sample}"_R2.fastq.gz "${!sample}"/"${!sample}"_R2.fq.gz
        fi

        fastqc -t 2 \
            "${!sample}"/*_R1.fastq.gz \
            "${!sample}"/*_R2.fastq.gz &
            
    done

    wait

    echo ""

    echo "FastQC HTMLs now available. Please take a moment to review."

    read -p "Press 'Enter' to continue. "
    echo ""

fi

# uses cutadapt to trim off adapters and poly-G tails from sequencing data
if compgen -G "${sample_1}/trimmed_${sample_1}_R1.fastq.gz" > /dev/null; then
    echo "Skipping cutadapt trimming."
else
    echo -e "Now trimming with cutadapt.\n"

    for ((i=1;i<=num_samples; i++)); do
        sample="sample_${i}"

        cutadapt \
            --trim-n \
            -m 20 \
            -q 20 \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -a "G{20}" \
            -A "G{20}" \
            -o "${!sample}"/trimmed_"${!sample}"_R1.fastq.gz \
            -p "${!sample}"/trimmed_"${!sample}"_R2.fastq.gz \
            "${!sample}"/*_R1.fastq.gz \
            "${!sample}"/*_R2.fastq.gz \
            -j 12 &

    done
    
    wait

    echo -e "Trimming now complete.\n"

fi

# if quality before trimming was subpar (or had overrepresented sequences or high adapter content),
# users may recheck FastQC HTMLs
read -p "Would you like to re-run FastQC HTML generation? <yes/no>: "
if [[ ${REPLY} = "yes" ]]; then
    for ((i=1;i<=num_samples; i++)); do
        sample="sample_${i}"

        fastqc -t 2 \
            "${!sample}"/trimmed_"${!sample}"_R1.fastq.gz \
            "${!sample}"/trimmed_"${!sample}"_R2.fastq.gz &

    done

    wait

    echo ""

    echo "Re-run FastQC HTMLs now available. Please take a moment to review."
    read -p "Press 'Enter' to continue. "
fi

echo ""

# aligning reads with reference genomes using Bowtie2
if [[ ! -f "${index}.1.bt2" ]]; then
    echo "Warning: index files matching "${index}" not found. Check the path."
    exit 1
else
    echo -e "Index files matching "${index}" found. Continuing to Bowtie2.\n"
fi

if compgen -G "${sample_1}/mapped_${sample_1}.bam" > /dev/null; then
    echo -e "Skipping Bowtie2 alignment.\n"
else
    echo -e "Beginning Bowtie2 alignment.\n"

    for ((i=1;i<=num_samples; i++)); do
        sample="sample_${i}"

        bowtie2 \
            -p 12 \
            -x ${index} \
            -1 "${!sample}"/trimmed_"${!sample}"_R1.fastq.gz \
            -2 "${!sample}"/trimmed_"${!sample}"_R2.fastq.gz \
            | samtools view -bS -> "${!sample}"/mapped_"${!sample}".bam &
    done

    wait

    echo -e "Alignment complete.\n"
    
fi

# sorting & indexing aligned files
if compgen -G "${sample_1}/sorted_mapped_${sample_1}.bam" > /dev/null; then
    echo -e "Skipping sorting & indexing.\n"
else
    if compgen -G "${sample_1}/mapped_${sample_1}.bam" > /dev/null; then
        echo -e "Now sorting .bam files.\n"

        for ((i=1;i<=num_samples; i++)); do
            sample="sample_${i}"

            samtools sort \
                -@ 50 \
                "${!sample}"/mapped_"${!sample}".bam \
                -o "${!sample}"/sorted_mapped_"${!sample}".bam &
        done

        wait

        echo -e "Now creating indices for .bam files.\n"
        for ((i=1;i<=num_samples; i++)); do
            sample="sample_${i}"

            samtools index \
                -b "${!sample}"/sorted_mapped_"${!sample}".bam &
        done

        wait

        echo -e "Sorting complete. Moving forward to normalization.\n"
    else
        echo "Error: Bowtie2 output does not seem to exist. Exiting now."
        exit 1
    fi
fi

# performs normalization of sequence reads
if [[ -f "${sample_1}/sorted_mapped_${sample_1}.bw" ]]; then
    echo "Skipping normalization."
else
    if [[ "${norm}" == "RPKM" || "${norm}" == "rpkm" ]]; then
        echo -e "Normalizing using RPKM with a bin size of 20bp.\n"

        eval "$(conda shell.bash hook)"
        conda activate deeptools_env

        for ((i=1;i<=num_samples; i++)); do
            sample="sample_${i}"

            bamCoverage \
                -b "${!sample}"/sorted_mapped_"${!sample}".bam \
                -o "${!sample}"/sorted_mapped_"${!sample}".bw \
                --normalizeUsing RPKM \
                --binSize 20 \
                -p 50 \
                --minMappingQuality 10 &
        done

        wait

    elif [[ "${norm}" == "CPM" || "${norm}" == "cpm" ]]; then
        echo -e "Normalizing using CPM with a bin size of 20 and a smooth length of 60bp.\n"

        eval "$(conda shell.bash hook)"
        conda activate deeptools_env

        for ((i=1;i<=num_samples; i++)); do
            sample="sample_${i}"

            bamCoverage \
                -b "${!sample}"/sorted_mapped_"${!sample}".bam \
                -o "${!sample}"/sorted_mapped_"${!sample}".bw \
                --normalizeUsing CPM \
                --binSize 20 \
                -p 12 \
                --minMappingQuality 10 \
                --smoothLength 60 &
        done

        wait

    elif [[ "${norm}" == "spike-in" ]]; then
        echo -e "Normalizing using spike-in control.\n"

        if [[ ! -f spike_index ]]; then
            bowtie2-build \
                spike.fa spike_index
        fi

        for ((i=1;i<=num_samples; i++)); do
            sample="sample_${i}"

            bowtie2 \
                -p 24 \
                -x spike_index \
                -1 "${!sample}"/trimmed_"${!sample}"_R1.fastq.gz \
                -2 "${!sample}"/trimmed_"${!sample}"_R2.fastq.gz \
                -S "${!sample}"_spike.sam 2> "${!sample}"_spike_stats.txt &
        
            wait

            declare -A spike_counts
            concordant_1=$(grep "concordantly exactly 1 time" "${!sample}"_spike_stats.txt | awk '{print $1}')
            concordant_multi=$(grep "concordantly >1 times" "${!sample}"_spike_stats.txt | head -1 | awk '{print $1}')
            total_align=$((concordant_1 + concordant_multi))
            spike_counts["${!sample}"]=${total_align}
            echo ""${!sample}" spike-in aligned reads: ${total_align}"
            
            min_count=""
            if [[ -z "$min_count" || ${spike_counts["${!sample}"]} -lt $min_count ]]; then
                min_count=${spike_counts["${!sample}"]}
            fi

            declare -A sf
            sf["${!sample}"]=$(awk -v c=${min_count} -v n=${spike_counts["${!sample}"]} 'BEGIN{printf "%.6f", c/n}')
            echo ""${!sample}" scale factor: ${sf["${!sample}"]}"
        done

        eval "$(conda shell.bash hook)"
        conda activate deeptools_env

        for ((i=1;i<=num_samples; i++)); do
            sample="sample_${i}"

            bamCoverage \
                -b "${!sample}"/sorted_mapped_"${!sample}".bam \
                -o "${!sample}"/sorted_mapped_"${!sample}".bw \
                --scaleFactor ${sf["${!sample}"]} \
                --binSize 20 \
                -p 24 \
                --minMappingQuality 10 &
        done
        
        wait
    fi

echo ""

echo -e "Normalization complete. \n"

echo -e "Now creating matrix for plotting profiles and heatmaps with average signal.\n"

# computing a matrix for plotting profiles & heatmaps
if [[ -f "${anno}" ]]; then

    eval "$(conda shell.bash hook)"
    conda activate deeptools_env

    for ((i=1;i<=num_samples; i++)); do
        name="sample_${i}"
        computing_samples_file+="${!name}/sorted_mapped_${!name}.bw "
        computing_samples+="${!name} "
    done

    computeMatrix scale-regions \
        -S ${computing_samples_file} \
        -R ${anno} \
        -o ${job_name}_matrix \
        -m 5000 \
        --startLabel TSS \
        --endLabel TES \
        -b 2000 \
        -a 2000 \
        --binSize 50 \
        --samplesLabel ${computing_samples} \
        -p 20 \
        --missingDataAsZero
else
    echo "The annotation file path is not correct. Exiting now."
    exit 1
fi

echo -e "Matrix creation complete. Now moving forward to plotting profiles and heatmaps.\n"

eval "$(conda shell.bash hook)"
conda activate deeptools_env

# plotting average signal profiles
plotProfile \
    -m ${job_name}_matrix \
    -o ${job_name}_profile.pdf \
    --averageType mean \
    --yAxisLabel Average_Signal \
    --plotTitle ${job_name} \
    --legendLocation best \
    --perGroup \
    --startLabel TSS \
    --endLabel TES 

eval "$(conda shell.bash hook)"
conda activate deeptools_env

# plotting average signal heatmap
plotHeatmap \
    -m ${job_name}_matrix \
    -o ${job_name}_heatmap.pdf \
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
    --yAxisLabel Log2Input \
    --samplesLabel ${computing_samples}

echo ""

echo -e "Profile & heatmap plotting complete. Upstream analysis now complete.\n"

exit 0
