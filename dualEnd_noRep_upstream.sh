#!/usr/bin/env bash
set -e

# arguments
usage() {
    echo "Help: This script takes in -s <sample_1_name>, -S <sample_2_name>,"
    echo "-i <Bowtie2_index_path>, -n <normalization method: CPM, RPKM, or spike-in>,"
    echo "-a <annotation path>, and -h for help."
}

while getopts "s:S:i:n:a:h" opt; do
    echo "SETTINGS:"
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

                declare -g "sample_${i}=${sample_name}" # sample_1=30min_1, sample_1_1=30min_1

                OPTIND=$((OPTIND + 1))
            done
            ;;
        S)
            exp_name="$OPTARG"
            echo "Experiment name: $exp_name"
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

echo -e "Beginning dual end sequencing upstream analysis.\n"

# FastQC HTML generation to check sequencing quality
if compgen -G "${sample_1}/*_2_fastqc.html" > /dev/null; then
    echo -e "Skipping FastQC HTML generation.\n"
else
    echo "Starting FastQC HTML generation."

    for ((i=1;i<=num_samples; i++)); do
        sample="sample_${i}"

        fastqc -t 2 \
            "${!sample}"/*_1.fastq.gz \
            "${!sample}"/*_2.fastq.gz &
    done

    wait

    echo ""

    echo "FastQC HTMLs now available. Please take a moment to review."

    read -p "Press 'Enter' to continue. "
    echo ""
fi

echo -e "Now trimming with cutadapt.\n"

# trimming adapters, poly-G tails
if compgen -G "${sample_1}/trimmed_*_2.fastq.gz" > /dev/null; then
    echo -e "Skipping cutadapt trimming.\n"
else
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
            -o "${!sample}"/trimmed_"${!sample}"_1.fastq.gz \
            -p "${!sample}"/trimmed_"${!sample}"_2.fastq.gz \
            "${!sample}"/*_1.fastq.gz \
            "${!sample}"/*_2.fastq.gz \
            -j 12 &

    done
    
    wait

    echo ""
fi

echo -e "Trimming now complete.\n"

# if quality before trimming was subpar (or had overrepresented sequences or high adapter content), recheck FastQC HTMLs
read -p "Would you like to re-run FastQC HTML generation? <yes/no>: "
if [[ ${REPLY} = "yes" ]]; then
    for ((i=1;i<=num_samples; i++)); do
        sample="sample_${i}"

        fastqc -t 2 \
            "${!sample}"/trimmed_"${!sample}"_1.fastq.gz \
            "${!sample}"/trimmed_"${!sample}"_2.fastq.gz &
        fastqc -t 2 \
            "${!}"/trimmed_"${!}"_1.fastq.gz \
            "${!}"/trimmed_"${!}"_2.fastq.gz &
    done

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

if compgen -G "${sample_1}/mapped_*.bam" > /dev/null; then
    echo -e "Skipping Bowtie2 alignment.\n"
else
    echo -e "Beginning Bowtie2 alignment.\n"

    for ((i=1;i<=num_samples; i++)); do
        sample="sample_${i}"

        bowtie2 \
            -p 12 \
            -x ${index} \
            -1 "${!sample}"/trimmed_"${!sample}"_1.fastq.gz \
            -2 "${!sample}"/trimmed_"${!sample}"_2.fastq.gz \
            | \
            samtools view -bS -> "${!sample}"/mapped_"${!sample}".bam &
    done

    wait

    echo ""
fi

echo -e "Alignment complete. Continuing to sorting.\n"

# sorting & indexing files
if compgen -G "${sample_1}/sorted_mapped_*.bam" > /dev/null; then
    echo -e "Skipping sorting.\n"
else
    if compgen -G "${sample_1}/mapped_*.bam" > /dev/null; then
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

# normalization of read numbers
if [[ -f "${sample_1}/sorted_mapped_${sample_1}.bw" ]]; then
    echo -e "Skipping normalization.\n"
else
    if [[ "${norm}" == "RPKM" || "${norm}" == "rpkm" ]]; then
        echo -e "Normalizing using RPKM.\n"

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
        echo -e "Normalizing using CPM.\n"

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
                -1 "${!sample}"/trimmed_"${!sample}"_1.fastq.gz \
                -2 "${!sample}"/trimmed_"${!sample}"_2.fastq.gz \
                -S "${!sample}"_spike.sam 2> "${!sample}"_spike_stats.txt &
        
            wait

            declare -A spike_counts
            for f in "${!sample}" "${!}"; do
                concordant_1=$(grep "concordantly exactly 1 time" ${f}_spike_stats.txt | awk '{print $1}')
                concordant_multi=$(grep "concordantly >1 times" ${f}_spike_stats.txt | head -1 | awk '{print $1}')
                total_align=$((concordant_1 + concordant_multi))
                spike_counts[${f}]=${total_align}
                echo "${f} spike-in aligned reads: ${total_align}"
            done
            
            min_count=""
            for f in "${!sample}" "${!}"}; do
                if [[ -z "$min_count" || ${spike_counts[${f}]} -lt $min_count ]]; then
                    min_count=${spike_counts[${f}]}
                fi
            done

            echo "min count is: ${min_count}"

            declare -A sf
            for f in "${!sample}" "${!}"}; do
                sf[${f}]=$(awk -v c=${min_count} -v n=${spike_counts[${f}]} 'BEGIN{printf "%.6f", c/n}')
                echo "${f} scale factor: ${sf[${f}]}"
            done
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
fi

echo ""

echo -e "Normalization complete. Now merging replicates.\n"

# merging replicates of sample_1 and of sample_2
if [[ -f "${sample_1}_merged.bw" ]]; then
    echo -e "Skipping replicate merging.\n"
else
    eval "$(conda shell.bash hook)"
    conda activate deeptools_env

    for ((i=1;i<=num_samples; i++)); do
        sample="sample_${i}"

        cp "${!sample}"/sorted_mapped_"${!sample}".bw ${!sample}_merged.bw

    done
        
    wait
fi

echo ""

echo -e "Replicate merging complete. Now creating matrix for plotting profiles and heatmaps.\n"

# computing a matrix for profile & heatmap
if [[ -f "${anno}" ]]; then

    eval "$(conda shell.bash hook)"
    conda activate deeptools_env

    for ((i=1;i<=num_samples; i++)); do
        name="sample_${i}"
        computing_samples_file+="${!name}_merged.bw "
        computing_samples+="${!name} "
    done

    computeMatrix scale-regions \
        -S ${computing_samples_file} \
        -R ${anno} \
        -o ${exp_name}_matrix \
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

# creating a profile
plotProfile \
    -m ${exp_name}_matrix \
    -o ${exp_name}_profile.pdf \
    --averageType mean \
    --yAxisLabel Average_Signal \
    --plotTitle ${exp_name} \
    --legendLocation best \
    --perGroup \
    --startLabel TSS \
    --endLabel TES 

eval "$(conda shell.bash hook)"
conda activate deeptools_env

# creating a heatmap
plotHeatmap \
    -m ${exp_name}_matrix \
    -o ${exp_name}_heatmap.pdf \
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

echo -e "Profile & heatmap plotting complete. Upstream analysis now complete."
exit 0
