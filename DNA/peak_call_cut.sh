#!/usr/bin/env bash
set -e

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
                declare -g "sample_${i}_1=${sample_name}_1"
                declare -g "sample_${i}_2=${sample_name}_2"

                OPTIND=$((OPTIND + 1))
            done
            ;;
        S)
            sample_name="$OPTARG"
            ;;
    esac
done

#macs3 callpeak
if [[ -f "macs3_broad/${sample_1_1}_broad_peaks.broadPeak" ]]; then
    echo "Skipping MACS3 peak calling."
else
    echo "Beginning MACS3 peak calling."

    eval "$(conda shell.bash hook)"
    conda activate macs3

    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        macs3 callpeak \
            -t "${!sample1}"/sorted_mapped_"${!sample1}".bam \
            -c controls/phos/sorted_mapped_IgG_1.bam controls/phos/sorted_mapped_IgG_2.bam \
            -n "${!sample1}"_broad \
            -f BAMPE \
            -g hs \
            --outdir macs3_broad \
            --broad \
            --broad-cutoff 0.1 &
        macs3 callpeak \
            -t "${!sample2}"/sorted_mapped_"${!sample2}".bam \
            -c controls/phos/sorted_mapped_IgG_1.bam controls/phos/sorted_mapped_IgG_2.bam \
            -n "${!sample2}"_broad \
            -f BAMPE \
            -g hs \
            --outdir macs3_broad \
            --broad \
            --broad-cutoff 0.1 &
    done

    wait
fi

# creating 
if [[ -f "${sample_1_1}_${sample_name}.tab" ]]; then
    echo "Skipping MACS3 peak calling."
else
    eval "$(conda shell.bash hook)"
    conda activate bigwig
    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        bigWigAverageOverBed \
            "${!sample1}"/sorted_mapped_"${!sample1}".bw \
            hg38_genes.bed \
            "${!sample1}"_${!sample_name}.tab &
        bigWigAverageOverBed \
            "${!sample2}"/sorted_mapped_"${!sample2}".bw \
            hg38_genes.bed \
            "${!sample2}"_${!sample_name}.tab &
    done

    wait
fi

if [[ -f "${sample_1_1}_90.txt" ]]; then
    echo "Skipping 10% cutoff thresholding."
else
    for ((i=1;i<=num_samples; i++)); do
        sample1="sample_${i}_1"
        sample2="sample_${i}_2"

        cutoff=$(awk 'NR > 1 {print $5}' ${!sample1}_${!sample_name}.tab | \
            sort -n | \
            awk '{
                values[NR] = $1
            }
            END {
                n = NR
                pos = int(0.10 * (n - 1)) + 1
                print values[pos]
            }')
            
        echo "10th percentile cutoff: $cutoff"

        awk -v cutoff="$cutoff" 'NR > 1 && $5 >= cutoff {
            print $1
        }' ${!sample1}_${!sample_name}.tab > ${!sample1}_90.txt &

        cutoff=$(awk 'NR > 1 {print $5}' ${!sample2}_${!sample_name}.tab | \
            sort -n | \
            awk '{
                values[NR] = $1
            }
            END {
                n = NR
                pos = int(0.10 * (n - 1)) + 1
                print values[pos]
            }')
            
        echo "10th percentile cutoff: $cutoff"

        awk -v cutoff="$cutoff" 'NR > 1 && $5 >= cutoff {
            print $1
        }' ${!sample2}_${!sample_name}.tab > ${!sample2}_90.txt &
    done
fi
