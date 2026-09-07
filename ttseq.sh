#!/usr/bin/env bash
set -euo pipefail
set -x

THREADS=16
BIN=50
UP=5000
DOWN=5000

mkdir -p matrices profiles logs

# -----------------------------
# bigWig files (replicates)
# plus genes -> forward strand
# minus genes -> reverse strand
# order must match LABELS exactly
# -----------------------------
F_BWS=(
 "WT_forward_log2ratio.bw"
 "T4V_forward_log2ratio.bw"
)

R_BWS=(
 "WT_reverse_log2ratio.bw"
 "T4V_reverse_log2ratio.bw"
)

LABELS=(
 "WT_log2"
 "T4V_log2"
)

# -----------------------------
# BED sets
# strand-specific, non-overlap, mutually exclusive length bins
# -----------------------------
declare -A PLUS_BEDS
declare -A MINUS_BEDS

PLUS_BEDS[nonoverlap_all]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_plus.nonoverlap.bed"
MINUS_BEDS[nonoverlap_all]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_minus.nonoverlap.bed"

PLUS_BEDS[len_lt5kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_plus.nonoverlap.lt5kb.bed"
MINUS_BEDS[len_lt5kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_minus.nonoverlap.lt5kb.bed"

PLUS_BEDS[len_ge5kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_plus.nonoverlap.ge5kb.bed"
MINUS_BEDS[len_ge5kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_minus.nonoverlap.ge5kb.bed"

PLUS_BEDS[len_5to50kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_plus.nonoverlap.5to50kb.bed"
MINUS_BEDS[len_5to50kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_minus.nonoverlap.5to50kb.bed"

PLUS_BEDS[len_50to120kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_plus.nonoverlap.50to120kb.bed"
MINUS_BEDS[len_50to120kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_minus.nonoverlap.50to120kb.bed"

PLUS_BEDS[len_gt120kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_plus.nonoverlap.gt120kb.bed"
MINUS_BEDS[len_gt120kb]="/stor/work/ZhangYJ/Qian_test/annotation/ttseq_bed/strand_specific_length_bins/hg38_PCG_minus.nonoverlap.gt120kb.bed"

BED_KEYS=(
 nonoverlap_all
 len_lt5kb
 len_ge5kb
 len_5to50kb
 len_50to120kb
 len_gt120kb
)

echo "[INFO] forward bigWigs: ${F_BWS[*]}"
echo "[INFO] reverse bigWigs: ${R_BWS[*]}"
echo "[INFO] labels: ${LABELS[*]}"

# -----------------------------
# loop over bed groups
# -----------------------------
for setname in "${BED_KEYS[@]}"; do
 PLUS_BED="${PLUS_BEDS[$setname]}"
 MINUS_BED="${MINUS_BEDS[$setname]}"

 [[ -f "$PLUS_BED" ]] || { echo "[ERROR] missing plus bed for ${setname}: ${PLUS_BED}"; exit 1; }
 [[ -f "$MINUS_BED" ]] || { echo "[ERROR] missing minus bed for ${setname}: ${MINUS_BED}"; exit 1; }

 # shorter body for short genes
 if [[ "$setname" == "len_lt5kb" ]]; then
  BODY=3000
 else
  BODY=15000
 fi

 echo "========================================"
 echo "[INFO] Processing bed set: ${setname}"
 echo "[INFO] PLUS bed: ${PLUS_BED}"
 echo "[INFO] MINUS bed: ${MINUS_BED}"
 echo "[INFO] regionBodyLength: ${BODY}"
 echo "========================================"

 # ----------------------------------
 # plus genes: use forward-strand bigWigs
 # ----------------------------------
 echo "[INFO] scale-regions matrix for plus genes ..."
 computeMatrix scale-regions \
  -S "${F_BWS[@]}" \
  -R "${PLUS_BED}" \
  -b ${UP} -a ${DOWN} \
  --regionBodyLength ${BODY} \
  --binSize ${BIN} \
  --skipZeros \
  --missingDataAsZero \
  --numberOfProcessors ${THREADS} \
  --samplesLabel "${LABELS[@]}" \
  -o matrices/${setname}.plus.scale.gz \
  2> logs/${setname}.plus.scale.log

 # ----------------------------------
 # minus genes: use reverse-strand bigWigs
 # ----------------------------------
 echo "[INFO] scale-regions matrix for minus genes ..."
 computeMatrix scale-regions \
  -S "${R_BWS[@]}" \
  -R "${MINUS_BED}" \
  -b ${UP} -a ${DOWN} \
  --regionBodyLength ${BODY} \
  --binSize ${BIN} \
  --skipZeros \
  --missingDataAsZero \
  --numberOfProcessors ${THREADS} \
  --samplesLabel "${LABELS[@]}" \
  -o matrices/${setname}.minus.scale.gz \
  2> logs/${setname}.minus.scale.log

 # ----------------------------------
 # merge plus/minus -> sense matrix
 # ----------------------------------
 echo "[INFO] Merge scale sense matrix ..."
 computeMatrixOperations rbind \
  -m matrices/${setname}.plus.scale.gz matrices/${setname}.minus.scale.gz \
  -o matrices/${setname}.sense.scale.gz

 # ----------------------------------
 # profile plot
 # ----------------------------------
 plotProfile \
  -m matrices/${setname}.sense.scale.gz \
  --plotTitle "TTseq_scale_Sense_${setname}" \
  --regionsLabel "genes" \
  --samplesLabel "${LABELS[@]}" \
  --perGroup \
  --legendLocation upper-right \
  -out profiles/TTseq_scale_Sense_${setname}.svg

plotHeatmap \
    -m matrices/${setname}.sense.scale.gz \
    -o heatmaps/TTseq_scale_Sense_${setname}.svg \
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
    --yAxisLabel genes \
    --samplesLabel "${LABELS[@]}"

done

echo "========================================"
echo "[DONE] All scale-region matrices and plots are generated."
echo "[DONE] Matrices: matrices/"
echo "[DONE] Plots:  profiles/"
echo "[DONE] Logs:   logs/"
echo "========================================"
set +x
