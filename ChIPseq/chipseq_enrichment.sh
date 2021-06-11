#!/bin/bash

mkdir deeptools

## TSS enrichment

awk '$1 == "X"' Mus_musculus.GRCm38.97.chr.gtf > GRCm38.97.chrX.gtf

parallel -j 8 \
computeMatrix reference-point \
-p 8 \
-S bigwig/*_{1}_*_{2}_*_{3}.bw \
-R GRCm38.97.chrX.gtf \
--upstream 5000 \
--downstream 5000 \
--skipZeros \
--nanAfterEnd \
--outFileName deeptools/{1}_{2}_{3}_TSS_enrichment.mat.gz \
::: H2AK119Ub H3K4me1 H3K4me3 H3K9ac H3K27ac H3K27me3 H4ac Input ::: 0h 4h 8h 12h 24h ::: C57 CAST

## Enhancer enrichment
mkdir peaks/consensus/enhancer

# extract TSS and extend 1kb each direction
awk '$3 == "gene" && $1 == "X" {print ($7=="+")?$1"\t"$4-1000"\t"$4+1000:$1"\t"$5-1000"\t"$5+1000}' GRCm38.97.chrX.gtf > GRCm38.97.chrX.TSS.bed

# intersect H3K4me1 & H3k27ac and exclude promoter regions
get_enhancer(){
	bedtools intersect \
	-a $1 \
	-b $2 \
	| awk '$1 == "X"{print $1"\t"$2"\t"$3}' \
	| bedtools subtract -A -a stdin -b $3
}
export -f get_enhancer

parallel --link get_enhancer peaks/consensus/ChIPseq_H3K4me1_TX1072_WT_{}_Dox_TOTAL_peaks_consensus_filtered.broadPeak peaks/consensus/ChIPseq_H3K27ac_TX1072_WT_{}_Dox_TOTAL_peaks_consensus_filtered.broadPeak GRCm38.97.chrX.TSS.bed  ">" peaks/consensus/enhancer/Enhancers_{}.bed ::: 0h 4h 8h 12h 24h
cat peaks/consensus/enhancer/Enhancers_* | bedtools sort -i stdin | bedtools merge -d 500 -i stdin > peaks/consensus/enhancer/Enhancers_merged.bed

parallel -j 8 \
computeMatrix reference-point \
-p 8 \
-S bigwig/*_{1}_*_{2}_*_{3}.bw \
-R peaks/consensus/enhancer/Enhancers_merged.bed \
--referencePoint center \
--upstream 5000 \
--downstream 5000 \
--skipZeros \
--outFileName deeptools/{1}_{2}_{3}_Enhancer_enrichment.mat.gz \
::: H2AK119Ub H3K4me1 H3K4me3 H3K9ac H3K27ac H3K27me3 H4ac Input ::: 0h 4h 8h 12h 24h ::: C57 CAST
