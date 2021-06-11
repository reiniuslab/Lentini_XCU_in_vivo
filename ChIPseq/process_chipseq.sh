#!/bin/bash
## download data

mkdir fastq

download_fastq(){
	echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${1:0:6}/00${1:9:1}/${1}/${1}_{1..2}.fastq.gz | xargs -n 1 -P 2 wget -q -P fastq
}
export -f download_fastq

nohup parallel -j 10 download_fastq ::: SRR75143{15..94} > wget.log &

## trim data
mkdir clean clean/fastp_reports
parallel mkdir clean/fastp_reports/{} ::: SRR75143{15..94}

nohup parallel -j 64 \
fastp \
--in1 fastq/{}_1.fastq.gz \
--in2 fastq/{}_2.fastq.gz \
--out1 clean/{}_1.fastq.gz \
--out2 clean/{}_2.fastq.gz \
--json clean/fastp_reports/{}/fastp.json --html clean/fastp_reports/{}/fastp.html -R {} \
::: SRR75143{15..94} > fastp.log & # 7138

## QC
multiqc -d -f -o . clean/fastp_reports

## align data
mkdir bam

map_bowtie(){
	echo $1
	bowtie2 \
	-N 1 --mm -p 16 \
	-x bowtie2/GRCm38.mgp.v5.CAST_EiJ_Nmasked/index \
	-1  clean/${1}_1.fastq.gz \
	-2  clean/${1}_2.fastq.gz \
	| bamsort sortthreads=16 readthreads=16 writethreads=16 blockmb=10000 inputformat=sam markduplicates=1 index=1 indexfilename=bam/${1}.bam.bai O=bam/${1}.bam
}
export -f map_bowtie

nohup parallel -j 6 map_bowtie ::: SRR75143{15..94} > bowtie2.log &

## split alleles
mkdir allelic

nohup parallel -j 12 SNPsplit -o allelic --snp_file all_SNPs_CAST_EiJ_GRCm38.txt.gz {} ::: bam/*.bam > snpsplit.log & # 31305
# sort output
parallel  -j 8 sambamba -q sort -m 10GB -t 12 ::: allelic/*.genome{1..2}.bam

## bigwig tracks
mkdir bigwig bigwig/relative

# get ENCODE blacklist and convert to ENSEMBL chromosome names
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
zcat mm10-blacklist.v2.bed.gz | awk '{gsub(/^chr/,""); print}' | gzip > mm10-blacklist.v2.ENSEMBL.bed.gz

# use metadata to input run accession and output run name
nohup parallel -j 8 --colsep '\t' --header : bamCoverage -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --skipNAs --ignoreDuplicates --centerReads --blackListFileName mm10-blacklist.v2.ENSEMBL.bed.gz --bam bam/{1}.bam --outFileName bigwig/{2}_TOTAL.bw :::: PRJNA480803.tsv > bamCoverageTOTAL.log &
nohup parallel -j 8 --colsep '\t' --header : bamCoverage -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --skipNAs --ignoreDuplicates --centerReads --blackListFileName mm10-blacklist.v2.ENSEMBL.bed.gz --bam allelic/{1}.genome1.sorted.bam --outFileName bigwig/{2}_C57.bw :::: PRJNA480803.tsv > bamCoverageC57.log &
nohup parallel -j 8 --colsep '\t' --header : bamCoverage -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --skipNAs --ignoreDuplicates --centerReads --blackListFileName mm10-blacklist.v2.ENSEMBL.bed.gz --bam allelic/{1}.genome2.sorted.bam --outFileName bigwig/{2}_CAST.bw :::: PRJNA480803.tsv > bamCoverageCAST.log &

## peak calling
mkdir peaks

# call vs input (7 modifications, repeat input 7 times)
nohup parallel -j 12 --link --colsep '\t' macs2 callpeak --broad --broad-cutoff 0.01 -f BAMPE -g 2652783500 --outdir peaks -n {2}_TOTAL -t bam/{1}.bam -c bam/{3}.bam :::: <(awk '!/Input|accession/' PRJNA480803.tsv) :::: <(awk '/Input/ {print; print; print; print; print; print; print;}' PRJNA480803.tsv) > callpeakTOTAL.log & # 124503
nohup parallel -j 12 --link --colsep '\t' macs2 callpeak --broad --broad-cutoff 0.01 -f BAMPE -g 2652783500 --outdir peaks -n {2}_C57 -t allelic/{1}.genome1.sorted.bam -c allelic/{3}.genome1.sorted.bam :::: <(awk '!/Input|accession/' PRJNA480803.tsv) :::: <(awk '/Input/ {print; print; print; print; print; print; print;}' PRJNA480803.tsv) > callpeakC57.log & # 125300
nohup parallel -j 12 --link --colsep '\t' macs2 callpeak --broad --broad-cutoff 0.01 -f BAMPE -g 2652783500 --outdir peaks -n {2}_CAST -t allelic/{1}.genome2.sorted.bam -c allelic/{3}.genome2.sorted.bam :::: <(awk '!/Input|accession/' PRJNA480803.tsv) :::: <(awk '/Input/ {print; print; print; print; print; print; print;}' PRJNA480803.tsv) > callpeakCAST.log & # 126097

# filter blacklisted peaks
mkdir peaks/filtered

parallel bedtools subtract -A -a {} -b mm10-blacklist.v2.ENSEMBL.bed.gz ">" peaks/filtered/{/.}_filtered.broadPeak ::: peaks/*.broadPeak

# consensus peaks
mkdir peaks/consensus

parallel --link --rpl '{%}  s:.*/::; s:Rep.*::;' bedtools intersect -f 0.25 -r -a  {1} -b {2} ">" peaks/consensus/{1%}TOTAL_peaks_consensus_filtered.broadPeak ::: peaks/filtered/*Rep1_TOTAL_peaks_filtered.broadPeak ::: peaks/filtered/*Rep2_TOTAL_peaks_filtered.broadPeak
parallel --link --rpl '{%}  s:.*/::; s:Rep.*::;' bedtools intersect -f 0.25 -r -a  {1} -b {2} ">" peaks/consensus/{1%}C57_peaks_consensus_filtered.broadPeak ::: peaks/filtered/*Rep1_C57_peaks_filtered.broadPeak ::: peaks/filtered/*Rep2_C57_peaks_filtered.broadPeak
parallel --link --rpl '{%}  s:.*/::; s:Rep.*::;' bedtools intersect -f 0.25 -r -a  {1} -b {2} ">" peaks/consensus/{1%}CAST_peaks_consensus_filtered.broadPeak ::: peaks/filtered/*Rep1_CAST_peaks_filtered.broadPeak ::: peaks/filtered/*Rep2_CAST_peaks_filtered.broadPeak

## reference peak ratios
mkdir allelic_quantification

parallel -j 12 \
multiBigwigSummary BED-file \
-p 8 \
-b bigwig/*_{1}_*_{2}_*C57.bw bigwig/*_{1}_*_{2}_*CAST.bw bigwig/*Input*_{2}_*C57.bw bigwig/*Input*_{2}_*CAST.bw \
--labels {1}_{2}_C57_rep1 {1}_{2}_C57_rep2 {1}_{2}_CAST_rep1 {1}_{2}_CAST_rep2 Input_{2}_C57_rep1 Input_{2}_C57_rep2 Input_{2}_CAST_rep1 Input_{2}_CAST_rep2 \
--outRawCounts allelic_quantification/{1}_{2}_TOTAL.tab \
-o allelic_quantification/{1}_{2}_TOTAL.npz \
--BED peaks/consensus/ChIPseq_{1}_TX1072_WT_{2}_Dox_TOTAL_peaks_consensus_filtered.broadPeak \
::: H2AK119Ub H3K27ac H3K27me3 H3K4me1 H3K4me3 H3K9ac H4ac ::: 0h 4h 8h 12h 24h

## cleanup
rm fastq/*
rm clean/*
