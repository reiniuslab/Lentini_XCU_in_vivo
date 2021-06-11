#!/bin/bash

SNPsplit_genome_preparation --vcf_file mgp.v5.merged.snps_all.dbSNP142.vcf.gz --dual_hybrid --strain 129S1_SvImJ --strain2 CAST_EiJ --reference_genome ref

cd 129S1_SvImJ_CAST_EiJ_dual_hybrid.based_on_GRCm38_N-masked/

## build bowtie index
bowtie2-build chr1.N-masked.fa,chr2.N-masked.fa,chr3.N-masked.fa,chr4.N-masked.fa,chr5.N-masked.fa,chr6.N-masked.fa,chr7.N-masked.fa,chr8.N-masked.fa,chr9.N-masked.fa,chr10.N-masked.fa,chr11.N-masked.fa,chr12.N-masked.fa,chr13.N-masked.fa,chr14.N-masked.fa,chr15.N-masked.fa,chr16.N-masked.fa,chr17.N-masked.fa,chr18.N-masked.fa,chr19.N-masked.fa,chrMT.N-masked.fa,chrX.N-masked.fa,chrY.N-masked.fa bowtie2/GRCm38.mgp.v5.129xCAST_Nmasked/index

## digest genome
hicup_digester --genome GRCm38_129xCAST --re1 ^GATC,MboI chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,MT,X,Y}.N-masked.fa

## get raw data
cd ../HiC/

mkdir fastq
cd fastq

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/005/SRR7473265/SRR7473265_{1..2}.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/006/SRR7473266/SRR7473266_{1..2}.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/007/SRR7473267/SRR7473267_{1..2}.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/008/SRR7473268/SRR7473268_{1..2}.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/005/SRR7473275/SRR7473275_{1..2}.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/006/SRR7473276/SRR7473276_{1..2}.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/007/SRR7473277/SRR7473277_{1..2}.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR747/008/SRR7473278/SRR7473278_{1..2}.fastq.gz

cd ..

## run HiCUP
mkdir hicup

nohup parallel -j 4 --link \
hicup \
--threads 16 \
--bowtie2 /usr/bin/bowtie2 \
--index bowtie2/GRCm38.mgp.v5.129xCAST_Nmasked/index \
--digest 129S1_SvImJ_CAST_EiJ_dual_hybrid.based_on_GRCm38_N-masked/Digest_GRCm38_129xCAST_MboI_None_16-50-59_16-12-2020.txt \
--shortest 50 \
--longest 700 \
--zip \
--outdir hicup \
{1} {2} \
::: fastq/*_1.fastq.gz ::: fastq/*_2.fastq.gz &

## split to alleles
mkdir snpsplit

nohup parallel -j 4 \
SNPsplit \
--hic \
-o snpsplit \
--snp_file all_CAST_EiJ_SNPs_129S1_SvImJ_reference.based_on_GRCm38.txt \
::: hicup/*_1_2.hicup.bam &

## convert to juicer format
HiCUP-0.8.0/Conversion/hicup2juicer --zip snpsplit/*.G1_G1.bam snpsplit/*.G1_UA.bam snpsplit/*.G2_G2.bam snpsplit/*.G2_UA.bam

## merge per genome
mkdir prejuicer

function merge_prejuicer(){
	# merge and sort by chromosome
	zcat ${1}_1_2.hicup.G1_G1.bam.prejuicer.gz ${1}_1_2.hicup.G1_UA.bam.prejuicer.gz ${2}_1_2.hicup.G1_G1.bam.prejuicer.gz ${2}_1_2.hicup.G1_UA.bam.prejuicer.gz | sort -k3,3d -k7,7d | gzip > ${3}_G1.txt.gz
	zcat ${1}_1_2.hicup.G2_G2.bam.prejuicer.gz ${1}_1_2.hicup.G2_UA.bam.prejuicer.gz ${2}_1_2.hicup.G2_G2.bam.prejuicer.gz ${2}_1_2.hicup.G2_UA.bam.prejuicer.gz | sort -k3,3d -k7,7d | gzip > ${3}_G2.txt.gz
}
export -f merge_prejuicer

merge_prejuicer snpsplit/SRR7473266 snpsplit/SRR7473275 prejuicer/D0
merge_prejuicer snpsplit/SRR7473267 snpsplit/SRR7473276 prejuicer/D3
merge_prejuicer snpsplit/SRR7473268 snpsplit/SRR7473277 prejuicer/D7
merge_prejuicer snpsplit/SRR7473265 snpsplit/SRR7473278 prejuicer/D10

## generate genome information for juicer
cat 129S1_SvImJ_CAST_EiJ_dual_hybrid.based_on_GRCm38_N-masked/chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,MT,X,Y}.N-masked.fa > 129xCAST.fa

./juicer/misc/generate_site_positions.py MboI mm10 129xCAST.fa
awk 'BEGIN{OFS="\t"}{print $1, $NF}' mm10_MboI.txt > mm10.chrom.sizes

## calculate contact maps
nohup parallel --plus -j 8 \
java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar pre \
-d \
-f mm10_MboI.txt \
-q 10 \
-r 500000,1000000,2500000 \
{} \
juicer/{/..}.hic \
mm10.chrom.sizes \
::: prejuicer/*.txt.gz &

## dump maps
mkdir dump
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar dump observed KR {} X X BP 1000000 dump/{/.}_1mb_chrX.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar dump observed KR {} 5 5 BP 1000000 dump/{/.}_1mb_chr5.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar dump observed KR {} 7 7 BP 1000000 dump/{/.}_1mb_chr7.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar dump observed KR {} 9 9 BP 1000000 dump/{/.}_1mb_chr9.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar dump observed KR {} 13 13 BP 1000000 dump/{/.}_1mb_chr13.txt ::: juicer/*.hic

# hi-res
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar dump observed KR {} X X BP 500000 dump/{/.}_500kb_chrX.txt ::: juicer/*.hic

## calculate pearson correlations
mkdir pearsons
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar pearsons KR {} X BP 1000000 pearsons/{/.}_1mb_chrX.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar pearsons KR {} 5 BP 1000000 pearsons/{/.}_1mb_chr5.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar pearsons KR {} 7 BP 1000000 pearsons/{/.}_1mb_chr7.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar pearsons KR {} 9 BP 1000000 pearsons/{/.}_1mb_chr9.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar pearsons KR {} 13 BP 1000000 pearsons/{/.}_1mb_chr13.txt ::: juicer/*.hic

##
mkdir eigen
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar eigenvector KR {} X BP 1000000 eigen/{/.}_1mb_chrX.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar eigenvector KR {} 5 BP 1000000 eigen/{/.}_1mb_chr5.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar eigenvector KR {} 7 BP 1000000 eigen/{/.}_1mb_chr7.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar eigenvector KR {} 9 BP 1000000 eigen/{/.}_1mb_chr9.txt ::: juicer/*.hic
parallel -j 4 java -Xmx8g -jar ./juicer/juicer_tools_1.22.01.jar eigenvector KR {} 13 BP 1000000 eigen/{/.}_1mb_chr13.txt ::: juicer/*.hic
