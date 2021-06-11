#!/bin/bash
mkdir -m 777 bam bam/bulk fastq fastq/clean/ fastq/clean/qc_reports

ulimit -n 4000 # increase open file limit

## connect read name with barcode info
mkdir fastq/tagged/
parallel --plus -j 32  --rpl '{%} s{^.*/|_R\d.*$}{}g;' "(zcat {} | sed -e 's/@/@{%}:/' | gzip > fastq/tagged/{/..}.tagged.fastq.gz)" ::: fastq/HCB*.fastq.gz

## quality & adapter trim
parallel --plus --link -j 12 --rpl '{%} s{^.*/|_R\d.*$}{}g;' fastp --in1 {1} --in2 {2} --out1 fastq/clean/{1/} --out2 fastq/clean/{2/} --json fastq/clean/qc_reports/{1%}.fastp.json --html fastq/clean/qc_reports/{1%}.fastp.html -R {1%} ::: fastq/tagged/HCB*R1*.tagged.fastq.gz ::: fastq/tagged/HCB*R2*.tagged.fastq.gz

multiqc fastq/clean/qc_reports/ -f -n qc_fastq -m fastp

## align
align_atac(){
	# $1 read1
	# $2 read2
	# $3 output
	(bowtie2 \
	-x bowtie2/GRCm38.mgp.v5.CAST_EiJ_Nmasked/index \
	--mm \
	--very-sensitive -N 1 -X 2000 -k 10 \
	-p 8 \
	-1 $1 \
	-2 $2 \
	| bamsort sortthreads=8 blockmb=6144 inputformat=sam M=fastq/clean/qc_reports/${3}.bamsort.log \
	| bammarkduplicates markthreads=8 M=fastq/clean/qc_reports/${3}.bammarkduplicates.log index=1 indexfilename=bam/${3}.bam.bai O=bam/${3}.bam \
	)2>fastq/clean/qc_reports/${3}.bowtie2.log
}
export -f align_atac

parallel -j 8 --link --rpl '{%} s{^.*/|_R\d.*$}{}g;'  align_atac {1} {2} {1%} ::: fastq/clean/HCB*R1*.fastq.gz ::: fastq/clean/HCB*R2*.fastq.gz

multiqc fastq/clean/qc_reports -f -n qc_mapping -m bowtie2 -m biobambam2

## merge bam files
samtools merge -@ 63 -O bam bam/merge/scATAC_all.sorted.bam bam/*.bam

## allelic expression
# split data
sambamba -q sort -N -m 100GB -t 64 -o bam/merge/scATAC_all.bam bam/merge/scATAC_all.sorted.bam
SNPsplit --paired  --no_sort -o bam/merge/ --snp_file all_SNPs_CAST_EiJ_GRCm38.txt.gz bam/merge/scATAC_all.bam

# sort
sambamba -q sort -m 20GB -t 64 bam/merge/scATAC_all.genome1.bam
sambamba -q sort -m 20GB -t 64 bam/merge/scATAC_all.genome2.bam
