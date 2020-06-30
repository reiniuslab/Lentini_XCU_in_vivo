#!/bin/bash

#############
## PRJNA130405 ##
#############

mkdir PRJNA130405 PRJNA130405/RNAseq PRJNA130405/RNAseq/quant

cd PRJNA130405
wget https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-23943/E-GEOD-23943.sdrf.txt

curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA130405&result=read_run&fields=study_accession,secondary_study_accession,experiment_accession,run_accession,fastq_ftp,sample_title&format=tsv&download=true" > PRJNA130405.txt

awk -F '\t' '/RNA-Seq/{print $40}' E-GEOD-23943.sdrf.txt | sort | uniq | wget -i - -P RNAseq

find RNAseq/ -maxdepth 1 -not -type d | parallel --rpl '{%} s{^.*/|\..*}{}g;' -j 5 salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r {} -o RNAseq/quant/{%}

cd ..

#############
## PRJNA354946 ##
#############

mkdir PRJNA354946 PRJNA354946/RNAseq PRJNA354946/RNAseq/quant

cd PRJNA354946
curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA354946&result=read_run&fields=study_accession,secondary_study_accession,experiment_accession,run_accession,fastq_ftp,sample_title&format=tsv&download=true" > PRJNA354946.txt

awk '{print $6}' GSE90516.txt | wget -i - -P RNAseq

# merge runs
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054337.fastq.gz RNAseq/SRR5054338.fastq.gz -o RNAseq/quant/SRX2375307
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054339.fastq.gz RNAseq/SRR5054340.fastq.gz -o RNAseq/quant/SRX2375308
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054341.fastq.gz RNAseq/SRR5054342.fastq.gz -o RNAseq/quant/SRX2375309
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054343.fastq.gz RNAseq/SRR5054344.fastq.gz -o RNAseq/quant/SRX2375310
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054345.fastq.gz RNAseq/SRR5054346.fastq.gz -o RNAseq/quant/SRX2375311
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054347.fastq.gz RNAseq/SRR5054348.fastq.gz -o RNAseq/quant/SRX2375312
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054349.fastq.gz RNAseq/SRR5054350.fastq.gz -o RNAseq/quant/SRX2375313
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054351.fastq.gz RNAseq/SRR5054352.fastq.gz -o RNAseq/quant/SRX2375314
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054353.fastq.gz RNAseq/SRR5054354.fastq.gz -o RNAseq/quant/SRX2375315
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054355.fastq.gz RNAseq/SRR5054356.fastq.gz -o RNAseq/quant/SRX2375316
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054357.fastq.gz RNAseq/SRR5054358.fastq.gz -o RNAseq/quant/SRX2375317
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054359.fastq.gz RNAseq/SRR5054360.fastq.gz -o RNAseq/quant/SRX2375318
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054361.fastq.gz RNAseq/SRR5054362.fastq.gz -o RNAseq/quant/SRX2375319
salmon quant --validateMappings -l A -p 12 -i /samus.ssd/genome_index/salmon/GRCm38.p6.gencode.vM22 -r RNAseq/SRR5054363.fastq.gz RNAseq/SRR5054364.fastq.gz -o RNAseq/quant/SRX2375320

cd ..
