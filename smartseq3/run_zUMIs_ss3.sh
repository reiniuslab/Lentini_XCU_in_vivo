#!/bin/bash

zUMIs/zUMIs-master.sh -y mESC_NEW.yaml
zUMIs/zUMIs-master.sh -y mESC_clone5_diff.yaml

# from https://github.com/sandberg-lab/Smart-seq3/tree/master/allele_level_expression
Rscript get_variant_overlap_CAST.R --yaml mESC_NEW.yaml --vcf /huxley.hdd/projects/smartseq3/CAST.SNPs.superset.ENSEMBL.vcf.gz
Rscript get_variant_overlap_CAST.R --yaml mESC_clone5_diff.yaml --vcf /huxley.hdd/projects/smartseq3/CAST.SNPs.superset.ENSEMBL.vcf.gz
