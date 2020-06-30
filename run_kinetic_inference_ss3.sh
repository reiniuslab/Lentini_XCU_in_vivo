#!/bin/bash

## tpm
parallel -j 1 ./txburstML.py ::: kinetics/input/tpm_ss3/*.csv
mv *.pkl kinetics/output/tpm_ss3
parallel -j 10 ./pkl_to_tsv.py ::: kinetics/output/tpm_ss3/*.pkl

./split_kinetics.R kinetics/output/tpm_ss3/ kinetics/output/tpm_ss3/split/tpm_ss3 2,3,4,1

## umi
parallel -j 1 ./txburstML.py ::: kinetics/input/umi_ss3/*.csv
mv *.pkl kinetics/output/umi_ss3
parallel -j 10 ./pkl_to_tsv.py ::: kinetics/output/umi_ss3/*.pkl

./split_kinetics.R kinetics/output/umi_ss3/ kinetics/output/umi_ss3/split/umi_ss3 3,4,5,2