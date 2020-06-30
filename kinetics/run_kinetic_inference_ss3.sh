#!/bin/bash

## tpm
parallel -j 1 ./txburstML.py ::: input/tpm_ss3/*.csv
mv *.pkl output/tpm_ss3
parallel -j 10 ./pkl_to_tsv.py ::: output/tpm_ss3/*.pkl

./split_kinetics.R output/tpm_ss3/ output/tpm_ss3/split/tpm_ss3 2,3,4,1

## umi
parallel -j 1 ./txburstML.py ::: input/umi_ss3/*.csv
mv *.pkl output/umi_ss3
parallel -j 10 ./pkl_to_tsv.py ::: output/umi_ss3/*.pkl

./split_kinetics.R output/umi_ss3/ output/umi_ss3/split/umi_ss3 3,4,5,2
