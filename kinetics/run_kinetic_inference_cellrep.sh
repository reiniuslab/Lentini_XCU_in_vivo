#!/bin/bash

## tpm
parallel -j 1 ./txburstML.py ::: input/tpm/*.csv
mv *.pkl kinetics/output/tpm
parallel -j 10 ./pkl_to_tsv.py ::: output/tpm/*.pkl

./split_kinetics.R output/tpm/ output/tpm/split/postimplantation 2,3,4,1
