#!/bin/bash

mkdir kinetics kinetics/output/ kinetics/output/rpkm kinetics/output/tpm

## tpm
parallel -j 1 ./txburstML.py ::: kinetics/input/tpm/*.csv
mv *.pkl kinetics/output/tpm
parallel -j 10 ./pkl_to_tsv.py ::: kinetics/output/tpm/*.pkl

./split_kinetics.R kinetics/output/tpm/ kinetics/output/tpm/split/postimplantation 2,3,4,1

#python3 transform_and_regress_out_burst_frequency.py -i kinetics/output/tpm/split/postimplantation_bf.tsv kinetics/output/tpm/split/postimplantation_bs.tsv -o kinetics/output/tpm/split/postimplantation_bf_regressed.tsv kinetics/output/tpm/split/postimplantation_bs_regressed.tsv --revert_transform
#python3 transform_and_regress_out_burst_frequency.py -i kinetics/output/tpm/split/postimplantation_bf.tsv kinetics/output/tpm/split/postimplantation_bs.tsv -o kinetics/output/tpm/split/postimplantation_bf_regressed_log.tsv kinetics/output/tpm/split/postimplantation_bs_regressed_log.tsv

# ./txburstPL.py --file split/rpkm.c57.M.1.XaXi.csv --MLFile rpkm.c57.M.1.XaXi_ML.pkl
#./txburstPL.py --file split/rpkm.cast.M.1.XaXi.csv --MLFile rpkm.cast.M.1.XaXi_ML.pkl
# ./txburstTEST.py --file1 split/rpkm.c57.M.1.XaXi.csv --file2 split/rpkm.cast.M.1.XaXi.csv --ML1 rpkm.c57.M.1.XaXi_ML.pkl --ML2 rpkm.cast.M.1.XaXi_ML.pkl