#!/usr/bin/python3
import sys
import re
import pandas as pd

in_file = sys.argv[1]
out_file = re.sub(".pkl",".tsv",in_file)

df = pd.read_pickle(in_file)
df[['k_on','k_off','k_syn']] = pd.DataFrame(df.iloc[:,0].values.tolist(), index= df.index)
df = df.drop(columns=[0])

df.to_csv(out_file,sep="\t")