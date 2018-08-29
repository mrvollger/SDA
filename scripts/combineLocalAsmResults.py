#!/usr/bin/env python
import os 
import pandas as pd
import glob
from collections import Counter

dfs = [] 
tables = glob.glob("*/*.table.tsv")
for df in tables:
	#print(df)
	dfs.append( pd.read_csv(df, sep = "\t" ))

merged = pd.concat(dfs)
print(len(tables))
print(Counter(merged.Status))

merged.to_csv("localAssemblyStats.tsv",sep="\t", index=False)

