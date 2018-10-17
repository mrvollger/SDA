#!/usr/bin/env python
import os 
import pandas as pd
import glob
from collections import Counter

dfs = [] 
tables = glob.glob("*/canu.sda.table.tsv")
for df in tables:
	#print(df)
	if(os.path.getsize(df) > 0):
		dfs.append( pd.read_csv(df, sep = "\t" ))

merged = pd.concat(dfs)
print(len(tables))
print(Counter(merged.Status))

merged.to_csv("localAssemblyStats.tsv",sep="\t", index=False)

