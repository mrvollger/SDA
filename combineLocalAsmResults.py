#!/usr/bin/env python
import os 
import pandas as pd
import glob
dfs = []
for df in glob.glob("*/abp.table.tsv"):
	print(df)
	dfs.append( pd.read_csv(df, sep = "\t" ))

merged = pd.concat(dfs)
merged.to_csv("localAssemblyStats.tsv",sep="\t", index=False)

