#!/usr/bin/env python
import os 
import pandas as pd
import glob
from collections import Counter


for asm in ["canu", "miniasm", "wtdbg"]:
	dfs = [] 
	tables = glob.glob("*/{}.sda.table.tsv".format(asm))
	for df in tables:
		#print(df)
		if(os.path.getsize(df) > 0):
			dfs.append( pd.read_csv(df, sep = "\t" ))

	merged = pd.concat(dfs, sort = True)
	#print(len(tables))
	print(asm, len(tables))
	print(Counter(merged.Status))

	merged.to_csv("{}.localAssemblyStats.tsv".format(asm), sep="\t", index=False)

