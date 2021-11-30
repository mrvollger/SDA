#!/usr/bin/env python
import argparse
import sys
import os

parser = argparse.ArgumentParser(description="")
parser.add_argument("infiles", help="collapse file(s) from SDA (bed).", nargs="+")
parser.add_argument("-c", "--coverage", help="expected wgs coverage", type=float)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
import pandas as pd

def count_collapse(df, fname = None):
	x = sum(df.end - df.start)
	y = sum( (df.end - df.start) * df["mean"]/args.coverage )
	#print("{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}".format(fname, x,y, x/10**6, y/10**6))
	return(x,y)

colnames ="File\tCollapsed_bp\tExpanded_bp\tCollapses_Mbp\tExpanded_Mbp".split("\t")

data = []

for f in args.infiles:
	df = pd.read_csv(f, sep="\t", names=["contig", "start", "end", "mean", "median", "cr_bp", "len"])
	fname=os.path.basename(f)
	x,y = count_collapse(df, fname)
	data.append( (fname, x, y, x/10**6, y/10**6 ) )
	#count_collapse(df[df.contig == "chrX_v0.7"])	
	

data = pd.DataFrame(data, columns=colnames) 
data.sort_values(by=["Expanded_bp", "Collapsed_bp"], inplace=True)
print(data, file=sys.stderr)


#pd.options.display.float_format = '{:.2f}'.format
#print(data)
data.to_csv(sys.stdout, sep="\t",index=False)

