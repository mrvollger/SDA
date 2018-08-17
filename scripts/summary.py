#!/usr/bin/env python
#import ABPUtils as abp
#import networkx as nx
#import numpy as np
from LocalAssembly3 import LocalAssembly
import pandas as pd
import os
import glob
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("--assembler", help="the name of the assembler used, canu, wtdbg, miniasm")
parser.add_argument("--summary", help="output location for the summary")
args = parser.parse_args()


# generate a summary of the assmebly 

cwd = os.getcwd()
LA = LocalAssembly(cwd, args.assembler)
df = LA.all
#data = glob.glob("abp.table.tsv")[0]
#df = pd.read_csv(data, sep = "\t")
row = df.iloc[0]
#print(row)
region = row[ "collapse" ]
copies = row[ "copiesInRef" ]
numAsm = row["numR"] + row["numPR"]
numCC =  row["numOfCCgroups"]
refRegions = row["refRegions"]
mysplit = refRegions.split(";")
refRegions = ""
for idx, reg in enumerate(mysplit):
    refRegions += reg + "  "
    if( (idx + 1) % 3 == 0 ):
        refRegions += "\n"

notHeader =[ "numR","numF", "numPR", "numMA", "collapse", "copiesInRef", "numOfCCgroups",
        "refRegions", "collapseLen", "totalReads", "totalPSVs", "numOfAssemblies", "totalPSVs",
		"bestStart", "bestEnd", "bestChr", "aveRefLength"]

header = ["CC_ID", "Status", "perID_by_matches", "bestMatch", "Length", "numPSVs", "query_name"]


# add in the psv matrix
nums = ( range(1000))
for col in list(df):
	if(col in nums):
		header.append(col)

print(header)

#df = df[headerList]
#df.drop(notHeader, axis=1, inplace=True)
#df = df[df["Status"] != "Failed"]
df = df[header]
df = df.sort_values(["CC_ID", "query_name"])
pd.set_option('display.width', 200)


#df=df.rename(columns = {'BestRegionInTheHumanReference':'BestInRef'})
#df=df.rename(columns = {"num_of_sam_reads":'reads'})
#df=df.rename(columns = {"averageLength":'aveLen'})
#df=df.rename(columns = {"averagePerID":'avePerID'})
#df=df.rename(columns = {"bestLength":'bestLen'})

text = '''region: {}
copies: {}
number of assemblies: {}
number of CC groups: {}
Regions in reference:
{}
Group Stats:
{}
'''.format(region, copies, numAsm, numCC, refRegions, df.to_string(index=False))

#print(text)
outfile = open(args.summary, "w+")
outfile.write(text)
outfile.close()


