#!/usr/bin/env python
import ABPUtils as abp
import networkx as nx
import numpy as np
import LocalAssembly
import pandas as pd
import os

cwd = os.getcwd()
print(cwd)


# old commands ot make a psvgraph, I do not use this anymore 
'''
myfile = "mi.cuts.gml"
out = "out.png"
g = abp.ReadGraph(myfile)

# get the labels
color = nx.get_node_attributes(g, "color")
degree = g.degree()
labels={}
for key in color:
    labels[key] = str(color[key]) + ":" + str(degree[key])


# set the node positions of a few values
alreadySeen = set()
fixed = []
for key in color:
    c = color[key]
    if(c not in alreadySeen):
        fixed.append(key)
        alreadySeen.add(c)

x = int(np.sqrt(len(fixed))) + 1
ind=list(np.ndindex(x,x))
clusterPos={} # color and then pos
for idx, c in enumerate(alreadySeen):
    clusterPos[c] = ind[idx]

print(clusterPos)

for n in g.node.keys():
    c = g.node[n]["color"]
    g.node[n]["pos"]=clusterPos[c]

pos = nx.get_node_attributes(g, "pos")
'''


# generate a summary of the assmebly 
# get text to add to the bottom of the graph. 
LA = LocalAssembly.LocalAssembly(cwd)
df = LA.asPD()
row = df.iloc[0]
region = row[ "region_in_falcon" ]
copies = row[ "copies_in_reference" ]
numAsm = row["number_of_psv_assemblies"]
numCC =  row["number_of_CC_groups"]
refRegions = row["refRegions"]
mysplit = refRegions.split(";")
refRegions = ""
for idx, reg in enumerate(mysplit):
    refRegions += reg + "  "
    if( (idx + 1) % 3 == 0 ):
        refRegions += "\n"

headerList =[ "CC_ID", "Status","numOfPSVs","length", 
                "num_of_sam_reads",
                "averagePerID", "averageLength", 
                "bestPerID", "bestLength", "BestRegionInTheHumanReference"]
df = df[headerList] 
#df = df[df["Status"] != "Failed"]
df = df.sort_values(["Status", "CC_ID"])
pd.set_option('display.width', 200)

text = '''region: {}
copies: {}
number of assemblies: {}
number of CC groups: {}
Regions in reference:
{}
Group Stats:
{}
'''.format(region, copies, numAsm, numCC, refRegions, df.to_string(index=False))

print(text)
outfile = open("summary.txt", "w+")
outfile.write(text)
outfile.close()


#abp.DrawGraph(g, out,  labels=labels, pos=pos, addText = text)




