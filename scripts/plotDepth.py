#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument("depth", help="tab seperated file, contig   pos   depth" )
parser.add_argument("png",help="output png, the name must have the .png ext",default=None)
parser.add_argument("--m5",help="an m5 output from blasr that indicates a mapping of some sequences onto the depth profile", default=None)
parser.add_argument("--sam", help="a sam or bam file of sequences to plot against the depth profile(contigs against a depth profile)", default=None)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()
DEBUG=args.d

import glob
import os
import sys
import re
import itertools 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
import pysam 
import runCmd
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


colNames = ["contig", "pos", "depth"]
df = pd.read_csv(args.depth, sep="\t", names=colNames)

#df.plot(x="pos", y="depth")
grouped = df.groupby('contig')
numGroups = len(grouped)
if(numGroups == 0):
    print("No reads map to this contig.")
    exit(0)
print(grouped)

plt.rc('font', family='serif')
fig, axsList = plt.subplots(numGroups, figsize=(16,9))
# convert to array
if(type(axsList) != type(np.array([0]))):
    axsList = np.array([axsList])
# convert to dictionary 
contigs = df.contig.unique()
axs={}
for idx, ax in enumerate(axsList):
	# associate a ax with a specific contig
	axs[ contigs[idx] ] = ax
print(axs)

# will plot mutliple contigs if appropriate
counter = 0
for key, group in grouped:
    print(key)
    ax = axs[key]
    #group.plot(ax=axs[counter], kind='scatter', x='pos', y='depth', color="black", label=key)
    ax.plot(group['pos'], group['depth'], 'o', color="black",  markeredgewidth=0.0, markersize=3, label=key)
    ax.set_xlim(xmin = 0, xmax = max(group["pos"]))
    ax.set_ylim(ymin = 0, ymax = max(group["depth"]))
    counter += 1
    
    ax.set_xlabel('BP Position, ' + key)
    ax.set_ylabel('Depth')

    ylabels = [format(label, ',.0f') for label in ax.get_yticks()]
    xlabels = [format(label, ',.0f') for label in ax.get_xticks()]
    ax.set_yticklabels(ylabels)
    ax.set_xticklabels(xlabels)

    # Hide the right and top spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
 
max_y = max(df["depth"])
min_y = min(df["depth"])


# add ref sequences if they are there
if(args.m5 is not None):
    contigMap = {}
    # each row will be x, y, width, height  
    pos = []
    lines = open(args.m5).readlines()
    num_aln = len(lines)
    for idx, line in enumerate(lines):
        tokens = line.split()
        pos.append([0,0,0,0])
        pos[idx][0] = int(tokens[7])
        pos[idx][1] = idx * max_y/(num_aln*3)
        pos[idx][2] = int(tokens[8]) - pos[idx][0] 
        pos[idx][3] = max_y/(num_aln*3)
        # template sequence, reads sequence name
        name = re.sub("group.","",tokens[0])
        name = re.sub("_quiver","",name)
        contigMap[idx] = {"target":tokens[5], "name":name}
        print(pos[idx]) 
    
    # add a rectangle for the alignemnt on each appropraite subplot
    counter = 0
    for key, group in grouped:
        ax = axs[key]    
        for idx in contigMap:
            contig = contigMap[idx]["target"]
            label = contigMap[idx]["name"]
            if(contig == key):
                ax.add_patch( patches.Rectangle((pos[idx][0], pos[idx][1]),   # (x,y)
                    pos[idx][2], pos[idx][3], alpha=0.75, linewidth=3, color = "darkred")) # width, height 
                # add an annotation
                ax.text(pos[idx][0], pos[idx][1], label )
        counter += 1


def addPatchToAx(ax, x, y, width, height, operation):
	col = {2:"gray", 7:"darkgreen", 8:"darkred"}
	ax.add_patch( patches.Rectangle((x,y), width, height, alpha=0.75, linewidth=0, color = col[operation] )) 


def squaresFromCigar(ax, cigar, start, y, height):
	x = start 
	for pair in cigar:
		operation = pair[0]
		width = pair[1]
		if(operation in [2,7,8]):
			addPatchToAx(ax, x, y, width, height, operation)
			x += width
	return(x)

# add contigs from sam file
if(args.sam is not None):
	import pysam
	samfile = pysam.AlignmentFile(args.sam)
	# a map from reference to names of queriers that mathc  
	ref={}
	reads={}
	for read in samfile.fetch():
		if(read.is_unmapped):
			continue 
		# add reads to a dictionary by query
		reads[read.query_name] = read
		# add queries to a dictionary by contig they map to
		if(read.reference_name not in ref):
			ref[read.reference_name] = []
		ref[read.reference_name].append(read.query_name)
		

	for contig in ref:
		ax = axs[contig]
		queries = ref[contig]
		group = grouped.get_group(contig)
		max_y = max(group["depth"])
		height = max_y/(len(queries)*3)
		# add rectangles for each alignment / cigar string tupple 
		for counter, query in enumerate(queries):
			y = counter * max_y/(len(queries)*3)
			read = reads[query]
			end = squaresFromCigar(ax, read.cigar, read.reference_start, y , height)
			#assert(end == read.reference_end)
			label="{}:{}-{};{};{}".format(query, read.query_alignment_start, read.query_alignment_end, 
					read.is_reverse, read.infer_query_length())
			ax.text(read.reference_start, y, label)



# this shoudl fix overlapping lables 
plt.tight_layout()
plt.show()
plt.savefig(args.png)

