#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse 


parser = argparse.ArgumentParser()
parser.add_argument("--nucfreq", help="assembly.consensus.nucfreq", default="assembly.consensus.nucfreq")
parser.add_argument("--plot", help="name of output plot", default="Coverage.png")
parser.add_argument("--psvsites", help="mi.gml.sites", default=None)
parser.add_argument("--automin", help="autoMinCoverage", default="autoMinCoverage")
parser.add_argument("--automax", help="autoMaxCoverage", default="autoMaxCoverage")
args = parser.parse_args()

nucfreq=args.nucfreq
autoMin=args.automin
autoMax=args.automax


colnames = ["contig", "pos", "A", "C", "G", "T", "deletion", "insertion"]

f = open(nucfreq)
first  = []
second = []
third = []
truepos= []
for line in f:
    line = line.split()
    truepos.append(int(line[1]))
    bases = []
    for basepair in line[2:6]:
        bases.append(int(basepair))
    bases = sorted(bases, reverse=True)
    first.append(bases[0])
    second.append(bases[1])
    third.append(bases[3])
pos = np.array( range(0,len(second)) ) 
first = np.array(first)
second = np.array(second)
third = np.array(third)
truepos = np.array(truepos)
#truepos = pos

#plt.rc('font', family='serif')
matplotlib.rcParams.update({'font.size': 18})
fig, ax = plt.subplots( figsize=(16,9) )
prime, = plt.plot(truepos, first, 'o', color="black", markeredgewidth=0.0, markersize=2, label = "most frequent base pair")
#plt.gca().set_ylim(top=300)
sec, = plt.plot(truepos, second,'o', color="red",   markeredgewidth=0.0, markersize=2, label = "second most frequent base pair")
#tri, = plt.plot(truepos, third,'o', color="green",   markeredgewidth=0.0, markersize=1, label = "forth most frequent base pair")

plotnum = 2
if(args.psvsites is not None):
		cuts = {}
		for idx, line in enumerate(open(args.psvsites).readlines()):
			try:	
				vals = line.strip().split()
				cuts[idx] = list(map(int, vals))
				#make plot
				plotnum += 1
				x = cuts[idx]
				idxs = (np.isin(truepos, x))
				y = second[idxs]
				plt.plot(x,y, label="group:{}".format(idx), alpha = 0.5)
			except Exception as e:
				print("Skipping because error: {}".format(e), file=sys.stderr)
				continue


ax.set_xlabel('Collapse Position (bp)')
ax.set_ylabel('Sequence Read Depth')

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


lgnd = plt.legend(loc="upper right", ncol = int(np.ceil(plotnum/2.0)), fontsize = 12)
for handle in lgnd.legendHandles:
	handle._sizes = ([50.0])

plt.savefig(args.plot, dip=900)




