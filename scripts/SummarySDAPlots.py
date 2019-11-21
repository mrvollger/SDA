#!/usr/bin/env python
import argparse
import sys

parser = argparse.ArgumentParser(description="")
parser.add_argument("summary", help=" input sda summary file")
parser.add_argument("--dir",  help="sda direcotry file", default=None)
parser.add_argument("--prefix",  help="sda direcotry file", default=None)
parser.add_argument("--length",  help="plot of sda assembly lengths", default=None)
parser.add_argument("--collapse",  help="length of collapse vs length of sda assemblies", default=None)
parser.add_argument("--pair",  help="pairwise plot of all numeric columns in the summary", default=None)
parser.add_argument("--bar",  help="bar plot of total mbp of collapse and total mbp of sda", default=None)
parser.add_argument('-d', action="store_true", default=False)
args = parser.parse_args()

import os
import pandas as pd
import numpy as np
from scipy import stats
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [16, 9]
plt.rc("axes.spines", top=False, right=False)

import seaborn as sns
sns.set_context("poster")


def qfilter(array, bot=.05, top=.95, log = False):
    top = array.quantile(top)
    bot = array.quantile(bot)
    rtn =  array[ (array>= bot) & (array <= top) ] 
    if(log):
        rtn = np.log10(rtn)
    return(rtn)

def corrfunc(x, y, **kws):
	r, _ = stats.pearsonr(x, y)
	ax = plt.gca()
	ax.annotate("r = {:.2f}".format(r), xy=(.1, .9), xycoords=ax.transAxes)

# In[185]:


PRE = args.prefix
DIR = args.dir

df = pd.read_csv(args.summary, sep="\t")


if("collapseLength" not in list(df)):
    import glob
    from Bio import SeqIO
    collen = []
    for prefix in df["prefix"]:
        ref = f"{DIR}/{PRE}.LocalAssemblies/{prefix}/ref.fasta"
        assert os.path.exists(ref)
        collen.append( len(list(SeqIO.parse(ref, "fasta"))[0].seq ) )
    df["collapseLength"] = collen


numeric_vars = list( df.select_dtypes(include=[np.number]) )
numeric_vars.remove("cut")
numeric_vars.remove("uid")

correlation = df[numeric_vars].corr()


## PLOTS ##
if(args.length):
	g3 = sns.distplot( np.log10(df["length"].dropna()) , kde=False, hist_kws=dict(edgecolor="k", linewidth=2))
	g3.set_xticklabels( [ "$10^{{{:.1f}}}$".format(item) for item in g3.get_xticks() ]  , rotation=0)

	x = df["length"].dropna()
	mean = x.mean()
	median = x.median()
	n50 = np.sort(x)[ np.where( np.cumsum(np.sort(x)) >= (sum(x)/2)  )[0][0] ]
	plt.axvline(x=np.log10(  mean ), linestyle="dashed", color="orange", label = "Mean: {:.0f}".format(mean))
	plt.axvline(x=np.log10(  median ), linestyle="dashed", color="darkred", label = "Median: {:.0f}".format(median))
	plt.axvline(x=np.log10(  n50 ), linestyle="dashed", color="darkgreen", label = "N50: {:.0f}".format(n50))

	g3.set_ylabel("# of assemblies")
	g3.set_xlabel("Length of assembly")
	plt.legend()
	plt.savefig(args.length)
	plt.clf()

# In[190]:
if(args.collapse):
	g4 = sns.jointplot(df["length"], df["collapseLength"], height=12)#, xlim=(smallest, largest), ylim=(smallest, largest))
	g4.ax_joint.set_xscale('log')
	g4.ax_joint.set_yscale('log')
	g4.set_axis_labels("Length of assembly", "Length of collapse")
	plt.savefig(args.collapse)
	plt.clf()


# In[192]:
if(args.bar):
	x=["collapsed sequence", "SDA assembly"]
	y=[df.drop_duplicates(subset="prefix")["collapseLength"].sum()/10**6, df["length"].sum()/10**6]
	g6 = sns.barplot(x=x, y=y, ci="sd")
	
	counter = 0
	for  name, mbp in zip(x,y):
		print(name, mbp)
		g6.text(counter, mbp, "{:.2f}".format(mbp), horizontalalignment='center' )
		counter += 1
	
	g6.set_ylabel("Mbp")
	plt.savefig(args.bar)
	plt.clf()

# In[191]:
if(args.pair):
	g5 = sns.pairplot(df[numeric_vars].dropna(axis=0), height = 6)
	g5.map_lower(corrfunc)
	g5.map_upper(corrfunc)
	plt.savefig(args.pair)
	plt.clf()




"""
# In[186]:


print(list(df))
print(correlation)


# In[187]:


g = sns.jointplot(df["length"], df["nReads"], height=12)
g.ax_joint.set_xscale('log')
g.ax_joint.set_yscale('log')


# In[188]:


g2 = sns.distplot( qfilter(11000*df["nReads"]/df["length"], bot= 0).dropna() , kde=False, hist_kws=dict(edgecolor="k", linewidth=2))


# In[189]:

"""
