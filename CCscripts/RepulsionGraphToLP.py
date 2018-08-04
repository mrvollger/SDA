#!/usr/bin/env python

import ABPUtils
import sys



import argparse

ap = argparse.ArgumentParser(description="Given a graph that has repulsion edges (color=red), write an LP")
ap.add_argument("--graph", help="Input file.", required=True)
ap.add_argument("--lp", help="Output linear program",required=True)
ap.add_argument("--maximize", help="Produce MaxAgree LP.", default=False, action='store_true')
ap.add_argument("--wn", help="Negative weight", default=1,type=int)
ap.add_argument("--wp", help="Positive weight", default=1, type=int)
ap.add_argument("--relaxed", help="Allow non-binary solution.", default=False, action='store_true')
args = ap.parse_args()


g = ABPUtils.ReadGraph(args.graph)

lp = open(args.lp, 'w')

def Var(i,j):
    return "x_" + str(i) + "_" + str(j)


def GetWeight(edge):
    # for now use unite weights.
    return 1

    # later full weight will be used
    if ('weight' in edge):
        return float(edge['weight'])
    else:
        return 1

attract = []
repulse = []
for i in g.nodes():
    for j in g[i]:
        if (j > i):
            #
            # Check for repulsion edge
            #

            if 'cost' in g[i][j] and g[i][j]['cost'] < 0 :
                repulse.append( " {:2.2f} ".format(abs(args.wn*g[i][j]['cost'])) + Var(i,j))
            else:
                attract.append( " {:2.2f} ".format(1.0) + Var(i,j))


if (args.maximize == False):
    lp.write("Minimize\n")    
    nVars = 0
    lp.write("repulsion: ")
                    
    lp.write(" + ".join(attract) + " - " + " - ".join(repulse))
    lp.write("\n")

else:
    lp.write("Maximize\n")
    lp.write("attraction:")
    lp.write(" + ".join(attract) + " - " + " - ".join(repulse))
    lp.write("\n")
    
                
lp.write("Subject To\n")
#
# Write triangle inequality constraint
#
cn = 0
for i in g.nodes():
    for k in g[i]:
        if (k > i):
            for j in g[i]:
                if (j > i and j < k and g.has_edge(j,k)):
                    lp.write("c" + str(cn) + ": " + Var(i,k) +  " - " + Var(j,k) + " - " + Var(i,j) + " <= 0\n")
                    cn+=1                    
                    lp.write("c" + str(cn) + ": " + Var(i,j) +  " - " + Var(i,k) + " - " + Var(j,k) + " <= 0\n")
                    cn+=1                    
                    lp.write("c" + str(cn) + ": " + Var(j,k) +  " - " + Var(i,k) + " - " + Var(i,j) + " <= 0\n")                    
                    cn+=1

print "LP has " + str(cn)+ " variables"
lp.write("\n")
lp.write("Bounds\n")
for i in g.nodes():
    for j in g[i]:
        if (j > i):
            lp.write(" 0 <= " + Var(i,j) + " <= 1\n")
if (args.relaxed == False):
    lp.write("binary\n")
    for i in g.nodes():
        for j in g[i]:
            if (j > i):
                lp.write(Var(i,j) + "\n")

lp.write("End\n")
lp.close()

                    

