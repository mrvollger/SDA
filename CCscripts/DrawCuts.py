#!/usr/bin/env python
import ABPUtils

import argparse

ap = argparse.ArgumentParser(description="Plot cuts individually")
ap.add_argument("graph", help="Original graph file.")
ap.add_argument("cuts", help="Cuts file.")
ap.add_argument("--out", help="Output file.", default="./")
args = ap.parse_args()


g = ABPUtils.ReadGraph(args.graph)
cuts = ABPUtils.ReadCuts(args.cuts)
cutIndex = 1
ABPUtils.ColorGraphByCut(g,cuts)
for cut in cuts:
    sub = g.subgraph(list(cut))
    sub.graph['NumVertexColors']=len(cuts)
    ABPUtils.DrawGraph(sub, "{}subgraph.{}.png".format(args.out,cutIndex))
    cutIndex+=1
