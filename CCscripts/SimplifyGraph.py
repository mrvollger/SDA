#!/usr/bin/env python

import argparse
import networkx as nx
import ABPUtils
import sys
import numpy as np
import itertools

ap = argparse.ArgumentParser(description="Prune non-solution edges from a graph")
ap.add_argument("--graph", help="Input graph", required=True)
ap.add_argument("--max-nodes", help="Maximum number of nodes in the graph.", type=int, default=None, dest='maxNodes')
ap.add_argument("--max-edges", help="Maximum number of edges in the graph.", type=int, default=None, dest='maxEdges')
ap.add_argument("--max-tri-diff",help="Maximum difference between positive and negative edges.",type=int, default=20, dest='maxTriDiff')
ap.add_argument("--pos", help="Positions of variant sites, and minor allele at each site", required=True)
ap.add_argument("--out", help="Output file.", required=True)
ap.add_argument("--snv", help="Fragment snv file.", required=True)
ap.add_argument("--clique", help="Condense cliques of this isze", type=int, default=None)
ap.add_argument("--min", help="Minimum count to include in a cluster", default=3, type=int)
ap.add_argument("--triangle-radius", help="Radius to search around points for triangles.", type=int, default=6, dest="triangleRadius")
ap.add_argument("--write-iterations", help="Write iterations to this file.", default=None, dest='iterations')
ap.add_argument("--iter-step", help="Step between writing iterations", default=1, dest='iterStep', type=int)
ap.add_argument("--plot-layout", help="Use the layout in this file when rendering the plot", default=None, dest="plotLayout")
args = ap.parse_args()
g = ABPUtils.ReadGraph(args.graph)

C = {}


def RankClusters(P, C):
    cs = [0]*len(C)
    start = [0]*len(C)
    end = [None]*len(C)
    ci = 0
    for c in C:
        s = None
        e = None
        for p in P:
            if p in c:
                cs[ci]+=1
                if s is None or p < s:
                    s = p
                if e is None or p > e:
                    e = p
        start[ci] = s
        end[ci] = e
                    
        ci+=1
    return (cs, start, end)

# read snv information


snvFile = open(args.snv)
snvs = [ABPUtils.ParseSNVLine(line) for line in snvFile ]
# Each row in snvPos are the positions for SNVs for read i
snvPos = [np.array(sorted(v[1].keys())) for v in snvs]
# Each row in snvVal are the values for SNVs for read i
snvVal = [np.array([v[1][i] for i in sorted(v[1].keys())]) for v in snvs]

posFile = open(args.pos)
pos = {}
for line in posFile:
    v = line.split()
    if (v[3] == 'N'):
        v[3] = 1
    pos[int(v[0])] = [int(v[1]), int(v[2]), int(v[3]), int(v[4])]

graphPos = {}


def Condense(g, c):
    #
    # If there are no nodes to condense, just return
    #
    if len(c) <= 1:
        return
    #
    # Pick the first node is the new target
    #
    t = c[0]
    for ci in c[1:]:
        for d in g[ci]:
            if d in g[t]:
                g[t][d]['cost']+=g[ci][d]['cost']
            else:
                g.add_edge(t,d, g[ci][d])
        g.remove_node(ci)


def ScorePairCut(g, a, b):
    cost = 0
    for na in g[a]:
        #  Check if neigbor of a is in neighbor of b
        for nb in g[b]:
            if nb in g[na]:
                cost += g[na][nb]['cost']
    return cost

def ScoreTriangleCut(g, t):
    costs = [ScorePairCut(g,t[0], t[1]),
             ScorePairCut(g,t[1], t[2]),
             ScorePairCut(g,t[0], t[2])]
    print costs
    return min(costs)


def GetNeighbors(g,n,neighbors):
    if neighbors is None:
        return sorted(g[n].keys())
    else:
        return neighbors[n]
    
def CountMissingNeighbors(g, l, i, neighbors):
    #
    # Consider a triangle with nodes a, b, and c.
    # The missing triangle neighbors are neighbors of a that are not a
    # neighbor of b nor c.  This is an estimate of how well connected
    # the nodes are.
    #
    diff = 0
    shared = GetNeighbors(g, l[i], neighbors)
    r = range(0,len(l))
    r.remove(i)
    for nodeIndex in r:
        node   = l[nodeIndex]
        shared = np.intersect1d(shared, GetNeighbors(g, node,neighbors))
    return len(GetNeighbors(g, l[i], neighbors)) - len(shared)

def ScoreNeighbors(g,t,i,phaseNeighbors, repulsionNeighbors):
    diff = 0
    neighbors = phaseNeighbors[t[i]] 
    r = range(0,len(t))
    r.remove(i)
    repulsion = 0
    for nodeIndex in r:
        node = t[nodeIndex]
        repulsion += len(np.intersect1d(neighbors, repulsionNeighbors[node]))

    return repulsion

def ScoreTriangleRepulsion(g, t, phaseNeighbors, repulsionNeighbors):
    return sum([ScoreNeighbors(g,t,i,phaseNeighbors, repulsionNeighbors) for i in range(0,len(t))])

def ScoreTrianglePhase(g, t, phaseNeighbors):
    return sum([ScoreNeighbors(g,t,i,phaseNeighbors, phaseNeighbors) for i in range(0,len(t))])
               
    
def ScoreTriangleDifference(g, t, neighbors):
    #
    # The triangle difference is the total missing neighbors.
    #
    
    return sum([CountMissingNeighbors(g,t,i, neighbors) for i in range(0,len(t))])

for n in g.nodes():
    g.node[n]['pos'] = int(g.node[n]['pos'])-1
    graphPos[g.node[n]['pos']] = n

def GetNodes(m, l):
    n = []
    for i in l:
        if i in m:
            n.append(m[i])
    return n
            
def Cost(edge):
    if 'cost' in edge:
        return edge['cost']
    else:
        return 0

def ScoreCluster(G, c):
    score = 0
    for i in c:
        for j in c:
            if i in G and j in G[i]:
                score += Cost(G[i][j])
    return score

                
def ScoreCut(G, a, b):
    score = 0
    for i in a:
        for j in b:
            if i in G and j in G[i]:
                score += Cost(G[i][j])
    return score

def ScoreClusterPair(G, a, b):
    return ScoreCut(G,a,b) + ScoreCluster(G,a) + ScoreCluster(G,b)
    
minorAlleleFrags = []
minorAlleleNode = []
minorAlleleNodes = []
# Retain minor allele sites for each fragment
for frag in snvs:
    minorAlleleFrag = []
    minorAlleleNode = []
    sites = sorted(frag[1].keys())

    for site in sites:
        if site in pos and frag[1][site] == pos[site][2] and site in graphPos:
            minorAlleleFrag.append(site)
            minorAlleleNode.append(graphPos[site])

    minorAlleleFrags.append(minorAlleleFrag)
    minorAlleleNodes.append(minorAlleleNode)
    



    
fragLengths = sorted([(len(minorAlleleFrags[f]), f) for f in range(0,len(minorAlleleFrags))], reverse=True)

fragIndex = [f[1] for f in fragLengths]
# process fragments sorted by length
radius = args.triangleRadius
triangles = {}
nTri = 0

for f in fragIndex:
    frag = minorAlleleNodes[f]

    if len(frag) < 3:
        continue
    for i in range(0,len(frag)-2):
        
        for j in range(i+1, min(i+radius-1, len(frag)-1)):
            
            if frag[j] not in g[frag[i]]:
                continue
            for k in range(j+1, min(i+radius, len(frag))):
                if frag[k] not in g[frag[i]] or frag[k] not in g[frag[j]]:
                    continue
                t = tuple(sorted([frag[i],frag[j],frag[k]]))
                if t not in triangles:
                    triangles[t] = 1
                triangles[t] += 1
                nTri+=1

    
def GetPhasedNeighbors(g, node, phased=True):
    neighbors = []
    for neighbor in g[node]:
        if phased is True and Cost(g[node][neighbor]) >= 0:
            neighbors.append(neighbor)            
        elif phased is False and Cost(g[node][neighbor]) < 0:
            neighbors.append(neighbor)
    return sorted(neighbors)


phasedNeighbors = { node : GetPhasedNeighbors(g, node) for node in g.nodes() }
repulsionNeighbors = { node : GetPhasedNeighbors(g, node, False) for node in g.nodes() }

triList = sorted([(v,k) for k,v in triangles.iteritems()], reverse=True)

triDiff = [ScoreTriangleDifference(g, triList[i][1], neighbors=phasedNeighbors) for i in range(0,len(triList))]
orderedTriDiff = sorted([[triDiff[i],i] for i in range(0,len(triDiff))])



condensedNodes = {}
# Determine the number of triangles that would have to be merged to reach the target size.
i = 0
while i < len(orderedTriDiff) and len(g.node) - len(condensedNodes)  > args.maxNodes:
    idx = orderedTriDiff[i][1]
    #
    #Record all of the nodes that will be condensed in some form or another.
    #
    for t in triList[idx][1][1:]:
        condensedNodes[t] = True
    i+=1
    
numCondensedTriangles = i

        
#from IPython import embed
#embed()

triPhase = [ScoreTrianglePhase(g, triList[i][1], phasedNeighbors) for i in range(0,len(triList))]
triRepulsion = [ScoreTriangleRepulsion(g, triList[i][1], phasedNeighbors, repulsionNeighbors) for i in range(0,len(triList))]

def GetRoot(tree, node):
    while node in tree:
        if tree[node] == node:
            return node
        node = tree[node]
    return node

def AddNode(tree, node, parent):
    root = GetRoot(tree,parent)
    tree[node] = root

if args.maxNodes is None:
    args.maxNodes = len(g.nodes())
if args.maxEdges is None:
    args.maxEdges = len(g.edges())

triIndex = 0
ct = {}
iter = 0

sys.stdout.write( "Number of candidate triangles: " + str(numCondensedTriangles) + "\n")


if args.plotLayout is not None:
    plotLayout = ABPUtils.ReadLayout(args.plotLayout)


allColors = { ABPUtils.GetColor(g[n]) : True for n in g.nodes()}
numColors = len(allColors) + 1
while triIndex < numCondensedTriangles and (len(g.nodes()) > args.maxNodes or len(g.edges()) > args.maxEdges):
    idx = orderedTriDiff[triIndex][1]
    tri = triList[idx][1]
    t = (GetRoot(ct, tri[0]), GetRoot(ct,tri[1]), GetRoot(ct,tri[2]))
    if t[0] != t[1] and t[1] != t[2] and t[0] != t[2]:
        print "tri " + str(tri)
        print "condensing" + str(t)
    
        if args.iterations is not None:
            if iter % args.iterStep == 0:
                iterFileName = args.iterations + "." + "{:0>4}".format(iter) + ".png"
                g2 = g.copy()
                g2.node[t[0]]['color'] = 99999
                g2.node[t[1]]['color'] = 99999
                g2.node[t[2]]['color'] = 99999
                print "drawing graph " + iterFileName
                ABPUtils.DrawGraph(g2, plotLayout, iterFileName, drawRepulsion=True, minColors=numColors)
                print "done"
            
        Condense(g, t)
        print "#Edges: " + str(len(g.edges())) + " #Nodes: " + str(len(g.nodes()))
            
        AddNode(ct, t[1], t[0])
        AddNode(ct, t[2], t[0])
        if args.iterations is not None:
            iterFileName = args.iterations + "." + "{:0>4}".format(iter+1) + ".png"
            if iter % args.iterStep == 0:
                print "drawing graph " + iterFileName                    
                #ABPUtils.DrawGraph(g, plotLayout, iterFileName, drawRepulsion=True, minColors=len(allColors) + 1)
                print "done"
            iter+=2
    triIndex+=1

print "End simplify graph " + str(len(g.nodes())) + "\t" + str(len(g.edges()))    
ABPUtils.WriteGraph(g, args.out)
