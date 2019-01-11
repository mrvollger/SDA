#!/usr/bin/env python
nucToIndex = {'A':0,'C':1,'G':2,'T':3, '-':4}
nucs = ['A', 'C','G','T']
import networkx as nx
import sys
import pickle
import argparse
import intervaltree
import bisect
import numpy as np
import scipy as sp
import scipy.stats
import time
import ast
import re
import math
from collections import *
import matplotlib
matplotlib.use('Agg') # this forces matplotlib to use the backend instead of X server, which will fail on the cluster. mrv.  
import matplotlib.pyplot as plt


preNRe = re.compile("(n+).*")
sufNRe = re.compile(".*[^n]([n]+)$")

def Boundaries(gt):
    boundaries = [ [0,0] for i in range(0,gt.shape[0]) ]
    for i in range(0,gt.shape[0]):
        try:
            m = np.where(gt[i][:] != 'n')[0]
        except IndexError:
            sys.exit(0)
        if (len(m) > 1):
            boundaries[i][0],boundaries[i][1] = min(m),max(m)
        else:
            boundaries[i][0],boundaries[i][1] = 0,0
    return boundaries

        
def EvaluateTruthPartitions(g, vg):
    colors = { g.node[v]['color']: True for v in g.nodes() }
    idx = 0
    for c in colors.keys():
        colors[c] =idx
        idx+=1
    partitions = {}
    for n in g.nodes():
        partitions[n] = colors[g.node[n]['color']]
    score = EvaluatePartitions(partitions, vg)
    return score


def EvaluatePartitions(partitions, vectorizedGraph):
    score = 0
    for s in range(0,len(vectorizedGraph)):
        if (s not in partitions):
            continue
        p = partitions[s]
        for t in vectorizedGraph[s]:
            if (len(t) == 0):
                continue
            if (partitions[t[0]] == p):
                score += t[1]
            else:
                score -= t[1]
    return score
    
def VectorizeGraph(g):
    nodes = g.nodes()
    vg = []
    for i in range(0,max(nodes)+1):
        vg.append([ (), ()])
    for n in nodes:
        vg[n] = [(k,v['capacity']) for (k,v) in g[n].iteritems() ]
    return vg


def RemoveIsolated(g):
    toRemove = []
    for n in g.nodes():
        if (len(g[n])== 0):
            toRemove.append(n)
    for n in toRemove:
        g.remove_node(n)

import random

def RandomSwap(partitions, nodes, maxNeighborPartitions):
    L = len(maxNeighborPartitions)
    import random
    for i in range(0,L):
        #
        # Get partition of this node.
        #
        n = nodes[i]
        curPartition = partitions[n]
        curMaxPartition =maxNeighborPartitions[i]
        if (curPartition != curMaxPartition):
            #
            # Swap with another node whose max neighbor partition is this one, at random.
            #
            swapStart = random.randint(0,L)
            foundSwap = False
            for j in range(0,L):
                idx = (swapStart + j)%L
                swapi = nodes[idx]
                swapMaxPartition = maxNeighborPartitions[idx]
                if (partitions[swapi] != swapMaxPartition and
                    swapMaxPartition == curPartition and
                    curMaxPartition == partitions[swapi]):
                    swapPartition = partitions[swapi]
                    partitions[swapi] = curPartition
                    partitions[i] = swapPartition
                    foundSwap = True
                    break
                    

def OptimizePartitions(partitions, vectorizedGraph, nodes, nIter=1000):
    prevIter = 0
    for i in range(0,nIter):
        nNodes = len(nodes)
        mnp = [ MaxNeighborPartition(partitions,vectorizedGraph,nodes[i]) for i in range(0,nNodes)]
        RandomSwap(partitions, nodes, mnp)
        score = EvaluatePartitions(partitions, vectorizedGraph)
        if (score == prevIter):
            return score
        print score
        prevIter = score
        
    
    

def ScoreNode(partitions, vectorGraph, node):
    score = 0
    adjList = vectorGraph[node]
    p = partitions[node]
    for adj in adjList:
        if (p == partitions[adj[0]]):
            score += adj[1]
        else:
            score -= adj[1]
    return score
    

def GetCost(g, a,b):
    if 'cost' in g[a][b]:
        return g[a][b]['cost']
    else:
        return 1
    
def MaxNeighborPartition(partitions, vectorGraph, node):
    np = len(partitions)
    pw = [0]*np
    adjList = vectorGraph[node]
    if (len(adjList) == 0):
        return -1
    for adj in adjList:
        p = partitions[adj[0]]
        pw[p] += adj[1]
    mp = 0
    mpi = -1
    for i in range(0,np):
        if (pw[i] > mp):
            mp = pw[i]
            mpi = i
    return mpi
    

def EvaluateSwap(partitions, vectorGraph, nodes):
    #
    # Evaluate current
    #
    
    np = [[partitions[nodes[0]], partitions[nodes[1]]],
          [partitions[nodes[1]], partitions[nodes[0]]]]

    scores = [0,0]
    for it in range(0,2):
        score = 0    
        for i in range(0,2):
            n = nodes[i]
            p = np[it][i]
            adjList = vectorGraph[n]
            for adj in adjList:
                if (p == partitions[adj[0]]):
                    score += adj[1]
                else:
                    score -= adj[1]

        scores[it] = score
    return scores
    


from sets import Set
        
def GetCutNeighbors(gAdjList, cut):
    neighbors = Set([])
    for node in cut:
        neighbors.update(gAdjList[node])
    neighbors.difference_update(cut)
    return neighbors

def GetAttrNeighbors(g,n):
    return GetNeighbors(g,n,1)

def GetReplNeighbors(g,n):
    return GetNeighbors(g,n,-1 )



def ParseFragTuple(fragTuple):
    v = fragTuple.split(",")
    return [int(v[0]), v[1],v[2],v[3],int(v[4]),int(v[5])]


def ParseFragLine(line, keepRef=False):
    v = line.split()
    name, target, snvs= v[0],v[2], [ParseFragTuple(t) for t in v[3:]]
    snvPos  = np.array([s[0] for s in snvs], dtype=int)
    snvRef  = np.array([s[1] for s in snvs], dtype='U')
    snvAlt  = np.array([s[2] for s in snvs], dtype='U')
    snvRead = np.array([s[3] for s in snvs], dtype='U')
    snvPre  = np.array([s[4] for s in snvs], dtype=np.int16)
    snvPost = np.array([s[5] for s in snvs], dtype=np.int16)

    if (keepRef == True):
        return name, target, snvPos, snvRef, snvAlt, snvRead, snvPre, snvPost
    else:
        i = (snvRead != snvRef) & (snvRead != '-')
        return name, target, snvPos[i], snvRef[i], snvAlt[i], snvRead[i], snvPre[i], snvPost[i]

def ParseSNVLine(line):
    ts = line.find("\t")
    name = line[0:ts]
    rem = line[ts+1:]
    snvs = ast.literal_eval(rem)
    return(name, snvs)


class Fragment:
    def __init__(self, line, keepRef=False):
        self.name, self.target, self.snvPos, self.snvRef, self.snvAlt, self.snvRead, self.snvPre, self.snvPost = ParseFragLine(line, keepRef)


    def Length(self):
        return len(self.snvPos)

    def Range(self):
        if (len(self.snvPos) == 0):
            return (0,0)
        else:
            return self.snvPos[0], self.snvPos[-1]

def FragmentOverlap(a,b,i=None):

    #
    # Find the indices in a an b where the fragments share snv's
    #
    ai = np.in1d(a.snvPos, b.snvPos, assume_unique=True)
    bi = np.in1d(b.snvPos, a.snvPos, assume_unique=True)
#    print len(ai)
#    print len(bi)
#    print str(ai)
#    print str(bi)
#    print str(a.snvRead[ai])
#    print str(b.snvRead[bi])
    return sum(ai), sum(a.snvRead[ai] == b.snvRead[bi])

def BuildOverlapGraph(frags, minOverlap, outFile=None):
    frags = sorted(frags, key=lambda f: f.snvPos[0])
    
    active = {}
    activeKeys = []
    print "Building overlap graph."

    prevTime = time.time()
    num=0
    overlaps = []
    for f in frags:
        r = f.Range()
        num+=1        

        if (r[1] <= r[0]):
            continue
        if (r is None):
            continue

        #
        # Firt remove any active set elements ending before the start of this element.
        #
        toRemove = []
        
        for p in activeKeys:
            if (p >= r[0]):
                break
            else:
                toRemove.append(p)
        if (len(toRemove) > 0):
            for p in toRemove:
                del active[p]
            del activeKeys[0:len(toRemove)]

        
        if (len(active) > 0):

            for p in activeKeys:
                for fo in active[p]:
                    ovp = FragmentOverlap(f, fo)
                    overlaps.append([f,fo, ovp])
        if (r[1] not in active):
            active[r[1]] = [f]
            activeKeys.append(r[1])
            activeKeys.sort()
        else:
            active[r[1]].append(f)
                

    return overlaps

def ReadCuts(cutsFileName):
    f = open(cutsFileName)
    return [Set([int(i) for i in l.split()]) for l in f]

def ReadGenotypeMatrix(matFile):
    groups = {}
    groupList = []
    readNames = []
    gtList = []
    index = 0
    for line in matFile:
        v = line.split()
        gtl = np.array(list(v[0]))
        gtList.append(gtl)
        #
        # If the name of the read is packed at the end, add this.
        #
        if (len(v) >= 1):
            readNames.append(v[1])
        #
        # If there is ground truth information for each group, add this
        #

        if (len(v) >= 3):
            if (v[2] not in groups):
                groups[v[2]] = []
            groups[v[2]].append(index)
            groupList.append(v[2])
        index +=1
            
    gt = np.array(gtList)
    if readNames == []:
        readNames = None
    if groupList == []:
        groupList = None

    groupListArray = np.array(groupList)
    return {'mat':gt, 'readNames':readNames, 'groupList':groupListArray, 'groups':groups}



def GetRefAltCounts(gt,alph='1.n'):
    (a,r,g) = list(alph)
    altList = []
    refList = []
    for i in range (0,gt.shape[1]):
        altList.append(len(np.where(gt[:,i] == a)[0]))
        refList.append(len(np.where(gt[:,i] == r)[0]))
            
            
    return (refList, altList)


def GetMinorAllele(ref, alt, minCov, maxCov):
    minor = None
    major = None
    if (ref > 0 and alt > 0):
        if ( ref >= minCov and ref <= maxCov and  alt > maxCov):
            minor = '.'
            major = '1'
        elif (alt >= minCov and alt <= maxCov and ref > maxCov):
            minor = '1'
            major = '.'

    return (minor, major)




from itertools import groupby

def GetMinorAlleles(gt, minCov, maxCov):

    (ref,alt) = GetRefAltCounts(gt)
    minorAllele= np.array(['.' for i in range(0,len(ref))])
    for i in range(0,len(ref)):
        (minor,major) = GetMinorAllele(ref[i], alt[i], minCov, maxCov)
        minorAllele[i] = minor
    return minorAllele


def GetReadIndices(gt,key):
    index = []
    
    for i in range(0,len(key)):
        if (key[i] is None):
            index.append([])
        else:
            index.append(np.where(gt[:,i] == key[i])[0])
    return index


def GetMaxName(idx, names):
    maxCount = 0
    maxName = "None"
    if (len(idx) > 0):
        nameCount = [ (k,len(list(v))) for k,v in groupby(sorted(names[idx])) ]

        for i in range(0,len(nameCount)):
            if (nameCount[i][1] > maxCount):
                maxCount = nameCount[i][1]
                maxName  = nameCount[i][0]

    return (maxName, maxCount)

def WriteMinorIndexTable(gt,minCov,maxCov,names, outFileName):
    outFile = open(outFileName, 'w')
    minorAllele = GetMinorAlleles(gt, minCov, maxCov)
    readIndices = GetReadIndices(gt, minorAllele)
    nameEnumeration = {}
    ni = 1
    for n in names:
        if (n not in nameEnumeration):
            nameEnumeration[n] = ni
            ni+=1

    for i in range(0,len(readIndices)):
        (name,count) = GetMaxName(readIndices[i],names)
        if (len(readIndices[i]) == 0):
            nameIndex = 0
        else:
            nameIndex = nameEnumeration[name]
        for j in range(0,len(readIndices[i])):
            outFile.write("{}\t{}\t{}\n".format(i,readIndices[i][j],nameIndex))
    
class MI:
    def __init__(self, nSharedMinorAllele, totalOverlapping, fracIMinor, fracJMinor, lrt, maxGroup):
        self.nSharedMinorAllele = nSharedMinorAllele
        self.totalOverlapping = totalOverlapping
        self.fracIMinor = fracIMinor
        self.fracJMinor = fracJMinor
        self.maxGroup = maxGroup
        self.lrt = lrt

def End(a):
    if (len(a) > 0):
        return a[-1]
    else:
        return 0
def Start(a):
    if (len(a) > 0):
        return a[0]
    else:
        return 0
    
def FindMutualInformation(gt, minCov, maxCov, names, filt=None, miOutFileName=None, minNShared=0, accuracy=0.8, minLRT=6, vcf=[], conflicts=False, indices=None, pos=None):
    boundaries = Boundaries(gt)

    i = 0
    
    if (miOutFileName is not None):
        outFile = open(miOutFileName, 'w')
        miOutStr = ""
        
    (ref,alt) = GetRefAltCounts(gt)
    minorIdx = []
    majorIdx = []
    allIdx   = []

    mi = [ {} for i in range(0,len(ref)) ]

    if (indices is None):
        indices = range(0,len(ref))
    for i in indices:
        
        (minor,major) = GetMinorAllele(ref[i], alt[i], minCov, maxCov)

        if (minor is None):
            minorIdx.append([])
        else:
            minorIdx.append(np.where(gt[:,i] == minor)[0])
        if (major is None):
            majorIdx.append([])
        else:
            majorIdx.append(np.where(gt[:,i] == major)[0])

#        print(str(minorIdx))
        allIdx.append(np.where(gt[:,i]!='n')[0])

    
    ends   = [End(minorIdx[i]) for i in range(0,len(minorIdx))]
    starts = [Start(minorIdx[i]) for i in range(0,len(minorIdx))]

    siteConflicts = Counter()#{i: 0 for i in range(0,len(gt[0]))}
    siteAgreements = Counter()# {i: 0 for i in range(0,len(gt[0]))}
        
    for i in range(0,len(indices)-1):

        nOverlapping = 0
        nTested = 0

        if (len(minorIdx[i]) > 0):
            miList = []
            for j in range(i+1,len(indices)):
                if (starts[j] > ends[i]):
                    continue
                nMinorShared = 0
                nMinorSharedNames = 0
                maxGroup = -1
                nAllShared = 0
                maxGroup = "None"

                #
                # Find and count all reads that span from i to j
                #
                allOverlapping   = np.intersect1d(allIdx[i],allIdx[j], assume_unique=True)
                totalOverlapping = len(allOverlapping)

                #
                # Of the reads spanning from i to j, how many in i/j are minor allele
                #

                iMinor = np.intersect1d(minorIdx[i], allOverlapping, assume_unique=True)
                jMinor = np.intersect1d(minorIdx[j], allOverlapping, assume_unique=True)
                

                nIMinor = len(iMinor)
                nJMinor = len(jMinor)

                #
                # Now how many of those two sets are overlapping
                #
                sharedMinor = np.intersect1d(iMinor, jMinor, assume_unique=True)
                nSharedMinor = len(sharedMinor)

                
                if (conflicts == True and nIMinor > 0 and nJMinor > 0):
                    iMajor = np.intersect1d(majorIdx[i], allOverlapping, assume_unique=True)
                    jMajor = np.intersect1d(majorIdx[j], allOverlapping, assume_unique=True)

                    iConflicts = np.intersect1d(iMinor, jMajor, assume_unique=True)
                    jConflicts = np.intersect1d(jMinor, iMajor, assume_unique=True)
                    
                    nIConflicts = len(iConflicts)
                    nJConflicts = len(jConflicts)
                    nConflicts = nIConflicts + nJConflicts
                        
                    #print str(pos[i]) + "\t" + str(pos[j]) + "\t" + str(totalOverlapping) + "\t" + str(nSharedMinor) + "\t" + str(nConflicts) + "\t" + str(nIConflicts) + "\t" + str(nJConflicts) + "\t{:2.2f}".format(nConflicts/(2.0*(nConflicts+nSharedMinor)))
                    if (conflicts == True and False):
                        if (pos is not None):
                            print iMinor
                            print jMinor
                            print sharedMinor
                            print iConflicts
                            print jConflicts
                        else:
                            print str(iMinor)
                            print str(jMinor)
                        
                

                #
                # Compute the likelihood of nSharedMinor conditional that they are linked.
                #


                

                #
                # Now comapre the fraction of iMinor/jMinor and the sharedMinor to the
                # total coverage in this region.
                #
                (maxGroup, maxVal) = GetMaxName(sharedMinor, names)

                nMinorSharedNames = maxVal

                #
                # Possibly write everything to a graph so it doesn't need to be passed back.
                #
                lrti = 0
                lrtj = 0
                if (nSharedMinor >= minNShared):
                    lrti = np.log10(sp.stats.binom_test(nSharedMinor, nIMinor, accuracy) / sp.stats.binom_test(nSharedMinor, nIMinor, 1-accuracy)+1)
                    lrtj = np.log10(sp.stats.binom_test(nSharedMinor, nJMinor, accuracy) / sp.stats.binom_test(nSharedMinor, nJMinor, 1-accuracy)+1)
                    lrt = min(lrti, lrtj)


                else:
                    lrt = 0
                #
                # Store the color or this group if the names are known.
                #
                
                if (nSharedMinor > 0 and nMinorSharedNames /float(nSharedMinor) < 0.60 ):
                    maxGroup = "None"

                #
                # Add this edge to the pairwise graph.
                #
                nTested +=1
                if (lrt >= minLRT):
                    mi[i][j] = MI(nSharedMinor,
                                  totalOverlapping,
                                  nSharedMinor/float(nIMinor),
                                  nSharedMinor/float(nJMinor),
                                  lrt,
                                  maxGroup)
                    nOverlapping+=1

                miList.append((vcf[j], lrt))
                if (miOutFileName is not None):
                    miOutStr += str(vcf[i]) + "\t" + str(vcf[j]) + "\t" + str(lrt) + "\t" + str(lrti) + "\t" + str(lrtj) + "\t" + str(nIMinor) + "\t" + str(nJMinor) + "\t" + str(nSharedMinor) + "\n"
                
    if (miOutFileName is not None):
        #
        # Write whole thing in one list.
        #
        outFile.write(miOutStr)


#                outFile.write(str(vcf[i]) + "\t" str(vcf[j]) + " ".join([str(v[0]) + "," + str(v[1]) for v in miList] ) + "\n")
#                if (nSharedMinor >= minNShared):
#                    if (i < len(vcf) and j < len(vcf)):
#                        iPos = vcf[i]
#                        jPos = vcf[j]
#                    else:
#                        iPos = -1
#                        jPos = -1
#                    outFile.write("\t".join(str(i) for i in [iPos, jPos,
#                                                             i,j,
#                                                             nIMinor, nJMinor,
#                                                             nSharedMinor, nMinorSharedNames,
#                                                             totalOverlapping,
#                                                             "{:2.2f}".format(lrt),
#                                                             "{:2.2f}".format(nSharedMinor/float(nIMinor)),
#                                                             "{:2.2f}".format(nSharedMinor/float(nJMinor)) ]) + "\n")
#


    return mi
            
        

def GetRefAltLists(gt, alph='1.n'):
    (a,r,g) = list(alph)
    altList = []
    refList = []
    for i in range (0,gt.shape[0]):
        altList.append(np.where(gt[i] == a)[0])
        refList.append(np.where(gt[i] == r)[0])
    return (altList, refList)


def GetMean(a):
    if (len(a) == 0):
        return 0
    else:
        return np.mean(a)

def GetStart(v, g):
    ng = np.where(v != g)[0]
    if (len(ng) == 0):
        return 0
    else:
        return ng[0]


def GetRange(v, g):
    ng = np.where(v != g)[0]
    if (len(ng) < 2):
        return (0,0)
    else:
        return (ng[0],ng[-1])

def Intersect(a,b):
    na = len(a)
    nb = len(b)
    ai = 0
    bi = 0
    u = 0
    while (ai < na and bi < nb):
        if (a[ai] == b[bi]):
            u+=1
            ai+=1
            bi+=1
        elif (a[ai] < b[bi]):
            ai+=1
        else:
            bi+=1
    return u

def StoreDistanceMatrix(gt,penalty=1):
    (altList, refList) = GetRefAltLists(gt)
    ngt = len(gt)
    dm = np.ndarray(shape=(ngt,ngt), dtype=int)
    for i in range(0,ngt):
        scores = np.array([0]*ngt)
        for j in range(i,ngt):
            nMatch = len(np.intersect1d(altList[i],altList[j], assume_unique=True)) +\
              len(np.intersect1d(refList[i],refList[j], assume_unique=True))
            
            nMis   = (len(np.intersect1d(altList[i],refList[j], assume_unique=True))+\
                len(np.intersect1d(refList[i],altList[j], assume_unique=True)))*penalty
            
            dm[i,j] = nMatch - nMis
            dm[j,i] = nMatch - nMis            

    return(dm)


def DistanceMatrixToAdjList(dm, scoreCutoff=10):
    adjList = [ [] for i in range(0,len(dm))]
    for i in range(0,len(dm)):
        for j in range(0,len(dm)):
            if (j == i):
                continue
            if (dm[i][j] >= scoreCutoff):
                adjList[i].append([j, dm[i][j]])
    return adjList


def DistanceMatrixToGraph(dm, nameLabels=None,colorLabels=None, scoreCutoff=10):
    g = nx.Graph()

    for i in range(0,dm.shape[1]):
        g.add_node(i)
    if (nameLabels is not None):
        nx.set_node_attributes(g, 'name', dict(zip(range(0,len(g)),nameLabels)))
    if (colorLabels is not None):
        nx.set_node_attributes(g, 'color', dict(zip(range(0,len(g)),colorLabels)))
    
    for i in range(0,len(dm)-1):
        for j in range(i+1,len(dm)):
            if (j == i):
                continue
            if (dm[i][j] >= scoreCutoff):
                g.add_edge(i,j,capacity=dm[i][j])
    return g


def TrimEdges(g, maxEdges):
    nRemoved =0
    for n in g:
        edgeTuples = sorted([ (e[2]['capacity'], (e[0],e[1])) for e in g.edges(n,data=True)], reverse=True)
        if (len(edgeTuples) > maxEdges):
            for ei in range(maxEdges, len(edgeTuples)):
                g.remove_edge(*edgeTuples[ei][1])
                nRemoved+=1
    print "Removed " + str(nRemoved)
    return g
        
def WriteAdjList(adjList, graphName, colorLabels=None):
    g = nx.Graph()

    for i in range(0,len(adjList)):
        if (colorLabels is not None):
            g.add_node(i, color=colorLabels[i], id=str(i))
        else:
            g.add_node(i)

    idx = 0
    for i in range(0,len(adjList)):
        for j in range(0,len(adjList[i])):
            g.add_edge(i, adjList[i][j][0], capacity=adjList[i][j][1])

    if (graphName.find("gml") >= 0):
        nx.write_gml(g, graphName)
    elif (graphName.find("gexf") >= 0):
        nx.write_gexf(g, graphName)
    elif (graphName.find("graphml") >= 0):
        nx.write_graphml(g, graphName)
    

def WriteGraph(g, graphName, colorLabels=None):

    if (graphName.find("gml") >= 0):
        nx.write_gml(g, graphName)
    elif (graphName.find("gexf") >= 0):
        nx.write_gexf(g, graphName)
    elif (graphName.find("graphml") >= 0):
        nx.write_graphml(g, graphName)


def ReadGraph(graphName):
    if (graphName.find("gml") >= 0):
        g= nx.read_gml(graphName)
    elif (graphName.find("gexf") >= 0):
        g= nx.read_gexf(graphName)
    elif (graphName.find("graphml") >= 0):
        g=nx.read_graphml(graphName)
    else:
        print "ERROR, could not determine graph format " + graphName
        sys.exit(0)
    return g
        
def PruneGraph(g, minShared=5, minOverlap=5):
    adjList = [[] for i in range(0,max(g.nodes())+1)]
    for i in g.nodes():
        adjList[i] = np.array(sorted(g[i].keys()))

    numPruned = 0
    toRemove = {}
    for i in range(0,len(adjList)-1):
        if (len(adjList[i]) > 0):
            adjListI = adjList[i]
        else:
            continue
        retain = []
        idx = 0
        for j in adjListI:
            if (len(adjList[j]) > 0):
                adjListJ = adjList[j]
                nShared = len(np.intersect1d(adjListI, adjListJ, assume_unique=True))
                if (nShared < minShared):
                    badEdge = (min(i,j), max(i,j))
                    toRemove[badEdge] = True
                    

    for e in toRemove.keys():
        g.remove_edge(*e)
            
    return g


def GetComponents(g,minSize=0):
    comps = [[c for c in comp] for comp in nx.connected_components(g)]
    comps = sorted(comps, key=len, reverse=True)

    filtComps = []
    for comp in comps:
        if (len(comp) >= minSize):
            filtComps.append(comp)
    return filtComps


def WriteGenotypeMatrix(mat, outFile):
    m = mat['mat']
    rn = mat['readNames']
    gl = mat['groupList']
    for i in range(0,len(m)):
        row = ''.join(m[i,])
        if (rn is not None):
            row += "\t" + rn[i]
        if (gl is not None):
            row += "\t" + gl[i]
        row += "\n"
        outFile.write(row)

def OpenOutput(name):
    if (name == "/dev/stdout"):
        outFile  = sys.stdout
    else:
        outFile = open(name,'w')
    return outFile
        
def WriteMETISFile(g, name):
    mf = open(name,'w')
    mf.write("{} {} {}\n".format(g.number_of_nodes(), g.number_of_edges(), 1))
    vertexEnumeration = {}
    for n in g.nodes():
        adjList = np.array(sorted(g[n].keys()))
        for ai in range(0,len(adjList)-1):
            a = adjList[ai]
            w = g[n][a]['capacity']
            mf.write("{} {} ".format(a+1,w))
        if (len(adjList) > 0):
            a = adjList[-1]
            w = g[n][a]['capacity']
            mf.write("{} {}".format(a+1,w))
        mf.write("\n")
             
    
def ParseGHFile(ghFileName):

    ghFile = open(ghFileName)
    adjList = {}
    ghTree = nx.Graph()
    nodes = {}
    for line in ghFile:
        vals = np.array([int(i) for i in line.split()])
        adjList[vals[0]] = (vals[1],vals[2])
        nodes[vals[0]]=True
        nodes[vals[1]]=True

    for n in nodes.keys():
        ghTree.add_node(n)
    for e in adjList.keys():
        ghTree.add_edge(e,adjList[e][0], cut=adjList[e][1])
    return ghTree

def PruneGHTree(gh, minCut, minComponentSize = 0):
    toRemove = []
    for e in gh.edges_iter(data=True):
        if (e[2]['cut'] < minCut):
            toRemove.append((e[0],e[1]))
    for e in toRemove:
        gh.remove_edge(*e)

    if (minComponentSize > 0):
        comps = [[ci for ci in c] for c in nx.connected_components(gh)]        
        for comp in comps:
            if (len(comp) < minComponentSize):
                for n in comp:
                    gh.remove_node(n)
    return gh

def ParseNode(v, i):
    
    if v[i] != "node":
        print "ParseNode error" + v[i] + "\t" + str(i)
        sys.exit(0)
    depth = 0
    x = "0"
    y = "0"
    i+=1
    if v[i] != "[":
        print "ParseNode error" + v[i] + "\t" + str(i)
        sys.exit(0)
    i+=1
    depth = 1
    
    while i < len(v) and depth > 0:
        if v[i] == "[":
            depth +=1
        if v[i] == "x":
            x = float(v[i+1])
        if v[i] == "y":
            y = float(v[i+1])
        if v[i] == "z":
            z = float(v[i+1])
            
        if v[i] == "]":
            depth -=1
        if v[i] == "label":
            try:
                nlabel = int(float(v[i+1].replace("\"","")))
            except ValueError:
                print "bad label value with " + v[i+1]
                nlabel = 0
                
        i+=1
    return (nlabel, x, y, z, i)



def ReadLayout(graphFileName):
    pf = open(graphFileName)
    l = pf.readlines()
    v = "\t".join(l).split()

    i = 0
    while i < len(v) and v[i] != "node":
        i+=1

    graphics = {}
    while i < len(v) and v[i] == "node":
        (nlabel,x,y,z,i) = ParseNode(v,i)
        graphics[nlabel] = { 'x': x, 'y': y, 'z': z}

    pos = { k: np.array([graphics[k]['x'], graphics[k]['y']]) for k in graphics.keys() }

    return pos



def BuildLayoutArray(g):
    if len(g.nodes()) == 0:
        return None
    if 'x' not in g.node[g.nodes()[0]]:
        return None
    return { n: np.array([g.node[n]['x'], g.node[n]['y']]) for n in g.nodes()}

def ApplyLayout(g, posGraphFileName):
    pos = ReadLayout(posGraphFileName)
    if pos is not None:
        nx.set_node_attributes(g, 'x', dict(zip(g.nodes(), [pos[n][0] for n in g.nodes()])))
        nx.set_node_attributes(g, 'y', dict(zip(g.nodes(), [pos[n][1] for n in g.nodes()])))


def GetEdgeColor(e):
    if 'color' in e:
        return e['color']
    else:
        return 0

def GetColor(n):
    if 'color' in n:
        return n['color']
    else:
        return 0

def GetEdgeColorByCost(g, ep):
    e = g[ep[0]][ep[1]]
    
    if 'cost' in e:
        if e['cost'] < 0:
            return 'red'
        else:
            return 'black'
    return 'black'

def GetEdgeWidth(g, ep):
    e = g[ep[0]][ep[1]]
    
    if 'cost' in e:
        if e['cost'] < 0:
            return 0.1
        else:
            return 1
    return 1

def DrawGraph(g, plotFileName, title=None, subset=None, labels=None):
    #
    # It is likely that the graph is a subset of pos, fix that here
    #
    pos = BuildLayoutArray(g)
    if pos is None and len(g.nodes()) > 0:
        pos = nx.spring_layout(g, k=4*1/math.sqrt(len(g.nodes())), weight=2, scale=10)

    # Try and enumerate the colors
    colors = { GetColor(g.node[i]) : True for i in g.nodes() }
    nColors = len(colors)
    ck = sorted(colors.keys())

    for i in range(0,len(ck)):
        colors[ck[i]] = i

    node_cm = [ GetColor(g.node[i]) for i in g.nodes() ]
    edge_cm = [ GetEdgeColor(g[e[0]][e[1]]) for e in g.edges()]
    edge_lw=[0.5]*len(g.edges())
    plt.plot()
    # mrv edit: I think plt figure must be before plt.title or else doesnt work, 
    # to revert switch current order
    plt.figure(figsize=(12,12))
    if title is not None:
        plt.title(title)


    numColors = len(colors)+1
    if 'NumVertexColors' in g.graph:
        numVertexColors= g.graph['NumVertexColors']

    node_sm = [ 100 for n in g.nodes()]
    nx.draw_networkx_nodes(g,pos,with_labels=False,
                           node_color=node_cm,
                           node_size=node_sm,
                           cmap=plt.cm.Set1, vmin=0,vmax=numVertexColors, alpha=0.8)
    nx.draw_networkx_edges(g,pos, edge_color=edge_cm, width=edge_lw)
    if labels is None:
        labels = { n: str(n) for n in g.nodes()}
    nx.draw_networkx_labels(g,pos,labels,font_size=10)

    plt.savefig(plotFileName)
    plt.clf()
    plt.close()



#    if drawRepulsion:
#        edge_cm = [GetEdgeColor(g, e) for e in g.edges()]
#        edge_lw = [GetEdgeWidth(g, e) for e in g.edges()]
#    else:
#        for ep in g.edges():
#            e = g2[ep[0]][ep[1]]
#            if 'cost' in e and e["cost"] < 0:
#                toRemove.append(ep)
#        g2.remove_edges_from(toRemove)
#        edge_cm=['black']*len(g2.edges())
#        edge_lw=[1]*len(g2.edges())
#        
def ClearGraphColors(g):
    for n in g.nodes():
        g.node[n]['color'] = 0
    for e in g.edges():
        g[e[0]][e[1]]['color'] = 0

def ColorGraphByCut(g, cuts):
    ClearGraphColors(g)
    cutIndex = 0
    for cut in cuts:
        for node in cut:
            g.node[node]['color'] = cutIndex
        cutIndex +=1

