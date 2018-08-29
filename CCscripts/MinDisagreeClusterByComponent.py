#!/usr/bin/env python


import argparse
ap = argparse.ArgumentParser(description="Heuristic for min-disagree clustering")
ap.add_argument("--graph", help="Graph file.", required=True)
ap.add_argument("--radius", help="Radius for repulsion", type=int, default=None)
ap.add_argument("--penalty", help="Repulsion edge penalty.", type=float, default=1.0)
ap.add_argument("--factor", help="Number of repulstion nodes to add relative to neighboring", type=float,default=1)
ap.add_argument("--swap", help="Perform this many iterations of swapping nodes at cut boundaries", type=int, default=0)
ap.add_argument("--layout", help="Apply layout from this graph.", default=None)
ap.add_argument("--out", help="Output graph.")
ap.add_argument("--plot", help="Plot iterations", default=None)
ap.add_argument("--plotRepulsion", help="Plot repulsion", default=False,action='store_true')
ap.add_argument("--repulsion", help="input replusion file instead of calcualting it internally", default=None)
ap.add_argument("--niter", help="Number of clustering iterations.", type=int,default=10)
ap.add_argument("--embed", help="Stop in ipython shell.", action='store_true', default=False)
ap.add_argument("--cuts", help="Write cuts to this file.", default=None)
ap.add_argument("--scores", help="Write the different CC scores to an output file", default="CC.scores.txt")
ap.add_argument("--starts", help="Number of times to rerun CC looking for a better score", type=int, default=15)
<<<<<<< HEAD
ap.add_argument("--minlen", help="minimum length of a psv cluster", type=int, default=9000)
=======
>>>>>>> c5f142477186657859a42b2bb5fa5f5213be97cf
ap.add_argument("--sites", help="Write sites of cuts to this file.", default=None)
ap.add_argument("--seed", help="Do not seed the random number generator with a fixed value. This will cause runs on the same dataset to have different results.", default=True, action='store_false')
# added by mrv 2/21/18
ap.add_argument("--exportEachIteration", 
		help="This will export each iteration of cuts to a gml file", default=False, action='store_true')

args = ap.parse_args()


import networkx as nx
import ABPUtils
import itertools
import numpy as np
import pdb
from sets import Set
import copy
import random
from multiprocessing import Pool
import operator
import math

if args.seed is True:
    random.seed(3261827) #RIP LVB (an inside joke)


g = ABPUtils.ReadGraph(args.graph)


sites = [ (int(g.node[n]['pos']), n) for n in g.node ]
sites = sorted(sites)

i = 0
nRepulsion = 0
nOK = 0
#mrv addition: i jsut need a golbal counter for the number of exports in the case of exportEachIteration
exportCounter = 0
cutCounter = -1


# mrv note: determines weather to add a repulsion nodes
def StoreRepulsionAdjacency(g, n, pos, radius=None):
    repulsion = Set()
    if radius is None:
        neighbors = g[n].keys()        
        start = min([pos[node] for node in neighbors])
        end   = max([pos[node] for node in neighbors])
        
    for i in g.nodes():
        addRepulsion = False
        if radius is None:
            if pos[i] > start and  pos[i] < end and i not in g[n] and i != n:
                addRepulsion = True
        else:

            dist = abs(int(pos[i]) - int( pos[n]))
            if dist < radius and i not in g[n]:
                if i != n:
                    addRepulsion = True
        if addRepulsion:
            repulsion.update([i])
    return repulsion


#
# read repulsion from file, this is an mrv addition  
#
def ReadInRepulsion(g):
	repfile = open(args.repulsion).readlines()
	repulsion = { node:Set() for node in g.nodes() }
	site = { int(g.node[i]['pos']):i for i in g.nodes()}
	for line in repfile:
		line = line.strip().split()
		i = int(line[0])
		j = int(line[1])
		if(i in site and j in site):
			i = site[i]
			j = site[j]
			repulsion[i].add(j)
			repulsion[j].add(i)
	return(repulsion)

# mrv note: addes repulsion edges for each node useing storerepulsionadj.. 
# and then does a check to make sure that repulsion edges are reflected in both target and source node
def StoreRepulsion(g, radius=None):
    # mrv addition, I think comaprision of strings of ints might be a problem here, so I added the int bit
    pos = {i:int(g.node[i]['pos']) for i in g.nodes()}
    if(args.repulsion is None): # marks original 
        repulsion = { i: StoreRepulsionAdjacency(g,i,pos,radius) for i in g.nodes() }
    else: # mrv addition
        print("READING IN REPULSION")
        repulsion = ReadInRepulsion(g)
    nUpdated = 0
    # Make sure the repulsion sets are reflexive
    for src in repulsion:
        for dest in repulsion[src]:
            if src not in repulsion[dest]:
                repulsion[dest].update([src])
                nUpdated+=1
    return repulsion

# mrv note: appears to get rid of repulsion edges if there are more 
# repulsion edges than positve edges times a factor. 
# how it picks which ones to remove I am not sure, 
# it appears to sort the repulsion list and then only pick the top x
def SubsampleRepulsion(gAdjList, repulsion, factor=1.0):
    nReduced = 0
    for n in gAdjList.keys():
        if len(gAdjList[n]) *factor < len(repulsion[n]):
            toRemove = int(len(repulsion[n]) - len(gAdjList[n])*factor)
            repList = sorted(list(repulsion[n]))
            repSub  = Set([])
            fraction= len(repulsion[n])/(len(gAdjList[n])*factor)
            for i in range(0,len(repulsion[n])):
                if len(repSub) * fraction < i:
                    repSub.update(Set([repList[i]]))
                else:
                    nReduced+=1
            repulsion[n] = repSub

    return nReduced

def NeighborSimilarity(gAdjList,src,dest):
    return len(gAdjList[src].intersection(gAdjList[dest]))


def GetCost(g, i, j, r, m):
    if j in g[i]:
        return 1
    else:
        if j in m[i]:
            return -1
        else:
            return 0
    
allNodeSet = Set(g.nodes())



def CountRepulsion(cut, nodes, repulsion):
    n=0
    for node in nodes:
        n+= len(repulsion[node].intersection(cut))
    return n


def CountEdgesOutside(gadj, cut):
    #
    # return the number of edges between nodes in a and b
    #
    
    n=0
    for src in cut:
        n += len(gadj[src].difference(cut))
    return n
        

def ScoreCut(gadj, cut, repulsion):
    #
    # Count repulsion edges inside the cut.
    #
    nNeg = CountRepulsion(cut, cut, repulsion)/2

    #
    # Count edges from the cut to other nodes.
    #
    nPos = CountEdgesOutside(gadj, cut)
    return nNeg*args.penalty + nPos

def Neighboring(g, subgraph, n):
    neighbors = Set(sorted(g[n].keys()))
    return list(neighbors.difference(subgraph))

def GrowSubgraph(g, subgraph, neighbors, n):
    subgraph.update(n)
    for i in n:
        neighbors.update(g[i].keys())
    neighbors.difference_update(subgraph)
    neighbors.difference_update(n)


def GetAdjacent(gadj, subgraph, exclude=None):
    adjacent = Set()
    for n in subgraph:
        #
        # Possibly faster method to check for neighbors.
        #
        adjacent.update(gadj[n])

    if exclude is not None:
        adjacent.difference_update(exclude)

    # Make sure adjacent is not self referencing
    adjacent.difference_update(subgraph)
    return adjacent

def CountUnion(a,b):
    if (len(a) < len(b)):
        smaller = a
        larger = b
    else:
        smaller = b
        larger = a
    n = 0
    for i in smaller:
        if i in larger:
            n+=1
    return n


def ScoreCutExpansion(g, gadj, cut, expansion, repulsion):

    #
    # Thrifty method to compute new score of cut when nodes from
    # expansion are added to cut.
    #
    numEdgesCrossingCut = 0
    # mrv note: I think this counts the num of neighbors to the cut and the expanshion that are not
    # in the cut or the expansion, this is the num of positive edges continubuting to the score I think
    for src in expansion:
        neighbors = Set(g[src].keys())
        neighbors.difference_update(cut)
        neighbors.difference_update(expansion)
        numEdgesCrossingCut += len(neighbors)
    # mrv bug: is this double counting? should this not be divided by two since you will count u,v and v,u 
    numEdgesInExpansion = 0
    for src in expansion:
        numEdgesInExpansion += len(gadj[src].intersection(expansion))

    # 
    # 'black' edges that were previously from the cut to the
    # expansion.  The score will decrease by this.
    #

    numEdgesBetweenExpansionAndCut = 0
    for src in expansion:
        for dest in gadj[src]:
            if dest in cut:
                numEdgesBetweenExpansionAndCut += 1

    #
    # 'red' edges that were outside the cut but are now inside.
    #
    numNewIncludedRepulsion = CountRepulsion(cut, expansion, repulsion)
    # mrv bug note: should we not also be counting the number of repulsion edges within the cut itself?
	# something like  numNewIncludedRepulsion += CountRepulsion(expansion, expansion, repulsion)/2

#    print "  Score cut expansion exp-tot {}\texp-repl {}\texp-cut {}\texp-exp{}".format(numEdgesCrossingCut , numNewIncludedRepulsion ,  numEdgesBetweenExpansionAndCut , numEdgesInExpansion)

    return numEdgesCrossingCut + numNewIncludedRepulsion*args.penalty - numEdgesBetweenExpansionAndCut - numEdgesInExpansion

def IsolatedCutScore(gadj, testCut):
    n=0
    for node in testCut:
        for dest in gadj[node]:
            if dest > node:
                n+=1
    return n


class GraphOutput:
    def __init__(self, fileName):
        self.fileName = fileName
        self.vertexColor = 1
        self.iter = 0

def TestCutExpansion(g, gAdjList, cut, repulsion, expandedCut):
    curScore = ScoreCut(gAdjList, expandedCut, repulsion)
    #
    # Find the cost of all of the neighbors in what will be expanded.
    #
    isolatedTestCutScore = ScoreCut(gAdjList, expandedCut, repulsion)

    #
    # Try and come up with a thrifty method to compute expansion.
    #
    cutExpansion = ScoreCutExpansion(g, gAdjList, cut, expandedCut, repulsion)
    newCutScore = curScore + cutExpansion
    return curScore, newCutScore

    

def GrowCut(g, cut, repulsion, minNeighborSimilarity=3):
    # mrv additions 
    global cutCounter, exportCounter 
    cutCounter += 1
    exportCounter = 0
	
    neighbors = GetAdjacent(g, cut)
    grewCut = True
    iter = 0
    nodes = g.nodes()
    adjList = g.adjacency_list()
    gAdjList = { nodes[i]: Set(adjList[i]) for i in range(0,len(nodes)) }
    it = 0

    curScore = ScoreCut(gAdjList, cut, repulsion)
    while grewCut:
        grewCut = False
        nSearched = 0
        searched = []
        neighbors = Set([])
        cutNodes = list(cut)
        random.shuffle(cutNodes)
        newCutScore = 0
        iter +=1
        for cutNode in cutNodes:
            cutNeighbors = gAdjList[cutNode].difference(cut)
            for n in cutNeighbors:
				
                if NeighborSimilarity(gAdjList, cutNode, n) < minNeighborSimilarity:
                    continue
                # 
                # Define the nodes that will be expanded
                #
                expandedCut = Set([n])# #O(n)
            
                #
                # When the expansion is a search, e.g. Set([n]+ gAdjList[n]), remove overlap with current cut
                #
                expandedCut.difference_update(cut) # O(n)
        
                #
                # Find the cost of all of the neighbors in what will be expanded.
                #
                isolatedTestCutScore =  ScoreCut(gAdjList, expandedCut, repulsion)

                #
				# Try and come up with a thrifty method to compute expansion.
                #
                cutExpansion = ScoreCutExpansion(g, gAdjList, cut, expandedCut, repulsion)

            
                #
                # Get the score of the entire new test cut.  Ideally this is equal to the cut expansion
                #testCut = expandedCut.union(cut)
                #newCutScore = ScoreCut(gAdjList, testCut, repulsion) 
                newCutScore = curScore + cutExpansion

                #
                # The benefit of adding this cut is equal to the new test cut minus
                # all the nodes in isolation. I'm not sure if this is a good thing --
                # maybe the cost of the cut should be equal to the assumption that 
                # all other nodes are in the same cut.
                # 
                # mrv note: maybe this is the progblem we are having in 000028F.13115500.13158799. 
                # this could add to very disjoint sets connected by one node that have a lot of repulsion
                # maybe instead check if more nodes are added than negative edges added within the cut.
                newCutScoreDiff = newCutScore - isolatedTestCutScore

                it += 1
#                print str((newCutScore, isolatedTestCutScore, newCutScore - isolatedTestCutScore , curScore))
                if  newCutScore - isolatedTestCutScore < curScore:
                    GrowSubgraph(g, cut, neighbors, expandedCut)
                    grewCut = True
                    curScore = newCutScore
                    break
                nSearched+=1
            # mrv bug, pretty sure this was meant to be indented one extra level
            # but I dont think searched is used for anything so it does not matter
                searched.append(str(n))
        #mrv addition, export each cut growth so I cam watch them grow 
        if(args.exportEachIteration):
            exportIteration(g, cut, repulsion)
    return cut

class Input:
    def __init__(self,g,it,repulsion):
        self.g = g
        self.it = it
        self.repulsion = repulsion



# mrv note: kind of the first step in CC
# calls Grow Cut after picking a random unclustered node to start the cut. 
def SampleCuts(g, nIter, repulsion):
    clustered = Set([])
    cuts = []
    it = 0
    g.graph['NumVertexColors'] = nIter+1
    unclustered = None
    while (len(clustered) < len(g.nodes()) and it < nIter):
        it +=1
        unclustered = list(Set(g.nodes()).difference(clustered))
        if len(unclustered) > 0:
            seed = unclustered[random.randint(0,len(unclustered)-1)]
            cut = Set([seed])
            GrowCut(g, cut, repulsion)
            clustered.update(cut)
            cuts.append(cut)
            #print "iter: " + str(it) + "\tseed: " + str(seed) + "\t" + str(len(cuts[-1])) + "\t" + str(len(clustered)) + "/" + str(len(g.nodes()))
        else:
            break
    return (cuts, unclustered)
        
def SampleCutsIt(cl):
    return SampleCuts(cl.g, cl.it, cl.repulsion)

def ReciprocalOverlap(a,b,f):
    if f == 0:
        return True
    ovp = len(a.intersection(b))
    if ovp == 0:
        return False
    else:
        return float(ovp) / len(a) > f and float(ovp) / len(b) > f
    
# mrv note: makes sense what this does, but I am not sure why it is done. 
# if it ever does something it will presumably be against CC? 
def MergeCuts(cuts, ratio=0.9):
    #
    # Do a greedy merge of cuts.
    #
    i = 0
    cutMerged = True
    while cutMerged:
        cutMerged = False
        while i < len(cuts) - 1:
            j = i + 1
            while j < len(cuts):
                if ReciprocalOverlap(cuts[i], cuts[j], ratio):
                    cuts[i].update(cuts[j])
                    del cuts[j]
                    cutMerged = True
                else:
                    j+=1
            i+=1


# Now try swapping nodes that border two cuts

def GetCut(n, cuts):
    for i in range(0,len(cuts)):
        if n in cuts[i]:
            return i
    return None

#
# Check overlap between a set of nodes and all other cuts
#
def CountCutMembership(nodeSet, cuts):
    membership = { i: len(nodeSet.intersection(cuts[i])) for i in range(0,len(cuts)) }
    for cutIndex in membership.keys():
        if membership[cutIndex] == 0:
            del membership[cutIndex]
    return membership



def GetCutIndices(cuts):
    return  { j : i for i in range(0,len(cuts)) for j in cuts[i] }

def TestCutDeletion(adjList, repulsion, cuts, cutIndex, node):
    cuts[cutIndex].difference_update(Set([node]))            
    score = ScoreCut(adjList, cuts[cutIndex], repulsion)
    cuts[cutIndex].update(Set([node]))
    return score

def PickCutDeletion(gAdjList, repulsion, cuts, cutIndicesDict, node):
    cutIndices = sorted(cutIndicesDict.keys())
    scores = [TestCutDeletion(gAdjList, repulsion, cuts, i, node) for i in cutIndices]
    lowestScore = min(scores)
    lowestScore = None
    for i in range(0,len(scores)):
        if lowestScore == None or scores[i] < lowestScore:
            lowestScore = scores[i]
            minScoreIndex = i

    for i in range(0,minScoreIndex):
        
        cuts[cutIndices[i]].difference_update(Set([node]))
    for i in range(minScoreIndex+1, len(cutIndices)):
        cuts[cutIndices[i]].difference_update(Set([node]))

# mrv note: for every node counts how many connections is mas to each cut with count cut membership (dict)
# then if it has membership to more than one group it does pickCutDeletion
# which finds the best scoring cut for the node to be part of and then removes it form all other cuts
def AssignNodesToUniqueClusters(gAdjList, repulsion, cuts):

    for node in gAdjList.keys():
        membership = CountCutMembership(Set([node]), cuts)
        if len(membership) > 1:
            PickCutDeletion(gAdjList, repulsion, cuts, membership, node)

def TestCutSwap(adjList, repulsion, cuts, cutScores, a, b, swap):
    cuts[a].difference_update(swap)
    cuts[b].update(swap)
    scoreA = ScoreCut(adjList, cuts[a], repulsion)
    scoreB = ScoreCut(adjList, cuts[b], repulsion)
    #
    # Put things back where they came from
    #
    cuts[a].update(swap)
    cuts[b].difference_update(swap)
    return (cutScores[a] + cutScores[b] - scoreA - scoreB, scoreA, scoreB)

def OptimizeBySwappingNodes(gAdjList, cuts, repulsion, nIter=0, minGain=10, neighborSimilarityCutoff=3):

    cutScores = [ScoreCut(gAdjList, cut, repulsion) for cut in cuts]    
    cutIndices = GetCutIndices(cuts)
    gain = 0
    swapMade = True
    iterIndex = 0
    while  (nIter == 0 or iterIndex < nIter) and swapMade:
        swapMade = False
        iterIndex+=1

        for cutIt in range(0,len(cuts)):
            cut = cuts[cutIt]
            neighbors = ABPUtils.GetCutNeighbors(gAdjList, cut)
            maxNeighbor = None
            maxGain     = 0
            maxNeighborCut = 0
        
            for node in neighbors:
                if node not in cutIndices:
                    continue
                maxNeighborSimilarity = 0
                for cutNode in cut:
                    if node in gAdjList[cutNode]:
                        neighborSimilarity= NeighborSimilarity(gAdjList,cutNode, node)
                        maxNeighborSimilarity=max( neighborSimilarity , maxNeighborSimilarity)
                if maxNeighborSimilarity < neighborSimilarityCutoff:
                    continue
                            
                neighborCut = cutIndices[node]
                
                (swapScore, score0, score1) = TestCutSwap(gAdjList,
                                                          repulsion,
                                                          cuts,
                                                          cutScores,
                                                          neighborCut,
                                                          cutIt,
                                                          Set([node]))
                if swapScore > minGain:
                    swapMade = True
                    cutScores[neighborCut]  = score0
                    cutScores[cutIt]        = score1
                    cuts[cutIt].update(Set([node]))
                    cuts[neighborCut].difference_update(Set([node]))
                    cutIndices[node] = cutIt
                    gain += swapScore

    return gain

def FilterConflictingNodes(gAdjList, repulsion, cuts, nodes, cutoff):
    cutScores = [ScoreCut(gAdjList, cut, repulsion) for cut in cuts]    
    cutIndices = GetCutIndices(cuts)
    gain = 0
    for node in nodes:
        if node not in cutIndices:
            continue
        neighborMembers = CountCutMembership(gAdjList[node], cuts)
        if len(neighborMembers) > 1:
            if min(neighborMembers.values()) > cutoff:
                cuts[cutIndices[node]].difference_update(Set([node]))

<<<<<<< HEAD
# this function removes small clusters (less than 3), or cuts that are really samll in length (9kbp)
# this is because it could be a LINE element, and it is also the minimum default collapse size
def RemoveSmallCuts(cuts, minCutSize=3):
	i = 0
	while i < len(cuts):
		sites = sorted( [int(g.node[n]['pos']) for n in cuts[i]] )
		minpos = min(sites)
		maxpos = max(sites)
		if( (len(cuts[i]) < minCutSize)  or ( (maxpos - minpos) < args.minlen) ):
			cuts.remove(cuts[i])
		else:
			i+=1
=======
def RemoveSmallCuts(cuts, minCutSize=3):
    i = 0
    while i < len(cuts):
        if len(cuts[i]) < minCutSize:
            cuts.remove(cuts[i])
        else:
            i+=1
>>>>>>> c5f142477186657859a42b2bb5fa5f5213be97cf
    

def AddRepulsionEdges(g, repulsion):
    for src in repulsion.keys():
        for dest in repulsion[src]:
            g.add_edge(src, dest, weight=0.1, color=1)

# mrv function: this is a function that will allow me to export each iteration of the graph
def exportIteration(doNotTouchGraph, cut, repulsion):
	global exportCounter
	# need to make a copy so I do not mess up the node coloring for the rest of the program
	g = copy.deepcopy(doNotTouchGraph) # still not working right
	# need adj list for scoring
	adjList = g.adjacency_list()
	nodes = g.nodes()
	gAdjList = { nodes[i]: Set(adjList[i]) for i in range(0,len(nodes)) }

	#need to add colors based on cut
	for n in g.nodes():
		g.node[n]['color'] = 'b'
	for n in cut:
		g.node[n]['color'] = 'r'
	print(len(nodes), len(cut))
	# need to rename edge colors 
	for (u, v) in g.edges():
		g[u][v]['color'] = 0
	
	
	nNeg = CountRepulsion(cut, cut, repulsion)/2
	nPos = CountEdgesOutside(gAdjList, cut)
	# this is the same score calculation from ScoreCut, I just want both the pos and neg
	score = nNeg * args.penalty + nPos 

	# cut score, see if it is imporving
	info = " Score: {}, PosScore: {}, NegScore: {}, Cut Size: {} ".format( score, nPos, nNeg, len(cut) )
	
	title = "iteration.{cut:03d}.{it:04d}".format(cut = cutCounter, it=exportCounter)
	print(title)
	ABPUtils.DrawGraph(g, "extraCCplots/" + title + ".png", title=title+info)
	
	# plot relevant repulsion edges 
	AddRepulsionEdges(g, repulsion)

	# remove irrelevant edges
	for (u,v) in g.edges():
		edge = g[u][v]
		isRep = ((u in cut) or (v in cut)) and (edge['color'] == 1)
		# if it is not a repulsion node touching our cut remove it
		if(not isRep):
			g.remove_edge(u,v)
	
	ABPUtils.DrawGraph(g, "extraCCplots/" + title + ".rep.png", title=title+info)

	exportCounter += 1
	return 



# mrv addition,
# this check all paris of cuts, and if a pair of cuts would have a better score if it existed as only one cut
# it collapese into one cut
def MergeCutsCC(g, cuts, gAdjList, repulsion):
	cutlen = len(cuts)
	if(cutlen <= 1):
		return
	diffs = {}
	#score all the cuts, this is some simple memoization
	scores = { key:ScoreCut(gAdjList, cut, repulsion) for key, cut in enumerate(cuts) }
	neighs = { key:GetAdjacent(g,cut) for key, cut in enumerate(cuts) }
	
	# check the pairs
	idxs = range(cutlen)
	for idx1, idx2 in itertools.combinations(idxs ,2):
		# check if these two clusters are touching eachother 
		neigh1 = neighs[idx2]
		neigh2 = neighs[idx1]
		inter1 = len(neigh1.intersection(cuts[idx2]))
		inter2 = len(neigh2.intersection(cuts[idx1]))
		if( True):#inter1 > 0 or inter2 > 0 ):
			newCutScore = ScoreCut(gAdjList, cuts[idx1].union(cuts[idx2]), repulsion)
			diffs[ newCutScore - scores[idx1] - scores[idx2] ] = (idx1, idx2)  
	
	if(len(diffs) > 0 ):
		best = min(diffs.keys())
		if(best < 0 ):
			idx1 = diffs[best][0]
			idx2 = diffs[best][1]
			cuts[idx1] = cuts[idx1].union(cuts[idx2])
			del cuts[idx2]
			# calculate overall score 
			totalScore = sum( [ScoreCut(gAdjList, cut, repulsion) for cut in cuts ])
			#print(idx1, idx2, best, len(cuts), totalScore)
	
	# if we have merged a cut, recursivly call on the new set. 
	if(len(cuts) != cutlen ):
		MergeCutsCC(g, cuts, gAdjList, repulsion)
	return


# mrv addition, add positions based on mi.gml 
def AddLayoutPreCC(g):
	pos = ABPUtils.BuildLayoutArray(g)
	if pos is None and len(g.nodes()) > 0:
		pos = nx.spring_layout(g, k=1/math.sqrt(len(g.nodes())), 
				weight=math.sqrt(len(g.nodes())), scale=100, iterations=1000)

	if pos is not None:
		nx.set_node_attributes(g, 'x', dict(zip(g.nodes(), [pos[n][0] for n in g.nodes()])))
		nx.set_node_attributes(g, 'y', dict(zip(g.nodes(), [pos[n][1] for n in g.nodes()])))
	


starts = args.starts
scores = {}
scoreVals = []
for idx in range(starts):
	# mrv note: the 4 means the compenent must have at least 4 nodes
	components = ABPUtils.GetComponents(g,4)
	repulsionTotal = {}
	totalScore = 0
	allCuts = []

	# mrv note: runs ABP for each unconnected compenent in the graph. 
	for comp in components:
		sub = g.subgraph(comp)

		adjList = sub.adjacency_list()
		nodes = list(sub.nodes())
		#
		# Store the adjacency list ofr every node.
		#
		gAdjList = { n: Set(sub[n].keys()) for n in sub.nodes() }

		repulsion = StoreRepulsion(sub, args.radius)
		#nReduced = SubsampleRepulsion(gAdjList, repulsion, factor=args.factor)
		#print nReduced
		repulsionTotal.update( repulsion )
		
		# mrv addition, adding the layout before running cc 
		if(args.exportEachIteration):
			AddLayoutPreCC(sub)

		step = Input(sub, args.niter, repulsion)
		(cuts, unclustered) = SampleCutsIt(step)
		if args.layout:
			ABPUtils.ApplyLayout(sub,args.layout)

		MergeCuts(cuts,0.75)
		AssignNodesToUniqueClusters(gAdjList, repulsion, cuts)

		cutLengths = [len(c) for c in cuts]
		# mrv addition, merge cuts that get better if they are merged
		print "before/after"
		cutLengths = [len(c) for c in cuts]; print cutLengths

		MergeCutsCC(sub, cuts, gAdjList, repulsion)
		cutLengths = [len(c) for c in cuts]; print cutLengths

		RemoveSmallCuts(cuts)
		cutLengths = [len(c) for c in cuts]; print cutLengths

		if args.swap > 0:
			OptimizeBySwappingNodes(gAdjList, cuts, repulsion, args.swap)
		cutLengths = [len(c) for c in cuts]; print cutLengths

		RemoveSmallCuts(cuts)            
		cutLengths = [len(c) for c in cuts]; print cutLengths

		allCuts += cuts
	
		if(len(cuts)>0):
			for cut in cuts:
				totalScore += ScoreCut(gAdjList, cut, repulsion)
	print(totalScore)
	scoreVals.append(totalScore)
	scores[ totalScore ] = allCuts

print( scores.keys() )

#### WRITE SCORES TO OUTPUT FILE ###
f = open(args.scores, "w+")
for idx, tmpscore in enumerate(scoreVals):
	f.write( "{}\t{}\t{}\n".format(idx+1, tmpscore, tmpscore - min(scoreVals) ) )
f.close()
######


#------------------------------------------------------------------------------------------#
cuts = scores[ min( scores.keys() ) ]
print("------------------")
print(len(cuts))
print("------------------")
ABPUtils.ColorGraphByCut(g,cuts)

if args.plotRepulsion:
    if(repulsionTotal): # checks to make sure it is not empty, and then adds the total of repulsion, mvr addition
        AddRepulsionEdges(g, repulsionTotal)
    

if args.plot is not None:
    pass

if args.cuts is not None:
    cutsFile = open(args.cuts,'w')
    for cut in cuts:
        cutsFile.write("\t".join([str(g.node[c]["index"]) for c in sorted(list(cut))]) + "\n")
    cutsFile.close()

if args.sites is not None:
<<<<<<< HEAD
	sitesFile = open(args.sites, 'w')
	for cut in cuts:
		sites = sorted([int(g.node[n]['pos']) for n in cut])
		sitesFile.write("\t".join(str(s) for s in sites ) + "\n")
=======
    sitesFile = open(args.sites, 'w')
    for cut in cuts:
	sites = sorted([int(g.node[n]['pos']) for n in cut])
        sitesFile.write("\t".join(str(s) for s in sites ) + "\n")
>>>>>>> c5f142477186657859a42b2bb5fa5f5213be97cf
    
if args.out is not None:
    allCuts = Set()
    for c in cuts:
        allCuts.update(c)
    notCut = Set(g.nodes()) - allCuts
    notCutList = sorted(list(notCut))
    g.remove_nodes_from(notCutList)
    for n in g.nodes():
        g.node[n]['color'] = -1
        
    for i in range(0,len(cuts)):
        for c in cuts[i]: 
            g.node[c]['color']= i
	# mrv note: getting weird error looking for source: 
	# nx.write_gml(g, graphName)raise NetworkXError('%r is not a string' % (value,))
    # only happens sometimes, when using the export each iteration option. Not a big problem. 
    print(args.out)
    for cut in cuts: print(cut)
    ABPUtils.WriteGraph(g, args.out)


if args.embed:
    IPython.embed()

#AddRepulsionEdges(g,repulsion)
#ABPUtils.DrawGraph(g, "test.repl.png")

