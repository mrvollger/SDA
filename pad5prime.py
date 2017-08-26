#!/bin/env python
from falcon_kit.FastaReader import FastaReader
import networkx as nx

RCMAP = dict(zip("ACGTacgt","TGCAtgca"))
def rc_seq(seq):
    return "".join([RCMAP[c] for c in seq[::-1]])

def main(*argv):
    ctg_g = nx.DiGraph()
    ctg_path = {}
    with open("p_ctg_tiling_path") as f:
        for row in f:
            row = row.strip().split()
            ctg_id, v, w, edge_rid, b, e  = row[:6]
            ctg_path.setdefault( ctg_id, [] )
            ctg_path[ctg_id].append( (v, w) )
            ctg_g.add_edge(v, w)

    padding_read_ids = set()
    for ctg_id in ctg_path:
        left_end = ctg_path[ ctg_id ][0][0]
        if ctg_g.in_degree(left_end) == 0:
            left_read = left_end.split(":")[0]
            padding_read_ids.add(left_read)


    f = FastaReader("preads4falcon.fasta")
    padding_reads = {}
    for r in f:
        if r.name not in padding_read_ids:
            continue
        else:
            padding_reads[r.name] = r.sequence


    p_ctg_seq = {}
    f = FastaReader("p_ctg.fa")
    for r in f:
        p_id = r.name.split()[0]
        p_ctg_seq[p_id] = r.sequence
        left_end = ctg_path[ p_id ][0][0]
        left_read, end = left_end.split(":")
        if left_read in padding_reads:
            seq = padding_reads[left_read]
            if end == "B":
                seq = rc_seq(seq)
            print ">"+p_id+"_p"
            print seq + r.sequence
        else:
            print ">"+p_id
            print r.sequence
main()
