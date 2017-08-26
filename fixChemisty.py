#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('inputFile', nargs='?', help="autodetects plus does regex, but the regex requires quotes around it", default=None)
parser.add_argument('inputFile2', nargs='?', help="autodetects plus does regex, but the regex requires quotes around it", default=None)
#parser.add_argument('--selfMap', help="map the bac assmbly to itsefl", action='store_true')
#parser.set_defaults(selfMap=False)
args = parser.parse_args()

h51=args.inputFile
h52=args.inputFile2

# other imports, a little slow 
import glob
import subprocess
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import h5py
import numpy as np
# global varibles

def h5Read(h5): 
    h5file = h5py.File(h5, 'r')
    # top level keys: PulseData ScanData
    # Scan data does not have anything useful for us
    # next level BaseCalls Regions 
    # pretty sure we want the dir PulseData/Regions
    #for key, value in h5file.iteritems():
    #    print(key)
    for key, value in h5file["ScanData"].iteritems():
        print(key, value)
    
    for key, value in h5file["ScanData"]["DyeSet"].iteritems():
        print(key,value)
    print(h5file["ScanData"]["RunInfo"])
    
    exit()
   
    data = h5file["PulseData"]["Regions"] 
    # creates a copy
    data = data[:] 
    h5file.close()
    trimmedInserts = holeGroup(data)
    names = []
    lengths = []
    start = 2; end = 3
    for insert in trimmedInserts:
        length = insert[end] - insert[start]
        if( length > 0 ):
            lengths.append(length)
            names.append(insert[0])
    printLengths(names, lengths)
    
    '''  # the slow lame way
    tempFasta = "ForReadLengthTemp.fasta"
    runPls2Fasta(h5, tempFasta)
    names, lengths = fastaLengths(tempFasta)
    printLengths(names, lengths)
    os.system("rm ForReadLengthTemp.fasta")
    '''

def main():
    h5Read(h51)

main()
