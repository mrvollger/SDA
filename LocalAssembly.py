#!/usr/bin/env python
import glob
import subprocess
import argparse
import os
import sys
import re
import itertools 
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO


def runCmd(cmd):
    #print("starting", cmd)
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    err = err.decode("utf-8")
    out = out.decode("utf-8")
    if(len(err) > 1):
        print(cmd)
        print(err)
    #print("ending cmd", err, out)
    return(out)

def existNotZero(myfile, prefix=""):
    my_path = os.path.join(prefix, myfile)
    mybool = os.path.exists(my_path) and os.path.getsize(my_path) > 0
    return mybool

class LocalAssembly:

    def __init__(self, mydir, psvGraphLoc=None, psvURLPATH=None):
        self.mydir = os.path.abspath(mydir.strip())
        self.psvGraphLoc = psvGraphLoc
        self.psvURLPATH = psvURLPATH 
        self.basename()
        self.mi_gml_png()
        self.readRef()
        self.readAssemblies()
        self.numOfPSVs()
        self.samReads()
        self.readIdentity()

    def basename(self):
        self.basename = self.mydir.split("/")[-1].strip()
        print(self.basename)

    def mi_gml_png(self):
        # skip if not specified 
        if(self.psvGraphLoc == None or self.psvURLPATH == None):
            self.psvURL = "None"
            return

        self.psvGraph= os.path.join(self.psvGraphLoc, self.basename + ".pdf")
        self.orginPNG = os.path.abspath(os.path.join(self.mydir, "mi.cuts.gml.pdf"))
        cmd = "ln -sf {} {}".format(self.orginPNG, self.psvGraph)
        self.psvURL= self.psvURLPATH + "/" + self.basename + ".pdf"
        #print(cmd)
        runCmd(cmd)
    

    def samReads(self):
        self.num_sam_reads = []
        for samFile in glob.glob( os.path.join(self.mydir,"group.*/H2.WH.sam") ):
            if(os.path.exists(samFile)):
                cmd = 'grep -c "^[^@]" ' + samFile
                self.num_sam_reads.append( int(runCmd(cmd)) )

    def numOfPSVs(self):
        self.PSVs = []
        mi_gml_cuts = os.path.join(self.mydir, "mi.gml.cuts")
        if(not os.path.exists(mi_gml_cuts)):
            return 

        f = open(mi_gml_cuts)
        for line in f:
            temp=len(line.strip().split())
            self.PSVs.append(temp)
    

    def readIntersect(self, isc_path):
        f = open(isc_path)
        rtn = {}
        for line in f:
            line = line.strip().split("\t")
            if(len(line)==0):
                continue
            key = "{}:{}-{}".format(line[0], line[1], line[2])
            if(key not in rtn):
                rtn[key] = []
            rtn[key].append(line[6])
        return rtn 

    def readRef(self):
        #
        my_max = os.path.join(self.mydir, "dup_max.intersect")
        my_mean = os.path.join(self.mydir, "dup_mean.intersect")
        self.max = {}
        self.mean = {}
        if(os.path.exists(my_max) and os.path.exists(my_mean) ):
            self.max  = self.readIntersect(my_max) 
            self.mean = self.readIntersect(my_mean)

        # read in fasta files for ref regions
        dups = os.path.join(self.mydir, "duplications.fasta")
        if(os.path.exists(dups)): 
            self.refs = list(SeqIO.parse(dups, "fasta"))
        else:
            self.refs = []
        

        # turn ref regions and seg dup identities into text
        self.dupNumber = len(self.refs)
        self.refRegions = ""
        for ref in self.refs:
            key = ref.description.strip() 
            des = key
            if(key in self.max):
                des += "[" + ",".join(self.max[key]) + "]"
            if(key in self.mean):
                des += "[" + ",".join(self.mean[key])+ "]"
            
            self.refRegions += des + ";"
        
        self.refRegions = self.refRegions[:-1]
 

    def readIdentity(self):
        self.bestID = []
        self.bestMatch = []
        self.averageID = []
        self.bestLength = []
        self.averageLength = []
        for group in glob.glob(os.path.join(self.mydir, "group.*/")):
            my_path1 = os.path.join(group, "WH.average.m5")
            my_path2 = os.path.join(group, "WH.best.m5")
            if(os.path.exists(os.path.join(self.mydir, "duplications.fasta")) and
                    os.path.exists(my_path1) and os.path.getsize(my_path1) > 0 and
                    os.path.exists(my_path2) and os.path.getsize(my_path2) > 0 ):
                #print("PCT ID Files exist")
                #average = pd.read_csv( my_path1, sep=" ", header = 0)
                average = pd.read_csv( my_path1, sep=" ", header=None)
                best =    pd.read_csv( my_path2, sep=" ", header=None)
                print(best)
                # skip if there are no alingments 
                if( len(average) < 1 or len(best) < 1):
                    #print("too short")
                    self.bestID.append("NA")
                    self.averageID.append("NA")
                    self.bestMatch.append("NA")
                    self.bestLength.append("NA")
                    self.averageLength.append("NA")
                else:
                    #print("alignments exist")
                    average["pctsimilarity"] = 100*average[11]/(average[11] + average[12])
                    self.averageID.append(round(average["pctsimilarity"].mean(), 2))
                    aveLen = int(abs(average[2] - average[3]).mean())
                    self.averageLength.append(aveLen) 
                    
                    best["pctsimilarity"] = 100*best[11]/(best[11] + best[12])
                    pctID = round(best["pctsimilarity"].iloc[0],2)
                    group = best[0].iloc[0]
                    region = best[5].iloc[0]
                    start = best[7].iloc[0]
                    end = best[8].iloc[0]
                    length = abs(end - start)
                    self.bestID.append(pctID)
                    self.bestMatch.append(region)
                    self.bestLength.append(length)
            else:
                self.bestID.append("NA")
                self.averageID.append("NA")
                self.bestMatch.append("NA")
                self.bestLength.append("NA")
                self.averageLength.append("NA")

    def readAssemblies(self):
        self.ASMs = []
        self.lens = []
        self.reads = []
        self.groups = []
        self.status = []
        self.asmNumber=0
        for group in glob.glob( os.path.join(self.mydir,"group.*/") ):
            asmFile = os.path.join(group, "WH.assembly.fasta")
            groupID = group.split("/")[-2]
            groupID = int(groupID.split(".")[1])
            self.groups.append(groupID)
            if(os.path.exists(asmFile)):
                asm = list(SeqIO.parse(asmFile, "fasta"))
                if(len(asm) == 1):
                    asm = asm[0]
                    self.status.append("Success")
                    self.ASMs.append(asm.description)
                    self.lens.append(asm.description.split(" ")[1] )
                    self.reads.append(asm.description.split(" ")[2] )
                    self.asmNumber += 1
                else:
                    self.status.append("MultipleAsm")
                    self.ASMs.append("NA")
                    self.lens.append("NA")
                    self.reads.append("NA")
            else:
                self.status.append("Failed")
                self.ASMs.append("NA")
                self.lens.append("NA")
                self.reads.append("NA")

        self.CCNumber = len(self.ASMs)
        return

    # output stuff
    #tableFormat = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
 
    def getHeader(self):
        headerList =["region_in_falcon", "copies_in_reference","number_of_psv_assemblies",
                "number_of_CC_groups", 
                "CC_ID", "Status","numOfPSVs","length", "num_of_asm_reads",
                "num_of_sam_reads","psvURL",
                "averagePerID", "averageLength", 
                "bestPerID", "bestLength", "BestRegionInTheHumanReference", 
                "refRegions"]
        header = ""
        for item in headerList:
            header += item + "\t"
        header = header[:-1] + "\n"
        return(header)

    def __str__(self):
        rtn = ""
        debug = False
        if(debug):
            for key in vars(self):
                rtn += str(key) + ":"  + str(vars(self)[key]) + "\n"
        else:
            for idx, asm in enumerate(self.ASMs):
                psvs = self.PSVs[self.groups[idx]]
                rtnList = [self.basename, self.dupNumber, self.asmNumber, 
                        self.CCNumber, 
                        self.groups[idx], self.status[idx], psvs, self.lens[idx], self.reads[idx], 
                        self.num_sam_reads[idx], self.psvURL,
                        self.averageID[idx], self.averageLength[idx], 
                        self.bestID[idx], self.bestLength[idx], self.bestMatch[idx], 
                        self.refRegions] 
                for item in rtnList:
                    rtn += str(item) + "\t"
                rtn = rtn[:-1] + "\n"
        return(rtn)

    def asPD(self):
        text = self.getHeader()
        text += str(self)
        text = StringIO(text)         
        df = pd.read_csv(text, sep = "\t", header = 0)
        return(df)






