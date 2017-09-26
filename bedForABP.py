#!/usr/bin/env python
import re

class BedRec:
    def __init__(self, chr=None, start=None, end=None, 
            name="name", score="0", strand="+", 
            thickS=None, thinkE=None, rgb="0,0,0"):
        self.chr = chr
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        if(thickS is None):
            self.thickS = self.start
            self.thickE = self.end
        else:
            self.thickS = thickS
            self.thickE = thickE
        self.rgb = rgb

    def __str__(self):
        rtn = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.chr, self.start, self.end,
            self.name, self.score, self.strand,
            self.thickS, self.thickE, self.rgb)
        return(rtn)

rtn = ""

colapsePos = open("ref.fasta.bed") 
for line in colapsePos:
    token = line.strip().split("\t")
    chr = token[0]
    start = token[1]
    end = token[2]
    rgb="255,0,0"
    rtn += str(BedRec(chr=chr, start=start, end=end, rgb=rgb)) +"\n"

m5 = open("WH_dup.m5")
for idx, line in enumerate(m5):
    token = line.strip().split(" ")
    key = token[5]
    temp = re.split(":|-", key)
    chr = temp[0]
    start = int(temp[1])
    end = int(temp[2])
    contig = token[0]
    mstart = int(token[7])
    mend = int(token[8])
    strand = token[9]
    diff = mend - mstart 
    perID = float(token[11])/(float(token[11]) + float(token[12]))
    realStart = start + mstart
    realEnd = realStart + diff
    

    rgb = "0,0,0"
    if(perID >= 0.998):
        rgb = "0,255,0"
    elif(perID >= 0.95):
        rgb="128,128,128"
    else:
        rgb="255,255,0"

    if(idx == 0):
        print("browser position " + "{}:{}-{}".format(chr, start,end) )
        print('track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"')
    
    B = BedRec(chr=chr, start=realStart, end=realEnd, name=contig, strand=strand, rgb=rgb)
    rtn += str(B) + "\n"


print(rtn)
f=open("asm.bed", "w+")
f.write(rtn)
f.close()

