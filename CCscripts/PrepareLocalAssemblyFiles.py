#!/usr/bin/env python
import sys
import argparse
ap = argparse.ArgumentParser(description="Given a high coverage region file, prepare local directories with read pileups")
ap.add_argument("cov", help="Coverage file.")
ap.add_argument("fofn", help="Read fofn.")
ap.add_argument("ref", help="Reference file.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")

args = ap.parse_args()

outFile = open(args.out, 'w')
covFile = open(args.cov)
outFile.write( "#!/usr/bin/env bash\n")
outFile.write("set v\n")
for line in covFile:
    vals = line.split()
    regionDir = ".".join(vals[0:3])
    region = vals[0] + ":" + vals[1] + "-" + vals[2]
    outFile.write("echo " + region + "\n")
#    outFile.write("mkdir " + regionDir + "\n")
#    outFile.write("/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/CollectRegion.sh {} {} {}/reads.orig.bam".format(args.fofn, region, regionDir) + "\n")
#    outFile.write("samtools index {}/reads.orig.bam".format(regionDir) + "\n")
    outFile.write("bounds=`bedtools bamtobed -i {}/reads.orig.bam | /net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/local_asm/GetRegionBounds.py`".format(regionDir) + "\n")
    # remember, start is now 0-based
    outFile.write("fastart=`echo $bounds | awk '{print $1+1;}'`" + "\n")    
    outFile.write("start=`echo $bounds | awk '{print $1;}'`" + "\n")
    outFile.write("endp=`echo $bounds | awk '{print $2;}'`" + "\n")
    outFile.write("samtools faidx {} {}:$fastart-$endp | /net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/local_asm/FixRef.py > {}/ref.fasta".format(args.ref, vals[0], regionDir) + "\n") 
    outFile.write("samtools view -h {}/reads.orig.bam | /net/eichler/vol5/home/mchaisso/projects/AssemblyByPhasing/scripts/local_asm/FixPos.py {} $start | samtools view -bS - | samtools sort -T tmp -o {}/reads.bam".format(regionDir, vals[0], regionDir) + "\n")
    outFile.write("samtools index {}/reads.bam".format(regionDir) + "\n")

outFile.close()
