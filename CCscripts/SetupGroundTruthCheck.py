#!/usr/bin/env python
import sys
import argparse
ap = argparse.ArgumentParser(description="Given a high coverage region file, prepare local directories with read pileups")
ap.add_argument("cov", help="Coverage file.")
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
    outFile.write("blasr {}/ref.fasta /var/tmp/mchaisso/GRCh38.fasta -sam -bestn 30 -out {}/ref.fasta.sam\n".format(regionDir, regionDir))
    outFile.write("samToBed {}/ref.fasta.sam --reportIdentity | awk '{{ if ($3-$2 >10000 && $9 > 0.85) print $1\":\"$2\"-\"$3 }}' > {}/ref.fasta.rgn\n".format(regionDir, regionDir))
    outFile.write("samtools faidx /var/tmp/mchaisso/GRCh38.fasta `cat {}/ref.fasta.rgn` > {}/ref.dups.fasta\n".format(regionDir, regionDir))
    outFile.write("samtools view {}/reads.bam | awk '{{ print \">\"$1;print $10;}}' > {}/reads.fasta\n".format(regionDir, regionDir))
    outFile.write("blasr {}/reads.fasta {}/ref.dups.fasta -bestn 1 -m 4 -nproc 8 -out {}/reads.dups.m4\n".format(regionDir, regionDir,regionDir))
outFile.close()
