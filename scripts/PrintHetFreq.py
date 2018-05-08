#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Filter heterozygous positions.\nOutput starts at 1.")
ap.add_argument("minCount", help="Minimum count. Remove rows below this.", type=int)
ap.add_argument("--minTotal", help="Only consider rows where the minimum coverage is above this value.",type=int, default=0)
ap.add_argument("--maxCount", help="Maximum count, remove rows when there are more than two values above this.", type=int,default=0)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()


minCount = int(args.minCount)
sys.stdin.readline()
for line in sys.stdin:
    vals = line.split()
    v = [int(i) > minCount for i  in vals[2:6]]
    
    if args.maxCount > 0:
        m = [int(i) > args.maxCount for i in vals[2:6]]
    else:
        m = [0]
    nv = sum(v)
    nm = sum(m)
    intVals = [int(i) for i in vals[2:6]]
    if sum(intVals) < args.minTotal:
        continue
    if (nv > 1 and nm <= 2):
        sys.stdout.write(line)
         
