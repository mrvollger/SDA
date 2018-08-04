#!/usr/bin/env python

import argparse
import networkx
ap = argparse.ArgumentParser(description="Print components of file")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

