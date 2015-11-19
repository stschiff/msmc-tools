#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(description="parse output from ms and scrm simulators.")
parser.add_argument("chr", metavar="chrom", type=str, help="chromosome label to use")
parser.add_argument("L", metavar="length", type=int, help="length of chromosome")

args = parser.parse_args()

positions = []
alleles = []
for line in sys.stdin:
    if line.startswith("positions:"):
        fields = line.strip().split()
        for f in fields[1:]:
            positions.append(float(f))
            alleles.append([])
        continue
    
    if len(positions) > 0:
        for i, letter in enumerate(line.strip()):
            alleles[i].append(letter)

lastPos = 0
for i in range(len(positions)):
    realPos = int(positions[i] * args.L)
    if realPos > lastPos:
        print(args.chr, realPos, realPos - lastPos, "".join(alleles[i]), sep="\t")
    lastPos = realPos
