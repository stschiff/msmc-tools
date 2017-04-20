#!/usr/bin/env python

import sys
import gzip
import re

class MaskGenerator:
  def __init__(self, filename, chr):
    self.lastCalledPos = -1
    self.lastStartPos = -1
    self.file = gzip.open(filename, "w")
    self.chr = chr
  
  # assume 1-based coordinate, output in bed format
  def addCalledPosition(self, pos):
    if self.lastCalledPos == -1:
      self.lastCalledPos = pos
      self.lastStartPos = pos
    elif pos == self.lastCalledPos + 1:
      self.lastCalledPos = pos
    else:
      self.file.write("{}\t{}\t{}\n".format(self.chr, self.lastStartPos - 1, self.lastCalledPos))
      self.lastStartPos = pos
      self.lastCalledPos = pos


if len(sys.argv) < 3:
    print("too few arguments:")
    print("Usage: ./vcfCaller.py <chrom> <mask_out> > <vcf_out>")
    print("Reading VCF with all called sites, including hom-ref from stdin")
    sys.exit(1)

chr = sys.argv[1]
mask_filename = sys.argv[2]

mask = MaskGenerator(mask_filename, chr)

lastPos = 0
line_cnt = 0
for line in sys.stdin:
  if line[0] == '#':
    print line,
    continue
  fields = line.strip().split('\t')
  pos = int(fields[1])
  refAllele = fields[3]
  altAllele = fields[4]
  info = fields[7]
  genotypes = fields[9]
  if line_cnt % 10000 == 0:
    sys.stderr.write("parsing position {}\n".format(pos))
  line_cnt += 1
  if re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg\.]$", altAllele):
    if genotypes[0] != '.' and genotypes[2] != '.':
      mask.addCalledPosition(pos)
      if genotypes[0] == '1' or genotypes[2] == '1':
        print line,
  

