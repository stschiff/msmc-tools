#!/usr/bin/env python

import sys
import gzip
import re
from collections import defaultdict

class MaskGenerator:
    def __init__(self, filename, chrom):
        self.chrom = chrom
        self.lastCalledPos = -1
        self.lastStartPos = -1
        self.file = gzip.open(filename, "wt")
    def addCalledPosition(self, pos):
        if self.lastCalledPos == -1:
            self.lastCalledPos = pos
            self.lastStartPos = pos
        elif pos == self.lastCalledPos + 1:
            self.lastCalledPos = pos
        else:
            self.file.write(f"{self.chrom}\t{self.lastStartPos - 1}\t{self.lastCalledPos}\n")
            self.lastStartPos = pos
            self.lastCalledPos = pos
    def close(self):
        if self.lastCalledPos != -1:
            self.file.write(f"{self.chrom}\t{self.lastStartPos - 1}\t{self.lastCalledPos}\n")
        self.file.close()

if len(sys.argv) < 3:
    print("Usage: ./vcfCaller.py <input_vcf> <output_prefix>")
    print("Reads a VCF file, processes all chromosomes, and outputs separate mask BED files and filtered VCFs for each chromosome.")
    sys.exit(1)

input_vcf = sys.argv[1]
output_prefix = sys.argv[2]

prev_chrom = None

# Dictionaries to store VCF file handles and MaskGenerators for each chromosome
vcf_files = {}
mask_generators = {}
header_lines = []  # Store header lines

with gzip.open(input_vcf, "rt") if input_vcf.endswith(".gz") else open(input_vcf, "r") as vcf_file:
    line_cnt = 0
    for line in vcf_file:
        line = line.rstrip()
        if line.startswith('#'):
            # Write header to all chromosome VCF files when they are created
            header_lines.append(line)
            continue
        fields = line.split('\t')
        chrom, pos, refAllele, altAllele, info, genotypes = (
            fields[0], int(fields[1]), fields[3], fields[4], fields[7], fields[9]
        )
        if chrom != prev_chrom:
            if prev_chrom is not None:
                # Close previous chromosome files
                vcf_files[prev_chrom].close()
                mask_generators[prev_chrom].close()
            # Open new files for the new chromosome
            vcf_files[chrom] = gzip.open(f"{output_prefix}_{chrom}.vcf.gz", "wt")
            mask_generators[chrom] = MaskGenerator(f"{output_prefix}_{chrom}_mask.bed.gz", chrom)
            # Write headers for the new chromosome
            for header in header_lines:
                vcf_files[chrom].write(header)
            prev_chrom = chrom  # Update previous chromosome tracker
        if line_cnt % 10000 == 0:
            sys.stderr.write(f"Parsing position {pos} on {chrom}\n")
        line_cnt += 1
        if re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg\.]$", altAllele):
            if genotypes[0] != '.' and genotypes[2] != '.':
                mask_generators[chrom].addCalledPosition(pos)
                if genotypes[0] == '1' or genotypes[2] == '1':
                    vcf_files[chrom].write(line + "\n")

if prev_chrom:
    vcf_files[prev_chrom].close()
    mask_generators[prev_chrom].close()

print("Processing complete. Separate output files per chromosome are compressed with bgzip.")
