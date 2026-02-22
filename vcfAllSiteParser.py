#!/usr/bin/env python3

import sys
import gzip
import re
import argparse

class MaskGenerator:
    def __init__(self, filename, chrom):
        self.lastCalledPos = -1
        self.lastStartPos = -1
        self.file = gzip.open(filename, "wt")  # text mode
        self.chrom = chrom

    # assume 1-based coordinate, output in bed format
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


def parse_args():
    parser = argparse.ArgumentParser(
        description="vcfAllSiteParser.py\n\n"
                    "Reads a VCF from stdin and outputs a mask file (BED) and a filtered VCF.\n"
                    "Original mode: single chromosome (use positional <chrom> <mask_out> [<vcf_out>]).\n"
                    "Split mode: multiple chromosomes using --splitChromosomes."
    )
    parser.add_argument("--splitChromosomes", action="store_true",
                        help="Enable split mode for multiple chromosomes (VCF must be sorted by chromosome)")
    # Positional arguments
    parser.add_argument("arg1",
                        help="Chromosome name (<chrom> original mode), or prefix for mask files (<mask_out> split mode)")
    parser.add_argument("arg2",
                        help="Mask output file (<mask_out> original mode) or prefix for VCF files (<vcf_out> split mode)")
    parser.add_argument("arg3", nargs='?', default=None,
                        help="VCF output file (<vcf_out> original mode)")
    return parser.parse_args()


def split_mode(mask_prefix, vcf_prefix):
    vcf_files = {}
    mask_generators = {}
    header_lines = []
    prev_chrom = None
    line_cnt = 0

    for line in sys.stdin:
        line = line.rstrip()
        if line.startswith('#'):
            header_lines.append(line + "\n")
            continue

        fields = line.split('\t')
        chrom = fields[0]
        pos = int(fields[1])
        refAllele = fields[3]
        altAllele = fields[4]
        info = fields[7]
        genotypes = fields[9]

        if chrom != prev_chrom:
            if prev_chrom is not None:
                vcf_files[prev_chrom].close()
                mask_generators[prev_chrom].close()

            vcf_filename = f"{vcf_prefix}_{chrom}.vcf.gz"
            mask_filename = f"{mask_prefix}_{chrom}_mask.bed.gz"
            vcf_files[chrom] = gzip.open(vcf_filename, "wt")
            mask_generators[chrom] = MaskGenerator(mask_filename, chrom)

            for header in header_lines:
                vcf_files[chrom].write(header)

            prev_chrom = chrom

        if line_cnt % 10000 == 0:
            sys.stderr.write(f"Parsing position {pos} on {chrom}\n")
        line_cnt += 1

        if re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg\.]$", altAllele):
            if genotypes[0] != '.' and genotypes[2] != '.':
                mask_generators[chrom].addCalledPosition(pos)
                if genotypes[0] == '1' or genotypes[2] == '1':
                    vcf_files[chrom].write(line + "\n")

    if prev_chrom is not None:
        vcf_files[prev_chrom].close()
        mask_generators[prev_chrom].close()


def original_mode(chrom, mask_out, vcf_out=None):
    mask = MaskGenerator(mask_out, chrom)
    line_cnt = 0

    # Determine output VCF handle
    if vcf_out:
        if vcf_out.endswith(".gz"):
            vcf_handle = gzip.open(vcf_out, "wt")  # text mode
            to_close = True
        else:
            vcf_handle = open(vcf_out, "a")
            to_close = True
    else:
        vcf_handle = sys.stdout
        to_close = False  # don't close stdout

    for line in sys.stdin:
        if line.startswith('#'):
            vcf_handle.write(line)
            continue

        fields = line.strip().split('\t')
        pos = int(fields[1])
        refAllele = fields[3]
        altAllele = fields[4]
        info = fields[7]
        genotypes = fields[9]

        if line_cnt % 10000 == 0:
            sys.stderr.write(f"Parsing position {pos}\n")
        line_cnt += 1

        if re.match("^[ACTGactg]$", refAllele) and re.match("^[ACTGactg\.]$", altAllele):
            if genotypes[0] != '.' and genotypes[2] != '.':
                mask.addCalledPosition(pos)
                if genotypes[0] == '1' or genotypes[2] == '1':
                    vcf_handle.write(line + "\n")

    if to_close:
        vcf_handle.close()
    mask.close()


def main():
    args = parse_args()

    if args.splitChromosomes:
        # In split mode, 'chrom' is used as mask prefix
        mask_prefix = args.arg1
        vcf_prefix = args.arg2 if args.arg2 else sys.exit(
            "Error: split mode requires both mask and VCF prefixes"
        )
        print("VCF with multiple chromosome input must be sorted by chromosome", file=sys.stderr)
        split_mode(mask_prefix, vcf_prefix)
    else:
        # Original mode: chrom = chromosome name
        chrom = args.arg1
        mask_out = args.arg2
        vcf_out = args.arg3
        original_mode(chrom, mask_out, vcf_out)


if __name__ == "__main__":
    main()
