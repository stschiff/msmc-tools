#!/usr/bin/env python3

import sys
import gzip
import string
import copy
import argparse
import io

class MaskIterator:
    def __init__(self, filename, negative=False):
        if filename[-3:] == ".gz":
            self.file = io.TextIOWrapper(gzip.open(filename, "r"))
        else:
            self.file = open(filename, "r") #io.TextIOWrapper(open(filename, "r"))
        self.eof = False
        self.lastPos = 1
        self.negative = negative
        self.readLine()
  
    def readLine(self):
        try:
            line = next(self.file)
            fields = line.strip().split()
            if len(fields) == 2:
                self.start = int(fields[0])
                self.end = int(fields[1])
            else:
                self.start = int(fields[1]) + 1
                self.end = int(fields[2])
        except StopIteration:
            self.eof = True
  
    def getVal(self, pos):
        assert pos >= self.lastPos
        self.lastPos = pos
        while not self.eof and pos > self.end:
            self.readLine()
        if self.eof:
            return None
        if pos >= self.start and pos <= self.end:
            return True if not self.negative else False
        else:
            return False if not self.negative else True

class MergedMask:
    def __init__(self, mask_iterators):
        self.maskIterators = mask_iterators
  
    def getVal(self, pos):
        return all((m.getVal(pos) for m in self.maskIterators))

class VcfIterator:
    def __init__(self, filename):
        self.file = io.TextIOWrapper(gzip.open(filename, "r"))
    
    def __iter__(self):
        return self
    
    def __next__(self):
        line = next(self.file)
        while line[0] == "#":
            try:
                line = next(self.file)
            except StopIteration:
                return None
        fields = line.strip().split()
        chrom = fields[0]
        pos = int(fields[1])
        alleles = [fields[3]]
        for alt_a in fields[4].split(","):
            alleles.append(alt_a)
        geno = fields[9][:3]
        if len(geno) != 3 :
                print ("Non-diploid SNP found and considered as unphased data: %s" % geno, file=sys.stderr)
                phased = False
                geno = "%s/%s" % (geno[0], geno[0])
        else :
                phased = geno[1] == "|"
        return (chrom, pos, tuple(alleles), (int(geno[0]), int(geno[2])), phased)

class OrderedAlleles:
    def __init__(self):
        self.ordered_alleles = []
    
    def addGenotype(self, a1, a2, phasing):
        if len(self.ordered_alleles) == 0:
            self.ordered_alleles = [[a1, a2]]
            if not phasing and a1 != a2:
                self.ordered_alleles.append([a2, a1])
        else:
            new = []
            for o in self.ordered_alleles:
                new.append(o + [a1, a2])
                if not phasing and a1 != a2:
                    new.append(o + [a2, a1])
            self.ordered_alleles = new
  
    def phase(self, trio):
        child, father, mother = trio
        new = [] 
        for o in self.ordered_alleles:
            child_1, child_2 = o[2 * child : 2 * (child + 1)]
            pat = o[2 * father]
            mat = o[2 * mother]
            if child_1 == pat and child_2 == mat or child_2 == pat and child_1 == mat:
                new.append(o)
        if len(new) > 0:
            self.ordered_alleles = new
        self.ordered_alleles = unique(self.ordered_alleles)
    
    def getPrint(self, trios):
        child_indices = []
        for i in range(len(self.ordered_alleles[0])):
            for child, father, mother in trios:
                if i == 2 * child or i == 2 * child + 1:
                    child_indices.append(i)
        print_indices = [i for i in range(len(self.ordered_alleles[0])) if i not in child_indices]
        stripped_alleles = unique([[o[i] for i in print_indices] for o in self.ordered_alleles])
        if len(stripped_alleles[0]) == 2:
            return ''.join([self.ordered_alleles[0][i] for i in print_indices])
        else:
            return ','.join([''.join(o) for o in stripped_alleles])

def unique(list_of_lists):
    return list(set([tuple(l) for l in list_of_lists]))

class JoinedVcfIterator:
    def __init__(self, filenames, trios):
        self.vcfIterators = [VcfIterator(f) for f in filenames]
        self.current_lines = [next(v) for v in self.vcfIterators]
        self.trios = trios
    
    def __iter__(self):
        return self
    
    def __next__(self):
        minIndices = self.getMinIndices()
        chrom = self.current_lines[minIndices[0]][0]
        pos = self.current_lines[minIndices[0]][1]
        ref = self.current_lines[minIndices[0]][2][0]
      
        ordered_alleles = OrderedAlleles()
        
        for i, l in enumerate(self.current_lines):
            if i not in minIndices:
                ordered_alleles.addGenotype(ref, ref, True)
            else:
                alleles, geno, phased = l[2:5]
                ordered_alleles.addGenotype(alleles[geno[0]], alleles[geno[1]], phased)
                try:
                    self.current_lines[i] = next(self.vcfIterators[i])
                except StopIteration:
                    self.current_lines[i] = None
        for trio in self.trios:
            ordered_alleles.phase(trio)
        return (chrom, pos, ordered_alleles.getPrint(self.trios))
    
    def getMinIndices(self):
        activeLines = [(i, l) for i, l in enumerate(self.current_lines) if l]
        if len(activeLines) == 0:
            raise StopIteration
        if len(activeLines) == 1:
            return [activeLines[0][0]]
        else:
            minIndices = [activeLines[0][0]]
            minPos = activeLines[0][1][1]
            for a in activeLines[1:]:
                if a[1][1] == minPos:
                    minIndices.append(a[0])
                if a[1][1] < minPos:
                    minPos = a[1][1]
                    minIndices = [a[0]]
            return minIndices
    

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+", help="Input VCF files")
parser.add_argument("--mask", action="append", help="apply masks in bed format, should be given once for the calling mask from each individual, and in addition can be given for e.g. mappability or admixture masks. Mask can be gzipped, if indicated by .gz file ending.")
parser.add_argument("--negative_mask", action="append", help="same as mask, but interpreted as negative mask, so places where sites should be excluded")
parser.add_argument("--trio", action="append", help="declare trio-relationships. This should be a string with a format <child_index>,<father_index>,<mother_index>, where the three fields are the indices of the samples in the trio. This option will automatically phase parental and maternal haplotypes where possible and remove the child VCF file from the resulting file. Can be given multiple times if you have multiple trios.")
parser.add_argument("--chr", help="overwrite chromosomes in input files. Useful if chromosome names differ, such as chr1 vs. 1")

args = parser.parse_args()

trios = []
if args.trio:
    trios = [tuple(map(int, t.split(","))) for t in args.trio]

nrIndidividuals = len(args.files)
nrHaplotypes = 2 * (nrIndidividuals - len(trios))

sys.stderr.write("generating msmc input file with {} haplotypes\n".format(nrHaplotypes))

joinedVcfIterator = JoinedVcfIterator(args.files, trios)
maskIterators = []
if args.mask:
    for f in args.mask:
        sys.stderr.write("adding mask: {}\n".format(f))
        maskIterators.append(MaskIterator(f))
if args.negative_mask:
    for nm in args.negative_mask:
        sys.stderr.write("adding negative mask: {}\n".format(nm))
        maskIterators.append(MaskIterator(nm, True))

mergedMask = MergedMask(maskIterators)

def is_segregating(alleles):
    orders = alleles.split(",")
    for o in orders:
        for a in o[1:]:
            if a != o[0]:
                return True
    return False

pos = 0
nr_called = 0
for chrom, snp_pos, alleles in joinedVcfIterator:
    # sys.stderr.write("{}\t{}\t{}\n".format(chrom, snp_pos, alleles))
    while pos < snp_pos:
        pos += 1
        if mergedMask.getVal(pos):
            nr_called += 1
        if pos % 1000000 == 0:
            print("processing pos {}".format(pos), file=sys.stderr)
    if mergedMask.getVal(snp_pos):
        if is_segregating(alleles):
            c = chrom if not args.chr else args.chr
            print(c, snp_pos, nr_called, alleles, sep="\t")
            nr_called = 0
  
  
  
