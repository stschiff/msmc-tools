#!/usr/bin/env python3

import argparse
import plot_utils

parser = argparse.ArgumentParser()
parser.add_argument("crossCoalFile", help="Input file: Cross coalescence rate")
parser.add_argument("withinCoalFile1", help="Input file: coalescence rate within population 1")
parser.add_argument("withinCoalFile2", help="Input file: coalescence rate within population 2")

args = parser.parse_args()

cc = plot_utils.MSMCresult(args.crossCoalFile)
within1 = plot_utils.MSMCresult(args.withinCoalFile1)
within2 = plot_utils.MSMCresult(args.withinCoalFile2)

interp1_min, interp1_max, interp1 = within1.getInterp()
interp2_min, interp2_max, interp2 = within2.getInterp()

print("time_index", "left_time_boundary", "right_time_boundary", "lambda_00", "lambda_01", "lambda_11", sep="\t")

res = 10
for i, (tl, tr, l) in enumerate(zip(cc.times_left, cc.times_right, cc.lambdas[0])):
    lambda00 = 0.0
    lambda11 = 0.0
    for j in range(res):
        t = tl + j / float(res) * (tr - tl)
        if t <= interp1_min:
            l00 = within1.lambdas[0][0]
        elif t >= interp1_max:
            l00 = within1.lambdas[0][-1]
        else:
            l00 = interp1(t)

        if t <= interp2_min:
            l11 = within2.lambdas[0][0]
        elif t >= interp2_max:
            l11 = within2.lambdas[0][-1]
        else:
            l11 = interp2(t)
        
        lambda00 += l00 / float(res)
        lambda11 += l11 / float(res)
    print(i, tl, tr, lambda00, l, lambda11, sep="\t")
