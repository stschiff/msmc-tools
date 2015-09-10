#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fn")
parser.add_argument("--row", type=int, default=-1, help="loop number, starting with 0. Default is -1, indicating the last row")
args = parser.parse_args()

f = open(args.fn, "rt")
lines = f.readlines()
line = lines[args.row]
[timesStr, lambdaStr] = line.strip().split()[2:]
times = list(map(float, timesStr.split(",")))[:-1]
times.append(float("inf"))
lambdaVec = list(map(float, lambdaStr.split(",")))
assert len(lambdaVec) == len(times) - 1 or len(lambdaVec) == 3 * (len(times) - 1), "lambdaVec and times not of compatible length"

if len(lambdaVec) == 3 * (len(times) - 1):
    lambdaChunks = [lambdaVec[i:i+3] for i in range(0, len(lambdaVec), 3)] 
else:
    lambdaChunks = [[lambdaVec[i]] for i in range(len(lambdaVec))]

if len(lambdaChunks[0]) == 1:
    print("time_index", "left_time_boundary", "right_time_boundary", "lambda_00", sep="\t")
else:
    print("time_index", "left_time_boundary", "right_time_boundary", "lambda_00", "lambda_01", "lambda_11", sep="\t")

for i in range(len(times) - 1):
    print(i, abs(times[i]), times[i+1], *lambdaChunks[i], sep="\t") # The abs only fixes the -0.0 issue
        
