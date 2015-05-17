#!/usr/bin/env python3
import math
import bisect
from scipy.interpolate import interp1d 

class MSMCresult:
    """
    Simple class to read a MSMC result file. Constructor takes filename of MSMC result file
    """
    def __init__(self, filename):
        f = open(filename, "r")
        self.times_left = []
        self.times_right = []
        self.lambdas = []
        next(f) # skip header
        for line in f:
            fields = line.strip().split()
            nr_lambdas = len(fields) - 3
            if len(self.lambdas) == 0:
                self.lambdas = [[] for i in range(nr_lambdas)]
            time_left = float(fields[1])
            time_right = float(fields[2])
            self.times_left.append(time_left)
            self.times_right.append(time_right)
            for i in range(nr_lambdas):
                l = float(fields[3 + i])
                self.lambdas[i].append(l)
        self.T = len(self.times_left)
        self.times_left[0] = self.times_left[1] / 4.0
        self.times_right[-1] = self.times_right[-2] * 4.0
    
    def getInterp(self):
        x = [0.0] + [0.5 * (tl + tr) for tl, tr in zip(self.times_left, self.times_right)]
        y = [0.0] + self.lambdas[0]
        return (x[0], x[-1], interp1d(x, y))

def popSizeStepPlot(filename, mu=1.25e-8, gen=30.0):
    """
    to be used with a step-plot function, e.g. matplotlib.pyplot.steps.
    returns (x, y), where x contains the left point of each step-segment in years, and y contains the effective population size. Note that there are two ways to make a step-plot. You should make sure that your step-plot routine moves right and then up/down instead of the other way around.
    If plotted on a logarithmic x-axis, you should adjust x[0] = x[1] / 4.0, otherwise the leftmost segment will start at 0 and won't be plotted on a log-axis.
    
    Options:
        mu: Mutation rate per generation per basepair (default=1.25e-8)
        gen: generation time in years (default=30)
    """
    M = MSMCresult(filename)
    x = [t * gen / mu for t in M.times_left]
    y = [(1.0 / l) / (2.0 * mu) for l in M.lambdas[0]]
    return (x, y)

def coalRatePlot(filename, mu=1.25e-8, gen=30.0):
    M = MSMCresult(filename)
    x = [t * gen / mu for t in M.times_left]
    y = [l for l in M.lambdas[0]]
    return (x, y)

def crossCoalPlot(filename, mu=1.25e-8, gen=30.0):
    """
    returns (x, y) where x is the time in years and y is the relative cross coalescence rate.
    Check also the doc-string in popSizeStepPlot for Options and hints on how to plot.
    """
    M = MSMCresult(filename)
    x = [t * gen / mu for t in M.times_left]
    y = [2.0 * M.lambdas[1][i] / (M.lambdas[0][i] + M.lambdas[2][i]) for i in range(M.T)]
    return (x, y)

def crossCoalPlotCombined(filename1, filename2, filename12, mu=1.25e-8, gen=30.0):
    """
    returns (x, y) where x is the time in years and y is the relative cross coalescence rate.
    This function should be used with msmc2, where the three files correspond to the three different
    cases within-population1, within-population2 and across populations.
    Check also the doc-string in popSizeStepPlot for Options and hints on how to plot.
    """
    M1 = MSMCresult(filename1)
    M2 = MSMCresult(filename2)
    M12 = MSMCresult(filename12)
    I1 = M1.getInterp()
    I2 = M2.getInterp()
    x = []
    y = []
    resolution = 10
    for i in range(M12.T - 1):
        tLeft = M12.times_left[i]
        tRight = M12.times_right[i]
        avgWithinRate = 0.0
        for j in range(resolution):
            t = tLeft + j / float(resolution) * (tRight - tLeft)
            lambda1 = I1(t)
            lambda2 = I2(t)
            avgWithinRate += 0.5 * (lambda1 + lambda2) / float(resolution)
        x.append(tLeft * gen / mu)
        y.append(M12.lambdas[0][i] / avgWithinRate)
    return (x, y)     

def tmrcaDistribution(filename, resolution=10, lambda_index=0, mu=1.25e-8, gen=30, cdf=False):
    """
    returns (x, y) where x is the time in years, and y is the probability density for the tMRCA distribution.
    Options:
        resolution: sets the number of time points in each time-segment, default=10.
        lambda_index: sets which of the lambda-columns in the msmc result file should be used, default=0, i.e. lambda_00
        Check popSizeStepPlot for Options mu and gen.
    """
    fprob = get_tmrca_cumprob if cdf else get_tmrca_prob
    M = MSMCresult(filename)
    x = []
    y = []
    for i in range(M.T - 1):
        tLeft = M.times_left[i]
        tRight = M.times_right[i]
        for j in range(resolution):
            t = tLeft + j / float(resolution) * (tRight - tLeft)
            p = fprob(t, i, M.times_left, M.lambdas[lambda_index])
            x.append(t * gen / mu)
            y.append(p)
    return (x, y)

def get_tmrca_prob(t, left_index, time_boundaries, lambda_vals):
    """
    Helper function to compute the tMRCA probability density at time t
    """
    deltas = [time_boundaries[i + 1] - time_boundaries[i] for i in range(len(time_boundaries) - 1)]
    tleft = time_boundaries[left_index]
    lambda_ = lambda_vals[left_index]
    integ = sum(delta * lambda_prime for delta, lambda_prime in zip(deltas[:left_index], lambda_vals[:left_index])) + (t - tleft) * lambda_
    return lambda_ * math.exp(-integ)

def get_tmrca_cumprob(t, left_index, time_boundaries, lambda_vals):
    """
    Helper function to compute the tMRCA cumulative probability function at time t
    """
    deltas = [time_boundaries[i + 1] - time_boundaries[i] for i in range(len(time_boundaries) - 1)]
    tleft = time_boundaries[left_index]
    lambda_ = lambda_vals[left_index]
    integ = sum(delta * lambda_prime for delta, lambda_prime in zip(deltas[:left_index], lambda_vals[:left_index])) + (t - tleft) * lambda_
    return 1.0 - math.exp(-integ)
