import multiprocessing
import dadi
import pandas as pd
import numpy as np
import random
import argparse
from itertools import product

#Example usage:

#python3 dadi_sps.py -start 1 -inFile test.fs -outFile results.txt

#for i in $(seq 1 100); do python3 dadi_two_epoch.py -Na_lb -2 -Na_ub 2 -t_lb -2 -t_ub 3 -nu_reps 15 -t_reps 8 \
#-inFile dadi_inputFiles/sfs_files/demog_only/stationary/rr_fixed_mu_fixed/1Mb_rep"$i".fs \
#-outFile dadi_inputFiles/sfs_files/demog_only/stationary/rr_fixed_mu_fixed/results/"$i".txt; done

#np.exp(-6), np.exp(1)
#np.exp(-1), np.exp(6)

#Parse arguments
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-Na_lb', dest = 'Na_lb', action='store', nargs = 1, type = float, help = 'Na lower bound')
parser.add_argument('-Na_ub', dest = 'Na_ub', action='store', nargs = 1, type = float, help = 'Na upper bound')
parser.add_argument('-t_lb', dest = 't_lb', action='store', nargs = 1, type = float, help = 't lower bound')
parser.add_argument('-t_ub', dest = 't_ub', action='store', nargs = 1, type = float, help = 't upper bound')
parser.add_argument('-nu_reps', dest = 'nu_reps', action='store', nargs = 1, type = int, help = 'no. of different starting values for nu')
parser.add_argument('-t_reps', dest = 't_reps', action='store', nargs = 1, type = int, help = 'no. of different starting values for t')
parser.add_argument('-inFile', dest = 'inFile', action='store', nargs = 1, type = str, help = 'path to input file')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')

#read input parameters
args = parser.parse_args()
Na_lb = args.Na_lb[0]
Na_ub = args.Na_ub[0]
t_lb = args.t_lb[0]
t_ub = args.t_ub[0]
nu_reps = args.nu_reps[0]
t_reps = args.t_reps[0]
inFile = args.inFile[0]
outFile = args.outFile[0]


def two_epoch (params, ns, pts):
    """
    Instantaneous size change some time ago.

    params = (nu,T)
    ns = (n1,)

    nu: Ratio of contemporary to ancient population size
    T: Time in the past at which size change happened (in units of 2*Na 
       generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu , T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi , xx , T, nu)
    fs = dadi.Spectrum.from_phi(phi , ns , (xx ,))
    return fs

def sps (params, ns, pts):
    nu = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop (phi , xx , nu)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

fs = dadi.Spectrum.from_file(inFile)

def run_optimization(p0, lb, ub, outFile):
    # These are the grid point settings will use for extrapolation. Sample size is 100, so this is slightly higher
    pts_l = [210, 230, 260]
    func = two_epoch
    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(func)
    popt = dadi.Inference.optimize_log_fmin(p0, fs, func_ex, pts_l, lower_bound=lb, 
    upper_bound=ub, verbose=False, full_output=True, maxiter=100)
    #popt.append(p0)
    with open(outFile, 'a') as f:
        f.write(str(p0[0]) + '\t' + str(p0[1]) + '\t' + str(popt[0][0]) + '\t' + str(popt[0][1]) 
        + '\t' + str(popt[1]) + '\t' + str(popt[2]) + '\t' + str(popt[3]) + '\t' + str(popt[4]) + '\n')
    #return(popt)

#Create lists of initial guesses, evenly distributed acrss logspace.
p0_nu = np.logspace(Na_lb, Na_ub, base=10.0, num=nu_reps)
p0_T = np.logspace(t_lb, np.log10(t_ub), base=10, num=t_reps)
p0_list = [np.array(x) for x in product(p0_nu, p0_T)]

#Upper and lower parameter bounds
lb = [p0_nu[0], p0_T[0]]
ub = [p0_nu[-1], p0_T[-1]]

print ("Outputting results to " + outFile)

process_list = []
for i in p0_list:
    # Perturb our parameters before optimization. This does so by taking each
    # parameter a up to a factor of two up or down.
    p0 = dadi.Misc.perturb_params(i, fold=1, upper_bound=ub, lower_bound=lb)
    #p0=i
    p =  multiprocessing.Process(target=run_optimization, args=[p0, lb, ub, outFile])
    p.start()
    process_list.append(p)

