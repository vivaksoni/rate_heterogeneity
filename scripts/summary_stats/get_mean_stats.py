#script takes .stats file and .ms file and gets mean of summary statistics across windows for each replicate

#Example usage:
#python3  get_mean_stats.py \
#-statsPath "rr_mu_demog_inference/summary_stats/demog_only/stationary/rr_fixed_mu_fixed/" \
#-msPath "rr_mu_demog_inference/results/demog_only/stationary/rr_fixed_mu_fixed/" \
#-outFile "rr_mu_demog_inference/summary_stats/demog_only/stationary/rr_fixed_mu_fixed.stats" \
#-regionLen "1Mb" -minRep 1 -maxRep 150 -samples 100

from __future__ import print_function
import sys
import pandas as pd
import argparse

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-statsPath', dest = 'statsPath', action='store', nargs = 1, type = str, help = 'path to stats files')
parser.add_argument('-msPath', dest = 'msPath', action='store', nargs = 1, type = str, help = 'path to ms files (.ms format)')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = str, help = 'size of region (eg 1Mb, or 100kb')
parser.add_argument('-minRep', dest = 'minRep', action='store', nargs = 1, type = int, help = 'first replicate to begin loop')
parser.add_argument('-maxRep', dest = 'maxRep', action='store', nargs = 1, type = int, help = 'last replicate of loop')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'no. of samples in .ms file')

args = parser.parse_args()
statsPath =  args.statsPath[0]
msPath = args.msPath[0]
outFile = args.outFile[0]
regionLen = args.regionLen[0]
minRep = args.minRep[0]
maxRep = args.maxRep[0]
samples = args.samples[0]


#Get number of segregating sites from .ms file
def get_S(f_ms):
    #Return no. of segregating sites (line 1 of .ms file)
    pos_lines = [1]
    for position, line in enumerate(f_ms):
        if position in pos_lines:
            S = line.split()[1]
    #Set file pointer to start of file
    f_ms.seek(0)
    return int(S)


#Function to get list of allele frequencies from .ms file
def get_af_list(f_ms, samples):
    d_af = {}
    linecount = 0
    for line in f_ms:
        #Remove newline character
        line1 = line.strip('\n')
        #Only include lines of allele frequencies
        if "//" not in line1 and "segsites" not in line1 and "positions" not in line1:
            linecount += 1
            #Only include up to total number of individuals
            if linecount <= int(samples):
                col = 1
                #Loop through characters (ie sites) in line (ie in individual)
                for x in line1:
                    try:
                        d_af[col] = int(d_af[col]) + int(x)
                    except:
                        d_af[col] = int(x)
                    col += 1
    f_ms.close()
    #Convert dictionary into a list
    l_af = [d_af[x] for x in d_af.keys()]
    return(l_af)

#Loop through reps
for rep in range(minRep, maxRep+1):

    #Read summary stats file as df
    df = pd.read_csv(statsPath + regionLen + "_rep" + str(rep) + ".stats", sep='\t', header=0)
    #Drop columns of stats that involve counts (as windows overlap)
    df = df.drop(columns=['window', 'numexternalmutations', 'numpoly', 'numsingletons', 'rm', 'nhaps'])
    #Get means and transpose into columns
    rdf = pd.DataFrame(df.mean()).T

    #Open ms file and get singletons and segregating site counts
    with open(msPath + regionLen + "_rep" + str(rep) + ".ms") as f:
        segSites = get_S(f)
        singletons = get_af_list(f, samples).count(1)
    #Add to results df
    rdf['rep'] = rep
    rdf['segSites'] = segSites
    rdf['singletons'] = singletons

    #Reorder columns so that last three become first three
    cols = rdf.columns.tolist()
    cols = cols[-3:] + cols[:-3]
    rdf = rdf[cols]
    #Only include header for first rep
    if(rep==minRep):
        rdf.to_csv(outFile, header=True, index=False, sep='\t', mode='a')
    else:
        rdf.to_csv(outFile, header=False, index=False, sep='\t', mode='a')

