#script takes .ms file and calculates summary statistics across sliding windows using libsequence

#python3  stats_sliding_window_masked.py -msFile "rr_mu_demog_inference/results/demog_only/stationary/rr_fixed_mu_fixed/1Mb_rep1.ms" \
#-fixedFile "rr_mu_demog_inference/results/demog_only/stationary/rr_fixed_mu_fixed/1Mb_rep1.fixed" \
#-outFile  "rr_mu_demog_inference/results/demog_only/stationary/rr_fixed_mu_fixed/1Mb_rep1.stats" \
#-winSize 5000 -stepSize 2500 -regionLen 1000000 -samples 100 -N 7000 -burnIn 16

from __future__ import print_function
import libsequence
import sys
import pandas as pd
import numpy as np
import math
import argparse


#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-msFile', dest = 'msFile', action='store', nargs = 1, type = str, help = 'path to ms file (.ms format)')
parser.add_argument('-fixedFile', dest = 'fixedFile', action='store', nargs = 1, type = str, help = 'path to file containing fixed mutation data (.fixed format)')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, type = int, help = 'size of sliding windows')
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, type = int, help = 'step size for sliding windows')
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'size of simulated region')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'no. of samples in .ms file')
parser.add_argument('-N', dest = 'N', action='store', nargs = 1, type = int, help = 'Ancestral population size')
parser.add_argument('-burnIn', dest = 'burnIn', action='store', nargs = 1, type = int, help = 'Length of burn in (in N generations)')


args = parser.parse_args()
chr_len =  args.regionLen[0]
win_size = args.winSize[0]/float(chr_len)
step_size = args.stepSize[0]/float(chr_len)
samples = args.samples[0]
N = args.N[0]
burnIn = args.burnIn[0]
ms_file = args.msFile[0]
f_file = args.fixedFile[0]
out_file = args.outFile[0]



#function reads through .fixed file, 
def read_fixed_mutations(f_fixed, chr_len, N, burnIn):
    #Create empty dict to store substitutions
    d_subs = {}
    #Loop through lined in.fixed file
    for line in f_fixed:
        #strip newline character
        line1 = line.strip('\n')
        #split line into list
        lst = line1.split()
        #Skip first two lines to move straight to data lines
        if line1[0]!="#" and lst[0]!="Mutations:" and int(lst[-1])>(N*burnIn):
            #Estimate position using decimal notation
            posn = float(lst[3])/float(chr_len)
            #Assign base position as key
            #With each line, if position occurs again, value is incremented
            #(To account for repeat mutation in the same position)
            d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs 




def avg_divergence_win(d_subs, start, end):
    s_sum = 0
    for posn in d_subs.keys():
        if float(posn) <= end and float(posn) > start:
            s_sum = s_sum + 1
    return s_sum



#Get number of segregating sites from .ms file
def get_S(f_ms):
    #Return no. of segregating sites (line 1 of .ms file)
    pos_lines = [1]
    for position, line in enumerate(f_ms):
        if position in pos_lines:
            S = line.split()[1]
    #Set file pointer to start of file
    f_ms.seek(0)
    return S



#Function to parse .ms file data into nested list of positions and genotypes
#Function takes as input open reading file object and no. of samples
def get_nested_data_list(f_ms, samples):
    l_Pos = [] #list of positions of SNPs
    l_Genos = [] #list of alleles
    d_tmp = {} #dict to store individual allele info for each individual (values) at each site (keys)

    #positions on line 2
    pos_lines = [2]
    for position, line in enumerate(f_ms):
        if position in pos_lines:
            #store positions in list
            pos_list  = line.split()
    #Set file pointer to start of file
    f_ms.seek(0)

    i = 0
    #Loop through positions, storing in list
    for pos in pos_list[1:]:
        #Append position to l_Pos (after converting to float)
        l_Pos.append(float(pos))
        #Add dictionary key for each position, with empty value
        d_tmp[str(i)] = ""
        i += 1  

    #genotypes on line 3 onwards (use samples argument to determine length of file)
    g_lines = [x for x in range(3, samples + 4)]
    #Loop through lines (ie individuals)
    for position, line in enumerate(f_ms):
        if position in g_lines:
            #Remove newline character
            line1 = line.strip('\n')
            i = 0
            #For each individual, loop through each site, appending allele information for that individual to 
            #the site number in dict
            while i < len(line1):
                d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                i = i + 1

    f_ms.seek(0)

    #Create nested list of positions and genotypes
    l_data = [[j, d_tmp[str(i)]] for i,j in enumerate(l_Pos)]
    return(l_data)



#Function to create dictionary of polySIM summary statistics. Returns dictionary of summary stats
def get_polySIM_stats(sd):
    #Create polysim object
    ps = libsequence.PolySIM(sd)
    #Create list of methods (ie polySIM summaryStats)
    a = [method for method in dir(ps) if callable(getattr(ps, method)) if not method.startswith('_')]
    #Loop through methods, storing names as keys, and estimates as values in dict
    ss_dict = {}
    for method in a:
        ss_dict[method] = getattr(ps, method)()
        
    return(ss_dict)



#Function to create dictionary of LD stats. Returns dictionary.
def get_LD_stats(sd):
    ld = libsequence.ld(sd)
    df = pd.DataFrame(ld)
    ss_dict = {}
    ss_dict['meanrsq'] = sum(df['rsq'])/len(df['rsq'])
    ss_dict['meanD'] = sum(df['D'])/len(df['D'])
    ss_dict['meanDprime'] = sum(df['Dprime'])/len(df['Dprime'])
    return(ss_dict)


#Chromosome structure EDIT THESE!!!
totalGenes = 127
intergenicLength = 3811
exonLength = 588
exonsPerGene = 4
intronLength = 563
intronsPerGene = 3

exons = []
introns = []
intergenic = []
divergence = []

#Calculate length of gene
geneLength = (exonLength * exonsPerGene) + (intronLength * intronsPerGene)

#Loop through total number of genes to calculate coordinates for each chromosomal element type
for gene in range(0, totalGenes):
    geneStart = gene * (geneLength + intergenicLength) + 1
    #print(geneStart)
    for element in range(1, exonsPerGene+1):
        exonStart = geneStart + (exonLength * (element-1)) + (intronLength * (element-1))
        exonEnd = exonStart+(exonLength-1)
        exons.append([exonStart/chr_len, exonEnd/chr_len])

        if (element < exonsPerGene):
            intronStart = exonStart + exonLength
            intronEnd = intronStart+(intronLength-1)
            introns.append([intronStart/chr_len, intronEnd/chr_len])

    intergenicStart = exonEnd + 1
    intergenicEnd = intergenicStart + (intergenicLength-1)
    intergenic.append([intergenicStart/chr_len, intergenicEnd/chr_len])

#Get total non-coding region length
nc_region_length = chr_len -  (exonsPerGene * exonLength * totalGenes)



#Read in fixed mutations
f_subs = open(f_file, 'r')
d_subs = read_fixed_mutations(f_subs, chr_len, N, burnIn)
#Convert dictionary to df
df = pd.DataFrame([list(d_subs.keys()), list(d_subs.values())]).T
df.columns = ['position', 'fixations']
#Create empty list to store masking status
res = []
#Loop through positions in dataframe
for pos in df.position:
    #Set tmp to 0 (ie unmasked site)
    tmp = 0
    #Loop through exons, checking if site is exonic
    for i in exons:
        if((pos>i[0]) & (pos<i[1])):
        #If site is exonic, set tmp to 1 (ie mask site) and break out of loop (move to next site)
            tmp = 1
            break
    #Append mask status to results list
    res.append(tmp)
#Add column of masking status to df
df['masked'] = res
# #Remove masked sites
df = df[df.masked==0]
df = df[['position', 'fixations']]
#Convert df backinto dict
keys = list(df.position)
values = list(df.fixations)
d_subs = dict(zip(keys, values))



#Read in .ms file, create sd object for libsequence
f_ms = open(ms_file, 'r')
l_data = get_nested_data_list(f_ms, 100)

#Convert l_data to dataframe (in order to mask sites)
df = pd.DataFrame(l_data)
df.columns = ['position', 'genotype']

#Create empty list to store masking status
res = []
#Loop through positions in dataframe
for pos in df.position:
    #Set tmp to 0 (ie unmasked site)
    tmp = 0
    #Loop through exons, checking if site is exonic
    for i in exons:
        if((pos>i[0]) & (pos<i[1])):
        #If site is exonic, set tmp to 1 (ie mask site) and break out of loop (move to next site)
            tmp = 1
            break
    #Append mask status to results list
    res.append(tmp)
#Add column of masking status to df
df['masked'] = res
# #Remove masked sites
df = df[df.masked==0]
df = df[['position', 'genotype']]

#Convert back into nested lists
l_data = df.values.tolist()

#No. of segregating sites is length of df
S = len(df)
sd = libsequence.SimData(l_data)


#define sliding windows
wins = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
num_wins = len(wins)

#Loop through windows
for i, win in enumerate(wins):
	#Create empty dictionary, add replicate number
	d = {}
	d['window'] = i
	#Add divergence data to dict
	s_start = (i)*(1.0/float(num_wins))
	s_end = s_start + win_size
	d['divergence'] = avg_divergence_win(d_subs, s_start, s_end)
	#Combine dict output of functions with global dict
	d = {**d, **get_polySIM_stats(win)}
	if len(win.pos()) >= 5: #LD stats are pairwise. If only 1 site exists, it'll show an error.
		d = {**d, **get_LD_stats(win)}
	else:
		d['meanrsq'] = 'NA'
		d['meanD'] = 'NA'
		d['meanDprime'] = 'NA'

	#Convert dict to dataframe and append to file
	df = pd.DataFrame.from_dict(d, orient='index').T
	#Only include header for first window
	if(i==0):
		df.to_csv(out_file, sep='\t', index=False, header=True, mode='a')
	else:
		df.to_csv(out_file, sep='\t', index=False, header=False, mode='a')

	print ("Stats for window " + str(i) + " output to " + out_file)




