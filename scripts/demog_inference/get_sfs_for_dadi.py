import sys
import argparse
import dadi

#Example usage:
#python3 get_sfs_for_dadi.py \
#-inFile '/home/vivak/rr_mu_demography_DFE/results/rr_fixed_mu_fixed/1Mb_rep1.ms' \
#-outFile '/home/vivak/rr_mu_demography_DFE/results/1.fs' \
#-region_length 1000000 -samples 100 


#Parse arguments
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-inFile', dest = 'inFile', action='store', nargs = 1, type = str, help = 'path to input file')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')
parser.add_argument('-region_length', dest = 'region_length', action='store', nargs = 1, type = int, help = 'simulated region length in bps')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'number  of individuals for which the SFS will be made')

#read input parameters
args = parser.parse_args()
inFile = args.inFile[0]
outFile = args.outFile[0]
region_length = args.region_length[0]
samples = args.samples[0]


#Function to get list of allele frequencies from .ms file
def get_af_list(inFile, samples):
    f_ms = open(inFile, 'r')
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


#Function to extract sfs from dictionary of allele frequencies
def get_sfs(l_af, samples):
    d_sfs = {}
    s_seg = 0 #total number of truly segregating sites
    s_not_anc = 0 #required to know the d0_0 class
    #Loop through list of allele frequency categories
    for x in l_af:
        try:
            #if the category already exists, incrememnt by one
            d_sfs[x] = d_sfs[x] + 1
        except:
            #otherwise create the category
            d_sfs[x] = 1
        #add to count of segregating sites if not fixed or lost
        if int(x) > 0 and int(x) < int(samples):
            s_seg += 1
        if int(x) > 0:
            s_not_anc += 1
    return(d_sfs, s_seg, s_not_anc)

#Master function to produce spectrum and output to file
def output_sfs(inFile, samples, region_length, outFile=''):
    l_af = get_af_list(inFile, samples)
    sfs_list = get_sfs(l_af, samples)
    sfs = sfs_list[0]
    sfs = dict(sorted(sfs.items()))
    s_seg = sfs_list[1]#Number of polymorphic sites, excludes sites with AF=1
    s_not_anc = sfs_list[2] #Number of polymorphic sites + number of fixed sites in the sample      
    sfs_arr = list(sfs.values())
    sfs_arr = [0] + sfs_arr
    dadi.Spectrum(sfs_arr).to_file(outFile)   


output_sfs(inFile, samples, region_length, outFile)

print ("results output to " + outFile)


