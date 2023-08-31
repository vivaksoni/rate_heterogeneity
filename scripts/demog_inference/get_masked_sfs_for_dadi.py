import sys
import argparse
import numpy as np
import pandas as pd
import dadi

#Example usage:
#python3 get_masked_sfs_for_fsc.py -inFile '/home/vivak/rr_mu_demography_DFE/results/rr_fixed_mu_fixed/sim1_gene1_rep1.ms' \
#-outFile '/home/vivak/rr_mu_demography_DFE/results/sim1_gene1_rep1_DAFpop0.obs' \
#-region_length 96876 -samples 100

#for i in "${demog[@]}"; do for j in "${model[@]}"; do for k in $(seq 1 100); do python3 get_sfs_for_fsc.py -inFile ../results/demog_only/"$i"/"$j"/sim"$k"_gene1_rep"$k".ms \
#-outFile fsc_inputFiles/sfs_files/"$i"/"$j"/sim"$k"_gene1_rep"$k"_DAFpop0.obs -region_length 100000 -samples 100; done; done; done

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

#Function to parse .ms file data into df of positions and DAC
#Function takes as input open reading file object and no. of samples
def get_DAC(f_ms, samples, chr_len):
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
        x = int(np.round(float(pos)*chr_len,0))
        l_Pos.append(int(x))
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
    
    #Sum genotype information to get DACs
    for i,j in enumerate(l_data):
        l_data[i].append(sum([int(x) for x in l_data[i][1]]))
    
    #Convert dictionary to dataframe, removing genotype column
    df = pd.DataFrame(l_data)[[0,2]]
    #Rename columns
    df.columns = ['position', 'DAC']
    df['position'] = df.position / chr_len
    return(df)


def get_masked_af_list(inFile, samples, chr_len):
    #Open ms file
    f_ms = open(inFile, 'r')
    #Get df of derived allele counts by position
    df = get_DAC(f_ms, samples, chr_len)

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
    #Remove masked sites, and convert derived allele counts column into list
    l_af_masked = list(df[df['masked']==0].DAC)
    return(l_af_masked)



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
def output_sfs(inFile, samples, chr_len, outFile=''):
    l_af = get_masked_af_list(inFile, samples, chr_len)
    sfs_list = get_sfs(l_af, samples)
    sfs = sfs_list[0]
    sfs = dict(sorted(sfs.items()))
    s_seg = sfs_list[1]#Number of polymorphic sites, excludes sites with AF=1
    s_not_anc = sfs_list[2] #Number of polymorphic sites + number of fixed sites in the sample      
    sfs_arr = list(sfs.values())
    sfs_arr = [0] + sfs_arr
    dadi.Spectrum(sfs_arr).to_file(outFile) 


#Chromosome structure EDIT THESE!!!
chr_len = region_length
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
nc_region_length = chr_len -  (4 * 588 * 127)


output_sfs(inFile, samples, chr_len, outFile)

print ("results output to " + outFile)
