import sys
import argparse
import operator
import pandas as pd
import numpy as np

#Example usage:
#python3 get_input_for_Grampes.py -inPath '/home/vivak/rr_mu_demography_DFE/results/rr_fixed_mu_fixed/sim1_gene1_rep1.ms' \
#-outFile '/home/vivak/rr_mu_demography_DFE/results/sim1_gene1_rep1_DAFpop0.obs' \
#-region_length 96876 -samples 100 -minRep 1 -maxRep 100

#for i in "${demog[@]}"; do for j in "${model[@]}"; do for k in $(seq 1 100); do python3 get_sfs_for_fsc.py -inFile ../results/demog_only/"$i"/"$j"/sim"$k"_gene1_rep"$k".ms \
#-outFile fsc_inputFiles/sfs_files/"$i"/"$j"/sim"$k"_gene1_rep"$k"_DAFpop0.obs -region_length 100000 -samples 100; done; done; done

#Parse arguments
parser = argparse.ArgumentParser(description='Information about creating Grapes input file from SLiM output')
parser.add_argument('-inPath', dest = 'inPath', action='store', nargs = 1, type = str, help = 'path to input files (such that _rep"$replicate".ms is added)')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')
parser.add_argument('-region_length', dest = 'region_length', action='store', nargs = 1, type = int, help = 'simulated region length in bps')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'number  of individuals for which the SFS will be made')
parser.add_argument('-minRep', dest = 'minRep', action='store', nargs = 1, type = int, help = 'first replicate to begin loop')
parser.add_argument('-maxRep', dest = 'maxRep', action='store', nargs = 1, type = int, help = 'last replicate of loop')

#read input parameters
args = parser.parse_args()
inPath = args.inPath[0]
outFile = args.outFile[0]
chr_len = args.region_length[0]
samples = args.samples[0]
minRep = args.minRep[0]
maxRep = args.maxRep[0]


#Function to parse .ms file data into nested list of positions and genotypes
#Function takes as input open reading file object and no. of samples
def get_nested_data_list(f_ms, samples, chr_len):
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
    return(df)


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


#Function to combine dictionaries
def combine_dicts(a, b, op=operator.add):
    return {**a, **b, **{k: op(a[k], b[k]) for k in a.keys() & b}}


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
        exons.append([exonStart, exonEnd])

        if (element < exonsPerGene):
            intronStart = exonStart + exonLength
            intronEnd = intronStart+(intronLength-1)
            introns.append([intronStart, intronEnd])

    intergenicStart = exonEnd + 1
    intergenicEnd = intergenicStart + (intergenicLength-1)
    intergenic.append([intergenicStart, intergenicEnd])

#Get total exonic region length
exonic_region_length = exonsPerGene * exonLength * totalGenes



#Site counters
s_sites = 0
ns_sites = 0
#Create empty dicts to store NS and S sfs
d_ns = {}
d_s = {}
#Loop through replicates
for i in range(minRep, maxRep+1):
    ms_file = inPath + '_rep' + str(i) + '.ms'
    f_ms = open(ms_file, 'r')
    #Use function to get list of positions and genotypes
    df = get_nested_data_list(f_ms, samples, chr_len)

    #Create empty list to store masking status
    res = []
    for pos in df.position:
        #Set tmp to 0 (ie unmasked site)
        tmp = 0
        #Loop through exons, checking if site is exonic
        for exon in exons:
            if((pos>=exon[0]) & (pos<=exon[1])):
            #If site is exonic, set tmp to 1 (ie exonic site) and break out of loop (move to next site)
                tmp = exon[0]
                break
        #Append mask status to results list
        res.append(tmp)
    #Add column of masking status to df
    df['masked'] = res
    #Subset df to include only exonic sites
    df = df[df.masked>0]
    #Determine whether snp is NS or S
    df['pos_in_exon'] = df.position - df.masked + 1
    df['csq'] = np.where(df.pos_in_exon%3==0, 'S', 'NS')

    ns = list(df[df.csq=='NS'].DAC)
    s = list(df[df.csq=='S'].DAC)
    #Get NS and S sfs
    d_ns2 = get_sfs(ns, samples)
    d_s2 = get_sfs(s, samples)
    #Combine current sfs with main sfs
    d_ns = combine_dicts(d_ns, d_ns2[0])
    d_s = combine_dicts(d_s, d_s2[0])
    #Get site counts
    ns_sites = ns_sites + (exonic_region_length - (exonic_region_length/3))
    s_sites = s_sites + (exonic_region_length/3)

#Write to file in Grapes-friendly format    
with open(outFile, 'w') as f:
    f.write(outFile + '\n#unfolded\n')
    f.write('all_genes')
    f.write('\t')
    f.write(str(samples))
    f.write('\t')
    f.write(str(int(ns_sites)))
    f.write('\t')
    for i in range(1, samples):
        f.write(str(d_ns.get(i, 0)))
        f.write('\t')
    f.write(str(int(s_sites)))
    f.write('\t')
    for i in range(1, samples):
        f.write(str(d_s.get(i, 0)))
        f.write('\t')        