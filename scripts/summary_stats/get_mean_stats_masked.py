#script takes .stats file and .ms file and gets mean of summary statistics across windows for each replicate

#Example usage:
#python3  get_mean_stats_masked.py \
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


#Chromosome structure EDIT THESE!!!
chr_len = 997204
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
        #Read in .ms file, create sd object for libsequence
        #f_ms = open(ms_file, 'r')
        l_data = get_nested_data_list(f, 100)

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
        segSites = len(df)
        #Get no. of singletons by counting how many sites have only a single change ('1')
        singletons = 0
        for site in l_data:
            if (site[1].count('1')==1):
                singletons+=1

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

