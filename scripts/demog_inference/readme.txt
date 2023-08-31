See manuscript for fastsimcoal2 command line.
get_sfs_for_fsc.py and get_sfs_for_dadi.py will generate input files from slim output.
Use masked versions to mask exonic regions.


#For fastimcoal2, get parameter estimate with smallest difference between the observed and expected likelihoods
awk '{print $1,$2,$3,$3-$2}' OFS='\t' sim1_rep1.bestlhoods | sort -nk 4 | head -n1

#For population size change
awk '{print $1,$2,$3,$4,$5,$6,$6-$5}' OFS='\t' 1Mb_rep1.bestlhoods | sort -nk 7 | head -n1