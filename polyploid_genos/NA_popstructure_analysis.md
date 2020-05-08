# Info about using NA's to analyze population structure

## Overview
* Convert genotypes to 1=genotype; 0=NA
* Generate distance matrices based on 1's and 0's
* Look at PCoA and NJ trees
* Look at --max-missing 0 (no filter), 0.2 (80% missing data), 0.5 (50% 
missing data) and 0.8 (20% missing data)

## Generate NA Distance Matrices
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/no_filtering

sbatch gen_NA_dists_test.sh
sbatch gen_Chr01K_NA_dists.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/few_missing
sbatch gen_Chr01K_few_miss_NA_dists.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/half_missing
sbatch gen_Chr01K_half_miss_NA_dists.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/high_missing
sbatch gen_Chr01K_high_miss_NA_dists.sh
```




