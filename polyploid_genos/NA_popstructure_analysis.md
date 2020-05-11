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
## Combine distance Matrices
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

Rscript \
/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/combine_NA_dists.r \
/global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/no_filtering \
NA_DistMat.rds Chr01K.polyploid.CDS.allsamps.NA_DistMat.total.rds 

Rscript \
/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/combine_NA_dists.r \
/global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/few_missing \
NA_DistMat.rds Chr01K.polyploid.CDS.allsamps.few_miss.NA_DistMat.total.rds

Rscript \
/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/combine_NA_dists.r \
/global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/half_missing \
NA_DistMat.rds Chr01K.polyploid.CDS.allsamps.half_miss.NA_DistMat.total.rds

Rscript \
/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/combine_NA_dists.r \
/global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/high_missing \
NA_DistMat.rds Chr01K.polyploid.CDS.allsamps.high_miss.NA_DistMat.total.rds
```

## Generate PCoA figure(s) from NA distance Matrices
* R script
  * `/GIT_PATH/sg_ha_effort/polyploid_genos/NA_PCoA_figs.r`
* PCo 1 vs 2 figure
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/Chr01K.CDS.NA_PCo_1v2.pdf`


