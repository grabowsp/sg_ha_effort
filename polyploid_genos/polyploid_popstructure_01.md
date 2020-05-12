# Preliminary Population Structure Analysis of Polyploid Genotypes

## Overview
* Compare ploidy-appropriate, disomic, and tetrasomic genotypes
* Calculate distance matrices
* Analyses:
  * PCoA
  * NJ Trees
  * Heterozygosity
  * Pi

## Generate Chr01K CDS Distance Matrix
* generate matrix for three genotype approaches
  * "polyploid" = 4X get disomic genotypes, 8X get tetrasomic genotypes
  * "disomic" = all samples get disomic genotypes
  * "tetrasomic" = all samples get tetrasomic genotypes
* Distance matrix from genotype dosages ranging from 2 (homozygous Ref) to 
0 (homozygous Alt) 
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/polyploid_dists

sbatch gen_poly_dists_test.sh
sbatch gen_poly_dists_Chr01K.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/disomic_dists

sbatch gen_disomic_dists_test.sh
sbatch gen_disomic_dists_Chr01K.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/tetrasomic_dists

sbatch gen_tetrasomic_dists_test.sh
sbatch gen_tetrasomic_dists_Chr01K.sh
```

## Combine Chr01K distance matrices
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

# polyploid distances

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/polyploid_dists

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/polyploid_dists/
FILE_SUF=ploidy_DistMat.rds
OUT_FILE=Chr01K.polyploid.CDS.allsamps.few_miss_ploidy_DistMat.total.rds

Rscript $SCRIPT_DIR'/combine_dists.r' $DATA_DIR $FILE_SUF $OUT_FILE

# disomic distances

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/disomic_dists

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/disomic_dists/
FILE_SUF=disomic_DistMat.rds
OUT_FILE=Chr01K.polyploid.CDS.allsamps.few_miss_disomic_DistMat.total.rds

Rscript $SCRIPT_DIR'/combine_dists.r' $DATA_DIR $FILE_SUF $OUT_FILE

# tetrasomic distances

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/tetrasomic_dists

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/tetrasomic_dists/
FILE_SUF=tetrasomic_DistMat.rds
OUT_FILE=Chr01K.polyploid.CDS.allsamps.few_miss_tetrasomic_DistMat.total.rds

Rscript $SCRIPT_DIR'/combine_dists.r' $DATA_DIR $FILE_SUF $OUT_FILE
```

## Generate General PCoA Figures
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

# polyploid distances

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/polyploid_dists

Rscript $SCRIPT_DIR'/combine_dists.r' \
/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/polyploid_dists/Chr01K.polyploid.CDS.allsamps.few_miss_ploidy_DistMat.total.rds



```

## Next Steps
* `combine_dists.r` to combine the distance matrices from each data type
* `make_general_PCoA_figs.r` for each data type
  * `R_analysis`
* `make_standard_NJ_trees.r` for each data type
  * `r_phylo_tools`
