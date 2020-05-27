# Pop structure analysis using geographic samples for subset of data

## `geo` Chr01K CDS Analysis
* Only includes "Natural Collection" 8X samples
* Using polyploid genotype
  * 4X get disomic genotypes, 8X get tetrasomic genotypes
### Generate distance matrices
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists

sbatch gen_geo_poly_dists_test.sh
sbatch gen_geo_poly_dists_Chr01K.sh
```
### Combine matrices
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists/
FILE_SUF=ploidy_DistMat.rds
OUT_FILE=Chr01K.polyploid.CDS.geosamps.ploidy_DistMat.total.rds

Rscript $SCRIPT_DIR'/combine_dists.r' $DATA_DIR $FILE_SUF $OUT_FILE
```
### General PCoA Figures
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists

TOT_DIST_MAT=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists/Chr01K.polyploid.CDS.geosamps.ploidy_DistMat.total.rds

Rscript $SCRIPT_DIR'/make_general_PCoA_figs.r' \
$TOT_DIST_MAT \
manhattan 826_geo_samps_polyploid_genotypes

Rscript $SCRIPT_DIR'/make_general_PCoA_figs.r' \
$TOT_DIST_MAT \
euclidean 826_geo_samps_polyploid_genotypes

```
### NJ Trees
```
module load python3/3.7-anaconda-2019.10
source activate r_phylo_tools

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists

TOT_DIST_MAT=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists/Chr01K.polyploid.CDS.geosamps.ploidy_DistMat.total.rds

Rscript $SCRIPT_DIR'/make_standard_NJ_trees.r' \
$TOT_DIST_MAT \
manhattan

Rscript $SCRIPT_DIR'/make_standard_NJ_trees.r' \
$TOT_DIST_MAT \
euclidean
```

## `expand_geo` Chr01K CDS Analysis
* Includes "Natural Collection" and "Cultivar" 8X samples
* Using polyploid genotypes
### Generate distance matrices
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists

sbatch gen_expandgeo_poly_dists_test.sh
sbatch gen_expandgeo_poly_dists_Chr01K.sh
```
### Combine matrices
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists/
FILE_SUF=ploidy_DistMat.rds
OUT_FILE=Chr01K.polyploid.CDS.expandgeosamps.ploidy_DistMat.total.rds

Rscript $SCRIPT_DIR'/combine_dists.r' $DATA_DIR $FILE_SUF $OUT_FILE
```
### General PCoA Figures
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists

TOT_DIST_MAT=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists/Chr01K.polyploid.CDS.expandgeosamps.ploidy_DistMat.total.rds

Rscript $SCRIPT_DIR'/make_general_PCoA_figs.r' \
$TOT_DIST_MAT \
manhattan 839_expanded_geo_samps_polyploid_genotypes

Rscript $SCRIPT_DIR'/make_general_PCoA_figs.r' \
$TOT_DIST_MAT \
euclidean 839_expanded_geo_samps_polyploid_genotypes

```
### NJ Trees
```
module load python3/3.7-anaconda-2019.10
source activate r_phylo_tools

SCRIPT_DIR=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists

TOT_DIST_MAT=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists/Chr01K.polyploid.CDS.expandgeosamps.ploidy_DistMat.total.rds

Rscript $SCRIPT_DIR'/make_standard_NJ_trees.r' \
$TOT_DIST_MAT \
manhattan

Rscript $SCRIPT_DIR'/make_standard_NJ_trees.r' \
$TOT_DIST_MAT \
euclidean
```





