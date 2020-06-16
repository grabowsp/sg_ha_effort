# Notes about PCA-based population structure analysis using adegenet

## Overview
* using `adegenet` package in R to run PCA for 4X/8X population structure
analysis
* As of June 1 2020, the PCA function in `adegenet` was not working properly
on Cori (there was an error with the eigen() function), so current workflow is:
  * Generate genlight object on Cori (documented in other .md file)
  * Transfer genlight object to HA
  * Run PCA on HA
  * Transfer PCA results to Cori
  * PCA-based figures and analysis on Cori

## PCA on all geo samples
* uses the `adegenet_pca.r` script
### Location of genlight object
* on Cori
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.rds`
* on HA
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.595K.polyploid.CDS.geosamps.genlight.rds`
### Run PCA on HA
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/

qsub combo595K_PCA_submit.sh
```
### Location of PCA results
* on Cori
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds`
* on HA
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds`
### Generate PCA figures
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

DATA_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds

PLOT_PRE=Geographic_samples_595K_SNPs

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE

```
### Divide Samples into genetic ecotype groups
* upland = PC1 < -25
* lowland = PC1 > 0

## PCA on Upland Geo samples
### Location of genlight object
* on Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/Combo.sub.polyploid.CDS.upgeosamps.genlight.rds`
* on HA:
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.sub.polyploid.CDS.upgeosamps.genlight.rds`
### Run PCA on HA
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/

qsub combo_upgeo_PCA_submit.sh
```
### Location of PCA results
* On Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/Combo.sub.polyploid.CDS.upgeosamps.genlight.PCAresults.rds`
* On HA:
  * '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.sub.polyploid.CDS.upgeosamps.genlight.PCAresults.rds'
### Generate PCA Figures
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/

DATA_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/Combo.sub.polyploid.CDS.upgeosamps.genlight.PCAresults.rds

PLOT_PRE=Upland_geo_samps

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE
```
### Split samples into preliminary groups for testing r^2
* Note: these are NOT necessarily true groupings, they are just groups that
    are similar samples
up_tet_1 <- PC1>0, PC2 > 12.5
up_tet_2 <- PC1 > 0, PC2 < 15
up_oct_1 <- PC1<0, PC3 > 15
up_oct_2 <- PC1 < 0, PC2 > 15
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)

pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/Combo.sub.polyploid.CDS.upgeosamps.genlight.PCAresults.rds'

up_pca_tot <- readRDS(pca_res_file)

up_pca <- up_pca_tot$scores

up_tet_1 <- intersect(which(up_pca[,1] > 0), which(up_pca[,2] > 12.5))
up_tet_2 <- intersect(which(up_pca[,1] > 0), which(up_pca[,2] < -15))

up_oct_1 <- intersect(which(up_pca[,1] < 0), which(up_pca[,3] > 15))
up_oct_2 <- intersect(which(up_pca[,1] < 0), which(up_pca[,2] > 15))

ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'ploidy_calling/sg_ploidy_results_v3.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/Reseq_Metadata_Sept_2019_Edited_for_R.tsv'

meta <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F,
  quote = "", comment.char = '$')

up_tet_1_meta_inds <- c()
for(i in up_tet_1){
  tmp_ind <- which(meta$LIBRARY == rownames(up_pca)[i])
  up_tet_1_meta_inds <- c(up_tet_1_meta_inds, tmp_ind)
}
# 60 - North Dakota

up_tet_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_1_libs.txt', sep = '')
write.table(meta$LIBRARY[up_tet_1_meta_inds], file = up_tet_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)


up_tet_2_meta_inds <- c()
for(i in up_tet_2){
  tmp_ind <- which(meta$LIBRARY == rownames(up_pca)[i])
  up_tet_2_meta_inds <- c(up_tet_2_meta_inds, tmp_ind)
}
# 67 - all over
up_tet_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_2_libs.txt', sep = '')
write.table(meta$LIBRARY[up_tet_2_meta_inds], file = up_tet_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

up_oct_1_meta_inds <- c()
for(i in up_oct_1){
  tmp_ind <- which(meta$LIBRARY == rownames(up_pca)[i])
  up_oct_1_meta_inds <- c(up_oct_1_meta_inds, tmp_ind)
}
#16 - east
up_oct_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_oct_1_libs.txt', sep = '')
write.table(meta$LIBRARY[up_oct_1_meta_inds], file = up_oct_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

up_oct_2_meta_inds <- c()
for(i in up_oct_2){
  tmp_ind <- which(meta$LIBRARY == rownames(up_pca)[i])
  up_oct_2_meta_inds <- c(up_oct_2_meta_inds, tmp_ind)
}
#18 - west/all over
up_oct_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_oct_2_libs.txt', sep = '')
write.table(meta$LIBRARY[up_oct_2_meta_inds], file = up_oct_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)
```


## PCA on Lowland Geo Samples
### Location of genlight object
* on Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/Combo.sub.polyploid.CDS.lowgeosamps.genlight.rds`
* on HA:
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.sub.polyploid.CDS.lowgeosamps.genlight.rds`
### Run PCA on HA
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/

qsub combo_lowgeo_PCA_submit.sh
```
### Location of PCA results
* On Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/Combo.sub.polyploid.CDS.lowgeosamps.genlight.PCAresults.rds`
* On HA:
  * '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.sub.polyploid.CDS.lowgeosamps.genlight.PCAresults.rds'
### Generate PCA Figures
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/

DATA_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/Combo.sub.polyploid.CDS.lowgeosamps.genlight.PCAresults.rds

PLOT_PRE=Lowland_geo_samps

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE
```
### Split samples into preliminary groups for testing r^2
* Note: these are NOT necessarily true groupings, they are just groups that
    are similar samples
low_tx_1 <- PC1 < -40, PC2 < -15
low_tx_2 <- PC1 < -40, PC2 > -15, PC2 < 5 
low_gc_1 <- PC2 > 35
low_ec_1 <- PC1 > 20, PC3 < -22
low_ec_2 <- PC1 > 20, PC3 > -3, PC3 < 6
* Some additional notes:
  * there are several lines with "NA" in meta$STATE that were used here
    * ex: the NAM lines from Noble, Timber, etc
  * These are synthetic/breeding lines and I don't think should be used in
climate or "natural" population analyses
  * I think the issue came because I used "LIB_CLIMATE" but I should
have used "LIB_BIOCLIM"

```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)

pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/Combo.sub.polyploid.CDS.lowgeosamps.genlight.PCAresults.rds'

low_pca_tot <- readRDS(pca_res_file)

low_pca <- low_pca_tot$scores

low_tx_1 <- intersect(which(low_pca[,1] < -40), which(low_pca[,2] < -15))
low_tx_2 <- intersect(which(low_pca[,1] < -40), 
  intersect(which(low_pca[,2] > -15), which(low_pca[,2] < 5)))

low_gc_1 <- which(low_pca[,2] > 35)

low_ec_1 <- intersect(which(low_pca[,1] > 20), which(low_pca[,3] < -22))
low_ec_2 <- intersect(which(low_pca[,1] > 20), 
  intersect(which(low_pca[,3] > -3), which(low_pca[,3] < 6)))

sum(duplicated(c(low_tx_1, low_tx_2, low_gc_1, low_ec_1, low_ec_2)))
#0

ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'ploidy_calling/sg_ploidy_results_v3.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/Reseq_Metadata_Sept_2019_Edited_for_R.tsv'

meta <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F,
  quote = "", comment.char = '$')

geno_meta_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/genotype.metadata.May2020.rds'

geno_meta <- readRDS(geno_meta_file)

low_tx_1_meta_inds <- c()
for(i in low_tx_1){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_tx_1_meta_inds <- c(low_tx_1_meta_inds, tmp_ind)
}
# 43 
table(meta$STATE[low_tx_1_meta_inds]), rest are 6 AK, 3 KS, 7 OK, 5 TX, 1 CO
# these are 1/2 NAM samples from NOBLE
low_tx_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_1_libs.txt', sep = '')
write.table(meta$LIBRARY[low_tx_1_meta_inds], file = low_tx_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_tx_2_meta_inds <- c()
for(i in low_tx_2){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_tx_2_meta_inds <- c(low_tx_2_meta_inds, tmp_ind)
}
# 73
table(meta$STATE[low_tx_2_meta_inds])
sum(is.na(meta$STATE[low_tx_2_meta_inds]))
# mainly TX and MX
low_tx_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_2_libs.txt', sep = '')
write.table(meta$LIBRARY[low_tx_2_meta_inds], file = low_tx_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_gc_1_meta_inds <- c()
for(i in low_gc_1){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_gc_1_meta_inds <- c(low_gc_1_meta_inds, tmp_ind)
}
#33
table(meta$STATE[low_gc_1_meta_inds])
# 16 MS, 8 LA, 8 LF, 1 AK
low_gc_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_gc_1_libs.txt', sep = '')
write.table(meta$LIBRARY[low_gc_1_meta_inds], file = low_gc_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_ec_1_meta_inds <- c()
for(i in low_ec_1){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_ec_1_meta_inds <- c(low_ec_1_meta_inds, tmp_ind)
}
#46
table(meta$STATE[low_ec_1_meta_inds])
# 14 RI, 15 NY, 7 CT, 3 ME, 3 MA, 3 NH
low_ec_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_1_libs.txt', sep = '')
write.table(meta$LIBRARY[low_ec_1_meta_inds], file = low_ec_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_ec_2_meta_inds <- c()
for(i in low_ec_2){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_ec_2_meta_inds <- c(low_ec_2_meta_inds, tmp_ind)
}
#52
table(meta$STATE[low_ec_2_meta_inds])
# 26 NY, 15 NJ, 6 DE, 2 ML, 1 NC, 2 VI
low_ec_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_2_libs.txt', sep = '')
write.table(meta$LIBRARY[low_ec_2_meta_inds], file = low_ec_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)


```


