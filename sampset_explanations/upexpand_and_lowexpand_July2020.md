# Steps for generating lists of upland and lowland samples in `expandgeo`

## Overview
* Divided `expandgeo` PCA along PC1 to separate upland and lowland samples
  * upland = PC1 < -20
  * lowland = PC1 > -3
* `expand_geo` PCA results
  * Result file: `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/combo.sub.polyploid.CDS.expandgeosamps.genlight.PCAresults.rds`
* Library lists
  * `expand_upland`: `/global/cscratch1/sd/grabowsp/sg_ploidy/expandupland_libs_July2020.txt`
    * 287 samples
  * `expand_lowland`: `/global/cscratch1/sd/grabowsp/sg_ploidy/expandlowland_libs_June2020.txt`
    * 486 samples
  * 12 samples in neither upland or lowland 

## Generate Library List
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)

res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/combo.sub.polyploid.CDS.expandgeosamps.genlight.PCAresults.rds'

pca_res <- readRDS(res_file)

pca_mat <- pca_res$scores

up_samps <- rownames(pca_mat)[which(pca_mat[,1] < -20)]
# 287 - 10 more than geo_samps
low_samps <- rownames(pca_mat)[which(pca_mat[,1] > -3)]
# 486 - 4 more than geo_samps

nrow(pca_mat)-(length(up_samps)+length(low_samps))
# 12 of 785 samples are in neither group; 1 fewer than for geo_samps

up_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/expandupland_libs_July2020.txt'
write.table(up_samps, file = up_out_file, quote = F, sep = '\t', row.names = F,
  col.names = F)

low_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/expandlowland_libs_June2020.txt'
write.table(low_samps, file = low_out_file, quote = F, sep = '\t',
  row.names = F, col.names = F)

```
