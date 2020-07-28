# Selecting Lowland and Upland Samples from the 'geo_samp' sample set

## Overview
* Divided `geo_samp` PCA along PC1 to seprate upland and lowland samples
  * upland = PC1 > 20
  * lowland = PC1 < 0

* `geo_samp` PCA results
  * Result file: `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/combo.sub.polyploid.CDS.geosamps.genlight.PCAresults.rds`
  * Figure:

## Generate Library List
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)

res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/combo.sub.polyploid.CDS.geosamps.genlight.PCAresults.rds'

pca_res <- readRDS(res_file)

pca_mat <- pca_res$scores

up_samps <- rownames(pca_mat)[which(pca_mat[,1] > 20)]
# 277
low_samps <- rownames(pca_mat)[which(pca_mat[,1] < 0)]
# 482

# 13 of 772 samples are in neither group

up_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/upland_libs_June2020.txt'
write.table(up_samps, file = up_out_file, quote = F, sep = '\t', row.names = F,
  col.names = F)

low_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/lowland_libs_June2020.txt'
write.table(low_samps, file = low_out_file, quote = F, sep = '\t',
  row.names = F, col.names = F)
```
