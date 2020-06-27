# Notes about selecting genetically-similar subgroups based on PCA results

## Overview
* Use plots from PCA results to visually set cutoffs for selecting genetically
similar samples
* Note: these are NOT necessarily true groupings, they are just groups that
    are similar samples

## Genetically-Upland samples
### Criteria for samples in each upland sub-group
up_tet_1 <- PC1>5, PC1 < 10, PC2 > 17.5, PC2 < 22, PC4 > 1, PC4 < 6
up_tet_2 <- PC1 > 17, PC2 < -17
up_oct_1 <- PC1<-25, PC3 > 10
up_oct_2 <- PC1 < -15, PC2 < -15
### Generate lists of libraries for each sub-gproup
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)

pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/combo.sub.polyploid.CDS.upgeosamps.genlight.PCAresults.rds'

up_pca_tot <- readRDS(pca_res_file)

up_pca <- up_pca_tot$scores

up_tet_1 <- intersect(which(up_pca[,1] > 5), intersect(which(up_pca[,1] < 10),
  intersect(which(up_pca[,2] > 17.5), intersect(which(up_pca[,2] < 22),
  intersect(which(up_pca[,4] > 1), which(up_pca[,4] < 6))))))
  
up_tet_2 <- intersect(which(up_pca[,1] > 17), which(up_pca[,2] < -17))

up_oct_1 <- intersect(which(up_pca[,1] < -25), which(up_pca[,3] > 10))
up_oct_2 <- intersect(which(up_pca[,1] < -15), which(up_pca[,2] < -15))

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
# 38 - Midwest: 15 IL, 3 IN, 15 MI, 5 WI

up_tet_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_1_libs.txt', sep = '')
write.table(meta$LIBRARY[up_tet_1_meta_inds], file = up_tet_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

up_tet_2_meta_inds <- c()
for(i in up_tet_2){
  tmp_ind <- which(meta$LIBRARY == rownames(up_pca)[i])
  up_tet_2_meta_inds <- c(up_tet_2_meta_inds, tmp_ind)
}
# 33 - North Dacotah
up_tet_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_2_libs.txt', sep = '')
write.table(meta$LIBRARY[up_tet_2_meta_inds], file = up_tet_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

up_oct_1_meta_inds <- c()
for(i in up_oct_1){
  tmp_ind <- which(meta$LIBRARY == rownames(up_pca)[i])
  up_oct_1_meta_inds <- c(up_oct_1_meta_inds, tmp_ind)
}
#33 - east, 19 PA, 2 WV, 9 NY, 1 NJ, 2 MA

up_oct_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_oct_1_libs.txt', sep = '')
write.table(meta$LIBRARY[up_oct_1_meta_inds], file = up_oct_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

up_oct_2_meta_inds <- c()
for(i in up_oct_2){
  tmp_ind <- which(meta$LIBRARY == rownames(up_pca)[i])
  up_oct_2_meta_inds <- c(up_oct_2_meta_inds, tmp_ind)
}
#13 - west/all over [this might be the Blackwell clade?]
up_oct_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_oct_2_libs.txt', sep = '')
write.table(meta$LIBRARY[up_oct_2_meta_inds], file = up_oct_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)
```

## Genetically-lowland samples
### Criteria for samples in each lowland sub-group
* Note: these are NOT necessarily true groupings, they are just groups that
    are similar samples
low_tx_1 <- PC1 < -40, PC2 > -19, PC2 < -11
low_tx_2 <- PC1 < -53, PC2 < -23 
low_gc_1 <- PC2 > 40
low_gc_2 <- PC1 < -20, PC2 > 20, PC2 < 27
low_ec_1 <- PC1 > 20, PC3 > 19
low_ec_2 <- PC1 > 21, PC3 < -15
### Generate lists of Libraries in each subgroup
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)

pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/combo.sub.polyploid.CDS.lowgeosamps.genlight.PCAresults.rds'

low_pca_tot <- readRDS(pca_res_file)

low_pca <- low_pca_tot$scores

low_tx_1 <- intersect(which(low_pca[,1] < -40), 
  intersect(which(low_pca[,2] > -19), which(low_pca[,2] < -11)))

low_tx_2 <- intersect(which(low_pca[,1] < -53), which(low_pca[,2] < -23))

low_gc_1 <- which(low_pca[,2] > 40)

low_gc_2 <- intersect(which(low_pca[,1] < -20), 
  intersect(which(low_pca[,2] > 20), which(low_pca[,2] < 27)))

low_ec_1 <- intersect(which(low_pca[,1] > 20), which(low_pca[,3] > 19))

low_ec_2 <- intersect(which(low_pca[,1] > 21), which(low_pca[,3] < -15))  

sum(duplicated(c(low_tx_1, low_tx_2, low_gc_1, low_gc_2, low_ec_1, low_ec_2)))
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
# 47
table(meta$STATE[low_tx_1_meta_inds])
# 25 TX 2 AR, 10 MX, 1 KS, 1 LA, 1 MS, 1 NC, 1 SC

low_tx_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_1_libs.txt', sep = '')
write.table(meta$LIBRARY[low_tx_1_meta_inds], file = low_tx_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_tx_2_meta_inds <- c()
for(i in low_tx_2){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_tx_2_meta_inds <- c(low_tx_2_meta_inds, tmp_ind)
}
# 17
table(meta$STATE[low_tx_2_meta_inds])
sum(is.na(meta$STATE[low_tx_2_meta_inds]))
# OK 7, TX 2, KS 3, CO 1, AR 1, 2 Kanlow
low_tx_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_2_libs.txt', sep = '')
write.table(meta$LIBRARY[low_tx_2_meta_inds], file = low_tx_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_gc_1_meta_inds <- c()
for(i in low_gc_1){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_gc_1_meta_inds <- c(low_gc_1_meta_inds, tmp_ind)
}
#31
table(meta$STATE[low_gc_1_meta_inds])
# 16 MS, 6 LA, 8 LF, 1 AR
low_gc_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_gc_1_libs.txt', sep = '')
write.table(meta$LIBRARY[low_gc_1_meta_inds], file = low_gc_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_gc_2_meta_inds <- c()
for(i in low_gc_2){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_gc_2_meta_inds <- c(low_gc_2_meta_inds, tmp_ind)
}
# 18
table(meta$STATE[low_gc_2_meta_inds])
# 11 AR, 3 FL, 3 LA, 1 TX
low_gc_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_gc_2_libs.txt', sep = '')
write.table(meta$LIBRARY[low_gc_2_meta_inds], file = low_gc_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_ec_1_meta_inds <- c()
for(i in low_ec_1){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_ec_1_meta_inds <- c(low_ec_1_meta_inds, tmp_ind)
}
#63
table(meta$STATE[low_ec_1_meta_inds])
# 15 RI, 32 NY, 4 CT, 3 ME, 7 MA, 2 NH
low_ec_1_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_1_libs.txt', sep = '')
write.table(meta$LIBRARY[low_ec_1_meta_inds], file = low_ec_1_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)

low_ec_2_meta_inds <- c()
for(i in low_ec_2){
  tmp_ind <- which(meta$LIBRARY == rownames(low_pca)[i])
  low_ec_2_meta_inds <- c(low_ec_2_meta_inds, tmp_ind)
}
# 59
table(meta$STATE[low_ec_2_meta_inds])
# 3 GA, 22 ML, 15 NC, 15 SC, 4 VI
low_ec_2_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_2_libs.txt', sep = '')
write.table(meta$LIBRARY[low_ec_2_meta_inds], file = low_ec_2_out_file,
  quote = F, sep = '\t', row.names = F, col.names = F)
```


