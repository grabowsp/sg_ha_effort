# Analysis of inheritance patterns in Octoploid (8X) switchgrass

## Overview
### Analysis ideas
* measure LD via r^2 in different groups
* Look at levels of heterozygosity within individuals
  * compare levels of 1:3, 2:2, 3:1 genotypes
* Look for fixed HET SNPs and/or high numbers of HET loci across population

## LD analysis
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/Chr01K.polyploid.CDS.up_oct_1.vcf_00

HEAD_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/CDS.up_oct_1.vcf.header.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/

MAX_DIST=10000

MAF_CUT=0.05

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
$IN_FILE $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT
 
```
### Chr01K analysis
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups

TOP_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/

HEAD_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/CDS.up_oct_1.vcf.header.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/

MAX_DIST=10000

MAF_CUT=0.05

for OCT_1_F in $TOP_DIR'Chr01K.polyploid.CDS.up_oct_1.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $OCT_1_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

for OCT_2_F in $TOP_DIR'Chr01K.polyploid.CDS.up_oct_2.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $OCT_2_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

for TET_1_F in $TOP_DIR'Chr01K.polyploid.CDS.up_tet_1.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $TET_1_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

for TET_2_F in $TOP_DIR'Chr01K.polyploid.CDS.up_tet_2.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $TET_2_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

```
#### Consolidate results
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups

# in R

file_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/'

group_name_vec <- c('up_tet_1', 'up_tet_2', 'up_oct_1', 'up_oct_2')

for(j in seq(length(group_name_vec))){
  res_file_pre <- paste('Chr01K.polyploid.CDS.', group_name_vec[j], '_', 
    sep = '')
  r2_file_vec <- system(paste('ls ', file_dir, res_file_pre, '*', 'r2.rds', 
    sep = ''), intern = T)
  r2_tot <- readRDS(r2_file_vec[1])
  for(i in c(2:length(r2_file_vec))){
    tmp_res <- readRDS(r2_file_vec[i])
    r2_tot <- rbind(r2_tot, tmp_res)
  }
  out_file <- paste(file_dir, res_file_pre, 'tot_r2.rds', sep = '')
  saveRDS(r2_tot, out_file)
}
```
#### Generate Plot
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

# in R
library(ggplot2)
library(patchwork)

oct_1_in_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/Chr01K.polyploid.CDS.up_oct_1_tot_r2.rds'

oct_1_res <- readRDS(oct_1_in_file)

gen_window_intervals <- function()

wind_size <- 10

wind_dists <- lapply(seq(10000/10), function(x) 
  ((x-1)*wind_size)+c(1:wind_size))

wind_inds <- lapply(wind_dists, function(x) 
  which(oct_1_res$comp_dist %in% x))

num_wind_inds <- unlist(lapply(wind_inds, length))

wind_mean <- unlist(lapply(wind_inds, function(x) 
  mean(oct_1_res$r2[x], na.rm = T)))

gen_r2_windows <- function(r2_df, window_size){
  max_dist <- max(r2_df$comp_dist)
}

test_mean <- generate_dist_window_df(dist_vec = oct_1_res$comp_dist, 
  value_vec = oct_1_res$r2, window_size = 10)

per_high_ld <- function(data_df, window_size = 10, max_ld = 0.9){
  max_dist <- max(data_df$comp_dist)
  dist_mult_vec <- seq(ceiling(max_dist/window_size))
    wind_dists <- lapply(dist_mult_vec, function(x) 
    ((x-1)*window_size) + c(1:window_size))
  window_inds <- lapply(wind_dists, function(x)
    which(data_df$comp_dist %in% x))
  num_window_inds <- unlist(lapply(window_inds, length))
  n_na_inds <- unlist(lapply(window_inds, function(x) 
    sum(is.na(data_df$r2[x]))))
  n_high_ld_inds <- unlist(lapply(window_inds, function(x) 
    sum(data_df$r2[x] >= max_ld, na.rm = T)))
  per_high_ld <- n_high_ld_inds/(num_window_inds-n_na_inds)
  tot_df <- data.frame(window_pos = dist_mult_vec * window_size,
    per_high_ld = per_high_ld, n_comps = num_window_inds, n_nas = n_na_inds,
    stringsAsFactors = F)
  return(tot_df)
}

test_hi_ld <- per_high_ld(data_df = oct_1_res, window_size = 10, max_ld = 0.9)

gg_oct_1 <- ggplot(oct_1_res, aes(x = comp_dist, y = r2)) +
  geom_point() + xlab('distance between SNPs (bp)') +
  ylab('r^2') +
  ggtitle('Upland 8X Group 1, Chr01K')

oct_2_in_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/Chr01K.polyploid.CDS.up_oct_2_tot_r2.rds'

oct_2_res <- readRDS(oct_2_in_file)

gg_oct_2 <- ggplot(oct_2_res, aes(x = comp_dist, y = r2)) +
  geom_point() + xlab('distance between SNPs (bp)') +
  ylab('r^2') +
  ggtitle('Upland 8X Group 2, Chr01K')

tet_1_in_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/Chr01K.polyploid.CDS.up_tet_1_tot_r2.rds'

tet_1_res <- readRDS(tet_1_in_file)

gg_tet_1 <- ggplot(tet_1_res, aes(x = comp_dist, y = r2)) +
  geom_point() + xlab('distance between SNPs (bp)') +
  ylab('r^2') +
  ggtitle('Upland 4X Group 1, Chr01K')

tet_2_in_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/Chr01K.polyploid.CDS.up_tet_2_tot_r2.rds'

tet_2_res <- readRDS(tet_2_in_file)

gg_tet_2 <- ggplot(tet_2_res, aes(x = comp_dist, y = r2)) +
  geom_point() + xlab('distance between SNPs (bp)') +
  ylab('r^2') +
  ggtitle('Upland 4X Group 2, Chr01K')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/Chr01K.polyploid.CDS.upgroups_r2.png'

png(out_file, width = 480*2, height = 400*2)
(gg_oct_1 + gg_oct_2) / (gg_tet_1 + gg_tet_2)
dev.off()




```
 
