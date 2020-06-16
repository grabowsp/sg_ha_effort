# Script for analyzing r^2 patterns

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)
library(patchwork)

general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)

### LOAD DATA ###
oct_1_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/', 
  'Chr01K.polyploid.CDS.up_oct_1_tot_r2.rds', sep = '')
oct_1_res <- readRDS(oct_1_res_file)

oct_2_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/',
  'Chr01K.polyploid.CDS.up_oct_2_tot_r2.rds', sep = '')
oct_2_res <- readRDS(oct_2_res_file)

tet_1_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/',
  'Chr01K.polyploid.CDS.up_tet_1_tot_r2.rds', sep = '')
tet_1_res <- readRDS(tet_1_res_file)

tet_2_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/',
  'Chr01K.polyploid.CDS.up_tet_2_tot_r2.rds', sep = '')
tet_2_res <- readRDS(tet_2_res_file)

### SET OUTPUTS ###
r2_out_fig_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/', 
  'Chr01K.polyploid.CDS.upgroups_r2.pdf', sep = '')

hi_ld_out_fig_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/', 
  'Chr01K.polyploid.CDS.upgroups_hiLD.pdf', sep = '')

### SET VARIABLES ###


######
### FUNCTIONS ###
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

###########

oct_1_ld_wind <- generate_dist_window_df(dist_vec = oct_1_res$comp_dist,
  value_vec = oct_1_res$r2, window_size = 10)

oct_1_hi_ld <- per_high_ld(data_df = oct_1_res, window_size = 10, max_ld = 0.9)

oct_2_ld_wind <- generate_dist_window_df(dist_vec = oct_2_res$comp_dist,
  value_vec = oct_2_res$r2, window_size = 10)

oct_2_hi_ld <- per_high_ld(data_df = oct_2_res, window_size = 10, max_ld = 0.9)

tet_1_ld_wind <- generate_dist_window_df(dist_vec = tet_1_res$comp_dist,
  value_vec = tet_1_res$r2, window_size = 10)

tet_1_hi_ld <- per_high_ld(data_df = tet_1_res, window_size = 10, max_ld = 0.9)

tet_2_ld_wind <- generate_dist_window_df(dist_vec = tet_2_res$comp_dist,
  value_vec = tet_2_res$r2, window_size = 10)

tet_2_hi_ld <- per_high_ld(data_df = tet_2_res, window_size = 10, max_ld = 0.9)

###

gg_oct1_r2 <- ggplot(oct_1_ld_wind, aes(x = window_pos, y = window_avg)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('r^2') +
  ylim(0,0.5) +
  ggtitle('r^2 in Upland 8X Group 1\nChr01K, 10bp windows')

gg_oct2_r2 <- ggplot(oct_2_ld_wind, aes(x = window_pos, y = window_avg)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('r^2') +
  ylim(0,0.5) +
  ggtitle('r^2 in Upland 8X Group 2\nChr01K, 10bp windows')

gg_tet1_r2 <- ggplot(tet_1_ld_wind, aes(x = window_pos, y = window_avg)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('r^2') +
  ylim(0,0.5) +
  ggtitle('r^2 in Upland 4X Group 1\nChr01K, 10bp windows')

gg_tet2_r2 <- ggplot(tet_2_ld_wind, aes(x = window_pos, y = window_avg)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('r^2') +
  ylim(0,0.5) +
  ggtitle('r^2 in Upland 4X Group 2\nChr01K, 10bp windows')

pdf(r2_out_fig_file, width = 12, height = 10)
(gg_oct1_r2 + gg_oct2_r2) / (gg_tet1_r2 + gg_tet2_r2)
dev.off()

###

gg_oct1_hi_ld <- ggplot(oct_1_hi_ld, aes(x = window_pos, y = per_high_ld)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('% r^2 > 0.9') +
  ylim(0, 0.3) +
  ggtitle('Percentage of r^2 values > 0.9\nin Upland 8X Group 1, Chr01K, 10bp windows')

gg_oct2_hi_ld <- ggplot(oct_2_hi_ld, aes(x = window_pos, y = per_high_ld)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('% r^2 > 0.9') +
  ylim(0, 0.3) +
  ggtitle('Percentage of r^2 values > 0.9\nin Upland 8X Group 2, Chr01K, 10bp windows')

gg_tet1_hi_ld <- ggplot(tet_1_hi_ld, aes(x = window_pos, y = per_high_ld)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('% r^2 > 0.9') +
  ylim(0, 0.3) +
  ggtitle('Percentage of r^2 values > 0.9\nin Upland 4X Group 1, Chr01K, 10bp windows')

gg_tet2_hi_ld <- ggplot(tet_2_hi_ld, aes(x = window_pos, y = per_high_ld)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('% r^2 > 0.9') +
  ylim(0, 0.3) +
  ggtitle('Percentage of r^2 values > 0.9\nin Upland 4X Group 2, Chr01K, 10bp windows')

pdf(hi_ld_out_fig_file, width = 12, height = 10)
(gg_oct1_hi_ld + gg_oct2_hi_ld) / (gg_tet1_hi_ld + gg_tet2_hi_ld)
dev.off()



