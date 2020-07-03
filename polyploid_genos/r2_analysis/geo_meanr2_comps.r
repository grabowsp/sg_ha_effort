# Code for generating r^2 figures for geo sampls

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)
library(patchwork)

general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)

### LOAD DATA ###

res_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results'

res_dir <- add_slash(res_dir)

combo_short <- 'combined.geo_subgroup.mean_r2.10bpwindow.rds_wrong'
combo_file <- paste(res_dir, combo_short, sep = '')

combo_r2 <- readRDS(combo_file)

chrom_short <- 'Chromosome.geo_subgroup.mean_r2.10bpwindow.rds_wrong'
chrom_file <- paste(res_dir, chrom_short)

# Goal: sort the r^2 values for each window and quantify the number
#  of windows where 8X values are above 4X (or quantify the rank of the 8X
#  groups

wind_r2_df <- data.frame(window_pos = combo_r2[[1]]$window_pos, 
  stringsAsFactors = F)

for(i in seq(length(combo_r2))){
  wind_r2_df[, names(combo_r2)[i]] <- combo_r2[[i]]$window_avg
}

r2_ord <- apply(wind_r2_df[, c(2:ncol(wind_r2_df))], 1, 
  function(x) order(x, decreasing = T))

r2_ord_mat <- matrix(r2_ord, ncol = 10, byrow = T)

up_8_1_ord <- apply(r2_ord_mat, 1, function(x) which(x == 3))

summary(up_8_1_ord)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  4.000   5.000   5.000   6.325   7.000  10.000

up_8_2_ord <- apply(r2_ord_mat, 1, function(x) which(x == 4))

summary(up_8_2_ord)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   2.00    3.00    3.00    3.11    3.00    5.00

r2_up_ord <- apply(wind_r2_df[, c(2:5)], 1,    
  function(x) order(x, decreasing = T))
r2_up_ord_mat <- matrix(r2_up_ord, ncol = 4, byrow = T)

up8x1_up_ord <- apply(r2_up_ord_mat, 1, function(x) which(x == 3))
summary(up8x1_up_ord)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  2.000   3.000   3.000   3.405   4.000   4.000

up8x2_up_ord <- apply(r2_up_ord_mat, 1, function(x) which(x == 4))
summary(up8x2_up_ord)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   1.081   1.000   2.000


