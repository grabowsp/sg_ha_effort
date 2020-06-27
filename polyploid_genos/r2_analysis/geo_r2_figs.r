# Code for generating r^2 figures for geo sampls

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)
library(patchwork)

general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)

### LOAD DATA ###













gg_oct1_r2 <- ggplot(oct_1_ld_wind, aes(x = window_pos, y = window_avg)) +
  geom_point() +
  xlab('distance between SNPs (bp), 10bp windows') +
  ylab('r^2') +
  ylim(0,0.5) +
  ggtitle('r^2 in Upland 8X Group 1\nChr01K, 10bp windows')




