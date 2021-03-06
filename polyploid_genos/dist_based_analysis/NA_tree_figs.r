# Script for generating trees from the pseudohaploid distances

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

# bash
# source activate r_phylo

#args = commandArgs(trailingOnly = TRUE)

# LOAD PACKAGES #
library(phangorn)

### LOAD DATA ###
#data_file <- args[1]
# File locations at HA
data_file_no_filt <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'polyploid_genos/dist_mats/', 
  'Chr01K.polyploid.CDS.allsamps.NA_DistMat.total.rds', sep = '')
data_file_high_miss <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'polyploid_genos/dist_mats/', 
  'Chr01K.polyploid.CDS.allsamps.high_miss.NA_DistMat.total.rds', 
   sep = '')
data_file_half_miss <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'polyploid_genos/dist_mats/',
  'Chr01K.polyploid.CDS.allsamps.half_miss.NA_DistMat.total.rds', 
   sep = '')
data_file_few_miss <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'polyploid_genos/dist_mats/',
  'Chr01K.polyploid.CDS.allsamps.few_miss.NA_DistMat.total.rds', 
   sep = '')

data_list <- list()
data_list[[1]] <- readRDS(data_file_no_filt)
data_list[[2]] <- readRDS(data_file_high_miss)
data_list[[3]] <- readRDS(data_file_half_miss)
data_list[[4]] <- readRDS(data_file_few_miss)

#ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
#  'ploidy_calling/sg_ploidy_results_v2.0.txt', sep = '')
ploidy_info_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/sg_ploidy_results_v2.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUTS ###
out_dir <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'polyploid_genos/dist_mats/', sep = '')

out_file_pre <- 'Chr01K.polyploid.'

out_file_suf_vec <- c('no_filt', 'high_miss', 'half_miss', 'few_miss')
# SET VARIABLES #

############
### Set Colors
#table(ploidy_info$SUBPOP_SNP)
#table(ploidy_info$ECOTYPE_SNP_CHLR)
#table(ploidy_info$total_ploidy)

subpop_col_vec <- rep(NA, times = length(table(ploidy_info$SUBPOP_SNP)))
names(subpop_col_vec) <- names(table(ploidy_info$SUBPOP_SNP))
subpop_col_vec['?Subpop'] <- 'black'
subpop_col_vec['Eastcoast'] <- 'orange2'
subpop_col_vec['Eastcoast_Admixed'] <- 'orangered3'
subpop_col_vec['Gulfcoast'] <- 'springgreen3'
subpop_col_vec['Gulfcoast_Admixed'] <- 'darkgreen'
subpop_col_vec['Midwest'] <- 'royalblue2'
subpop_col_vec['Midwest_Admixed'] <- 'blue2'
subpop_col_vec['Texas'] <- 'violetred2'
subpop_col_vec['Texas_Admixed'] <- 'darkorchid2'

#subpop_palette <- scale_colour_manual(name = 'Subpop', values = subpop_col_vec)

totPloid_col_vec <- rep(NA, times = length(table(ploidy_info$total_ploidy)))
names(totPloid_col_vec) <- names(table(ploidy_info$total_ploidy))
totPloid_col_vec['4X'] <- 'red2'
totPloid_col_vec['8X'] <- 'blue2'

#ploidy_palette <- scale_colour_manual(name = 'Ploidy', 
#  values = totPloid_col_vec)

#data_lab_vec <- c('No NA Filter', 'Max 80% NAs', 'Max 50% NAs', 'Max 20% NAs')

# Make Euclidean-distance PCoA data.frames

for(i in seq(length(data_list))){
  euc_dist_mat <- data_list[[i]][[2]]

  meta_ord <- c()
  for(j in seq(nrow(euc_dist_mat))){
    tmp_ind <- which(ploidy_info$lib == rownames(euc_dist_mat)[j])
    meta_ord <- c(meta_ord, tmp_ind)
  }

  subpop_plot_cols <- ploidy_info$SUBPOP_SNP[meta_ord]
  for(SP in names(subpop_col_vec)){
    tmp_inds <- which(subpop_plot_cols == SP)
    subpop_plot_cols[tmp_inds] <- subpop_col_vec[SP]
  }

  totPloid_plot_cols <- ploidy_info$total_ploidy[meta_ord]
  for(SP in names(totPloid_col_vec)){
    tmp_inds <- which(totPloid_plot_cols == SP)
    totPloid_plot_cols[tmp_inds] <- totPloid_col_vec[SP]
  }

  samp_labs <- paste(ploidy_info$PLANT_ID, ploidy_info$lib, sep = '_')[meta_ord]
  rownames(euc_dist_mat) <- colnames(euc_dist_mat) <- samp_labs

  euc_nj <- NJ(euc_dist_mat)

  tmp_out_file_sp <- paste(out_dir, out_file_pre, out_file_suf_vec[i], 
    '.subpopColors.NJtree.pdf', sep = '')

  pdf(file = tmp_out_file_sp, width = 10, height = 10)
  plot(euc_nj, type = 'unrooted', use.edge.length = T, no.margin = T,
    cex = 0.3, tip.color = subpop_plot_cols)
  dev.off()

  tmp_out_file_pl <- paste(out_dir, out_file_pre, out_file_suf_vec[i],
    '.ploidyColors.NJtree.pdf', sep = '')
  
  pdf(file = tmp_out_file_pl, width = 10, height = 10)
  plot(euc_nj, type = 'unrooted', use.edge.length = T, no.margin = T,
    cex = 0.3, tip.color = totPloid_plot_cols)
  dev.off()

}

quit(save = 'no')

