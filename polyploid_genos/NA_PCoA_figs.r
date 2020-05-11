# Script for generating trees from the pseudohaploid distances

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

#args = commandArgs(trailingOnly = TRUE)

# LOAD PACKAGES #
# library(phangorn)
library(ggplot2)
library(patchwork)

### LOAD DATA ###
#data_file <- args[1]
data_file_no_filt <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'na_popstructure/all_samps/', 
  'no_filtering/Chr01K.polyploid.CDS.allsamps.NA_DistMat.total.rds', sep = '')
data_file_high_miss <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'na_popstructure/all_samps/', 
  'high_missing/Chr01K.polyploid.CDS.allsamps.high_miss.NA_DistMat.total.rds', 
   sep = '')
data_file_half_miss <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'na_popstructure/all_samps/',
  'half_missing/Chr01K.polyploid.CDS.allsamps.half_miss.NA_DistMat.total.rds', 
   sep = '')
data_file_few_miss <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'na_popstructure/all_samps/',
  'few_missing/Chr01K.polyploid.CDS.allsamps.few_miss.NA_DistMat.total.rds', 
   sep = '')

data_list <- list()
data_list[[1]] <- readRDS(data_file_no_filt)
data_list[[2]] <- readRDS(data_file_high_miss)
data_list[[3]] <- readRDS(data_file_half_miss)
data_list[[4]] <- readRDS(data_file_few_miss)

ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'ploidy_calling/sg_ploidy_results_v2.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUTS ###
pco1v2_fig_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'na_popstructure/all_samps/',
  'Chr01K.CDS.NA_PCo_1v2.pdf', sep = '')

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

subpop_palette <- scale_colour_manual(name = 'Subpop', values = subpop_col_vec)

totPloid_col_vec <- rep(NA, times = length(table(ploidy_info$total_ploidy)))
names(totPloid_col_vec) <- names(table(ploidy_info$total_ploidy))
totPloid_col_vec['4X'] <- 'red2'
totPloid_col_vec['8X'] <- 'blue2'

ploidy_palette <- scale_colour_manual(name = 'Ploidy', 
  values = totPloid_col_vec)

data_lab_vec <- c('No NA Filter', 'Max 80% NAs', 'Max 50% NAs', 'Max 20% NAs')

# Make Euclidean-distance PCoA data.frames

cmd_list <- list()
cmd_var_list <- list()
for(i in seq(length(data_list))){
  euc_dist_mat <- data_list[[i]][[2]]
  euc_cmd <- cmdscale(euc_dist_mat, k = 200)
  euc_tot_var <- sum(apply(euc_cmd, 2, var))
  euc_per_var <- (apply(euc_cmd, 2, var)/euc_tot_var)*100

  euc_df <- data.frame(euc_cmd, stringsAsFactors = F)
  colnames(euc_df) <- paste('PCo_', seq(ncol(euc_df)), sep = '')

  meta_ord <- c()
  for(j in seq(nrow(euc_dist_mat))){
    tmp_ind <- which(ploidy_info$lib == rownames(euc_dist_mat)[j])
    meta_ord <- c(meta_ord, tmp_ind)
  }

  euc_df$samp <- paste(ploidy_info$lib, 
    ploidy_info$PLANT_ID, sep = '_')[meta_ord]
  euc_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
  euc_df$totPloid <- ploidy_info$total_ploidy[meta_ord]

  cmd_list[[i]] <- euc_df
  cmd_var_list[[i]] <- euc_per_var
}

pcoX <- 1
pcoY <- 2

gg_nf_s <- ggplot(cmd_list[[1]], aes(x = cmd_list[[1]][, 1], 
    y = cmd_list[[1]][, 2])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[1]][pcoX], 2), '%)', 
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[1]][pcoY], 2), '%)', 
    sep = '')) +
  ggtitle(paste(data_lab_vec[1], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_nf_p <- ggplot(cmd_list[[1]], aes(x = cmd_list[[1]][, 1],
    y = cmd_list[[1]][, 2])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[1]][pcoX], 2), '%)', 
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[1]][pcoY], 2), '%)', 
    sep = '')) +
  ggtitle(paste(data_lab_vec[1], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_hi_s <- ggplot(cmd_list[[2]], aes(x = cmd_list[[2]][, 1],
    y = cmd_list[[2]][, 2])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[2]][pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[2]][pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(data_lab_vec[2], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_hi_p <- ggplot(cmd_list[[2]], aes(x = cmd_list[[2]][, 1],
    y = cmd_list[[2]][, 2])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[2]][pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[2]][pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(data_lab_vec[2], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_half_s <- ggplot(cmd_list[[3]], aes(x = cmd_list[[3]][, 1],
    y = cmd_list[[3]][, 2])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[3]][pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[3]][pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(data_lab_vec[3], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_half_p <- ggplot(cmd_list[[3]], aes(x = cmd_list[[3]][, 1],
    y = cmd_list[[3]][, 2])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[3]][pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[3]][pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(data_lab_vec[3], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_few_s <- ggplot(cmd_list[[4]], aes(x = cmd_list[[4]][, 1],
    y = cmd_list[[4]][, 2])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[4]][pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[4]][pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(data_lab_vec[4], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_few_p <- ggplot(cmd_list[[4]], aes(x = cmd_list[[4]][, 1],
    y = cmd_list[[4]][, 2])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(cmd_var_list[[4]][pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(cmd_var_list[[4]][pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(data_lab_vec[4], ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(pco1v2_fig_file, width = 12, height = 20)
(gg_nf_s + gg_nf_p) / (gg_hi_s + gg_hi_p) / (gg_half_s + gg_half_p) / 
  (gg_few_s + gg_few_p)
dev.off()

quit(save = 'no')

