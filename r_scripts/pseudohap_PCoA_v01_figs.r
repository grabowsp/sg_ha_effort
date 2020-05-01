# Script for generating trees from the pseudohaploid distances

# source activate r_phylo

# LOAD PACKAGES #
# library(phangorn)
library(ggplot2)
library(patchwork)
### LOAD DATA ###
data_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/',
  'filtered_vcfs/', 'pseudohapdists_v01_total.rds', sep = '')
data <- readRDS(data_file)

ploidy_info_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/sg_ploidy_results_v2.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUTS ###
# manhattan-distance plots
man_PCo1v2_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_manhattan_PCo1vPCo2.pdf', sep = '')
man_PCo2v3_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_manhattan_PCo2vPCo3.pdf', sep = '')
man_PCo1v4_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_manhattan_PCo1vPCo4.pdf', sep = '')
man_PCo4v5_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_manhattan_PCo4vPCo5.pdf', sep = '')

man_PCo_tot_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_manhattan_PCo_figs.pdf', sep = '')


euc_PCo1v2_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_euclidean_PCo1vPCo2.pdf', sep = '')
euc_PCo2v3_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_euclidean_PCo2vPCo3.pdf', sep = '')
euc_PCo1v4_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_euclidean_PCo1vPCo4.pdf', sep = '')
euc_PCo4v5_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_euclidean_PCo4vPCo5.pdf', sep = '')

euc_PCo_tot_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/filtered_vcfs/',
  'pseudohapdists_v01_euclidean_PCo_figs.pdf', sep = '')



# SET VARIABLES #

############

# Set Colors
table(ploidy_info$SUBPOP_SNP)
table(ploidy_info$ECOTYPE_SNP_CHLR)
table(ploidy_info$total_ploidy)

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

ploidy_palette <- scale_colour_manual(name = 'Ploidy', values = totPloid_col_vec)

# Make Manhattan-distance PCoA plots 

man_dist_mat <- data[[2]][-902,-902]

man_cmd <- cmdscale(man_dist_mat, k = 200)
man_tot_var <- sum(apply(man_cmd, 2, var))
man_per_var <- (apply(man_cmd, 2, var)/man_tot_var)*100

man_dist_label_names <- gsub('^X', '', colnames(man_dist_mat))
man_dist_label_names[which(man_dist_label_names == '9001.3.BN389.69S')] <- '9001-3 BN389-69S'
man_dist_label_names[which(man_dist_label_names == '9020.5')] <- '9020-5'

meta_ord <- c()
for(DL in man_dist_label_names){
  meta_ord <- c(meta_ord, which(ploidy_info$PLANT_ID == DL))
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

man_df <- data.frame(man_cmd, stringsAsFactors = F)
colnames(man_df) <- paste('PCo_', seq(ncol(man_df)), sep = '')

man_df$samp <- man_dist_label_names
man_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
man_df$totPloid <- ploidy_info$total_ploidy[meta_ord]

pcoX <- 1
pcoY <- 2

gg_man_1v2_s <- ggplot(man_df, aes(x = man_df[, 1], y = man_df[, 2])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_man_1v2_p <- ggplot(man_df, aes(x = man_df[, 1], y = man_df[, 2])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(man_PCo1v2_fig, width = 12, height = 5)
gg_man_1v2_s + gg_man_1v2_p
dev.off()

pcoX <- 2
pcoY <- 3

gg_man_PCo3_s <- ggplot(man_df, aes(x = man_df[, 2], y = man_df[, 3])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_man_PCo3_p <- ggplot(man_df, aes(x = man_df[, 2], y = man_df[, 3])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(man_PCo2v3_fig, width = 12, height = 5)
gg_man_PCo3_s + gg_man_PCo3_p
dev.off()

pcoX <- 1
pcoY <- 4

gg_man_PCo4_s <- ggplot(man_df, aes(x = man_df[, 1], y = man_df[, 4])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_man_PCo4_p <- ggplot(man_df, aes(x = man_df[, 1], y = man_df[, 4])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(man_PCo1v4_fig, width = 12, height = 5)
gg_man_PCo4_s + gg_man_PCo4_p
dev.off()

pcoX <- 4
pcoY <- 5

gg_man_PCo5_s <- ggplot(man_df, aes(x = man_df[, 4], y = man_df[, 5])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_man_PCo5_p <- ggplot(man_df, aes(x = man_df[, 4], y = man_df[, 5])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(man_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(man_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(man_PCo4v5_fig, width = 12, height = 5)
gg_man_PCo5_s + gg_man_PCo5_p
dev.off()

pdf(man_PCo_tot_fig, width = 12, height = 20)
(gg_man_1v2_s + gg_man_1v2_p) / (gg_man_PCo3_s + gg_man_PCo3_p) / 
  (gg_man_PCo4_s + gg_man_PCo4_p) / (gg_man_PCo5_s + gg_man_PCo5_p)
dev.off()


###################
# Euclidean distances

euc_dist_mat <- data[[3]][-902,-902]

euc_cmd <- cmdscale(euc_dist_mat, k = 200)
euc_tot_var <- sum(apply(euc_cmd, 2, var))
euc_per_var <- (apply(euc_cmd, 2, var)/euc_tot_var)*100

euc_dist_label_names <- gsub('^X', '', colnames(euc_dist_mat))
euc_dist_label_names[which(euc_dist_label_names == '9001.3.BN389.69S')] <- '9001-3 BN389-69S'
euc_dist_label_names[which(euc_dist_label_names == '9020.5')] <- '9020-5'

euc_df <- data.frame(euc_cmd, stringsAsFactors = F)
colnames(euc_df) <- paste('PCo_', seq(ncol(euc_df)), sep = '')

euc_df$samp <- euc_dist_label_names
euc_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
euc_df$totPloid <- ploidy_info$total_ploidy[meta_ord]

pcoX <- 1
pcoY <- 2

gg_euc_1v2_s <- ggplot(euc_df, aes(x = euc_df[, 1], y = euc_df[, 2])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_euc_1v2_p <- ggplot(euc_df, aes(x = euc_df[, 1], y = euc_df[, 2])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(euc_PCo1v2_fig, width = 12, height = 5)
(gg_euc_1v2_s + gg_euc_1v2_p)
dev.off()

pcoX <- 2
pcoY <- 3

gg_euc_PCo3_s <- ggplot(euc_df, aes(x = euc_df[, 2], y = euc_df[, 3])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_euc_PCo3_p <- ggplot(euc_df, aes(x = euc_df[, 2], y = euc_df[, 3])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(euc_PCo2v3_fig, width = 12, height = 5)
gg_euc_PCo3_s + gg_euc_PCo3_p
dev.off()

pcoX <- 1
pcoY <- 4

gg_euc_PCo4_s <- ggplot(euc_df, aes(x = euc_df[, 1], y = euc_df[, 4])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_euc_PCo4_p <- ggplot(euc_df, aes(x = euc_df[, 1], y = euc_df[, 4])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(euc_PCo1v4_fig, width = 12, height = 5)
gg_euc_PCo4_s + gg_euc_PCo4_p
dev.off()

pcoX <- 4
pcoY <- 5

gg_euc_PCo5_s <- ggplot(euc_df, aes(x = euc_df[, 4], y = euc_df[, 5])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_euc_PCo5_p <- ggplot(euc_df, aes(x = euc_df[, 4], y = euc_df[, 5])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(euc_per_var[pcoX], 2), '%)', sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(euc_per_var[pcoY], 2), '%)', sep = '')) +
  ggtitle(paste('PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

pdf(euc_PCo4v5_fig, width = 12, height = 5)
gg_euc_PCo5_s + gg_euc_PCo5_p
dev.off()

quit(save = 'no')

