# Script for generating trees from the pseudohaploid distances

# source activate r_phylo

# LOAD PACKAGES #
library(phangorn)

### LOAD DATA ###
data_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/',
  'filtered_vcfs/', 'pseudohapdists_v01_total.rds', sep = '')
data <- readRDS(data_file)

ploidy_info_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/sg_ploidy_results_v2.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUTS ###
# output tree objects
man_nj_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_manhattan_nj.tre', 
  sep = '')
man_upgma_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_manhattan_upgma.tre', 
  sep = '')

euc_nj_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_euclidean_nj.tre',
  sep = '')
euc_upgma_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_euclidean_upgma.tre',
  sep = '')

# output tree figures
man_nj_fig_pre <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_manhattan_nj.',
  sep = '')
man_upgma_fig_pre <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_manhattan_upgma.',
  sep = '')

euc_nj_fig_pre <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_euclidean_nj.',
  sep = '')
euc_upgma_fig_pre <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_euclidean_upgma.',
  sep = '')


# SET VARIABLES #

############

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

totPloid_col_vec <- rep(NA, times = length(table(ploidy_info$total_ploidy)))
names(totPloid_col_vec) <- names(table(ploidy_info$total_ploidy))
totPloid_col_vec['4X'] <- 'red2'
totPloid_col_vec['8X'] <- 'blue2'

man_dist <- as.dist(data[[2]][-902, -902])

man_dist_label_names <- gsub('^X', '', labels(man_dist))
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

man_upgma <- upgma(man_dist)

man_nj <- NJ(man_dist)

write.tree(man_upgma, file = man_upgma_tree_file)
write.tree(man_nj, file = man_nj_tree_file)

man_nj_subpop_file <- paste(man_nj_fig_pre, 'subpop.pdf', sep = '')
pdf(file = man_nj_subpop_file, width = 10, height = 10)
plot(man_nj, type = 'unrooted', use.edge.length = T, no.margin = T, 
  cex = 0.3, tip.color = subpop_plot_cols)
dev.off()

man_nj_totploid_file <- paste(man_nj_fig_pre, 'TotalPloidy.pdf', sep = '')
pdf(file = man_nj_totploid_file, width = 10, height = 10)
plot(man_nj, type = 'unrooted', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = totPloid_plot_cols)
dev.off()

###

euc_dist <- as.dist(data[[3]][-902,-902])

euc_upgma <- upgma(euc_dist)
euc_nj <- NJ(euc_dist)

write.tree(euc_upgma, file = euc_upgma_tree_file)
write.tree(euc_nj, file = euc_nj_tree_file)

euc_dist_label_names <- gsub('^X', '', labels(euc_dist))
euc_dist_label_names[which(euc_dist_label_names == '9001.3.BN389.69S')] <- '9001-3 BN389-69S'
euc_dist_label_names[which(euc_dist_label_names == '9020.5')] <- '9020-5'

meta_ord <- c()
for(DL in euc_dist_label_names){
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

euc_nj_subpop_file <- paste(euc_nj_fig_pre, 'subpop.pdf', sep = '')
pdf(file = euc_nj_subpop_file, width = 10, height = 10)
plot(euc_nj, type = 'unrooted', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = subpop_plot_cols)
dev.off()

euc_nj_totploid_file <- paste(euc_nj_fig_pre, 'TotalPloidy.pdf', sep = '')
pdf(file = euc_nj_totploid_file, width = 10, height = 10)
plot(euc_nj, type = 'unrooted', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = totPloid_plot_cols)
dev.off()




##############



quit(save = 'no')

