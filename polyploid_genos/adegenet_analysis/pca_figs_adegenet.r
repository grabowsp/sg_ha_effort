# Script for generating PCA figures from adegenet results 

#module load python/3.7-anaconda-2019.07
#source activate r_adegenet_env

args = commandArgs(trailingOnly = TRUE)

### LOAD PACKAGES ###

library(adegenet)
library(parallel)
library(ggplot2)
library(patchwork)

### LOAD DATA ###

pca_res_file <- args[1]
#pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds'

pca_res <- readRDS(pca_res_file)

ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'ploidy_calling/sg_ploidy_results_v3.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUT ###
out_file <- gsub('.rds', '.figs.pdf', pca_res_file, fixed = T)

### SET VARIABLES ###
plot_title_pre <- args[2]
#plot_title_pre <- 'Geographic samples, 595K SNPs\n'

##########
pca_eig <- pca_res$eig
per_var_vec <- pca_eig / sum(pca_eig) * 100

pca_mat <- pca_res$scores
pca_df <- data.frame(lib = rownames(pca_mat), pca_mat, stringsAsFactors = F)

meta_ord <- c()
for(j in seq(nrow(pca_df))){
  tmp_ind <- which(ploidy_info$lib == pca_df$lib[j])
  meta_ord <- c(meta_ord, tmp_ind)
}

pca_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
pca_df$totPloid <- ploidy_info$total_ploidy_2[meta_ord]

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

# Need to adjust this to include 6X
totPloid_col_vec <- rep(NA, times = length(table(ploidy_info$total_ploidy_2)))
names(totPloid_col_vec) <- names(table(ploidy_info$total_ploidy_2))
totPloid_col_vec['4X'] <- 'red2'
totPloid_col_vec['8X'] <- 'blue2'
totPloid_col_vec['6X'] <- 'gray40'

ploidy_palette <- scale_colour_manual(name = 'Ploidy',
  values = totPloid_col_vec)

# Need to add % variance in eigenvectors to axis labels
pcX <- 1
pcY <- 2

gg_1_2_s <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = SUBPOP)) + subpop_palette + 
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) + 
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) + 
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

gg_1_2_p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = totPloid)) + ploidy_palette + 
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) + 
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) + 
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

pcX <- 1
pcY <- 3

gg_1_3_s <- ggplot(pca_df, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = SUBPOP)) + subpop_palette + 
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) + 
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) + 
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

gg_1_3_p <- ggplot(pca_df, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) +
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) +
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

pcX <- 1
pcY <- 4

gg_1_4_s <- ggplot(pca_df, aes(x = PC1, y = PC4)) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) +
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) +
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

gg_1_4_p <- ggplot(pca_df, aes(x = PC1, y = PC4)) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) +
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) +
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

pcX <- 1
pcY <- 5

gg_1_5_s <- ggplot(pca_df, aes(x = PC1, y = PC5)) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) +
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) +
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

gg_1_5_p <- ggplot(pca_df, aes(x = PC1, y = PC5)) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) +
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) +
  ggtitle(paste(plot_title_pre, ' PC', pcX, ' vs PC', pcY, sep = ''))

pdf(out_file, width = 12, height = 20)
(gg_1_2_s + gg_1_2_p) / (gg_1_3_s + gg_1_3_p) / (gg_1_4_s + gg_1_4_p) / 
  (gg_1_5_s + gg_1_5_p)
dev.off()

quit(save = 'no')

