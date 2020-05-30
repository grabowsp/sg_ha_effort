# Script for generating PCoA figs from either manhattan or euclidean distances

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

# LOAD PACKAGES #
# library(phangorn)
library(ggplot2)
library(patchwork)

### LOAD DATA ###
data_file <- args[1]
#data_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
#  'polyploid_genos_popstructure/polyploid_dists/', 
#  'Chr01K.polyploid.CDS.allsamps.few_miss_00_ploidy_DistMat.rds', sep = '')

data <- readRDS(data_file)

ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'ploidy_calling/sg_ploidy_results_v2.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

dist_type_in <- args[2]
#dist_type_in <- 'manhattan'

dt_first_letter <- unlist(strsplit(dist_type_in, split =''))[1]

if(length(intersect(dt_first_letter, c('m', 'M'))) > 0){
  dist_type <- 'manhattan_dist'
  out_file_dist_suf <- 'manDist'}
if(length(intersect(dt_first_letter, c('e', 'E'))) > 0){
  dist_type <- 'euclidean_dist'
  out_file_dist_suf <- 'eucDist'}
if(length(intersect(dt_first_letter, c('m', 'M', 'e', 'E'))) == 0){
  print('Must indicate if using manhattan or euclidean distance')
  quit()
}

### SET OUTPUTS ###
out_file_pre <- gsub('.rds', '', data_file)

out_file <- paste(out_file_pre, out_file_dist_suf,
  'general_PCoAs.pdf', sep = '.')

# SET VARIABLES #

plot_title_pre <- args[3]
#plot_title_pre <- 'Polyploid_Genotypes'

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

# Make PCoA data.frames

dist_mat <- as.matrix(data[[dist_type]])

dist_cmd <- cmdscale(dist_mat, k = 200)
dist_tot_var <- sum(apply(dist_cmd, 2, var))
dist_per_var <- (apply(dist_cmd, 2, var)/dist_tot_var)*100

dist_df <- data.frame(dist_cmd, stringsAsFactors = F)
colnames(dist_df) <- paste('PCo_', seq(ncol(dist_df)), sep = '')

meta_ord <- c()
for(j in seq(nrow(dist_mat))){
  tmp_ind <- which(ploidy_info$lib == rownames(dist_mat)[j])
  meta_ord <- c(meta_ord, tmp_ind)
}

dist_df$samp <- paste(ploidy_info$lib, 
  ploidy_info$PLANT_ID, sep = '_')[meta_ord]
dist_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
dist_df$totPloid <- ploidy_info$total_ploidy[meta_ord]

pcoX <- 1
pcoY <- 2
# using the variable to call columns for x and y doesn't work; so hard-coding
## it for now
# seems to be because the last pcoX and pcoY variables are called during
## plotting
# Can try pasting the column names and then using column name characters

gg_1_2_s <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 2])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)', 
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)', 
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_1_2_p <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 2])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

####

pcoX <- 1
pcoY <- 3

gg_1_3_s <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 3])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_1_3_p <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 3])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

####

pcoX <- 1
pcoY <- 4

gg_1_4_s <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 4])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_1_4_p <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 4])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

######

pcoX <- 1
pcoY <- 5

gg_1_5_s <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 5])) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

gg_1_5_p <- ggplot(dist_df, aes(x = dist_df[, 1], y = dist_df[, 5])) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab(paste('PCo_', pcoX, ' (', round(dist_per_var[pcoX], 2), '%)',
    sep = '')) +
  ylab(paste('PCo_', pcoY, ' (', round(dist_per_var[pcoY], 2), '%)',
    sep = '')) +
  ggtitle(paste(plot_title_pre, ' PCo_', pcoX, ' vs PCo_', pcoY, sep = ''))

#####

pdf(out_file, width = 12, height = 20)
(gg_1_2_s + gg_1_2_p) / (gg_1_3_s + gg_1_3_p) / (gg_1_4_s + gg_1_4_p) / 
  (gg_1_5_s + gg_1_5_p) 
dev.off()

quit(save = 'no')

