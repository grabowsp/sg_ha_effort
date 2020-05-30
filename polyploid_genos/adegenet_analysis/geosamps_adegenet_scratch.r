# Code for quick geosamps analysis using adegenet

bash
source activate r_adegenet_env

library(adegenet)
library(parallel)

data_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.595K.polyploid.CDS.geosamps.genlight.rds'

tot_gl <- readRDS(data_in)

tot_pca <- glPca(tot_gl, nf = 100, loadings = F, alleleAsUnit = F,
  useC = F)

sub_inds <- sort(sample(seq(nLoc(tot_gl)), size = 5000))

sub_gl <- tot_gl[ , sub_inds]

sub_pca <- glPca(sub_gl, nf = 100, loadings = F, alleleAsUnit = F, 
  useC = F)

sub_res_out <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/sub5k.polyploid.CDS.geosamps.PCAresults.rds'

saveRDS(sub_pca, sub_res_out)

#######

module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)
library(parallel)
library(ggplot2)
#library(patchwork)

pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sub5k.polyploid.CDS.geosamps.PCAresults.rds'

pca_res <- readRDS(pca_res_file)

pca_eig <- pca_res$eig

pca_mat <- pca_res$scores

pca_df <- data.frame(lib = rownames(pca_mat), pca_mat, stringsAsFactors = F)


ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'ploidy_calling/sg_ploidy_results_v3.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

meta_ord <- c()
for(j in seq(nrow(pca_df))){
  tmp_ind <- which(ploidy_info$lib == pca_df$lib[j])
  meta_ord <- c(meta_ord, tmp_ind)
}

pca_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
pca_df$totPloid <- ploidy_info$total_ploidy[meta_ord]

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
totPloid_col_vec <- rep(NA, times = length(table(ploidy_info$total_ploidy)))
names(totPloid_col_vec) <- names(table(ploidy_info$total_ploidy))
totPloid_col_vec['4X'] <- 'red2'
totPloid_col_vec['8X'] <- 'blue2'

ploidy_palette <- scale_colour_manual(name = 'Ploidy',
  values = totPloid_col_vec)

# Need to add % variance in eigenvectors to axis labels
pcX <- 'PC1'
pcY <- 'PC2'

gg_1_2_s <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = SUBPOP)) + subpop_palette + 
  xlab('PC1') + ylab('PC2') + 
  ggtitle('PC1 vs PC2 for 5k subsampled SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sub5k.polyploid.CDS.geosamps.PC1vPC2_subpop.pdf'

pdf(out_file, width = 6, height = 5)
gg_1_2_s
dev.off()

gg_1_2_p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab('PC1') + ylab('PC2') +
  ggtitle('PC1 vs PC2 for 5k subsampled SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sub5k.polyploid.CDS.geosamps.PC1vPC2_ploidy.pdf'

pdf(out_file, width = 6, height = 5)
gg_1_2_p
dev.off()

gg_1_3_s <- ggplot(pca_df, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab('PC1') + ylab('PC3') +
  ggtitle('PC1 vs PC3 for 5k subsampled SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sub5k.polyploid.CDS.geosamps.PC1vPC3_subpop.pdf'

pdf(out_file, width = 6, height = 5)
gg_1_3_s
dev.off()

gg_1_3_p <- ggplot(pca_df, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab('PC1') + ylab('PC3') +
  ggtitle('PC1 vs PC3 for 5k subsampled SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sub5k.polyploid.CDS.geosamps.PC1vPC3_ploidy.pdf'

pdf(out_file, width = 6, height = 5)
gg_1_3_p
dev.off()

gg_1_4_s <- ggplot(pca_df, aes(x = PC1, y = PC4)) +
  geom_point(aes(color = SUBPOP)) + subpop_palette +
  xlab('PC1') + ylab('PC4') +
  ggtitle('PC1 vs PC4 for 5k subsampled SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sub5k.polyploid.CDS.geosamps.PC1vPC4_subpop.pdf'

pdf(out_file, width = 6, height = 5)
gg_1_4_s
dev.off()

gg_1_4_p <- ggplot(pca_df, aes(x = PC1, y = PC4)) +
  geom_point(aes(color = totPloid)) + ploidy_palette +
  xlab('PC1') + ylab('PC4') +
  ggtitle('PC1 vs PC4 for 5k subsampled SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sub5k.polyploid.CDS.geosamps.PC1vPC4_ploidy.pdf'

pdf(out_file, width = 6, height = 5)
gg_1_4_p
dev.off()











