# Script used for generating plots for MNP counts for difference depths
## Adjust [seq_depth] and [genome_region] to generate differt sets of plots

### LOAD PACKAGES ###
# Need to load the R_analysis conda environment
# module load python/3.7-anaconda-2019.07
# source activate R_analysis

#args <- commandArgs(trainlingOnly=T)

library(ggplot2)

### SET INPUTS ###
seq_depth <- 3
#seq_depth <- as.numeric(args[1])
seq_depth_char <- as.character(seq_depth)

#genome_region <- 'genic'
#genome_region <- as.character(args[2])
# 'CDS', 'genic'

# Counts of CDS and genic MNPs
cds_mnp_rds <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_0', seq_depth_char, 
  '/CDS_MNP_count_list_0', seq_depth_char, 'depth.rds', 
  sep = '')
cds_mnp_list <- readRDS(cds_mnp_rds)

genic_mnp_rds <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_0', seq_depth_char,
  '/genic_MNP_count_list_0', seq_depth_char, 'depth.rds',
  sep = '')
genic_mnp_list <- readRDS(genic_mnp_rds)

# Get chromosome names
chrom_name_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/', 
  'chr_names.txt', sep = '')
chrom_names <- as.vector(
  read.table(chrom_name_file, header = F, stringsAsFactors = F)[,1])

# Seq depth mat
depth_mat_rds <- '/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth/sample_seq_depth_mat.rds'
depth_mat <- readRDS(depth_mat_rds)

# BAD samples from metadata
meta_lib_remove_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/sg_libs_to_remove_from_meta_Dec2019.txt'
meta_remove <- as.vector(read.table(meta_lib_remove_file, header = T,
  stringsAsFactors = F)[,1])

# reseq metadata
meta_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/Reseq_Metadata_Sept_2019_Edited_for_R.tsv'
meta <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F,
  quote = "", comment.char = '$')

# reseq cultivar info
cultivar_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/switchgrass_reseq_cultivars_corrected.tsv'
cult_info <- read.table(cultivar_file, sep = '\t', header = T,
  stringsAsFactors = F)

### SET OUTPUTS ###
fig_out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_0', seq_depth_char, '/', sep = '')
plot_height <- 5
plot_width <- 5.5

# results table
mnp_res_table_file <- paste(fig_out_dir, 'pos_0', seq_depth_char, 
  '_MNP_results_with_meta_info.txt', sep = '')

### SET VARIABLES ###

### SET CONSTANTS ###

##########

#sum(names(mnp_list) %in% meta$LIBRARY)
# [1] 921

miss_inds <- which(names(cds_mnp_list) %in% meta$LIBRARY == F)
miss_meta <- which(meta$LIBRARY %in% names(cds_mnp_list) == F)

mnp_overlap <- which(names(cds_mnp_list) %in% meta$LIBRARY)

# CONTINUE FROM HERE - NEED TO INCORPORATE BOTH CDS AND GENIC COUNTS

cds_lib_mnp_mat <- matrix(data = unlist(cds_mnp_list), byrow = T,
  nrow = length(cds_mnp_list))
rownames(cds_lib_mnp_mat) <- names(cds_mnp_list)
colnames(cds_lib_mnp_mat) <- chrom_names

cds_tot_mnp <- apply(cds_lib_mnp_mat, 1, sum)

genic_lib_mnp_mat <- matrix(data = unlist(genic_mnp_list), byrow = T,
  nrow = length(genic_mnp_list))
rownames(genic_lib_mnp_mat) <- names(genic_mnp_list)
colnames(genic_lib_mnp_mat) <- chrom_names

genic_tot_mnp <- apply(genic_lib_mnp_mat, 1, sum)

# Get the sequencing output for each sample
tot_coverage <- apply(depth_mat, 1, sum)

cds_mnp_corrected <- (cds_tot_mnp / tot_coverage) * min(tot_coverage)
genic_mnp_corrected <- (genic_tot_mnp / tot_coverage) * min(tot_coverage)

lib_info_0 <- data.frame(lib = names(cds_tot_mnp), cds_mnp_raw = cds_tot_mnp,
  cds_mnp_stand = cds_mnp_corrected, genic_mnp_raw = genic_tot_mnp,
  genic_mnp_stand = genic_mnp_corrected, seq_cov = tot_coverage,
  stringsAsFactors = F)

lib_info_df <- lib_info_0[mnp_overlap,]

lib_info_df$PLANT_ID <- NA
lib_info_df$PLOIDY <- NA
lib_info_df$PLD_OLD_CALL <- NA
lib_info_df$PLD_FLOW <- NA
lib_info_df$PLD_SNP <- NA
lib_info_df$ECOTYPE_SNP_CHLR <- NA
lib_info_df$SUBPOP_SNP <- NA

for(i in seq(nrow(lib_info_df))){
  tmp_meta_ind <- which(meta$LIBRARY == lib_info_df$lib[i])
  #print(i)
  #print(tmp_meta_ind)
  lib_info_df$PLANT_ID[i] <- meta$PLANT_ID[tmp_meta_ind]
  lib_info_df$PLOIDY[i] <- meta$PLOIDY[tmp_meta_ind]
  lib_info_df$PLD_OLD_CALL[i] <- meta$PLD_OLD_CALL[tmp_meta_ind]
  lib_info_df$PLD_FLOW[i] <- meta$PLD_FLOW[tmp_meta_ind]
  lib_info_df$PLD_SNP[i] <- meta$PLD_SNP[tmp_meta_ind]
  lib_info_df$ECOTYPE_SNP_CHLR[i] <- meta$ECOTYPE_SNP_CHLR[tmp_meta_ind]
  lib_info_df$SUBPOP_SNP[i] <- meta$SUBPOP_SNP[tmp_meta_ind]
}

lib_info_df$PLD_OLD_CALL[grep(4, lib_info_df$PLD_OLD_CALL)] <- '4X'
lib_info_df$PLD_OLD_CALL[grep(8, lib_info_df$PLD_OLD_CALL)] <- '8X'
lib_info_df$PLD_OLD_CALL[is.na(lib_info_df$PLD_OLD_CALL)] <- '?X'

lib_info_df$PLD_FLOW[is.na(lib_info_df$PLD_FLOW)] <- '?X'

lib_info_df$PLD_SNP[is.na(lib_info_df$PLD_SNP)] <- '?X'

lib_info_df$ECOTYPE_SNP_CHLR[is.na(lib_info_df$ECOTYPE_SNP_CHLR)] <- '?Ecotype'

lib_info_df$SUBPOP_SNP[is.na(lib_info_df$SUBPOP_SNP)] <- '?Subpop'

rm_inds <- c()
for(mr in meta_remove){
  tmp_ind <- which(lib_info_df$lib == mr)
  rm_inds <- c(rm_inds, tmp_ind)
}

filt_lib_df <- lib_info_df[-rm_inds, ]

cult_info_2 <- cult_info[which(cult_info$PLANT_ID %in% filt_lib_df$PLANT_ID),]

filt_lib_df$Cultivar_name <- '?'
filt_lib_df$Cultivar_ecotype <- '?'
filt_lib_df$Cultivar_ploidy <- '?'

for(ci in seq(nrow(cult_info_2))){
  tmp_meta_ind <- which(filt_lib_df$PLANT_ID == cult_info_2$PLANT_ID[ci])
  filt_lib_df$Cultivar_name[tmp_meta_ind] <- cult_info_2$Cultivar[ci]
  filt_lib_df$Cultivar_ecotype[tmp_meta_ind] <- cult_info_2$Ecotype[ci]
  filt_lib_df$Cultivar_ploidy[tmp_meta_ind] <- cult_info_2$Typical.Ploidy[ci]
}

# Save results table
write.table(filt_lib_df, file = mnp_res_table_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

# MNP vs depth color coded by different ploidy data

tmp_plot_param <- c('CDS', 'Flow')

cds_flow_dot <- ggplot(data = filt_lib_df) + 
  geom_point(aes(x = seq_cov, y = cds_mnp_stand, color = PLD_FLOW)) + 
  ggtitle(paste('Standardized ', tmp_plot_param[1], ' MNPs with ',
    seq_depth_char, '+ Depth vs Seq Depth;\nColors According to ', 
    tmp_plot_param[2], sep = ''))

tmp_plot_name <- paste(fig_out_dir, tmp_plot_param[1], '_MNP_0', seq_depth_char,
  '_by_', tmp_plot_param[2], '.pdf', sep = '')

pdf(tmp_plot_name, height = plot_height, width = plot_width)
cds_flow_dot
dev.off()

##########

tmp_plot_param <- c('Genic', 'Flow')

genic_flow_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = genic_mnp_stand, color = PLD_FLOW)) +
  ggtitle(paste('Standardized ', tmp_plot_param[1], ' MNPs with ',
    seq_depth_char, '+ Depth vs Seq Depth;\nColors According to ',
    tmp_plot_param[2], sep = ''))

tmp_plot_name <- paste(fig_out_dir, tmp_plot_param[1], '_MNP_0', seq_depth_char,
  '_by_', tmp_plot_param[2], '.pdf', sep = '')

pdf(tmp_plot_name, height = plot_height, width = plot_width)
genic_flow_dot
dev.off()

#########

tmp_plot_param <- c('CDS', 'SNPs')

cds_snp_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = cds_mnp_stand, color = PLD_SNP)) +
  ggtitle(paste('Standardized ', tmp_plot_param[1], ' MNPs with ',
    seq_depth_char, '+ Depth vs Seq Depth;\nColors According to ',
    tmp_plot_param[2], sep = ''))

tmp_plot_name <- paste(fig_out_dir, tmp_plot_param[1], '_MNP_0', seq_depth_char,
  '_by_', tmp_plot_param[2], '.pdf', sep = '')

pdf(tmp_plot_name, height = plot_height, width = plot_width)
cds_snp_dot
dev.off()

########

tmp_plot_param <- c('Genic', 'SNPs')

genic_snp_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = genic_mnp_stand, color = PLD_SNP)) +
  ggtitle(paste('Standardized ', tmp_plot_param[1], ' MNPs with ',
    seq_depth_char, '+ Depth vs Seq Depth;\nColors According to ',
    tmp_plot_param[2], sep = ''))

tmp_plot_name <- paste(fig_out_dir, tmp_plot_param[1], '_MNP_0', seq_depth_char,
  '_by_', tmp_plot_param[2], '.pdf', sep = '')

pdf(tmp_plot_name, height = plot_height, width = plot_width)
genic_snp_dot
dev.off()

#########

tmp_plot_param <- c('CDS', 'Cultivar')

cds_cultivar_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = cds_mnp_stand, color = Cultivar_ploidy)) +
  ggtitle(paste('Standardized ', tmp_plot_param[1], ' MNPs with ',
    seq_depth_char, '+ Depth vs Seq Depth;\nColors According to ',
    tmp_plot_param[2], sep = ''))

tmp_plot_name <- paste(fig_out_dir, tmp_plot_param[1], '_MNP_0', seq_depth_char,
  '_by_', tmp_plot_param[2], '.pdf', sep = '')

pdf(tmp_plot_name, height = plot_height, width = plot_width)
cds_cultivar_dot
dev.off()

#######

tmp_plot_param <- c('Genic', 'Cultivar')

genic_cultivar_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = genic_mnp_stand, color = Cultivar_ploidy)) +
  ggtitle(paste('Standardized ', tmp_plot_param[1], ' MNPs with ',
    seq_depth_char, '+ Depth vs Seq Depth;\nColors According to ',
    tmp_plot_param[2], sep = ''))

tmp_plot_name <- paste(fig_out_dir, tmp_plot_param[1], '_MNP_0', seq_depth_char,
  '_by_', tmp_plot_param[2], '.pdf', sep = '')
    
pdf(tmp_plot_name, height = plot_height, width = plot_width)
genic_cultivar_dot
dev.off()









metaploidy_2D_dot <- ggplot(data = filt_lib_df) + 
  geom_point(aes(x = seq_cov, y = mnp_stand, color = PLOIDY,
    shape = PLD_OLD_CALL)) + 
  ggtitle(paste(title_2D_prefix, title_metaploidy_suf, sep = ''))

fig_2D_metaploidy <- paste(fig_out_dir, fig_2D_prefix, fig_metaploidy_suf,
  sep = '')

pdf(fig_2D_metaploidy)
metaploidy_2D_dot
dev.off()

snpflow_2D_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = mnp_stand, color = PLD_SNP,
    shape = PLD_FLOW)) +
  ggtitle(paste(title_2D_prefix, title_snpflow_suf, sep = ''))

fig_2D_snpflow <- paste(fig_out_dir, fig_2D_prefix, fig_snpflow_suf,
  sep = '')

pdf(fig_2D_snpflow)
snpflow_2D_dot
dev.off()

ecotype_2D_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = mnp_stand, color = ECOTYPE_SNP_CHLR,
    shape = PLOIDY)) +
  ggtitle(paste(title_2D_prefix, title_ecotype_suf, sep = ''))

fig_2D_ecotype <- paste(fig_out_dir, fig_2D_prefix, fig_ecotype_suf,
  sep = '')

pdf(fig_2D_ecotype)
ecotype_2D_dot
dev.off()

subpop_2D_plot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = seq_cov, y = mnp_stand, color = SUBPOP_SNP,
    shape = PLOIDY)) +
  ggtitle(paste(title_2D_prefix, title_subpop_suf, sep = ''))

fig_2D_subpop <- paste(fig_out_dir, fig_2D_prefix, fig_subpop_suf,
  sep = '')

pdf(fig_2D_subpop)
subpop_2D_plot
dev.off()



