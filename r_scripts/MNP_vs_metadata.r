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

genome_region <- 'genic'
#genome_region <- as.character(args[2])
# 'CDS', 'genic'

# Counts of genic MNPSs with 10+ depth in each allele
mnp_rds <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_0', seq_depth_char, 
  '/', genome_region, '_MNP_count_list_0', seq_depth_char, 'depth.rds', 
  sep = '')
mnp_list <- readRDS(mnp_rds)

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

### SET OUTPUTS ###
fig_out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/', 
  'pos_0', seq_depth_char, '/', sep = '')

fig_1D_prefix <- paste(genome_region, '_MNP_0', seq_depth_char, 
  'depth_dotplot_by_', sep = '')

fig_2D_prefix <- paste(genome_region, '_MNP_0', seq_depth_char, 
  'depth_dot2D_by_', sep = '')

fig_metaploidy_suf <- 'MetaPloidy.pdf'
fig_snpflow_suf <- 'SNPandFlow.pdf'
fig_ecotype_suf <- 'EcotypeandMetaPloidy.pdf'
fig_subpop_suf <- 'SubpopandMetaPloidy.pdf'

### SET VARIABLES ###

title_1D_prefix <- paste('Standardized ', genome_region, ' MNPs with ', 
  seq_depth_char, '+ Depth; Colors and Shapes\nAccording to ', sep = '')

title_2D_prefix <- paste('Standardized ', genome_region, ' MNPs with ',
  seq_depth_char, '+ Depth vs Seq Depth;\nColors and Shapes According to ')

title_metaploidy_suf <- 'Current and Old Metadata Ploidy Designation'
title_snpflow_suf <- 'SNP and Flow Ploidy Designation'
title_ecotype_suf <- 'SNP Ecotype and Meta Ploidy Designation'
title_subpop_suf <- 'SNP Subpopulation and Meta Ploidy Designation'

### SET CONSTANTS ###

##########

#sum(names(mnp_list) %in% meta$LIBRARY)
# [1] 921

miss_inds <- which(names(mnp_list) %in% meta$LIBRARY == F)
miss_meta <- which(meta$LIBRARY %in% names(mnp_list) == F)

mnp_overlap <- which(names(mnp_list) %in% meta$LIBRARY)

lib_mnp_mat <- matrix(data = unlist(mnp_list), byrow = T,
  nrow = length(mnp_list))
rownames(lib_mnp_mat) <- names(mnp_list)
colnames(lib_mnp_mat) <- chrom_names

tot_mnp <- apply(lib_mnp_mat, 1, sum)

# Get the sequencing output for each sample
tot_coverage <- apply(depth_mat, 1, sum)

tot_mnp_corrected <- (tot_mnp / tot_coverage) * min(tot_coverage)

lib_info_0 <- data.frame(lib = names(tot_mnp), mnp = tot_mnp,
  mnp_stand = tot_mnp_corrected, seq_cov = tot_coverage,
  stringsAsFactors = F)

lib_info_df <- lib_info_0[mnp_overlap,]

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
  lib_info_df$PLOIDY[i] <- meta$PLOIDY[tmp_meta_ind]
  lib_info_df$PLD_OLD_CALL[i] <- meta$PLD_OLD_CALL[tmp_meta_ind]
  lib_info_df$PLD_FLOW[i] <- meta$PLD_FLOW[tmp_meta_ind]
  lib_info_df$PLD_SNP[i] <- meta$PLD_SNP[tmp_meta_ind]
  lib_info_df$ECOTYPE_SNP_CHLR[i] <- meta$ECOTYPE_SNP_CHLR[tmp_meta_ind]
  lib_info_df$SUBPOP_SNP[i] <- meta$SUBPOP_SNP[tmp_meta_ind]
}

lib_info_df$x_loc <- 1

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

metaploidy_1D_dot <- ggplot(data = filt_lib_df) + 
  geom_point(aes(x = x_loc, y = mnp_stand, color = PLOIDY, 
    shape = PLD_OLD_CALL),  
    position = position_jitter(w = 0.1, h = 0))  +
  ggtitle(paste(title_1D_prefix, title_metaploidy_suf, sep = ''))

fig_1D_metaploidy <- paste(fig_out_dir, fig_1D_prefix, fig_metaploidy_suf, 
  sep = '')

pdf(fig_1D_metaploidy)
metaploidy_1D_dot
dev.off()

snpflow_1D_dot <- ggplot(data = filt_lib_df) + 
  geom_point(aes(x = x_loc, y = mnp_stand, color = PLD_SNP, 
    shape = PLD_FLOW), 
    position = position_jitter(w = 0.1, h = 0))  +
  ggtitle(paste(title_1D_prefix, title_snpflow_suf, sep = ''))

fig_1D_snpflow <- paste(fig_out_dir, fig_1D_prefix, fig_snpflow_suf,
  sep = '')

pdf(fig_1D_snpflow)
snpflow_1D_dot
dev.off()

ecotype_1D_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = x_loc, y = mnp_stand, color = ECOTYPE_SNP_CHLR,
    shape = PLOIDY), 
    position = position_jitter(w = 0.1, h = 0))  +
  ggtitle(paste(title_1D_prefix, title_ecotype_suf, sep = ''))

fig_1D_ecotype <- paste(fig_out_dir, fig_1D_prefix, fig_ecotype_suf,
  sep = '')

pdf(fig_1D_ecotype)
ecotype_1D_dot
dev.off()

subpop_1D_dot <- ggplot(data = filt_lib_df) +
  geom_point(aes(x = x_loc, y = mnp_stand, color = SUBPOP_SNP,
    shape = PLOIDY), 
    position = position_jitter(w = 0.1, h = 0))  +
  ggtitle(paste(title_1D_prefix, title_subpop_suf, sep = ''))

fig_1D_subpop <- paste(fig_out_dir, fig_1D_prefix, fig_subpop_suf,
  sep = '')

pdf(fig_1D_subpop)
subpop_1D_dot
dev.off()

########
# Use # MNP vs seq cov to better discriminate the groups

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



