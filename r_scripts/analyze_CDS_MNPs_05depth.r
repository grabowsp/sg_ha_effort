# Script used for analyzing Genic MNPs

### LOAD PACKAGES ###
# Need to load the R_analysis conda environment
# module load python/3.7-anaconda-2019.07
# source activate R_analysis

library(ggplot2)

### SET INPUTS ###
# Counts of genic MNPSs with 10+ depth in each allele
mnp_rds <- '/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_05/CDS_MNP_count_list_05depth.rds'
mnp_list <- readRDS(mnp_rds)

# Get chromosome names
chrom_name_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/chr_names.txt'
chrom_names <- as.vector(
  read.table(chrom_name_file, header = F, stringsAsFactors = F)[,1])

# Seq depth mat
depth_mat_rds <- '/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth/sample_seq_depth_mat.rds'
depth_mat <- readRDS(depth_mat_rds)

# BAD samples from metadata
meta_lib_remove_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/sg_libs_to_remove_from_meta_Dec2019.txt'
meta_remove <- as.vector(read.table(meta_lib_remove_file, header = T,
  stringsAsFactors = F)[,1])

### SET OUTPUTS ###
fig_out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_05/'

fig_prefix <- 'CDS_MNP_05depth_'

title_prefix <- 'CDS MNPs with 5+ reads per allele'

### SET VARIABLES ###

### SET CONSTANTS ###

##########

lib_mnp_mat <- matrix(data = unlist(mnp_list), byrow = T,
  nrow = length(mnp_list))
rownames(lib_mnp_mat) <- names(mnp_list)
colnames(lib_mnp_mat) <- chrom_names

tot_mnp <- apply(lib_mnp_mat, 1, sum)
#summary(tot_mnp)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   196    3432    4656    4630    5360   14830  

# Get the sequencing output for each sample
tot_coverage <- apply(depth_mat, 1, sum)

# check that samples are in the same order
#sum(rownames(depth_mat) != rownames(lib_mnp_mat))
#[1] 0

tot_mnp_corrected <- (tot_mnp / tot_coverage) * min(tot_coverage)
#summary(tot_mnp_corrected)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  8.095   35.390   42.687   51.382   47.866 2013.156

lib_info_df <- data.frame(lib = names(tot_mnp), mnp = tot_mnp,
  mnp_stand = tot_mnp_corrected, seq_cov = tot_coverage,
  stringsAsFactors = F)

rm_inds <- c()
for(mr in meta_remove){
  tmp_ind <- which(lib_info_df$lib == mr)
  rm_inds <- c(rm_inds, tmp_ind)
}

filt_lib_df <- lib_info_df[-rm_inds, ]

tot_hist <- ggplot(data = lib_info_df) + geom_density(aes(x = mnp),
  fill = 'grey50') + ggtitle(title_prefix)

tot_hist_file <- paste(fig_out_dir, fig_prefix,'hist.pdf', sep = '')

pdf(tot_hist_file)
tot_hist
dev.off()

########
stand_hist <- ggplot(data = lib_info_df) + 
  geom_density(aes(x = mnp_stand),
  fill = 'grey50') + 
  ggtitle(paste(title_prefix, '\nstandardized by sequencing output', sep = ''))

stand_hist_file <- paste(fig_out_dir, fig_prefix, 'stand_hist.pdf', 
  sep = '')

pdf(stand_hist_file)
stand_hist
dev.off()
########

filt_hist <- ggplot(data = filt_lib_df) + 
  geom_density(aes(x = mnp), fill = 'grey50') + 
  ggtitle(paste(title_prefix, ' in Filtered Samples', sep = ''))

filt_hist_file <- paste(fig_out_dir, fig_prefix, 'hist_filt.pdf', 
  sep = '')

pdf(filt_hist_file)
filt_hist
dev.off()
##########

filt_stand_hist <- ggplot(data = filt_lib_df) +
  geom_density(aes(x = mnp_stand), fill = 'grey50') + 
  ggtitle(paste(title_prefix, 
    ' in Filtered Samples\nStandardized by Sequencing Output', sep = ''))

filt_stand_hist_file <- paste(fig_out_dir, 
  fig_prefix, 'stand_hist_filt.pdf', sep = '')

pdf(filt_stand_hist_file)
filt_stand_hist
dev.off()



