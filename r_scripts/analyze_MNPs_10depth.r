# Script used for analyzing Genic MNPs

### LOAD PACKAGES ###
# Need to load the R_analysis conda environment
# module load python/3.7-anaconda-2019.07
# source activate R_analysis

library(ggplot2)

### SET INPUTS ###
# Directory with result files
count_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10/'

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
fig_out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10/'

fig_prefix <- 'mnp_hist_10depth_'

title_prefix <- 'MNP with 10+ depth per allele'
### SET VARIABLES ###

### SET CONSTANTS ###

##########

tmp_libs <- rownames(depth_mat)

# Load the MNP counts
count_list <- list()

for(LIB in tmp_libs){
  tmp_count_file <- paste(count_dir, LIB, '_10_MNP_total_count.txt', sep = '')
  count_list[[LIB]] <- read.table(tmp_count_file, header = F, sep = ' ',
    stringsAsFactors = F)
}

tot_mnp <- unlist(lapply(count_list, function(x) sum(as.numeric(x[,2]))))

tot_coverage <- apply(depth_mat, 1, sum)

tot_mnp_corrected <- (tot_mnp / tot_coverage) * min(tot_coverage)

mnp_df <- data.frame(lib = tmp_libs, tot_mnp = tot_mnp,
  mnp_corrected = tot_mnp_corrected, seq_depth = tot_coverage,
  stringsAsFactors = F)

mnp_df_no_out <- mnp_df[-which(mnp_df$mnp_corrected > 1000), ]

meta_rm_inds <- c()
for(rm_lib in meta_remove){
  tmp_ind <- which(mnp_df$lib == rm_lib)
  meta_rm_inds <- c(meta_rm_inds, tmp_ind)
}

mnp_df_filt <- mnp_df[-meta_rm_inds, ]

tot_hist <- ggplot(data = mnp_df) + geom_density(aes(x = tot_mnp),
  fill = 'grey50') + 
  ggtitle(paste('All ', title_prefix, sep = ''))

tot_hist_file <- paste(count_dir, fig_prefix, 'tot.pdf', sep = '')

pdf(tot_hist_file)
tot_hist
dev.off()

tot_hist_scale <- ggplot(data = mnp_df) + geom_density(aes(x = mnp_corrected),
  fill = 'grey50') + ggtitle(paste('Corrected ', title_prefix, sep = ''))

tot_hist_scale_file <- paste(count_dir, fig_prefix, 'scaled_tot.pdf', sep = '')

pdf(tot_hist_scale_file)
tot_hist_scale
dev.off()
######

no_out_hist <- ggplot(data = mnp_df_no_out) + geom_density(aes(x = tot_mnp),
  fill = 'grey50') +
  ggtitle(paste('All ', title_prefix, '\nNo Outlier Samples', sep = ''))

no_out_hist_file <- paste(count_dir, fig_prefix, 'no_out.pdf', sep = '')

pdf(no_out_hist_file)
no_out_hist
dev.off()

no_out_hist_scale <- ggplot(data = mnp_df_no_out) +
  geom_density(aes(x = mnp_corrected),
  fill = 'grey50') +
  ggtitle(paste('Corrected ', title_prefix, '\nNo Outlier Samples', sep = ''))

no_out_hist_scale_file <- paste(count_dir, fig_prefix, 'scaled_no_out.pdf', 
  sep = '')

pdf(no_out_hist_scale_file)
no_out_hist_scale
dev.off()

#######

filt_hist <- ggplot(data = mnp_df_filt) + geom_density(aes(x = tot_mnp),
  fill = 'grey50') +
  ggtitle(paste('All ', title_prefix, '\nFiltered Samples',
    sep = ''))

filt_hist_file <- paste(count_dir, fig_prefix, 'filt.pdf', sep = '')

pdf(filt_hist_file)
filt_hist
dev.off()

filt_hist_scale <- ggplot(data = mnp_df_filt) +
  geom_density(aes(x = mnp_corrected),
  fill = 'grey50') +
  ggtitle(paste('Corrected ', title_prefix, '\nFiltered Samples', sep = ''))

filt_hist_scale_file <- paste(count_dir, fig_prefix, 'scaled_filt.pdf', 
  sep = '')

pdf(filt_hist_scale_file)
filt_hist_scale
dev.off()

###########

quit(save = 'no')


