# Script to find the bp window size that best fits the desired SNP window size

# module load python/3.7-anaconda-2019.07
# source activate local_PCA

args <- commandArgs(trailingOnly=T)
#rundir_args <- commandArgs(trailingOnly = F)

library(data.table)
library(lostruct)

repo_base_dir <- args[1]
# repo_base_dir <- '/global/homes/g/grabowsp/tools/'
if(rev(unlist(strsplit(repo_base_dir, split = '')))[1] != '/'){
  repo_base_dir <- paste(repo_base_dir, '/', sep = '')
}

locpca_func_file_short <- 'sg_ha_effort/polyploid_genos/local_PCA/localPCA_functions.r'
locpca_func_file <- paste(repo_base_dir, locpca_func_file_short, sep = '')
source(locpca_func_file)

### LOAD DATA ###

# VCF info
data_dir <- args[2]
# data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/'
if(rev(unlist(strsplit(data_dir, split = '')))[1] != '/'){
  data_dir <- paste(data_dir, '/', sep = '')
}

chr_name <- args[3]
# chr_name <- 'Chr01K'

# the text inbetween the chromosome name and the end of the name 
vcf_inbetween <- args[4]
#vcf_inbetween <- 'polyploid.CDS.expandv2'

vcf_pre <- paste(data_dir, paste(chr_name, vcf_inbetween, 'vcf_', sep = '.'), 
  sep = '')

head_in_short <- args[5]
#head_in_short <- 'CDS.expandv2.vcf.header.txt'
vcf_header_file <- paste(data_dir, head_in_short, sep = '')
vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

# Library info

# File that includes the file names of 4X, 8X, and bad library names
stand_input_file <- args[6]
# stand_input_file <- paste(repo_base_dir, 'sg_ha_effort/polyploid_genos/local_PCA/standard_input_files.r', sep = '')
source(stand_input_file)

#tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs <- intersect(tet_libs_0, vcf_header)

#oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs <- intersect(oct_libs_0, vcf_header)

#remove_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/local_pca_remove_samps_09152020.txt'
remove_libs <- as.vector(read.table(remove_lib_file, header = F,
  stringsAsFactors = F)[,1])

### SET VARIABLES ###

# set the interval of SNP windows to test
# can make these adjustable if need be
min_SNP_win_in <- 100
max_SNP_win_in <- 1000
snp_wind_interval <- 100
snp_win_size_vec <- seq(min_SNP_win_in, max_SNP_win_in, by = snp_wind_interval)

# set the interval of bp windows to test
wind_interval <- 1e4
max_wind_size <- 3e5
test_bp_window_sizes <- seq(from = wind_interval, to = max_wind_size,
  by = wind_interval)

### SET OUTPUTS ###

wind_test_out <- paste(data_dir, chr_name, '.', vcf_inbetween, '.', 
  'windowtests.txt', sep = '')

#######
#############

sys_com <- paste('ls ', vcf_pre, '*', sep = '')
vcf_files <- system(sys_com, intern = T)

#vcf_in <- vcf_files[1]

# Choose best window size

#snp_win_size <- 200
#snp_win_size_vec <- seq(100,1000, by = 100)

#wind_interval <- 1e4
#max_wind_size <- 3e5
#test_bp_window_sizes <- seq(from = wind_interval, to = max_wind_size,
#  by = wind_interval)

# vf <- 1

window_test_list <- list()

for(vf in seq(length(vcf_files))){
  print(paste('vcf subfile ', vf, sep = ''))
  vcf_in <- vcf_files[vf]
  vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')
  colnames(vcf) <- vcf_header
  #
  geno_mat_1 <- process_vcf(vcf = vcf, oct_libs = oct_libs, tet_libs = tet_libs,
    rm_libs = remove_libs)
  genomat_ngenos <- apply(geno_mat_1, 1,
    function(x) length(setdiff(unique(x), NA)))
  invar_loci <- which(genomat_ngenos == 1)
  geno_mat_1_filt <- geno_mat_1[-invar_loci, ]
  vcf_1_filt <- vcf[-invar_loci, ]
  #### Look at the bp sizes for the desired number of SNPs
#  n_good_vec <- c()
  tot_pos_wind_vec <- c()

  ####
  n_good_vec_list <- list()
#  tot_pos_wind_vec_list <- list()
  for(swzv in snp_win_size_vec){
    n_good_vec_list[[paste(swzv, '_SNPs', sep = '')]] <- c(NA)
#    tot_pos_wind_vec_list[[as.character(swzv)]] <- c(NA)
  }
  ######
  for(i in test_bp_window_sizes){
    test_wind_snps <- get_window_inds(vcf = vcf_1_filt, 
      min_pos = min(vcf_1_filt$POS), max_pos = max(vcf_1_filt$POS), 
      bp_win_size = i)
    tot_pos_wind_vec <- c(tot_pos_wind_vec, length(test_wind_snps[[2]]))    
    # find the windows with enough SNPs for each of the SNP window sizes
    for(j in snp_win_size_vec){
      n_good_winds <- length(which(unlist(
        lapply(test_wind_snps[[2]], length)) >= j))
      n_good_vec_list[[paste(j, '_SNPs', sep = '')]] <- c(
        n_good_vec_list[[paste(j, '_SNPs', sep = '')]], n_good_winds)
    }
  }
  n_good_vec_list <- lapply(n_good_vec_list, function(x) x[-which(is.na(x))])  
  max_good_wind_size_vec <- unlist(lapply(n_good_vec_list, function(x) 
    test_bp_window_sizes[which.max(x)]))
  good_wind_n_vec <- unlist(lapply(n_good_vec_list, max, na.rm = T))
  per_max_good_wind_vec <- unlist(lapply(n_good_vec_list, function(x)
    max(x, na.rm = T)/tot_pos_wind_vec[which.max(x)]))
  window_test_list[[vf]] <- list()
  window_test_list[[vf]][['best_window_vec']] <- max_good_wind_size_vec
  window_test_list[[vf]][['n_good_windows_vec']] <- good_wind_n_vec
  window_test_list[[vf]][['per_good_windows_vec']] <- per_max_good_wind_vec
}


best_window_mat <- matrix(data = unlist(
  lapply(window_test_list, function(x) x[['best_window_vec']])), 
  nrow = length(window_test_list), byrow = T)
colnames(best_window_mat) <- names(window_test_list[[1]][['best_window_vec']])

n_window_mat <- matrix(data = unlist(
  lapply(window_test_list, function(x) x[['n_good_windows_vec']])),
  nrow = length(window_test_list), byrow = T)
colnames(n_window_mat) <- names(window_test_list[[1]][['best_window_vec']])

per_good_window_mat <-  matrix(data = unlist(
  lapply(window_test_list, function(x) x[['per_good_windows_vec']])),
  nrow = length(window_test_list), byrow = T)
colnames(per_good_window_mat) <- names(
  window_test_list[[1]][['best_window_vec']])

# number of good windows summed across all VCFs
combo_n_windows <- apply(n_window_mat, 2, sum)

# mean, across VCFs, of bp window size that produces the most windows
mean_best_window <- apply(best_window_mat, 2, mean)

# median, across VCFs, of bp window size that produces the most windows
median_best_window <- apply(best_window_mat, 2, median)

# median, across VCFs, of percent good windows in best bp window size
median_per_good_window <- apply(per_good_window_mat, 2, median)

window_df <- data.frame(chr = chr_name, snp_window = names(combo_n_windows), 
  n_windows = combo_n_windows, mean_best_window, median_best_window, 
  median_percent_good_window = median_per_good_window)

write.table(window_df, wind_test_out, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')


