# Script to get info about the number and percentage of windows included
#  in a bp- and SNP-window combination

# module load python/3.7-anaconda-2019.07
# source activate local_PCA

args = commandArgs(trailingOnly=T)
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

head_in_short <- args[4]
#head_in_short <- 'CDS.expandv2.vcf.header.txt'
vcf_header_file <- paste(data_dir, head_in_short, sep = '')
vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

# Library info

# File that includes the file names of 4X, 8X, and bad library names
stand_input_file <- args[5]
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

snp_win_size <- as.numeric(args[6])
#snp_win_size <- 200

bp_window <- as.numeric(args[7])
bp_window <- 10000
chosen_window <- bp_window

### SET OUTPUTS ###

wind_info_out <- paste(data_dir, chr_name, '.', vcf_inbetween, '.', 
  SNP_window, 'SNP_', bp_window, 'bp_windows_info.txt', sep = '')

#######
#############

sys_com <- paste('ls ', vcf_pre, '*', sep = '')
vcf_files <- system(sys_com, intern = T)

chosen_window_list <- list()

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

  test_wind_snps <- get_window_inds(vcf = vcf_1_filt,
    min_pos = min(vcf_1_filt$POS), max_pos = max(vcf_1_filt$POS),
    bp_win_size = chosen_window)
  n_good_winds <- length(which(unlist(
    lapply(test_wind_snps[[2]], length)) >= snp_win_size))
  per_good_winds <- n_good_winds/length(test_wind_snps[[2]])

  chosen_window_list[[vf]] <- list()
  chosen_window_list[[vf]][['n_good_windows']] <- n_good_winds
  chosen_window_list[[vf]][['per_good_windows']] <- per_good_winds
}

chosen_good_window_vec <- unlist(lapply(chosen_window_list, 
  function(x) x[['n_good_windows']]))
total_good_windows <- sum(chosen_good_window_vec)
# for 1000 SNPs
# 403

per_good_window_vec <- unlist(lapply(chosen_window_list,
  function(x) x[['per_good_windows']]))

tot_window_vec <- chosen_good_window_vec/per_good_window_vec

total_potential_windows <- sum(tot_window_vec)
# For 1000 SNPs
#814

total_percent_good_windows <- total_good_windows/total_potential_windows
# For 1000 SNPs
# 0.495086
# When using 70k bp windows and 1000 SNPs, only include 50% of the windows 
#   across Chr01K

wind_info_df <- data.frame(chr = chr_name, SNP_win_size = snp_win_size, 
  bp_win_size = bp_window, n_good_windows = total_good_windows, 
  n_tot_windows = total_potential_windows, 
  per_good_windows = total_percent_good_windows, stringsAsFactors = F)

write.table(wind_info_df, wind_info_out, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')


out_text <- paste('For windows using ', snp_win_size, ' and ', bp_window, 
  'windows on ', chr_name, ', will include ', total_good_windows, ' of ',
  total_potential_windows, ' (', total_percent_good_windows, 
  '%) are included', sep = '')


write.table(out_text, wind_info_out, )








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

combo_n_windows <- apply(n_window_mat, 2, sum)
#100_SNPs  200_SNPs  300_SNPs  400_SNPs  500_SNPs  600_SNPs  700_SNPs  800_SNPs 
#   2520      1743      1179       924       753       648       572       505 
# 900_SNPs 1000_SNPs 
#      456       427 

mean_best_window <- apply(best_window_mat, 2, mean)
#100_SNPs  200_SNPs  300_SNPs  400_SNPs  500_SNPs  600_SNPs  700_SNPs  800_SNPs 
#10000.00  10000.00  18181.82  26363.64  34545.45  46363.64  54545.45  57272.73 
# 900_SNPs 1000_SNPs 
# 70000.00  74545.45 

median_best_window <- apply(best_window_mat, 2, median)
#100_SNPs  200_SNPs  300_SNPs  400_SNPs  500_SNPs  600_SNPs  700_SNPs  800_SNPs 
#   10000     10000     20000     20000     30000     40000     50000     50000 
# 900_SNPs 1000_SNPs 
#    70000     70000 

median_per_good_window <- apply(per_good_window_mat, 2, median)
#100_SNPs  200_SNPs  300_SNPs  400_SNPs  500_SNPs  600_SNPs  700_SNPs  800_SNPs 
#0.555294 0.3727679 0.4042553 0.4522613 0.4733333 0.5752212 0.6052632 0.5666667 
# 900_SNPs 1000_SNPs 
#0.6027397 0.6000000

window_df <- data.frame(snp_window = names(combo_n_windows), 
  n_windows = combo_n_windows, mean_best_window, median_best_window, 
  median_percent_good_window = median_per_good_window)

write.table(window_df, wind_test_out, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')


# How does chosen window size fit across genome

chosen_window <- median_best_window

chosen_window_list <- list()

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

  test_wind_snps <- get_window_inds(vcf = vcf_1_filt,
    min_pos = min(vcf_1_filt$POS), max_pos = max(vcf_1_filt$POS),
    bp_win_size = chosen_window)
  n_good_winds <- length(which(unlist(
    lapply(test_wind_snps[[2]], length)) >= snp_win_size))
  per_good_winds <- n_good_winds/length(test_wind_snps[[2]])

  chosen_window_list[[vf]] <- list()
  chosen_window_list[[vf]][['n_good_windows']] <- n_good_winds
  chosen_window_list[[vf]][['per_good_windows']] <- per_good_winds
}

chosen_good_window_vec <- unlist(lapply(chosen_window_list, 
  function(x) x[['n_good_windows']]))
total_good_windows <- sum(chosen_good_window_vec)
# for 1000 SNPs
# 403

per_good_window_vec <- unlist(lapply(chosen_window_list,
  function(x) x[['per_good_windows']]))

tot_window_vec <- chosen_good_window_vec/per_good_window_vec

total_potential_windows <- sum(tot_window_vec)
# For 1000 SNPs
#814

total_percent_good_windows <- total_good_windows/total_potential_windows
# For 1000 SNPs
# 0.495086
# When using 70k bp windows and 1000 SNPs, only include 50% of the windows 
#   across Chr01K


