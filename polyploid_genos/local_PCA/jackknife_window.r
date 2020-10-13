# Workflow
# 0) Choose which samples to omit
# A) Load file names for sub-vcf
# 1) Load first sub-vcf
# 2) Process genotypes
# 3) Select windows for sub-vcf
# 4) For all but final sub-vcf, remove final window to use with next vcf
# 4) Select 

# module load python/3.7-anaconda-2019.07
# source activate local_PCA

library(data.table)
library(lostruct)

repo_base_dir <- '/global/homes/g/grabowsp/tools/'
locpca_func_file_short <- 'sg_ha_effort/polyploid_genos/local_PCA/localPCA_functions.r'
locpca_func_file <- paste(repo_base_dir, locpca_func_file_short, sep = '')
source(locpca_func_file)

### LOAD DATA ###

# VCF info
#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/'
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/'

chr_name <- 'Chr01K'

vcf_inbetween <- 'polyploid.CDS.expandv2'

vcf_pre <- paste(data_dir, paste(chr_name, vcf_inbetween, 'vcf_', sep = '.'), 
  sep = '')

head_in_short <- 'CDS.expandv2.vcf.header.txt'
vcf_header_file <- paste(data_dir, head_in_short, sep = '')
vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

# Library info

tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs <- intersect(tet_libs_0, vcf_header)

oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs <- intersect(oct_libs_0, vcf_header)

remove_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/local_pca_remove_samps_09152020.txt'
remove_libs <- as.vector(read.table(remove_lib_file, header = F,
  stringsAsFactors = F)[,1])

#############

sys_com <- paste('ls ', vcf_pre, '*', sep = '')
vcf_files <- system(sys_com, intern = T)

vcf_in <- vcf_files[1]

bp_win_size <- 3e4
#snp_win_size <- 1000
snp_win_size <- 500
test_k <- 2
pcs_to_check <- c(1,2)

n_partitions <- 10

keep_wind_start_pos <- c()
keep_snp_pos <- c()

# vf <- 2

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
  if(vf != 1){
    if(length(remain_inds) > 0){
      geno_mat_1_filt <- rbind(remain_mat, geno_mat_1_filt)
      vcf_1_filt <- rbind(remain_vcf, vcf_1_filt)
    }
  }
  #####
  # Select windows to keep for analysis
  #min_pos <- min(vcf_filt$POS)
  if(vf == 1){
    min_pos <- 1
  } else{
    min_pos <- next_min_pos
  }
  max_pos <- max(vcf_1_filt$POS)
  #
  wind_list_out_1 <- get_window_inds(vcf = vcf_1_filt, min_pos = min_pos,
    max_pos = max_pos, bp_win_size = bp_win_size)
  # the starting position of the last window; is used for making windows
  #   in the next sub_vcf
  tmp_wind_start_pos <- wind_list_out_1[[1]]
  next_min_pos <- rev(tmp_wind_start_pos)[1]
  wind_list_1 <- wind_list_out_1[[2]]
  # remove the last window to combine with next sub_vcf except for the last
  #   sub_vcf
  if(vf == length(vcf_files)){
    wind_list <- wind_list_1
  } else{
    wind_list <- wind_list_1[-length(wind_list_1)]
    tmp_wind_start_pos <- tmp_wind_start_pos[-length(tmp_wind_start_pos)]
  }
  if(vf != length(vcf_files)){
    remain_inds <- unlist(wind_list_1[length(wind_list_1)])
    remain_vcf <- vcf_1_filt[remain_inds, ]
    remain_mat <- geno_mat_1_filt[remain_inds, ]
  }
  # Generate list of SNP indices to use for jackknife approach
  tmp_jack_inds <- get_jackknife_window_inds(bp_window_list = wind_list,
    snp_win_size = snp_win_size, n_part = n_partitions)
  # Select all jackknife indices to run using full SNP set
  tot_sub_inds <- sort(unique(unlist(tmp_jack_inds)))
  vcf_sub_filt <- vcf_1_filt[tot_sub_inds, ]
  geno_mat_sub <- geno_mat_1_filt[tot_sub_inds, ]
  # Run eigen_windows on full set of indices
  if(vf != length(vcf_files)){
    remain_inds <- unlist(wind_list_1[length(wind_list_1)])
    remain_vcf <- vcf_1_filt[remain_inds, ]
    remain_mat <- geno_mat_1_filt[remain_inds, ]
  }
  # run eigen_windows
  tmp_fn_ew <- eigen_window_fixNAs(test_mat = geno_mat_sub,
    win_size = snp_win_size, test_k = test_k, fna.verbose = T)
  if(vf == 1){
    tot_ew <- tmp_fn_ew
  } else{
    tot_ew <- rbind(tot_ew, tmp_fn_ew)
  }
  # 
  # Set partition and remaining window SNP sizes
  jack_part_size <- floor(snp_win_size/n_partitions)
  jack_win_size <- snp_win_size - jack_part_size
  # Calculate jackknife eigen-windows
  jack_ew_list <- list()
  for(tji in seq(length(tmp_jack_inds))){
    tmp_mat_sub <- geno_mat_1_filt[tmp_jack_inds[[tji]], ]
    tmp_ew_res <- eigen_window_fixNAs(test_mat = tmp_mat_sub,
      win_size = jack_win_size, test_k = test_k, fna.verbose = T)
    jack_ew_list[[tji]] <- tmp_ew_res
#    if(vf == 1){
#      jack_ew_list[[tji]] <- tmp_ew_res
#    } else{
#      jack_ew_list[[tji]] <- rbind(jack_ew_list[[tji]], tmp_ew_res)
#    }
  }
  if(vf == 1){
    tot_jack_ew <- jack_ew_list
  }  else {
    for(tje in seq(length(tot_jack_ew))){
      tot_jack_ew[[tje]] <- rbind(tot_jack_ew[[tje]], jack_ew_list[[tje]])
    }
  }
}

# Calculate variance-signal from full eigen_windows

pc1_tot_vec_cols <- grep('PC_1_', colnames(tot_ew))

tot_pc_mat <- tot_ew[, pc1_tot_vec_cols]
for(tpm in c(2:nrow(tot_pc_mat))){
  test_comp_1 <- sum(abs(tot_pc_mat[1, ] - tot_pc_mat[tpm, ]))
  test_comp_2 <- sum(abs(tot_pc_mat[1, ] + tot_pc_mat[tpm, ]))
  if(test_comp_2 < test_comp_1){
    tot_pc_mat[tpm, ] <- tot_pc_mat[tpm, ] * -1
  }
}
mean_pc_val <- apply(tot_pc_mat, 2, mean)
#
tot_pc_dif_sq <- matrix(NA, ncol = ncol(tot_pc_mat), nrow = nrow(tot_pc_mat))
for(ptds in seq(nrow(tot_pc_dif_sq))){
  tmp_dif_sq <- (tot_pc_mat[ptds, ] - mean_pc_val)^2
  tot_pc_dif_sq[ptds, ] <- tmp_dif_sq
}
#
mean_pc_dif <- apply(tot_pc_dif_sq, 1, mean)
tot_mean_variance <- mean(mean_pc_dif)
# this is the mean variance for signal
# 1000 SNPs
# [1] 0.0002419491
# 500 SNPs - prelim
# [1] 0.0004594831
# 250 SNPs
# [1] 0.0005042384
# 500 SNPs - 30kb windows, 3 sub-vcfs
# [1] 0.0006465006

# Calculate variance-noise from jackknive eigen_windows
pc1_vec_cols <- grep('PC_1_', colnames(tot_jack_ew[[1]]))

# wind1_pc1_jack1 <- jack_ew_list[[1]][1,pc1_vec_cols]
jack_pc_mat <- matrix(NA, nrow = length(tot_jack_ew), 
  ncol = length(pc1_vec_cols))
window_se_vec <- c()

for(wn in seq(nrow(tot_jack_ew[[1]]))){
  jack_pc_mat <- matrix(NA, nrow = length(tot_jack_ew),
    ncol = length(pc1_vec_cols))
  jack_pc_mat[1, ] <- tot_jack_ew[[1]][wn,pc1_vec_cols]
  for(jp in c(2:nrow(jack_pc_mat))){
    tmp_eigvec <- tot_jack_ew[[jp]][wn,pc1_vec_cols]
    test_dif_1 <- sum(abs(jack_pc_mat[1, ] - tmp_eigvec))
    test_dif_2 <- sum(abs(jack_pc_mat[1, ] + tmp_eigvec))
    if(test_dif_2 < test_dif_1){
      jack_pc_mat[jp, ] <- tmp_eigvec * -1
    } else{
      jack_pc_mat[jp, ] <- tmp_eigvec
    }
  }
  #
  mean_jack_vec <- apply(jack_pc_mat, 2, mean)
  #
  jack_samp_se <- c()
  for(mjv in seq(length(mean_jack_vec))){
    tmp_jack_dif_sq <- (jack_pc_mat[, mjv] - mean_jack_vec[mjv])^2
    tmp_samp_se <- sum(tmp_jack_dif_sq)*(9/10)
    jack_samp_se <- c(jack_samp_se, tmp_samp_se)
  }
  #
  jack_window_se <- mean(jack_samp_se)
  window_se_vec <- c(window_se_vec, jack_window_se)
}
  
noise_variance <- mean(window_se_vec)
# 1000 SNPs
#  [1] 0.0001921484
# 500 SNPs - prelim
# [1] 0.0003711886
# 250 SNPs
# [1] 0.0004971984
# 500 SNPs - 30kb windows, 3 sub-vcf
# [1] 0.0003891591






  # Set partition and remaining window SNP sizes
  jack_part_size <- floor(snp_win_size/n_partitions)
  jack_win_size <- snp_win_size - jack_part_size
  # Calculate jackknife eigen-windows
  jack_ew_list <- list()
  for(tji in seq(length(tmp_jack_inds))){
    tmp_mat_sub <- geno_mat_1_filt[tmp_jack_inds[[tji]], ]
    tmp_ew_res <- eigen_window_fixNAs(test_mat = tmp_mat_sub, 
      win_size = jack_win_size, test_k = test_k, fna.verbose = T)
    if(vf == 1){
      jack_ew_list[[tji]] <- tmp_ew_res
    } else{
      jack_ew_list[[tji]] <- rbind(jack_ew_list[[tji]], tmp_ew_res)
    }
  }
  # 
# then can re-incorporate starting from POS=1, etc


  ######
  # get sub-samples indices for SNPs in windows to keep and sub-select VCF
  #   and genotype matrix
#  good_inds_1 <- get_subsamp_window_inds(bp_window_list = wind_list,
#    snp_win_size = snp_win_size)
  tot_sub_inds <- sort(unique(unlist(tmp_jack_inds)))
  vcf_sub_filt <- vcf_1_filt[tot_sub_inds, ]
  geno_mat_sub <- geno_mat_1_filt[tot_sub_inds, ]
#  keep_snp_pos <- c(keep_snp_pos, vcf_sub_filt$POS)
  # get window starting positions for retained windows
#  keep_wind_start_pos <- c(keep_wind_start_pos, tmp_wind_start_pos[which(unlist(
#    lapply(wind_list, length)) >= snp_win_size)])
  #
  if(vf != length(vcf_files)){
    remain_inds <- unlist(wind_list_1[length(wind_list_1)])
    remain_vcf <- vcf_1_filt[remain_inds, ]
    remain_mat <- geno_mat_1_filt[remain_inds, ]
  }
  # run eigen_windows
  tmp_fn_ew <- eigen_window_fixNAs(test_mat = geno_mat_sub, 
    win_size = snp_win_size, test_k = test_k, fna.verbose = T)


  #
  # remove the last window to combine with next sub_vcf except for the last
  #   sub_vcf

  pc1_tot_vec_cols <- grep('PC_1_', colnames(tmp_fn_ew))


  tot_pc_mat <- tmp_fn_ew[, pc1_tot_vec_cols]
  for(tpm in c(2:nrow(tot_pc_mat))){
    test_comp_1 <- sum(abs(tot_pc_mat[1, ] - tot_pc_mat[tpm, ]))
    test_comp_2 <- sum(abs(tot_pc_mat[1, ] + tot_pc_mat[tpm, ]))
    if(test_comp_2 < test_comp_1){
      tot_pc_mat[tpm, ] <- tot_pc_mat[tpm, ] * -1
    }
  }


mean_pc_val <- apply(tot_pc_mat, 2, mean)

tot_pc_dif_sq <- matrix(NA, ncol = ncol(tot_pc_mat), nrow = nrow(tot_pc_mat))
for(ptds in seq(nrow(tot_pc_dif_sq))){
  tmp_dif_sq <- (tot_pc_mat[ptds, ] - mean_pc_val)^2
  tot_pc_dif_sq[ptds, ] <- tmp_dif_sq
}


mean_pc_dif <- apply(tot_pc_dif_sq, 1, mean)

tot_mean_variance <- mean(mean_pc_dif)
# this is the mean variance for signal
# 1000 SNPs
# [1] 0.0002419491
# 500 SNPs
# [1] 0.0004594831


tot_pc_mat <- matrix(NA, ncol = ncol(geno_mat_1_filt), 
  nrow = nrow(tmp_fn_ew))




  if(vf == 1){
    tot_ew <- tmp_fn_ew
  } else{
    tot_ew <- rbind(tot_ew, tmp_fn_ew)
  }
}

attr(tot_ew, 'npc') <- test_k
attr(tot_ew, 'w') <- 1

tot_dist <- pc_dist(tot_ew)
tot_cmd <- cmdscale(tot_dist, k = 100)

dist_tot_var <- sum(apply(tot_cmd, 2, var))
dist_per_var <- (apply(tot_cmd, 2, var)/dist_tot_var)*100

cmd_df <- data.frame(tot_cmd, stringsAsFactors = F)
colnames(cmd_df) <- paste('PCo_', seq(ncol(cmd_df)), sep = '')

library(ggplot2)

gg_c1 <- ggplot(cmd_df, aes(x = PCo_1, y = PCo_2)) +
  geom_point() +
  ggtitle('Chr01K expand_samps local PCA\n100kb, 1k SNP windows')

cmd_pdf_out_file <- paste(data_dir, 'localPCA_Chr01K_test1.pdf', sep = '')

pdf(cmd_pdf_out_file, width = 7, height = 7)
gg_c1
dev.off()






#####################



# NEXT - calculate distances

test_res <- tot_ew

attr(test_res, 'npc') <- 2
attr(test_res, 'w') <- 1

#tot_dist <- pc_dist(tot_ew)
tot_dist <- pc_dist(test_res)

tmp_cmd <- cmdscale(tmp_dist, k = 100)

library(ggplot2)

gg_c1 <- ggplot()


#############################









vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')
colnames(vcf) <- vcf_header

geno_mat_1 <- process_vcf(vcf = vcf, oct_libs = oct_libs, tet_libs = tet_libs, 
  rm_libs = remove_libs)

genomat_ngenos <- apply(geno_mat_1, 1, 
  function(x) length(setdiff(unique(x), NA)))

invar_loci <- which(genomat_ngenos == 1)

geno_mat_1_filt <- geno_mat_1[-invar_loci, ]

vcf_1_filt <- vcf[-invar_loci, ]





###

# Select SNPs to use with bp-based windows
#  window must have minimum number of SNPs
#  subsample SNPs so each window has same number of SNPs
# - This approach generates the matrix in the correct form to use
#     the standard eigen_windows() function

bp_win_size <- 1e5
snp_win_size <- 1000

#min_pos <- min(vcf_filt$POS)
min_pos <- 1
max_pos <- max(vcf_1_filt$POS)



wind_list_out_1 <- get_window_inds(vcf = vcf_1_filt, min_pos = min_pos, 
  max_pos = max_pos, bp_win_size = bp_win_size)

next_min_pos <- wind_list_out_1[[1]]

wind_list_1 <- wind_list_out_1[[2]]

wind_list <- wind_list_1[-length(wind_list_1)]

good_inds_1 <- get_subsamp_window_inds(bp_window_list = wind_list, 
  snp_win_size = snp_win_size)

vcf_sub_filt <- vcf_1_filt[good_inds_1, ]

geno_mat_sub <- geno_mat_1_filt[good_inds_1, ]

remain_inds <- unlist(wind_list_1[length(wind_list_1)])
remain_vcf <- vcf_1_filt[remain_inds, ]
remain_mat <- geno_mat_1_filt[remain_inds, ]

# run eigen_windows here

# Process remaining vcfs next

for(vf in c(2:(length(vcf_files)-1))){
  print(vf)
}



################

# Run eigen_windows
#  for each window with NA's
#    select overrepresented samples with NA's
#    substitute half of NA's with rowMean for the sample at that SNP
#    rerun eigen_windows
#    repeat if still have windows with NA's


win_size <- snp_win_size
test_mat <- geno_mat_sub
test_k <- 2

#############

# NEXT - calculate distances

tmp_dist <- pc_dist(tmp_ew)
tmp_cmd <- cmdscale(tmp_dist, k = 100)


#####################
# Generate bp-length-based windows

vcf_filt <- vcf[-invar_loci, ]

min_pos <- min(vcf_filt$POS)
max_pos <- max(vcf_filt$POS)
bp_win_size <- 1e4
snp_win_size <- 100

bp_start_pos <- seq(min_pos, max_pos, by = bp_win_size)
bp_end_pos <- c((bp_start_pos[-1] - 1), max_pos)

bp_window_list <- list()

for(i in seq(length(bp_start_pos))){
  tmp_inds <- which(vcf_filt$POS >= bp_start_pos[i] & 
    vcf_filt$POS < bp_end_pos[i])
  bp_window_list[[i]] <- tmp_inds
}

nsnps_per_wind <- unlist(lapply(bp_window_list, length))

bad_windows <- which(nsnps_per_wind < 50)

good_windows <- which(nsnps_per_wind >= snp_win_size)

bp_sub_window_list <- list()

for(j in seq(length(good_windows))){
  tmp_wind_inds <- bp_window_list[[good_windows[j]]]
  sub_inds <- sort(sample(tmp_wind_inds, size = snp_win_size))
  bp_sub_window_list[[j]] <- sub_inds
}

good_sub_snp_inds <- unlist(bp_sub_window_list)

vcf_filt_2 <- vcf_filt[good_sub_snp_inds, ]

geno_mat_sub <- geno_mat_filt[good_sub_snp_inds, ]

sub_ew <- eigen_windows(data = geno_mat_sub, k = 2, win = snp_win_size)
prob_windows <- which(is.na(sub_ew[,1]))



remove_snps <- unlist(bp_window_list[bad_windows])

summary(unlist(lapply(bp_window_list[-bad_windows], length)))

good_snp_inds <- unlist(bp_window_list[good_windows])

geno_mat_filt_2 <- geno_mat_filt[good_snp_inds, ]

bp_window_vec <- nsnps_per_wind[-bad_windows]

tmp_bp_ew <- eigen_windows(data = test_mat, k = test_k, 
  do.windows = bp_window_vec)

get_bp_window = function(geno_mat = geno_mat_filt_2, 
    window_snp_ind_list = bp_window_list, window_ind){
  window_snp_inds <- window_snp_ind_list[[window_ind]]
  out_mat <- geno_mat[window_snp_inds, ]
  return(out_mat)
}

tmp_bp_ew <- eigen_windows(data = get_bp_window, k = test_k,  
  do.windows = good_windows)

