# Functions for localPCA analysis

process_vcf <- function(vcf, oct_libs, tet_libs, rm_libs){
  # Processes VCF to generate matrix with dosage genotypes appropriate
  #   for each ploidy level
  # INPUTS
  # vcf = vcf data.frame file with header
  # oct_libs = vector of names of octoploid libraries
  # tet_libs = vector of names of tetraploid libraries
  # rm_libs = vector of names of libraries to remove
  # OUTPUT
  # matrix of ALT dosage genotypes ranging from 0 to 2; 
  #  rows = SNPS, columns = samples; 
  #######################
  oct_df <- vcf[, setdiff(oct_libs, rm_libs)]
  tet_df <- vcf[, setdiff(tet_libs, rm_libs)]
  #
  oct_df[oct_df == './.'] <- NA
  tet_df[tet_df == './.'] <- NA
  #
  geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
  oct_dosage_vec <- c('0', '0.5', '1', '1.5', '2')
  tet_dosage_vec <- c('0', '1', '1', '1', '2')
  #
  for(i in seq(length(geno_vec))){
    oct_df[oct_df == geno_vec[i]] <- oct_dosage_vec[i]
  }
  for(i in seq(length(geno_vec))){
    tet_df[tet_df == geno_vec[i]] <- tet_dosage_vec[i]
  }
  for(i in seq(ncol(oct_df))){
    oct_df[, i] <- as.numeric(oct_df[, i])
  }
  for(i in seq(ncol(tet_df))){
    tet_df[, i] <- as.numeric(tet_df[, i])
  }
  geno_df <- cbind(oct_df, tet_df)
  geno_mat <- matrix(unlist(geno_df), ncol = ncol(geno_df),
    byrow = F)
  rownames(geno_mat) <- seq(nrow(geno_mat))
  colnames(geno_mat) <- colnames(geno_df)
  return(geno_mat)
}

###

get_window_inds <- function(vcf, min_pos, max_pos, bp_win_size){
  # Get the indices for each SNP in a window
  # INPUTS
  # vcf = vcf data.frame with header;
  # min_pos = lowest position; where to start with making bp-based windows
  # max_pos = highest position; where to stop making windows
  # bp_win_size = the bp window size; ex: 5e5
  # OUTPUT
  # List with two elements
  # 'window_start_pos' = vector of the the starting position for each window
  # 'bp_window_list' <- list with each element containing the SNP indices
  #                       found in that corresponding window. For windows
  #                       that contain no SNPs, the element will be empty
  ################
  bp_start_pos <- seq(min_pos, max_pos, by = bp_win_size)
  bp_end_pos <- c((bp_start_pos[-1] - 1), max_pos)
  # Get snp indieces for each window
  bp_window_list <- list()
  #
  for(i in seq(length(bp_start_pos))){
    tmp_inds <- which(vcf$POS >= bp_start_pos[i] &
      vcf$POS < bp_end_pos[i])
    bp_window_list[[i]] <- tmp_inds
  }
  last_start_pos <- rev(bp_start_pos)[1]
  out_list <- list()
  out_list[['window_start_pos']] <- bp_start_pos
  out_list[['bp_window_list']] <- bp_window_list
  return(out_list)
}


###

get_subsamp_window_inds <- function(bp_window_list, snp_win_size){
  # Subsample the SNP indices within windows so that all windows have
  #   the same number of SNPs
  # INPUTS
  # bp_window_list = list of SNP indices in each window; is out_list[[2]]
  #                   generated by get_window_inds()
  # snp_win_size = the number of SNPs to include in each window; ex: 1000
  # OUTPUT
  # vector of indices for the snps that should be used for the analysis
  #   these indices are used to select from the VCF and geno_mat objects
  ################
  nsnps_per_wind <- unlist(lapply(bp_window_list, length))
  good_windows <- which(nsnps_per_wind >= snp_win_size)
  #
  bp_sub_window_list <- list()
  for(j in seq(length(good_windows))){
    tmp_wind_inds <- bp_window_list[[good_windows[j]]]
    sub_inds <- sort(sample(tmp_wind_inds, size = snp_win_size))
    bp_sub_window_list[[j]] <- sub_inds
  }
  #
  good_sub_snp_inds <- unlist(bp_sub_window_list)
  return(good_sub_snp_inds)
}

#############

eigen_window_fixNAs <- function(test_mat, win_size, test_k, fna.verbose = F){
  # run eigen_windows() and replace NA's in the genotype matrix if they
  #  are causing NA's in the output from eigenwindows
  # INPUTS
  # test_mat = genotype matrix used for eigen_windows()
  # win_size = number of SNPs in windows
  # test_k = k to test for PCA
  # fna.verbose = if T, then prints the problematic windows and windows that 
  #                  are getting fixed at the moment
  # OUTPUT
  # output from eigen_windows() with no NAs
  ############
  win_start <- seq(1, nrow(test_mat), win_size)
  win_end <- unique(c(seq(win_size, nrow(test_mat), win_size), nrow(test_mat)))
  tmp_ew <- eigen_windows(data = test_mat, k = test_k, win = win_size)
  prob_windows <- which(is.na(tmp_ew[,1]))
  if(fna.verbose){print(prob_windows)}
  if(length(prob_windows > 0)){
    while(length(prob_windows) > 0){
      for(pw in prob_windows){
        if(fna.verbose){print(pw)}
        test_cov <- cov(sweep(test_mat[c(win_start[pw]:win_end[pw]), ], 1,
          rowMeans(test_mat[c(win_start[pw]:win_end[pw]), ],
          na.rm = T), "-"), use = "pairwise")
        test_cov_nas <- which(is.na(test_cov), arr.ind = T)
        test_cov_prob_inds <- c(test_cov_nas[,1], test_cov_nas[,2])
        test_cov_prob_tab <- table(test_cov_prob_inds)
        test_prob_inds <- as.numeric(names(
          which(test_cov_prob_tab >= mean(test_cov_prob_tab))))
        #
        tmp_rM <- rowMeans(test_mat[c(win_start[pw]:win_end[pw]), ], na.rm = T)
        for(tpi in test_prob_inds){
          tmp_prob_snps <- which(is.na(
            test_mat[c(win_start[pw]:win_end[pw]), tpi]))
          tmp_snps_tochange <- sort(sample(tmp_prob_snps,
            length(tmp_prob_snps)/2))
          test_mat[c(win_start[pw]:win_end[pw])[tmp_snps_tochange],
            tpi] <- tmp_rM[tmp_snps_tochange]
        }
      }
      tmp_ew <- eigen_windows(data = test_mat, k = 2, win = win_size)
      prob_windows <- which(is.na(tmp_ew[,1]))
      if(fna.verbose){print(prob_windows)}
    }
  }
  return(tmp_ew)
}

#######

get_jackknife_window_inds <- function(bp_window_list, snp_win_size, 
  n_part = 10){
  # Subsample the SNP indices within windows so that all windows have
  #   the same number of SNPs
  # INPUTS
  # bp_window_list = list of SNP indices in each window; is out_list[[2]]
  #                   generated by get_window_inds()
  # snp_win_size = the number of SNPs to that would be used in each window; 
  #                 ex: 1000
  # n_part = the number of partitions for the jackknife approach; each window
  #            will be split into n_part parts, with one part removed and
  #            remaining indices retained
  # OUTPUT
  # vector of indices for the snps that should be used for the analysis
  #   these indices are used to select from the VCF and geno_mat objects
  ################
  nsnps_per_wind <- unlist(lapply(bp_window_list, length))
  good_windows <- which(nsnps_per_wind >= snp_win_size)
  #
  n_split_snps <- floor(snp_win_size/n_part)
  remain_snps <- snp_win_size - n_split_snps
  bp_sub_window_list <- list()
  for(j in seq(length(good_windows))){
    tmp_wind_inds <- bp_window_list[[good_windows[j]]]
    # select all subsample indices to use
    sub_inds <- sort(sample(tmp_wind_inds, size = snp_win_size))
    # remove partitions of indices for jackknive approach
    subpart_list <- list()
    for(spl in seq(n_part)){
      tmp_ind1 <- 1+((spl-1)*n_split_snps)
      tmp_ind2 <- (spl*n_split_snps)
      tmp_part_inds <- c(tmp_ind1:tmp_ind2)
      tmp_keep_inds <- sub_inds[-tmp_part_inds]
      subpart_list[[spl]] <- tmp_keep_inds
    }
    bp_sub_window_list[[j]] <- subpart_list
  }
  #
  jack_snp_ind_list <- list()
  for(bsw in seq(n_part)){
    tmp_jack_inds <- unlist(lapply(bp_sub_window_list, function(x) x[[bsw]]))
    jack_snp_ind_list[[bsw]] <- tmp_jack_inds
  }
  return(jack_snp_ind_list)
}

# NEED TO TEST THIS



