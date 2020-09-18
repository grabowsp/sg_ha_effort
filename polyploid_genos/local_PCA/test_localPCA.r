# cd /global/homes/g/grabowsp/tools
# module load python/3.7-anaconda-2019.07
# source activate local_PCA

library(data.table)
library(lostruct)

#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/'
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/'

#vcf_in_short <- 'Chr01K.polyploid.CDS.expandgeosamps.vcf_00'
vcf_in_short <- 'Chr01K.polyploid.CDS.expandgeosamps.vcf_00'

vcf_in <- paste(data_dir, vcf_in_short, sep = '')

vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

#head_in_short <- 'CDS.expandgeosamps.vcf.header.txt'
head_in_short <- 'CDS.expandgeosamps.vcf.header.txt'

vcf_header_file <- paste(data_dir, head_in_short, sep = '')

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf) <- vcf_header

tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs <- intersect(tet_libs_0, colnames(vcf))

oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs <- intersect(oct_libs_0, colnames(vcf))

oct_df <- vcf[, oct_libs]
tet_df <- vcf[, tet_libs]

oct_df[oct_df == './.'] <- NA
tet_df[tet_df == './.'] <- NA

geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
oct_dosage_vec <- c('0', '0.5', '1', '1.5', '2')
tet_dosage_vec <- c('0', '1', '1', '1', '2')

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

per_miss_samp <- apply(geno_mat, 2, function(x) sum(is.na(x)) / nrow(geno_mat))

bad_samps <- which(per_miss_samp > 0.05)

geno_mat_1 <- geno_mat[ , -bad_samps]
# genomat_ngenos <- apply(geno_mat, 1, function(x) length(table(x)))

genomat_ngenos <- apply(geno_mat_1, 1, 
  function(x) length(setdiff(unique(x), NA)))

invar_loci <- which(genomat_ngenos == 1)

geno_mat_filt <- geno_mat_1[-invar_loci, ]

vcf_filt <- vcf[-invar_loci, ]

###

# Select SNPs to use with bp-based windows
#  window must have minimum number of SNPs
#  subsample SNPs so each window has same number of SNPs
# - This approach generates the matrix in the correct form to use
#     the standard eigen_windows() function

bp_win_size <- 1e5
snp_win_size <- 1000

min_pos <- min(vcf_filt$POS)
max_pos <- max(vcf_filt$POS)

bp_start_pos <- seq(min_pos, max_pos, by = bp_win_size)
bp_end_pos <- c((bp_start_pos[-1] - 1), max_pos)

# Get snp indieces for each window
bp_window_list <- list()

for(i in seq(length(bp_start_pos))){
  tmp_inds <- which(vcf_filt$POS >= bp_start_pos[i] &
    vcf_filt$POS < bp_end_pos[i])
  bp_window_list[[i]] <- tmp_inds
}

nsnps_per_wind <- unlist(lapply(bp_window_list, length))

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

win_start <- seq(1, nrow(test_mat), win_size)
win_end <- unique(c(seq(win_size, nrow(test_mat), win_size), nrow(test_mat)))

tmp_ew <- eigen_windows(data = test_mat, k = test_k, win = win_size)

prob_windows <- which(is.na(tmp_ew[,1]))

while(length(prob_windows) > 0){
  for(pw in prob_windows){
    print(pw)
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
      tmp_snps_tochange <- sort(sample(tmp_prob_snps, length(tmp_prob_snps)/2))
      test_mat[c(win_start[pw]:win_end[pw])[tmp_snps_tochange], 
        tpi] <- tmp_rM[tmp_snps_tochange]
    }
  }
  tmp_ew <- eigen_windows(data = test_mat, k = 2, win = win_size)
  prob_windows <- which(is.na(tmp_ew[,1]))
  print(prob_windows)
}

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

