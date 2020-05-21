# Script for generating a distance matrix from a tetrasomic VCF using disomic 
#   genotype dosages for all samples
#
# # Distance matrices are based on disomic dosage genotypes:
#  RRRR (homozygous REF) = 2.0; AAAA (homozygous ALT) = 0.0
#  RRRA = RRAA = RAAA = 1.0
# 
###############
# Arguments
# [1]: (character) filename of the VCF; can be 'split' file because header is
#			not loaded from this file
# [2]: (character) filename of the header for the VCF;
# [3]: (character) directory where distance matrix will be saved
###############

###############
# OUTPUT
# List with 4 elements:
# [['nSNPS']] = the number of SNPs used to calculate the distance
# [['euclidean_dist']] = euclidean distance matrix based on dosage genotypes
# [['manhattan_dist']] = manhattan distance matrix based on dosage genotypes
# [['n_NAs']] = the number of NA genotypes in each sample; NA's are removed
#                       from distance calculations and normalized to nSNPs, so
#                       these values can help determine if/how missing data
#                       may be affecting results
###############

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

### LOAD INPUTS ###
vcf_in <- args[1]
#vcf_in <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/Chr01K.polyploid_100k.vcf'

vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

vcf_header_file <- args[2]
#vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/full_vcf_header.txt'

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

### SET OUTPUTS ###
out_dir <- args[3]
#out_dir <- 'test/'

# add slash to directory if not already there
out_dir_last_char <- rev(unlist(strsplit(out_dir, split = '')))[1]
if(out_dir_last_char != '/'){
  out_dir <- paste(out_dir, '/', sep = '')
}

out_suffix <- '_disomic_DistMat.rds'

file_pre <- gsub('.vcf', '', rev(unlist(strsplit(vcf_in, split = '/')))[1])

out_file <- paste(out_dir, file_pre, out_suffix, sep = '')

### SET VARIABLES ###


### LOAD PACKAGES ###


##################
#geno_df <- vcf[, c(10:ncol(vcf))]

oct_df <- vcf[, oct_libs]
tet_df <- vcf[, tet_libs]

oct_df[oct_df == './.'] <- NA
tet_df[tet_df == './.'] <- NA

geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
#oct_dosage_vec <- c('2', '1.5', '1', '0.5', '0')
tet_dosage_vec <- c('2', '1', '1', '1', '0')

for(i in seq(length(geno_vec))){
  oct_df[oct_df == geno_vec[i]] <- tet_dosage_vec[i]
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

geno_mat <- matrix(unlist(geno_df), ncol = nrow(geno_df), 
  byrow = T)

rownames(geno_mat) <- colnames(geno_df)

dist_euc <- dist(geno_mat, diag = T, upper = T, method = 'euclidean')
dist_man <- dist(geno_mat, diag = T, upper = T, method = 'manhattan')

nSNPs <- nrow(vcf)

n_NAs <- apply(geno_mat, 1, function(x) sum(is.na(x)))

tot_list <- list()

tot_list[['nSNPs']] <- nSNPs
tot_list[['euclidean_dist']] <- dist_euc
tot_list[['manhattan_dist']] <- dist_man
tot_list[['n_NAs']] <- n_NAs

saveRDS(tot_list, file = out_file)

quit(save = 'no')

