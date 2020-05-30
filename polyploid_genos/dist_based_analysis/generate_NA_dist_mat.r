# Script for generating an NA distance matrix from a VCF

###############
# Arguments
# [1]: (character) filename of the VCF; can be 'split' file because header is
#			not loaded from this file
# [2]: (character) filename of the header for the VCF;
# [3]: (character) directory where distance matrix will be saved
###############

###############
# OUTPUT
# List with 2 elements:
# [['nSNPS']] = the number of SNPs used to calculate the distance
# [['euclidean_dist']] = euclidean distance matrix based on NA's of the SNPs
#		NA's are given 0, present genotypes are given 1
#		Distance matrix is generated based on the matrix of 0's and 1's
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

### SET OUTPUTS ###
out_dir <- args[3]
#out_dir <- 'test/'

# add slash to directory if not already there
out_dir_last_char <- rev(unlist(strsplit(out_dir, split = '')))[1]
if(out_dir_last_char != '/'){
  out_dir <- paste(out_dir, '/', sep = '')
}

out_suffix <- '_NA_DistMat.rds'

file_pre <- gsub('.vcf', '', rev(unlist(strsplit(vcf_in, split = '/')))[1])

out_file <- paste(out_dir, file_pre, out_suffix, sep = '')

### SET VARIABLES ###


### LOAD PACKAGES ###


##################
geno_df <- vcf[, c(10:ncol(vcf))]

geno_df[geno_df == './.'] <- 0
geno_df[geno_df != '0'] <- 1

for(i in seq(ncol(geno_df))){
  geno_df[, i] <- as.numeric(geno_df[, i])
}

geno_mat <- matrix(unlist(geno_df), ncol = nrow(geno_df), 
  byrow = T)

rownames(geno_mat) <- colnames(geno_df)

dist_euc <- dist(geno_mat, diag = T, upper = T, method = 'euclidean')

nSNPs <- nrow(vcf)

tot_list <- list()

tot_list[['nSNPs']] <- nSNPs
tot_list[['euclidean_dist']] <- dist_euc

saveRDS(tot_list, file = out_file)

quit(save = 'no')

