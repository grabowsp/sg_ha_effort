# Script for calculating LD via r^2 using VCF data

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

### LOAD PACKAGES ###

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
# out_dir <- 'test/'

# add slash to directory if not already there
out_dir_last_char <- rev(unlist(strsplit(out_dir, split = '')))[1]
if(out_dir_last_char != '/'){
  out_dir <- paste(out_dir, '/', sep = '')
}

out_suffix <- '_r2.rds'

file_pre <- gsub('.vcf', '', rev(unlist(strsplit(vcf_in, split = '/')))[1])

out_file <- paste(out_dir, file_pre, out_suffix, sep = '')


### SET VARIABLES ###
# first column with sample genotypes
FIRST_SAMP = 4

# maximum distance in bp to measure LD
max_dist <- as.numeric(args[4])

### FUNCTIONS (to be saved in another file? ###
calc_snp_r2 <- function(vcf_df, max_dist, pos_column = 'pos', 
  name_column = 'rs', samp1_column = 4){
  #########
  r2_list <- list()
  for(tm in seq(nrow(vcf_df))){
####
#for(tm in seq(1000)){
    tmp_dist <- vcf_df[, pos_column] - vcf_df[tm, pos_column]
    tmp_inds <- which((tmp_dist > 0) & (tmp_dist <= max_dist))
    if(length(tmp_inds) == 0){next}
    tmp_cors <- cor(t(vcf_df[tm,c(samp1_column:ncol(vcf_df))]),
      t(vcf_df[tmp_inds, c(samp1_column:ncol(vcf_df))]),
      use = 'pairwise.complete.obs')
    tmp_r2 <- tmp_cors^2
    test_snp_vec <- rep(vcf_df[tm, name_column], times = length(tmp_r2))
    comp_snp_vec <- vcf_df[tmp_inds, name_column]
    comp_dist <- tmp_dist[tmp_inds]
    r2_list[[tm]] <- list()
    r2_list[[tm]][['test_snp_vec']] <- test_snp_vec
    r2_list[[tm]][['comp_snp_vec']] <- comp_snp_vec
    r2_list[[tm]][['comp_dist']] <- comp_dist
    r2_list[[tm]][['r2']] <- tmp_r2
####
  }
  r2_df <- data.frame(
    test_snp = unlist(lapply(r2_list, function(x) x[['test_snp_vec']])),
    comp_snp = unlist(lapply(r2_list, function(x) x[['comp_snp_vec']])),
    comp_dist = unlist(lapply(r2_list, function(x) x[['comp_dist']])),
    r2 = unlist(lapply(r2_list, function(x) x[['r2']])),
    stringsAsFactors = F)
  return(r2_df)
}

#test <- calc_snp_r2(vcf_df = toy_df, max_dist = 500)

#################

vcf[vcf == './.'] <- NA

geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
oct_dosage_vec <- c('2', '1.5', '1', '0.5', '0')
#tet_dosage_vec <- c('2', '1', '1', '1', '0')

for(i in seq(length(geno_vec))){
  vcf[vcf == geno_vec[i]] <- oct_dosage_vec[i]
}

for(i in c(FIRST_SAMP:ncol(vcf))){
  oct_df[, i] <- as.numeric(oct_df[, i])
}

vcf_r2 <- cal_snp_r2(vcf_df = vcf, max_dist = max_dist)

saveRDS(vcf_r2, out_file)

quit(save = 'no')


