# Script for calculating LD via r^2 using VCF data
#  v2 includes ability to calculate r^2 using subsets of samples and
#  simulated polyploid

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

### LOAD PACKAGES ###
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.basename <- dirname(script.name)
parent.dir <- dirname(script.basename)

poly_function_file <- file.path(parent.dir, 'polyploid_functions.r')
#poly_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/polyploid_functions.r'
source(poly_function_file)

general_function_file <- file.path(parent.dir, 'general_functions.r')
#general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)

### LOAD INPUTS ###
vcf_in <- args[1]
#vcf_in <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00'

vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

vcf_header_file <- args[2]
#vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt'

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf) <- vcf_header

# file showing how to group the samples in the VCF 
# if "all", then run all the samples together
# otherwise, this file has sample names in first column and population in
#   second column
samp_pop_in <- args[3]
if(samp_pop_in != 'all'){
  samp_pop_df <- read.table(samp_pop_in, header = T, stringsAsFactors = F,
    sep = '\t')
}

# whether the vcf needs to be converted from number of tetrasomic alleles to
#  dosage ranging from 0 to 2
vcf_type <- args[4]
#vcf_type <- 'allele_count' # or 'dosage'

tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs_1 <- intersect(tet_libs_0, vcf_header)

oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs_1 <- intersect(oct_libs_0, vcf_header)

### SET OUTPUTS ###
out_dir <- args[5]
# out_dir <- 'test/'

out_dir <- add_slash(out_dir)

out_suffix <- '_r2.rds'

#file_pre <- gsub('.vcf', '', rev(unlist(strsplit(vcf_in, split = '/')))[1])
file_pre <- args[6]

out_file <- paste(out_dir, file_pre, out_suffix, sep = '')

### SET VARIABLES ###
# first column with sample genotypes
FIRST_SAMP = 10

# maximum distance in bp to measure LD
max_dist <- as.numeric(args[7])
#max_dist <- 10000

# minimum minor allele frequency
maf_cut <- as.numeric(args[8])
#maf_cut <- 0.05

### FUNCTIONS (to be saved in another file? ###
calc_snp_r2 <- function(vcf_df, max_dist, pos_column = 'POS', 
  name_column = 'ID', samp1_column = FIRST_SAMP){
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

#vcf[vcf == './.'] <- NA

#geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
#oct_dosage_vec <- c('2', '1.5', '1', '0.5', '0')
#tet_dosage_vec <- c('2', '1', '1', '1', '0')

#for(i in seq(length(geno_vec))){
#  vcf[vcf == geno_vec[i]] <- oct_dosage_vec[i]
#}

#for(i in c(FIRST_SAMP:ncol(vcf))){
#  vcf[, i] <- as.numeric(vcf[, i])
#}

if(samp_pop_in == 'all'){
  samp_vec <- colnames(vcf)[c(FIRST_SAMP:ncol(vcf))]
  samp_pop <- rep(pop1, times = length(samp_vec))
  samp_pop_df <- data.frame(samp_vec, samp_pop, stringsAsFactors = F)
}

unique_pop_vec <- unique(samp_pop_df[,2])

r2_list <- list()

for(i in unique_pop_vec){
  print(i)
  tmp_inds <- which(samp_pop_df[,2] == i)
  tmp_samps <- samp_pop_df[tmp_inds,1]
  if(vcf_type == 'allele_count'){
    tmp_oct_libs <- intersect(oct_libs_1, tmp_samps)
    tmp_tet_libs <- intersect(tet_libs_1, tmp_samps)
    dosage_vcf <- cbind(vcf[ , c(1:(FIRST_SAMP-1))],
      generate_dosage_df(vcf = vcf[, tmp_samps], oct_libs = tmp_oct_libs, 
      tet_libs = tmp_tet_libs, R1 = F))
  } else{
    dosage_vcf <- cbind(vcf[ , c(1:(FIRST_SAMP-1))], vcf[, tmp_samps])
  }
  # Filter VCF based on maf
  n_nas <- apply( dosage_vcf[, c(FIRST_SAMP:ncol(dosage_vcf))], 1,
    function(x) sum(is.na(unlist(x))))
  n_gsamps <- length(c(FIRST_SAMP:ncol(dosage_vcf))) - n_nas
  sum_ref <- apply( dosage_vcf[ , c(FIRST_SAMP:ncol(dosage_vcf))], 1,
    function(x) sum(unlist(x), na.rm = T))
  #
  allele_freq <- sum_ref / (2*n_gsamps)
  allele_freq_1 <- 1-allele_freq
  #
  minor_af <- apply(cbind(allele_freq, allele_freq_1), 1, 
    function(x) min(x)[1])
  minor_af[which(is.na(minor_af))] <- 0
  #
  if(sum(minor_af < maf_cut) > 0){
    vcf_1 <- dosage_vcf[-which(minor_af < maf_cut), ]
  }else{vcf_1 <- dosage_vcf}
  # Filter out any non-variant SNPs, if they exist
  n_genos <- apply(vcf_1[, c(FIRST_SAMP:ncol(vcf_1))], 1, 
    function(x) length(table(unlist(x))))
  if(sum(n_genos==1) > 0) {
    vcf_2 <- vcf_1[-which(n_genos == 1), ]
  }else{vcf_2 <- vcf_1}
  vcf_r2 <- calc_snp_r2(vcf_df = vcf_2, max_dist = max_dist)
  r2_list[[i]] <- vcf_r2
}

saveRDS(r2_list, out_file)

quit(save = 'no')


