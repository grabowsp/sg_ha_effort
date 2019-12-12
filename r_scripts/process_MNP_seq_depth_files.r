# Script for analysing the sequencing depth files generated from MNP counts

### LOAD PACKAGES ###

### SET INPUTS ###
# directory with seq depth count files
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth/'
depth_files <- system(paste('ls ', data_dir, '*tot_seq_depth.txt', sep = ''), 
  intern = T)

# Libraries/samples flagged as BAD in metadata
meta_lib_remove_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/sg_libs_to_remove_from_meta_Dec2019.txt'
meta_remove <- as.vector(read.table(meta_lib_remove_file, header = T,
  stringsAsFactors = F)[,1])

### SET OUTPUT ###
coverage_mat_rds <- '/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth/sample_seq_depth_mat.rds'

### SET VARIABLES ###

### SET CONSTANTS ###

#####################
lib_names <- gsub(data_dir, '', gsub('_tot_seq_depth.txt', '', depth_files))

######
# check that all libraries are represented in the files
chr_libs <- read.table('/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/chr_file_libs.txt', header = F, stringsAsFactors = F)
full_libs <- read.table('/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/full_file_libs.txt', header = F, stringsAsFactors = F)

test_libs <- c(chr_libs[,1], full_libs[,1])
setdiff(test_libs, lib_names)
# character(0)
setdiff(lib_names, test_libs)
# character(0)
## all libraries are present
#########

#########
# Generate matrix containing seq depth/coverage/output for each chromosome
#   for each sample
#tmp_libs <- setdiff(lib_names, rm_libs)
tmp_libs <- lib_names
#tmp_files <- depth_files[-rm_file_inds]
tmp_files <- depth_files

depth_list <- list()

for(i in seq(length(tmp_libs))){
  depth_list[[tmp_libs[i]]] <- read.table(tmp_files[i], header = F,
    stringsAsFactors = F, sep = ' ')
}

#tot_coverage <- unlist(lapply(depth_list, function(x) sum(as.numeric(x[,2]))))

depth_mat <- matrix(
  data = unlist(lapply(depth_list, function(x) as.numeric(x[,2]))),
  nrow = length(depth_list), byrow = T
)
colnames(depth_mat) <- depth_list[[1]][,1]
rownames(depth_mat) <- names(depth_list)

saveRDS(depth_mat, coverage_mat_rds)
###########

###########
# Analyse Depth Results
tot_coverage <- apply(depth_mat, 1, sum)
summary(tot_coverage)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.390e+08 2.253e+10 2.628e+10 2.438e+10 2.846e+10 4.053e+10 

# Look the portion of genome-wide sequencing attributed to each chromosome
port_mat <- t(apply(depth_mat, 1, function(x) x/sum(x)))

port_mean <- apply(port_mat, 2, mean)

port_sd <- apply(port_mat, 2, sd)

# look for outliers - chromosomes that are 4 SD away from the mean
#  of the coverage proportion for that chromosome across all samples
cov_outlier_list <- list()
for(i in seq(ncol(port_mat))){
  cov_outlier_list[[i]] <- which(
    abs(port_mat[,i] - port_mean[i]) - (port_sd[i]*4) > 0)
}

cov_out_tab <- table(unlist(cov_outlier_list))

weird_cov_libs <- names(depth_list)[
  as.numeric(names(cov_out_tab)[cov_out_tab > 2])]
# [1] "IENM" "IIBU" "IICE" "IIDN" "INFT" "IRUH" "IRWF" "IRXF" "TCTS"
# These libraries have 3+ chromosomes that are 4+ SD away from the mean
#  sequencing portion across all the samples
# Bad coverage samples from metadata: IENM, IIBU, IICE, IRUH, IRWF,  
# Estimated 8X in metadata: IIDN
# No metadata info: INFT, TCTS
# Metadata says probable error: IRXF
#########

#########
# Repeat Depth Analysis after removing "Bad" Samples based on metadata

rm_inds <- c()
for(mr in meta_remove){
  tmp_ind <- which(rownames(depth_mat) == mr)
  rm_inds <- c(rm_inds, tmp_ind)
}

depth_mat_2 <- depth_mat[-rm_inds, ]

port_mat_2 <- t(apply(depth_mat_2, 1, function(x) x/sum(x)))
port_mean_2 <- apply(port_mat_2, 2, mean)
port_sd_2 <- apply(port_mat_2, 2, sd)

cov_outlier_list_2 <- list()
for(i in seq(ncol(port_mat_2))){
  cov_outlier_list_2[[i]] <- which(
    abs(port_mat_2[,i] - port_mean_2[i]) - (port_sd_2[i]*4) > 0)
}

cov_out_tab_2 <- table(unlist(cov_outlier_list_2))

weird_cov_libs_2 <- rownames(depth_mat_2)[
  as.numeric(names(cov_out_tab_2)[cov_out_tab_2 > 2])]
#[1] "IIDN" "INFT" "TCTS"
# The same 3 libraries that were flagged when including the "BAD" samples




