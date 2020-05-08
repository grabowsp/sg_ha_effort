# Script for combining all the NA-based distance matrices and nSNP info

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

# LOAD PACKAGES #

# LOAD DATA #
data_dir <- args[1]
#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/no_filtering'

data_dir_last_char <- rev(unlist(strsplit(data_dir, split = '')))[1]
if(data_dir_last_char != '/'){
  data_dir <- paste(data_dir, '/', sep = '')
}

in_file_suffix <- args[2]
#in_file_suffix <- 'NA_DistMat.rds'

# SET OUTPUTS #
out_file_short <- args[3]
#out_file_short <- 'Chr01K.polyploid.CDS.allsamps.NA_DistMat.total.rds'

out_file <- paste(data_dir, out_file_short, sep = '')

# SET VARIABLES #
#######
ls_com <- paste('ls ', data_dir, '*', in_file_suffix, sep = '')

dist_files <- system(ls_com, intern = T)

data_1 <- readRDS(dist_files[1])

tot_nSNPs <- data_1[[1]]

tot_euc_dist <- as.matrix(data_1[[2]])

for(DF in c(2:length(dist_files))){
  tmp_data <- readRDS(dist_files[DF])
  tot_nSNPs <- tot_nSNPs + tmp_data[[1]]
  tmp_euc_dist <- as.matrix(tmp_data[[2]])
  tot_euc_dist <- tot_euc_dist + tmp_euc_dist
  print(DF)
}

tot_list <- list()

tot_list[['nSNPs']] <- tot_nSNPs
tot_list[['euclidean_dist']] <- tot_euc_dist

saveRDS(tot_list, file = out_file)

quit(save = 'no')

