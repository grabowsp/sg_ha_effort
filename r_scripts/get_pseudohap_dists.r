# Generate allele ratio files from VCF

# LOAD PACKAGES #

# INPUTS #
args = commandArgs(trailingOnly=TRUE)

genos_in <- as.character(args[1])
genos <- readRDS(genos_in)
 
# OUTPUTS #
out_file <- paste(gsub('_pseudohapgenos_v01.rds', '', genos_in), 
  '_pseudohapdists_v01.rds', sep = '')

# SET VARIABLES #

################
full_mat <- matrix(unlist(genos[, c(6:ncol(genos))]), ncol = nrow(genos), 
  byrow = T)

rownames(full_mat) <- colnames(genos)[c(6:ncol(genos))]

full_dist_man <- dist(full_mat, diag = T, upper = T, method = 'manhattan')
full_dist_euc <- dist(full_mat, diag = T, upper = T, method = 'euclidean')

nSNPs <- nrow(genos)

tot_list <- list()

tot_list[['nSNPs']] <- nSNPs
tot_list[['manhattan_dist']] <- full_dist_man
tot_list[['euclidean_dist']] <- full_dist_euc

saveRDS(tot_list, file = out_file)

quit(save = 'no')


