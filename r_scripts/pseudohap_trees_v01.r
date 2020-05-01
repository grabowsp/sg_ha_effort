# Script for generating trees from the pseudohaploid distances

# source activate r_phylo

# LOAD PACKAGES #
library(phangorn)
library(Biostrings)
library(ggplot2)
library(devtools)
install_github('YuLab-SMU/ggtree')

library(ggtree, lib.loc = '/home/grabowsky/tools/r_tools/ggtree/appveyor.yml')
### LOAD DATA ###
data_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/',
  'filtered_vcfs/', 'pseudohapdists_v01_total.rds', sep = '')
data <- readRDS(data_file)


# test_genofile <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs/Chr09N_filt_split_07_pseudohapgenos_v01.rds'

# test_genos <- readRDS(test_genofile)

### SET OUTPUTS ###
man_nj_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_manhattan_nj.tre', 
  sep = '')
man_upgma_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_manhattan_upgma.tre', 
  sep = '')

euc_nj_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_euclidean_nj.tre',
  sep = '')
euc_upgma_tree_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
  'pseudohap/', 'filtered_vcfs/', 'pseudohapdists_v01_euclidean_upgma.tre',
  sep = '')

# SET VARIABLES #

options(gzip = '/bin/gzip')
Sys.setenv(TAR = '/bin/tar')

############


geno_mat <- matrix(unlist(test_genos[c(1:1000), c(6:ncol(test_genos))]),
  nrow = 1000, byrow = F)
colnames(geno_mat) <- colnames(test_genos)[c(6:ncol(test_genos))]

man_dist <- as.dist(data[[2]])

man_upgma <- upgma(man_dist)

man_nj <- NJ(man_dist)

write.tree(man_upgma, file = man_upgma_tree_file)
write.tree(man_nj, file = man_nj_tree_file)

###

euc_dist <- as.dist(data[[3]])

euc_upgma <- upgma(euc_dist)
euc_nj <- NJ(euc_dist)

write.tree(euc_upgma, file = euc_upgma_tree_file)
write.tree(euc_nj, file = euc_nj_tree_file)

##############



quit(save = 'no')

