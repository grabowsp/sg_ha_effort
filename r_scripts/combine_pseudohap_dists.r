# Script for combining all the pseudohap distances

# LOAD PACKAGES #

# LOAD DATA #

# SET OUTPUTS #
out_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/', 
  'filtered_vcfs/', 'pseudohapdists_v01_total.rds', sep = '')


# SET VARIABLES #
in_file_dir <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs/'
in_file_suffix <- 'pseudohapdists_v01.rds'

#######
ls_com <- paste('ls ', in_file_dir, '*', in_file_suffix, sep = '')

dist_files <- system(ls_com, intern = T)

data_1 <- readRDS(dist_files[1])

tot_nSNPs <- data_1[[1]]

tot_man_dist <- as.matrix(data_1[[2]])

tot_euc_dist <- as.matrix(data_1[[3]])

for(DF in c(2:length(dist_files))){
  tmp_data <- readRDS(dist_files[DF])
  tot_nSNPs <- tot_nSNPs + tmp_data[[1]]
  tmp_man_dist <- as.matrix(tmp_data[[2]])
  tot_man_dist <- tot_man_dist + tmp_man_dist
  tmp_euc_dist <- as.matrix(tmp_data[[3]])
  tot_euc_dist <- tot_euc_dist + tmp_euc_dist
  print(DF)
}

tot_list <- list()

tot_list[['nSNPs']] <- tot_nSNPs
tot_list[['manhattan_dist']] <- tot_man_dist
tot_list[['euclidean_dist']] <- tot_euc_dist

saveRDS(tot_list, file = out_file)

quit(save = 'no')

