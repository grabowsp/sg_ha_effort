# Script to combine the pick_window results from across chromosomes

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

args = commandArgs(trailingOnly=T)

### INPUTS ###

data_dir <- args[1]
# data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/'
if(rev(unlist(strsplit(data_dir, split = '')))[1] != '/'){
  data_dir <- paste(data_dir, '/', sep = '')
}

# the suffix of the individual result files to that only desired files
#  are chosen. Can be simply 'windowtests.txt' if directory only contains one
#  set of tests
res_suf <- args[2]
# res_suf <- 'polyploid.CDS.expandv2.windowtests.txt'

### SET OUTPUT ###

# the prefix of the output file
out_pre <- args[3]
#out_pre <- 'CDS.expandv2'

out_file <- paste(data_dir, out_pre, 'tot_best_bp_window_test.txt', sep = '')

###############

sys_com <- paste('ls ', data_dir, '*', res_suf, sep = '')

in_files <- system(sys_com, intern = T)

for(infi in in_files){
  res <- read.table(infi, header = T, sep = '\t', stringsAsFactors = F)
#  chr_name <- unlist(strsplit(sub(data_dir, '', infi), split = '.', 
#    fixed = T))[1]
#  res$chr <- chr_name
  if(infi == in_files[1]){
    tot_res <- res
  } else{
    tot_res <- rbind(tot_res, res)
  }
}

# median across chromosomes of median within-chrom best windows
tot_med_windows <- tapply(tot_res$median_best_window, tot_res$snp_window, 
  median)

# approximation of number of windows for best bp windows
## to get exact value, need to run `get_window_info.r` script
tot_n_windows <- tapply(tot_res$n_windows, tot_res$snp_window, sum)

# approximation of percent of good windows
## to get exact value, need to run `get_window_info.r` script
tot_per_windows <- tapply(tot_res$median_percent_good_window, 
  tot_res$snp_window, mean)

tot_wind_df <- data.frame(SNP_window = names(tot_med_windows),
  best_bp_window = tot_med_windows, stringsAsFactors = F)

tot_wind_df$n_good_windows <- tot_n_windows[rownames(tot_wind_df)]

tot_wind_df$per_good_windows <- tot_per_windows[rownames(tot_wind_df)]

tot_wind_df <- tot_wind_df[ order(as.numeric(
  sub('_SNPs', '', tot_wind_df$SNP_window))), ]

write.table(tot_wind_df, out_file, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')

