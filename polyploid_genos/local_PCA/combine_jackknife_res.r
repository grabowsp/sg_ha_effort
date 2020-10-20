# Script to combine the jackknife results from tests using multiple parameters

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
# res_suf <- 'jackknifetest.txt'

### SET OUTPUT ###

# the prefix of the output file
out_pre <- args[3]
#out_pre <- 'CDS.expandv2'

out_file <- paste(data_dir, out_pre, 'jackknife_results_combined.txt', sep = '')

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

tmp_res <- tot_res[order(tot_res$SNP_window), ]
tmp_res_2 <- tmp_res[order(tmp_res$PC), ]

tmp_res_2$sig_noise_dif <- (tmp_res_2$tot_mean_variance - 
  tmp_res_2$noise_variance)

write.table(tmp_res_2, out_file, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')

