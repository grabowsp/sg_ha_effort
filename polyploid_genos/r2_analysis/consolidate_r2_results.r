# Script for consolidating r^2 results run on sub-vcf files

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

### SET INPUTS ###
data_dir <- args[1]
#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/'

data_dir <- add_slash(data_dir)

res_files_suf <- args[2]
#res_files_suf <- 'geo_subgroup'

### SET OUTPUTS ###

chrom_out_file <- paste(data_dir, 'Chromosome.', res_files_suf, '.r2.rds', 
  sep = '')

combo_out_file <- paste(data_dir, 'combined.', res_files_suf, '.r2.rds', 
  sep = '')

### SET VARIABLES ### 

######

chr_vec <- paste(rep('Chr0', times = 18), rep(seq(9), each = 2),
  rep(c('K', 'N'), times = 9), sep = '')

tot_r2_ls <- list()

for(CHR_NAME in chr_vec){
  print(CHR_NAME)
  ls_com <- paste('ls ', data_dir, CHR_NAME, '*', res_files_suf, '_r2.rds', 
    sep = '')
  res_files <- system(ls_com, intern = T)
  tot_chr_r2 <- readRDS(res_files[1])
  if(length(res_files > 1)){
    for(RF in c(2:length(res_files))){
      tmp_r2 <- readRDS(res_files[RF])
      for(i in seq(length(tot_chr_r2))){
        tot_chr_r2[[i]] <- rbind(tot_chr_r2[[i]], tmp_r2[[i]])
      }
    }
  }
  tot_r2_ls[[CHR_NAME]] <- tot_chr_r2
}

saveRDS(tot_r2_ls, file = chrom_out_file)

combo_r2 <- tot_r2_ls[[1]]

for(j in c(2:length(tot_r2_ls))){
  print(j)
  for(k in seq(length(combo_r2))){
    combo_r2[[k]] <- rbind(combo_r2[[k]], tot_r2_ls[[j]][[k]])
  }
}

saveRDS(combo_r2, file = combo_out_file)

quit(save = 'no')

