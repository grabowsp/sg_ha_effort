# Script for generating window values for mean r^2

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

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
r2_data_in <- args[1]
#r2_data_in <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/Chromosome.geo_subgroup.r2.rds'

data_type <- args[2]
#data_type <- 'chromosomes' # or 'combo'

### SET VARIABLES ###
window_size <- as.numeric(args[3])
#window_size <- 10

hi_r2 <- as.numeric(args[4])
#hi_r2 <- 0.9 

### SET OUTPUTS ###
out_file_pre <- sub('r2.rds', '', r2_data_in, fixed = T)

mean_r2window_out <- paste(out_file_pre, 'mean_r2.', window_size, 
  'bpwindow.rds', sep = '')

########

if(data_type != 'chromosomes' & data_type != 'combo'){
  print('Need to indicate correct data_type')
  quit(save = 'no')
}

r2_res <- readRDS(r2_data_in)

if(data_type == 'combo'){
  perHighr2_windows <- lapply(r2_res, function(x)
    generate_perHighVal_window_df(dist_vec = x$comp_dist, value_vec = x$r2,
      window_size = window_size, max_val = hi_r2))
}

if(data_type == 'chromosomes'){
  perHighr2_windows <- lapply(r2_res, function(y)
    lapply(y, function(x)
      generate_perHighVal_window_df(dist_vec = x$comp_dist, value_vec = x$r2,
        window_size = window_size, max_val = hi_r2)
    )
  )
}

saveRDS(perHighr2_windows, file = mean_r2window_out)

quit(save = 'no')

