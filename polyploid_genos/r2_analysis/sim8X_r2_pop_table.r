# Generate Library and Population Table for simulated 8X Samples
#   Tables can be used for evaluating r^2 patterns


#module load python/3.7-anaconda-2019.07
#source activate R_analysis

### LOAD DATA ###
all_combo_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_samp_sim8X_lib_combos.txt'

all_combo <- read.table(all_combo_file, header = T, sep = '\t', 
  stringsAsFactors = F)

big_combo_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large_geo_samp_sim8X_lib_combos.txt'

big_combo <- read.table(big_combo_file, header = T, sep = '\t', 
  stringsAsFactors = F)

### SET OUTPUTS ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/'

# all sim8X
tot_df_out <- paste(out_dir, 'sim8X_libs_for_r2.txt', sep = '')

# large sim8X pops
large_df_out <- paste(out_dir, 'large_sim8X_libs_for_r2.txt', sep = '')

#################
tot_df <- all_combo[, c('sim_samp_name', 'sim_pop')]
colnames(tot_df) <- c('LIB', 'POP')
# 576 samples
write.table(tot_df, tot_df_out, quote = F, sep = '\t', row.names = F,
  col.names = T)

large_df <- big_combo[, c('sim_samp_name', 'sim_pop')]
colnames(large_df) <- c('LIB', 'POP')
# 630 samples
write.table(large_df, large_df_out, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')



