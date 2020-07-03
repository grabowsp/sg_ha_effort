# Script to generate combonation of 4X samples that will be used to make
#  simulated 8X genotypes

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

### LOAD PACKAGES ###
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.basename <- dirname(script.name)

general_function_file <- file.path(script.basename, 'general_functions.r')
#gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'

source(gen_function_file)

### INPUT FILES ###

sub_pop_file <- args[1]
#sub_pop_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt'

sub_pop_lib_df <- read.table(sub_pop_file, sep = '\t', header = T, 
  stringsAsFactors = F)

#up_tet_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_1_libs.txt', sep = '')
#up_tet_1 <- read.table(up_tet_1_in, header = F, stringsAsFactors = F)[,1]
#up_tet_2_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_2_libs.txt', sep = '')
#up_tet_2 <- read.table(up_tet_2_in, header = F, stringsAsFactors = F)[,1]
#low_tx_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_1_libs.txt', sep = '')
#low_tx_1 <- read.table(low_tx_1_in, header = F, stringsAsFactors = F)[,1]
#low_tx_2_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_2_libs.txt', sep = '')
#low_tx_2 <- read.table(low_tx_2_in, header = F, stringsAsFactors = F)[,1]
#low_gc_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_gc_1_libs.txt', sep = '')
#low_gc_1 <- read.table(low_gc_1_in, header = F, stringsAsFactors = F)[,1]
#low_ec_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_1_libs.txt', sep = '')
#low_ec_1 <- read.table(low_ec_1_in, header = F, stringsAsFactors = F)[,1]
#low_ec_2_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_2_libs.txt', sep = '')
#low_ec_2 <- read.table(low_ec_2_in, header = F, stringsAsFactors = F)[,1]

### SET OUTPUT ###

out_dir <- dirname(sub_pop_file)

combo_out_pre <- args[2]
#combo_out_pre <- 'geo_samp'
out_short <- paste(combo_out_pre, '_sim8X_lib_combos.txt')

combo_out_file <- file.path(out_dir, out_short)

#combo_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
#  'polyploid_vcfs/CDS_vcfs/', 'geo_samps/', 
#  'sim_8X_library_combinations.txt', sep = '')

### SET VARIABLES ###


###################################

# need to make list of sample combinations to use for all subfiles

#keep_list <- list(up_tet_1, up_tet_2, low_tx_1, low_tx_2, low_gc_1, low_ec_1,
#  low_ec_2)

#keep_pops <- c('up_tet_1', 'up_tet_2', 'low_tx_1', 'low_tx_2', 'low_gc_1',
#  'low_ec_1', 'low_ec_2')

#keep_samps <- unlist(keep_list)

#min_popsize <- min(unlist(lapply(keep_list, length)))
#33

keep_pops <- unique(sub_pop_lib_df$POP)

min_popsize <- min(tapply(sub_pop_lib_df$LIB, sub_pop_lib_df$POP, 
  function(x) length(x)))

sim_popsize <- min_popsize - 1

gen_samepop_sim_combo_libs <- function(lib_vec, n_combos, sim_popname){
  keep_1 <- c()
  keep_2 <- c()
  option_2 <- lib_vec
  for(i in sample(lib_vec, size = n_combos)){
    not_avail <- c()
    if(i %in% keep_2){
      tmp_ind <- which(keep_2 == i)
      not_avail <- keep_1[tmp_ind]
    }
    keep_1 <- c(keep_1, i)
    tmp_keep_2 <- sample(setdiff(option_2, c(i, not_avail)), size = 1)
    keep_2 <- c(keep_2, tmp_keep_2)
    option_2 <- setdiff(option_2, tmp_keep_2)
  }
  sim_sampnames <- paste(sim_popname, seq(n_combos), sep = '_')
  keep_df <- data.frame(keep_1, keep_2, sim_pop = sim_popname, 
    sim_samp_name = sim_sampnames, combo_type = 'intra_pop', 
    stringsAsFactors = F)
}

samepop_8X_list <- list()
for(j in keep_pops){
  tmp_libs <- sub_pop_lib_df$LIB[which(sub_pop_lib_df$POP == j)]
  samepop_8X_list[[j]] <- gen_samepop_sim_combo_libs(lib_vec = tmp_libs, 
    n_combos = sim_popsize, sim_popname = paste(j, '_8X', sep = ''))
}

samepop_8X_tot_df <- gen_unified_df(df_list = samepop_8X_list)

combopop_8X_list <- list()
comp_vec <- c(2:length(keep_pops))
for(j in c(1:(length(keep_pops)-1))){
  for(i in c((j+1):length(keep_pops))){
    tmp_libs_1 <- sub_pop_lib_df$LIB[which(sub_pop_lib_df$POP == keep_pops[j])]
    tmp_libs_2 <- sub_pop_lib_df$LIB[which(sub_pop_lib_df$POP == keep_pops[i])]
    keep_1 <- sample(tmp_libs_1, size = sim_popsize)
    keep_2 <- sample(tmp_libs_2, size = sim_popsize)
    sim_popname <- paste(keep_pops[j], keep_pops[i], '8X', sep = '_')
    sim_sampnames <- paste(sim_popname, seq(sim_popsize), sep = '_')
    keep_df <- data.frame(keep_1, keep_2, sim_pop = sim_popname, 
      sim_samp_name = sim_sampnames, combo_type = 'inter_pop', 
      stringsAsFactors = F)
    combopop_8X_list[[sim_popname]] <- keep_df
  }
}

combopop_8X_tot_df <- gen_unified_df(combopop_8X_list)

tot_8X_combo_df <- rbind(samepop_8X_tot_df, combopop_8X_tot_df)

write.table(tot_8X_combo_df, file = combo_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

quit(save = 'no')

