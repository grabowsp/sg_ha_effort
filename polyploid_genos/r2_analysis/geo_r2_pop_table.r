# Generate Library and Population Table for "Natural" Samples based
#   on previously-determined PCA results;
#  The tables can be used for evaluating r^2 patterns and generating
#    simulated 8X genotypes

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

### LOAD DATA ###
# Directories for loading library lists

up_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', sep = '')

low_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', sep = '')

# subpop names for inputing data
subpop_names <- c('up_tet_1', 'up_tet_2', 'up_oct_1', 'up_oct_2', 'low_tx_1',
  'low_tx_2', 'low_gc_1', 'low_gc_2', 'low_ec_1', 'low_ec_2')


### SET OUTPUTS ###

out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/geo_samps/', sep = '')

# all the subpopulations
tot_df_out <- paste(out_dir, 'subpop_libs_for_r2.txt', sep = '')

# only 4X subpopulations
tet_df_out <- paste(out_dir, '4X_subpop_libs_for_r2.txt', sep = '')

# only 4X subpops with N > min_size
big_tet_out <- paste(out_dir, 'large4X_subpop_libs_for_r2.txt', sep = '')

### SET VARIABLES ###
# minimum population size to be included in big_tet
min_size <- 25

###################

# Load the 
subset_list <- list()

for(spn in subpop_names){
  if(length(grep('^up_', spn) > 0)){
    tmp_dir <- up_dir
  } else{
    tmp_dir <- low_dir
  }
  tmp_in_file <- paste(tmp_dir, spn, '_libs.txt', sep = '')
  subset_list[[spn]] <- read.table(tmp_in_file, header = F, 
    stringsAsFactors = F)[,1]
}

lib_vec <- unlist(subset_list)
pop_vec <- c()
for(i in seq(length(subset_list))){
  tmp_pop_vec <- rep(names(subset_list)[i], times = length(subset_list[[i]]))
  pop_vec <- c(pop_vec, tmp_pop_vec)
}

tot_df <- data.frame(LIB = lib_vec, POP = pop_vec, stringsAsFactors = F)
# need to remove IEYM because is an 8X in a 4X population
tot_df <- tot_df[-which(tot_df$LIB == 'IEYM'), ]

write.table(tot_df, tot_df_out, quote = F, sep = '\t', row.names = F, 
  col.names = T)

tet_df <- tot_df[-grep('oct', tot_df$POP), ]
write.table(tet_df, tet_df_out, quote = F, sep = '\t', row.names = F, 
  col.names = T)

big_tet_pops <- names(table(tet_df$POP))[which(table(tet_df$POP) >= min_size)]

big_tet_df <- tet_df[which(tet_df$POP %in% big_tet_pops), ]
write.table(big_tet_df, big_tet_out, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')



