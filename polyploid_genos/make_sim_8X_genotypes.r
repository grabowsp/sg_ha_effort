#module load python/3.7-anaconda-2019.07
#source activate R_analysis

vcf_in <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00'
vcf_1 <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt'
vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf_1) <- vcf_header

tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs_1 <- intersect(tet_libs_0, vcf_header)

oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs_1 <- intersect(oct_libs_0, vcf_header)

up_tet_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_1_libs.txt', sep = '')
up_tet_1 <- read.table(up_tet_1_in, header = F, stringsAsFactors = F)[,1]
up_tet_2_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/up_geo_samps/', 'up_tet_2_libs.txt', sep = '')
up_tet_2 <- read.table(up_tet_2_in, header = F, stringsAsFactors = F)[,1]
low_tx_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_1_libs.txt', sep = '')
low_tx_1 <- read.table(low_tx_1_in, header = F, stringsAsFactors = F)[,1]
low_tx_2_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_tx_2_libs.txt', sep = '')
low_tx_2 <- read.table(low_tx_2_in, header = F, stringsAsFactors = F)[,1]
low_gc_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_gc_1_libs.txt', sep = '')
low_gc_1 <- read.table(low_gc_1_in, header = F, stringsAsFactors = F)[,1]
low_ec_1_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_1_libs.txt', sep = '')
low_ec_1 <- read.table(low_ec_1_in, header = F, stringsAsFactors = F)[,1]
low_ec_2_in <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/low_geo_samps/', 'low_ec_2_libs.txt', sep = '')
low_ec_2 <- read.table(low_ec_2_in, header = F, stringsAsFactors = F)[,1]

# need to make list of sample combinations to use for all subfiles

keep_list <- list(up_tet_1, up_tet_2, low_tx_1, low_tx_2, low_gc_1, low_ec_1,
  low_ec_2)

keep_pops <- c('up_tet_1', 'up_tet_2', 'low_tx_1', 'low_tx_2', 'low_gc_1',
  'low_ec_1', 'low_ec_2')

keep_samps <- unlist(keep_list)

min_popsize <- min(unlist(lapply(keep_list, length)))
#33

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
    sim_samp_name = sim_sampnames, stringsAsFactors = F)
}

samepop_8X_list <- list()
for(j in seq(length(keep_list))){
  samepop_8X_list[[j]] <- gen_samepop_sim_combo_libs(lib_vec = keep_list[[j]], 
    n_combos = sim_popsize, sim_popname = paste(keep_pops[j], '_8X', sep = ''))
}

gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(gen_function_file)

gen_unified_df <- function(df_list){
  n_cols <- ncol(df_list[[1]])
  df_colnames <- colnames(df_list[[1]])
  tot_list <- list()
  for(i in seq(n_cols)){
    tot_list[[i]] <- unlist(lapply(df_list, function(x) x[,i]))
  }
  tot_df <- data.frame(tot_list, stringsAsFactors = F)
  colnames(tot_df) <- colnames(df_list[[1]])
  return(tot_df) 
}

samepop_8X_tot_df <- gen_unified_df(df_list = samepop_8X_list)
# NEXT: combinations between populations

combopop_8X_list <- list()
comp_vec <- c(2:length(keep_list))




# Generate 4X dosage genotypes
# Select samples within and between groups to combine into 8X genotypes
# (dosage+dosage) / 2
# don't include SNPs with NA
# Save as rds object with informative names


