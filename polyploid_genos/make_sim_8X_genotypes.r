#module load python/3.7-anaconda-2019.07
#source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

### LOAD PACKAGES ###
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.basename <- dirname(script.name)

polyploid_function_file <- file.path(script.basename, 'polyploid_functions.r')
#polyploid_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/polyploid_functions.r'
source(polyploid_function_file)

gen_function_file <- file.path(script.basename, 'general_functions.r')
#gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(gen_function_file)

### LOAD INPUTS ###
vcf_in <- args[1]
#vcf_in <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00'
vcf_1 <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

vcf_header_file <- args[2]
#vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt'
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

sim_combo_libs_file <- args[3]
#sim_combo_libs_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_samp_sim8X_lib_combos.txt'

sim_combo_libs <- read.table(sim_combo_libs_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUT ###
out_dir <- args[4]
#out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs'
out_dir <- add_slash(out_dir)

out_pre <- args[5]
#out_pre <- 'geo_samp'
out_pre_long <- paste(out_pre, '_sim8X_AltDosage.vcf', sep = '')

in_pre <- basename(vcf_in)

out_short <- gsub('vcf', out_pre_long, in_pre)

out_file <- paste(out_dir, out_short, sep = '')

### SET VARIABLES ###

##################
tot_keep_libs <- unique(c(sim_combo_libs$keep_1, sim_combo_libs$keep_2))

#####
# These are old step that may need to use if the samples don't completely
#  overlap between the vcf and the chosen sample, for some reason
#miss_libs <- setdiff(tot_keep_libs_0, colnames(vcf_1))
#tot_keep_libs <- intersect(tot_keep_libs_0, colnames(vcf_1))
#scl_rm_inds <- unique(c(which(sim_combo_libs$keep_1 %in% miss_libs), 
#  which(sim_combo_libs$keep_2 %in% miss_libs)))
#sim_combo_libs_1 <- sim_combo_libs[-scl_rm_inds,]
######

pre_4X_dosage_df <- generate_dosage_df(vcf_1, oct_libs = c(), 
  tet_libs = tot_keep_libs, R1 = F)

sim_dosage_mat <- matrix(NA, ncol = nrow(sim_combo_libs), 
  nrow = nrow(pre_4X_dosage_df))
for(i in seq(nrow(sim_combo_libs))){
  tmp_samp_1 <- sim_combo_libs$keep_1[i]
  tmp_samp_2 <- sim_combo_libs$keep_2[i]
  sim_dosage_mat[,i] <- apply(pre_4X_dosage_df[, c(tmp_samp_1, tmp_samp_2)], 
    1, sum)
}

colnames(sim_dosage_mat) <- sim_combo_libs$sim_samp_name

sim_vcf <- data.frame(vcf_1[,c(1:9)], sim_dosage_mat, stringsAsFactors = F)

write.table(sim_vcf, file = out_file, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')

