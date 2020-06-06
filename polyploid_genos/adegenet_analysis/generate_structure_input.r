# Script to generate STRUCTURE input using adegenet input

# Steps
# * 2 lines for 4X, 4 lines for 8X
# * 1 column for each SNP
# top row is names of SNPs
# Col 1 = sample name
# For now, don't worry about PopData or PopFlag which could be cols 2 and 3
# For now, don't use LocData, Phenotype, or Extra Columns
# Col 2 = start with genotypes
# Use 1 = REF, 2 = ALT, -9 = NA

#module load python/3.7-anaconda-2019.07
#source activate r_adegenet_env

### LOAD PACKAGES ###
library(adegenet)
library(parallel)

# input arguments
input_args <- commandArgs(trailingOnly = T)

# get the r-script directory to load funciton files
rundir_args <- commandArgs(trailingOnly = F)

# get system path for helper files
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '', 
  rundir_args[grep(file.arg.name, rundir_args)])
script.basename <- dirname(script.name)
gen.basename <- sub('/adegenet_analysis', '', script.basename)

gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
#gen_function_file <- file.path(gen.basename, 'general_functions.r')
source(gen_function_file)

adegenet_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/adegenet_functions.r'
#adegenet_function_file <- file.path(script.basename, 'adegenet_functions.r')
source(adegenet_function_file)

### INPUT DATA ###

geno_in_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.rds'
#geno_in_file <- input_args[1]

gen_tot <- readRDS(geno_in_file)

### SET OUTPUTS ###
#out_short <- input_args[2]
out_short <- 'test_STRUC_genos.txt'

out_file <- file.path(dirname(geno_in_file), out_short)

### SET VARIABLES ###

n_snps <- 100
#n_snps <- as.numeric(input_args[2])

geno_type <- 'all_tet'
#geno_type <- input_args[3]

###############

keep_inds <- sort(sample(seq(nLoc(gen_tot)), size = n_snps))

keep_genos <- as.matrix(gen_tot[,keep_inds])

keep_snp_names <- colnames(keep_genos)

ploidy_vec <- gen_tot$ploidy

struc_list <- list()

for(i in seq(nrow(keep_genos))){
  struc_list[[i]] <- gen_struc_genos(genotypes = keep_genos[i, ], 
    lib_name = rownames(keep_genos)[i], ploidy = ploidy_vec[i], 
    geno_type = geno_type)
}

if(geno_type == 'all_dip'){
  struc_mat <- matrix(data = NA, nrow = (length(struc_list) * 2), 
    ncol = (ncol(struc_list[[1]])-1))
}else{
  struc_mat <- matrix(data = NA, nrow = (length(struc_list) * 4), 
    ncol = (ncol(struc_list[[1]])-1))
}

tmp_ind <- 1

for(j in seq(length(struc_list))){
  mat_inds <- c(tmp_ind:((tmp_ind-1) + nrow(struc_list[[j]])))
  struc_mat[mat_inds, ] <- matrix(
    data = unlist(struc_list[[j]][, -1]), nrow = nrow(struc_list[[j]]))
  tmp_ind <- tmp_ind + nrow(struc_list[[j]])
}

lib_vec <- unlist(lapply(struc_list, function(x) x[,1]))

rownames(struc_mat) <- lib_vec
colnames(struc_mat) <- keep_snp_names

write.table(struc_mat, file = out_file, quote = F, sep = '\t', row.names = T, 
  col.names = T)


