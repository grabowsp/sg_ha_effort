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

library(adegenet)
library(parallel)

gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(gen_function_file)

adegenet_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/adegenet_functions.r'
source(adegenet_function_file)

data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps'

data_dir <- add_slash(data_dir)

sub_in <- 'Combo.595K.polyploid.CDS.geosamps.genlight.rds'

in_tot <- paste(data_dir, sub_in, sep = '')

gen_tot <- readRDS(in_tot)

n_snps <- 10000

keep_inds <- sort(sample(seq(nLoc(gen_tot)), size = n_snps))

keep_genos <- as.matrix(gen_tot[,keep_inds])

keep_snp_names <- colnames(keep_genos)

ploidy_vec <- gen_tot$ploidy

struc_list <- list()

for(i in seq(nrow(keep_genos))){
  struc_list[[i]] <- gen_struc_genos(genotypes = keep_genos[i, ], 
    lib_name = rownames(keep_genos)[i], ploidy = ploidy_vec[i])
}




