# Call ploidy from MNP results

### LOAD PACKAGES ###
library(ggplot2)
library(reshape2)

### INPUT DATA ###
mnp_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_03/pos_03_MNP_results_with_meta_info.txt', sep = '')
mnp_res <- read.table(mnp_res_file, header = T, sep = '\t',
  stringsAsFactors = F)

nquire_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'sg_nquire/sg_CDS_nquire/sg_reseq_nQuire_results_summary_total_v2.0.txt',
  sep = '')
nquire_res <- read.table(nquire_res_file, header = T, sep = '\t',
  stringsAsFactors = F)


### SET OUTPUTS ###
fig_out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_03/', sep = '')
plot_height <- 5
plot_width <- 5.5

combined_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v1.0.txt'

### SET VARIABLES ###
cds_cut <- 100
genic_cut <- 400
thresh_per <- 0.1
cov_cut <- 1e10

##################
# Will only assign mnp ploidy beyond a 10% threshold from the "cutoff"

na_inds <- which(mnp_res$seq_cov < cov_cut)

cds_8X_inds <- which(mnp_res$cds_mnp_stand >= ((1+thresh_per)*cds_cut))
cds_4X_inds <- which(mnp_res$cds_mnp_stand <= ((1-thresh_per)*cds_cut))

mnp_res$cds_mnp_ploidy <- '?X'
mnp_res$cds_mnp_ploidy[cds_8X_inds] <- '8X'
mnp_res$cds_mnp_ploidy[cds_4X_inds] <- '4X'
mnp_res$cds_mnp_ploidy[na_inds] <- '?X'

###

genic_8X_inds <- which(mnp_res$genic_mnp_stand >= ((1+thresh_per)*genic_cut))
genic_4X_inds <- which(mnp_res$genic_mnp_stand <= ((1-thresh_per)*genic_cut))

mnp_res$genic_mnp_ploidy <- '?X'
mnp_res$genic_mnp_ploidy[genic_8X_inds] <- '8X'
mnp_res$genic_mnp_ploidy[genic_4X_inds] <- '4X'
mnp_res$genic_mnp_ploidy[na_inds] <- '?X'


###
# add nquire info
nquire_colnames <- paste('nquire_', colnames(nquire_res), sep = '')

mnp_res[, nquire_colnames] <- nquire_res

write.table(mnp_res, file = combined_res_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

quit(save = 'no')
