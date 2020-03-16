# Script for exploring the nQuire results for the switchgrass samples

# Ploidy designation
# c20_ploidy = ploidy with lowest proportion
# c20_clust_ploidy = ploidy calls based on (visually) clustering the
#   samples from the c20 proportions and for slope_vs_c50 proportions; these
#   clusters end up being identicle

### LOAD PACKAGES ###
library(ggplot2)
library(reshape2)

### LOAD DATA ###
nquire_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'sg_nquire/sg_CDS_nquire/sg_reseq_nQuire_results_summary_total.txt', 
  sep = '')
nquire_res_0 <- read.table(nquire_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

mnp_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/', 
  'pos_03/pos_03_MNP_results_with_meta_info.txt', sep = '')
mnp_res <- read.table(mnp_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

### SET OUTPUTS ###

#out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
#  'pos_03/', sep = '')

nquire_res_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'sg_nquire/sg_CDS_nquire/sg_reseq_nQuire_results_summary_total_v2.0.txt',
  sep = '')

######

nquire_good_libs <- which(nquire_res_0$samp %in% mnp_res$lib)

nquire_res <- nquire_res_0[nquire_good_libs, ]

# Ploidy from c20 results

c20_ploidy <- apply(nquire_res[,
  c('d_dip_portion_20', 'd_tri_portion_20', 'd_tet_portion_20')], 1,
  function(x) which(x == min(x)))
c20_ploidy <- paste((c20_ploidy + 1)*2, 'X', sep = '')

nquire_res$c20_ploidy <- c20_ploidy

######

dip_c20clust_ploidy <- rep(NA, times = length(nrow(nquire_res)))
dip_c20clust_ploidy[nquire_res$'d_dip_portion_20' < 0.4] <- '4X'
dip_c20clust_ploidy[nquire_res$'d_dip_portion_20' > 0.4] <- '8X'

nquire_res$c20cluster_ploidy <- dip_c20clust_ploidy

write.table(nquire_res, file = nquire_res_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

