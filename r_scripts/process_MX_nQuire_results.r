# Script for processing nQuire results

### LOAD PACKAGES ###

### LOAD DATA ###

args <- commandArgs(trailingOnly = T)

res_file <- as.character(args[1])
lrd_res <- read.table(file = res_file, header = T, stringsAsFactors = F,
  sep = '\t')

snp_file <- as.character(args[2])
snp_res <- read.table(snp_file, header = F, stringsAsFactors = F)

### SET OUTPUTS
# get sample name from file name
samp_name <- unlist(lapply(
    strsplit(lrd_res$file, split = '_'), function(x) x[1]))[1]

# Full matrix of results
samp_full_res_out <- paste(samp_name, '.lrd_results_processed.txt', sep = '')

# Summary with header
samp_summary_withhead <- paste(samp_name, '.nQuire_res_summary.txt', sep = '')

# Summary without header; for concatenating
samp_summary_nohead <- paste(samp_name, '.nQuire_res_summary.no_head.txt',
  sep = '')

### SET VARIABLES ###

##################
# get coverage and quality info
cov_vec <- as.numeric(gsub('c', '', 
  unlist(lapply(strsplit(lrd_res$file, split = '_'), function(x) x[2]))))
qual_vec <- as.numeric(gsub('-bedcc.bin', '', gsub('q', '', unlist(
  lapply(strsplit(lrd_res$file, split = '_'), function(x) x[3])))))

lrd_res$coverage <- cov_vec
lrd_res$qual <- qual_vec

# add SNP info
lrd_res$nSNPs <- snp_res[,1]

# calculate the proportion of d_loglikelihood
lrd_res$d_dip_portion <- (lrd_res$d_dip /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tri_portion <- (lrd_res$d_tri /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tet_portion <- (lrd_res$d_tet /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

######
# Generate summary stats

# Calculate lm slope
slope_dip <- summary(
  lm(lrd_res$d_dip_portion ~ lrd_res$coverage))$coefficients[2,1]
slope_tri <- summary(
  lm(lrd_res$d_tri_portion ~ lrd_res$coverage))$coefficients[2,1]
slope_tet <- summary(
  lm(lrd_res$d_tet_portion ~ lrd_res$coverage))$coefficients[2,1]

tot_c10_ind <- which(lrd_res$coverage == 10)
tot_c20_ind <- which(lrd_res$coverage == 20)
tot_c40_ind <- which(lrd_res$coverage == 40)

samp_res <- data.frame(samp = samp_name,
  #
  d_dip_portion_10 = lrd_res$d_dip_portion[tot_c10_ind],
  d_tri_portion_10 = lrd_res$d_tri_portion[tot_c10_ind],
  d_tet_portion_10 = lrd_res$d_tet_portion[tot_c10_ind],
  #
  d_dip_portion_20 = lrd_res$d_dip_portion[tot_c20_ind],
  d_tri_portion_20 = lrd_res$d_tri_portion[tot_c20_ind],
  d_tet_portion_20 = lrd_res$d_tet_portion[tot_c20_ind],
  #
  d_dip_portion_40 = lrd_res$d_dip_portion[tot_c40_ind],
  d_tri_portion_40 = lrd_res$d_tri_portion[tot_c40_ind],
  d_tet_portion_40 = lrd_res$d_tet_portion[tot_c40_ind],
  #
  nSNPS_10 = lrd_res$nSNPs[tot_c10_ind],
  nSNPS_20 = lrd_res$nSNPs[tot_c20_ind],
  nSNPS_40 = lrd_res$nSNPs[tot_c40_ind],
  #
  slope_dip = slope_dip,
  slope_tri = slope_tri,
  slope_tet = slope_tet,
  #
  stringsAsFactors = F)

#####
# Save results

write.table(lrd_res, file = samp_full_res_out, quote = F, sep = '\t',
  row.names = F, col.names = T)
write.table(samp_res, file = samp_summary_withhead, quote = F, sep = '\t',
  row.names = F, col.names = T)
write.table(samp_res, file = samp_summary_nohead, quote = F, sep = '\t',
  row.names = F, col.names = F)


quit(save = 'no')


