# Use nQuire to determine ploidy of samples

## Test workflow

```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_ABHM_workflow.sh

```



```
test_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/test_full_6_lrd_results.txt'

lrd_df <- read.table(file = test_file, header = T, stringsAsFactors = F, 
  sep = '\t')

snp_file <- ''
snp_res <- read.table(snp_file, header = F, stringsAsFactors = F)

cov_vec <- unlist(lapply(strsplit(lrd_df$file, split = '_'), function(x) x[2]))
qual_vec <- unlist(
             lapply(strsplit(lrd_df$file, split = '_'), function(x) x[3]))

samp_name <- unlist(lapply(
    strsplit(lrd_df$file, split = '_'), function(x) x[1]))[1]

#######

samp_full_res_out <- paste(samp_name, '.lrd_results_processed.txt', sep = '')
samp_summary_withhead <- paste(samp_name, '.nQuire_res_summary.txt', sep = '')
samp_summary_nohead <- paste(samp_name, '.nQuire_res_summary.no_head.txt', 
  sep = '')

######

lrd_df$coverage <- cov_vec
lrd_df$qual <- qual_vec

lrd_df$nSNPs <- snp_res[,1]

lrd_res$d_dip_portion <- (lrd_res$d_dip /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tri_portion <- (lrd_res$d_tri /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tet_portion <- (lrd_res$d_tet /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

####

slope_dip <- summary(
  lm(lrd_res$d_dip_portion ~ lrd_res$coverage))$coefficients[2,1]
slope_tri <- summary(
  lm(lrd_res$d_tri_portion ~ lrd_res$coverage))$coefficients[2,1]
slope_tet <- summary(
  lm(lrd_res$d_tet_portion ~ lrd_res$coverage))$coefficients[2,1]

tot_c20_ind <- which(lrd_res$coverage == 20)
tot_c50_ind <- which(lrd_res$coverage == 50)

samp_res <- data.frame(samp = samp_name, 
  d_dip_portion_20 = lrd_res$d_dip_portion[tot_c20_ind],
  d_tri_portion_20 = lrd_res$d_tri_portion[tot_c20_ind],
  d_tet_portion_20 = lrd_res$d_tet_portion[tot_c20_ind],
  #
  d_dip_portion_50 = lrd_res$d_dip_portion[tot_c50_ind],
  d_tri_portion_50 = lrd_res$d_tri_portion[tot_c50_ind],
  d_tet_portion_50 = lrd_res$d_tet_portion[tot_c50_ind],
  #
  nSNPS_20 = lrd_res$nSNPs[tot_c20_ind],
  nSNPS_50 = lrd_res$nSNPs[tot_c50_ind],
  #
  slope_dip = slope_dip,
  slope_tri = slope_tri,
  slope_tet = slope_tet,
  #
  stringsAsFactors = F)

#samp_full_res_out <- paste(samp_name, '.lrd_results_processed.txt', sep = '')
#samp_summary_withhead <- paste(samp_name, '.nQuire_res_summary.txt', sep = '')
#samp_summary_nohead <- paste(samp_name, '.nQuire_res_summary.no_head.txt', 
#  sep = '')

write.table(lrd_res, file = samp_full_res_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)
write.table(samp_res, file = samp_summary_withhead, quote = F, sep = '\t',
  row.names = F, col.names = T)
write.table(samp_res, file = samp_summary_nohead, quote = F, sep = '\t',
  row.names = F, col.names = F)


```
