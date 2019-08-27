# Script to correct the alleles for SNPs mapped to (-) strand in v5

# LOAD LIBRARIES AND FUNCTIONS
source(paste('/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/',
  'v4_to_v5_functions.r', sep = ''))

# LOAD DATA
args <- commandArgs(trailingOnly = T)

# snp_df_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files/sub00_try3_v4_and_v5info.rds'
snp_df_in <- args[1]
snp_df <- readRDS(snp_df_in)

# OUTPUTS
out_pre <- args[2]
# out_pre <- 'sub00'

new_snp_txt_out <- paste(out_pre, 'rev_corrected', 'v4_and_v5info.txt', 
  sep = '_')
new_snp_rds_out <- paste(out_pre, 'rev_corrected', 'v4_and_v5info.rds', 
  sep = '_')

#############3


snp_df$final_strand_flag <- NA

t1_good_inds <- which(snp_df$map_score == 60)
t2_good_inds <- which(snp_df$try2_map_score == 60)
t3_good_inds <- which(snp_df$try3_map_score== 60)

snp_df$final_strand_flag[t1_good_inds] <- snp_df$v5_strand_flag[t1_good_inds]
snp_df$final_strand_flag[t2_good_inds] <- snp_df$try2_v5_strand_flag[
  t2_good_inds]
snp_df$final_strand_flag[t3_good_inds] <- snp_df$try3_v5_strand_flag[
  t3_good_inds]

snp_df$v5_ref <- snp_df$ref
ref_alleles <- c('A', 'T', 'C', 'G')
rev_alleles <- c('T', 'A', 'G', 'C')

for(i in ref_alleles){
  tmp_rev_allele_inds <- intersect(which(snp_df$ref == i),  
    which(snp_df$final_strand_flag == '16'))
  if(length(tmp_rev_allele_inds) > 0){
    allele_ind <- which(ref_alleles == i)
    tmp_rev_allele <- rev_alleles[allele_ind]
    snp_df$v5_ref[tmp_rev_allele_inds] <- tmp_rev_allele
  }
}

tmp_allele_mat <- combn(c('A','T', 'C', 'G'), m = 2)
allele_vec_1 <- paste(tmp_allele_mat[1,], tmp_allele_mat[2,], sep = '')
allele_vec_2 <- paste(tmp_allele_mat[2,], tmp_allele_mat[1,], sep = '')
allele_vec <- c(allele_vec_1, allele_vec_2)
rev_allele_vec <- sapply(allele_vec, function(x) chartr('ATCG', 'TAGC', x))

snp_df$v5_alleles <- snp_df$alleles

for(i in allele_vec){
  tmp_inds <- intersect(which(snp_df$alleles == i),
    which(snp_df$final_strand_flag == 16))
  if(length(tmp_inds) > 0){
    snp_df$v5_alleles[tmp_inds] <- rev_allele_vec[[i]]
  }
}

write.table(snp_df, file = new_snp_txt_out, quote = F, sep = '\t',
  row.names = F, col.names = F)
saveRDS(snp_df, file = new_snp_rds_out)

quit(save = 'no')


