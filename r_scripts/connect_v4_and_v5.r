# Script for connecting exome capture v4 and v5 positions

# LOAD LIBRARIES AND FUNCTIONS
source(paste('/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/', 
  'v4_to_v5_functions.r', sep = ''))

# LOAD DATA
args <- commandArgs(trailingOnly = T)

# snp_df_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files/v4_exome_sub11_v4_snp_and_loci_info.rds'
snp_df_in <- args[1]
snp_df <- readRDS(snp_df_in)

# bwa_res_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files/v4_to_v5_aln_sub11.txt'
bwa_res_in <- args[2]
bwa_res_0 <- read.table(bwa_res_in, header = F, sep = '\n', 
  stringsAsFactors = F)

# OUTPUTS
out_pre <- args[3]
new_snp_txt_out <- paste(out_pre, '_v4_and_v5info.txt', sep = '')
new_snp_rds_out <- paste(out_pre, '_v4_and_v5info.rds', sep = '') 

##########

#re-format bwa results
bwa_res_list <- strsplit(bwa_res_0[,1], split = '\t')

bwa_res_1 <- data.frame(
  locus = unlist(lapply(bwa_res_list, function(x) x[1])),
  flag = as.numeric(unlist(lapply(bwa_res_list, function(x) x[2]))),
  chr_name = unlist(lapply(bwa_res_list, function(x) x[3])),
  pos = as.numeric(unlist(lapply(bwa_res_list, function(x) x[4]))),
  qual_score = as.numeric(unlist(lapply(bwa_res_list, function(x) x[5]))),
  cigar = unlist(lapply(bwa_res_list, function(x) x[6])),
  other_matches = unlist(lapply(bwa_res_list, function(x) x[length(x)])),
  stringsAsFactors = F)

#prim_hit_inds <- sort(
#  c(which(bwa_res_1$flag == 0), which(bwa_res_1$flag == 16)))

prim_hit_inds <- sort(
  c(which(bwa_res_1$flag == 0), which(bwa_res_1$flag == 16), 
  which(bwa_res_1$flag == 4)))

bwa_res <- bwa_res_1[prim_hit_inds, ]

miss_locs <- setdiff(unique(snp_df$locus), unique(bwa_res$locus))
if(length(miss_locs) > 0){
  n_locs <- length(miss_locs)
  miss_loc_df <- data.frame(locus = miss_locs, flag = 4, chr_name = '*',
  pos = NA, qual_score = 0, cigar = '*', other_matches = 'NA', 
  stringsAsFactors = F)
  bwa_res <- rbind(bwa_res, miss_loc_df)
}

snp_df$v5_chr <- NA
snp_df$v5_strand_flag <- NA
snp_df$v5_loc_start <- NA
snp_df$v5_cigar <- NA
snp_df$map_score <- NA

for(br in seq(nrow(bwa_res))){
  tmp_inds <- which(snp_df$locus == bwa_res[br,1])
  snp_df$v5_chr[tmp_inds] <- bwa_res[br,3]
  snp_df$v5_strand_flag[tmp_inds] <- bwa_res[br,2]
  snp_df$v5_loc_start[tmp_inds] <- bwa_res[br, 4]
  snp_df$v5_cigar[tmp_inds] <- bwa_res[br,6]
  snp_df$map_score[tmp_inds] <- bwa_res[br, 5]
}

snp_df$v5_adj <- unlist(sapply(bwa_res[,1], function(x)
  assign_v5_adj(x, snp_pos_df = snp_df)))

snp_df$v5_pos <- NA
same_strand_inds <- which(snp_df$v5_strand_flag == 0)
op_strand_inds <- which(snp_df$v5_strand_flag == 16)

snp_df$v5_pos[same_strand_inds] <- snp_df$v5_loc_start[same_strand_inds] + 
  (snp_df$dist_to_v4_start[same_strand_inds] + snp_df$v5_adj[same_strand_inds]
  -1)
snp_df$v5_pos[op_strand_inds] <- snp_df$v5_loc_start[op_strand_inds] + 
  (snp_df$dist_to_v4_end[op_strand_inds] + snp_df$v5_adj[op_strand_inds])

bad_map_inds <- which(snp_df$map_score < 60)

snp_df$v5_pos[bad_map_inds] <- NA

write.table(snp_df, file = new_snp_txt_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)
saveRDS(snp_df, file = new_snp_rds_out)

quit(save = 'no')


