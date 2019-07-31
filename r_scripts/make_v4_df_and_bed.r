# Script for making snp dataframe to connect v4 and v5 positions as well as
#   BED file to feed to BWA to map v4 loci to v5 genome

# LOAD LIBRARIES
v4_to_v5_func_file <- paste('/home/grabowsky/tools/workflows/sg_ha_effort/', 
  'r_scripts/v4_to_v5_functions.r', sep = '')
source(v4_to_v5_func_file)

# LOAD DATA
args <- commandArgs(trailingOnly = T)

snps_file_in <- args[1]
snps <- read.table(snps_file_in, header = F, sep = '\t', stringsAsFactors = F)

#snp_file_in <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
#  'exome/v4_snps/', 'combo_1168_v4_top1000.txt', sep = '')
#snps <- read.table(short_snps_file, header = T, sep = '\t',
#  stringsAsFactors = F)

# OUTPUTS
out_pre <- args[2]
bed_out <- paste(out_pre, '_v4_loci.bed', sep = '')
snp_df_rds_out <- paste(out_pre, '_v4_snp_and_loci_info.rds', sep = '')
snp_df_txt_out <- paste(out_pre, '_v4_snp_and_loci_info.txt', sep = '')

# CONSTANTS
# the maximum distance between SNPs to be in same locus
dif_cut <- as.numeric(args[3])
#dif_cut <- 500

# the distance from first and last SNP in locus to use for mapping
seq_edge <- as.numeric(args[4])
#seq_edge <- 100

# VARIABLES


######
snp_chrs <- unique(snps[ , 1])

bed_list <- list()
snp_df_list <- list()

for(sc in snp_chrs){
  sub_inds <- which(snps[ ,1] == sc)
  snp_sub <- snps[sub_inds, ]
  pos_1 <- snp_sub[1:(nrow(snp_sub)-1), 2]
  pos_2 <- snp_sub[2:(nrow(snp_sub)), 2]
  pos_dif <- pos_2 - pos_1
  dif_inds <- which(pos_dif > dif_cut)
  dif_start_inds <- c(1, dif_inds + 1)
  dif_end_inds <- c(dif_inds, nrow(snp_sub))
  #
  # make list of chr positions in each locus
  group_list <- list()
  for(i in seq(length(dif_start_inds))){
    group_list[[i]] <- snp_sub[c(dif_start_inds[i]:dif_end_inds[i]), 2]
  }
  seq_start_pos <- unlist(lapply(group_list, function(x) x[1] - seq_edge))
  seq_end_pos <- unlist(lapply(group_list, function(x) x[length(x)] + seq_edge))
  #
  locus_names <- paste('locus', sc, 
    unlist(lapply(group_list, function(x) x[1])), sep = '_')
  # make dataframe for BED file
  sc_bed_df <- data.frame(chr = snp_sub[1,1], start = seq_start_pos, 
    end = seq_end_pos, name = locus_names, score = 1, strand = '+', 
    stringsAsFactors = F)
  bed_list[[sc]] <- sc_bed_df
  # make snp info dataframe
  locus_df_list <- list()
  for(i in seq(length(group_list))){
    tmp_v4_pos <- group_list[[i]]
    tmp_loc_name <- locus_names[i]
    tmp_v4_start <- group_list[[i]][1] - seq_edge
    tmp_v4_end <- group_list[[i]][length(group_list[[i]])] + seq_edge
    tmp_dist_to_v4_start <- group_list[[i]] - tmp_v4_start
    tmp_dist_to_v4_end <- tmp_v4_end - group_list[[i]]
    locus_df_list[[i]] <- data.frame(v4_pos = tmp_v4_pos, locus = tmp_loc_name,
      v4_loc_start = tmp_v4_start, dist_to_v4_start = tmp_dist_to_v4_start,
      v4_loc_end = tmp_v4_end, dist_to_v4_end = tmp_dist_to_v4_end,
      stringsAsFactors = F)
  }
  locus_df_full <- data.frame(v4_chr = snp_sub[, 1], orig_pos = snp_sub[,2],
    v4_snp_name = paste(snp_sub[, 1], snp_sub[, 2], sep = '_'),
    ref = snp_sub[,3], alleles = snp_sub[,4],
    v4_pos = unlist(lapply(locus_df_list, function(x) x$v4_pos)),
    locus = unlist(lapply(locus_df_list, function(x) x$locus)),
    v4_loc_start = unlist(lapply(locus_df_list, function(x) x$v4_loc_start)),
    dist_to_v4_start = unlist(lapply(locus_df_list, function(x)
      x$dist_to_v4_start)),
    v4_loc_end = unlist(lapply(locus_df_list, function(x) x$v4_loc_end)),
    dist_to_v4_end = unlist(lapply(locus_df_list, function(x)
      x$dist_to_v4_end)),
    stringsAsFactors = F)
  #
  snp_df_list[[sc]] <- locus_df_full
}

full_bed_df <- bed_list[[1]]
full_loci_df <- snp_df_list[[1]]

if(length(snp_chrs) > 1){
  for(scl in c(2:length(snp_chrs))){
    full_bed_df <- rbind(full_bed_df, bed_list[[scl]])
    full_loci_df <- rbind(full_loci_df, snp_df_list[[scl]])
  }
}

write.table(full_bed_df, file = bed_out, quote = F, sep = '\t', row.names = F,
  col.names = F)
saveRDS(full_loci_df, file = snp_df_rds_out)
write.table(full_loci_df, file = snp_df_txt_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)

quite(save = 'no')


