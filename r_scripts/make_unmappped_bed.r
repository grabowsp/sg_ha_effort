# Script for making snp dataframe to connect v4 and v5 positions as well as
#   BED file to feed to BWA to map v4 loci to v5 genome

# LOAD LIBRARIES
v4_to_v5_func_file <- paste('/home/grabowsky/tools/workflows/sg_ha_effort/', 
  'r_scripts/v4_to_v5_functions.r', sep = '')
source(v4_to_v5_func_file)

# LOAD DATA
args <- commandArgs(trailingOnly = T)

# snps_file_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files/snps_sub15_v4_and_v5info.rds'
snps_file_in <- args[1]
snps <- readRDS(snps_file_in)

# OUTPUTS
out_pre <- args[2]

locus_pre <- args[3]
# locus_pre <- 'try2'

bed_out <- paste(out_pre, locus_pre, 'v4_loci.bed', sep = '_')
snp_df_rds_out <- paste(out_pre, locus_pre, 'v4_snp_and_loci_info.rds', 
  sep = '_')
snp_df_txt_out <- paste(out_pre, locus_pre, 'v4_snp_and_loci_info.txt', 
  sep = '_')

# CONSTANTS

# the maximum distance between SNPs to be in same locus
dif_cut <- as.numeric(args[4])
# dif_cut <- 50

# the distance from first and last SNP in locus to use for mapping
seq_edge <- as.numeric(args[5])
# seq_edge <- 100

# VARIABLES


######
new_loc_col <- paste(locus_pre, 'locus', sep = '_')
new_loc_start_col <- paste(locus_pre, 'v4_loc_start', sep = '_')
new_loc_start_dist_col <- paste(locus_pre, 'dist_to_v4_start', sep = '_')
new_loc_end_col <- paste(locus_pre, 'v4_loc_end', sep = '_')
new_loc_end_dist_col <- paste(locus_pre, 'dist_to_v4_end', sep = '_')

snps[, c(new_loc_col, new_loc_start_col, new_loc_start_dist_col, 
  new_loc_end_col, new_loc_end_dist_col)] <- NA

snp_chrs <- unique(snps$v4_chr)

bed_list <- list()

for(sc in snp_chrs){
  sub_inds <- intersect(which(snps$v4_chr == sc),which(is.na(snps$v5_pos)))
  if(length(sub_inds) == 0){
    next
  }
  snp_sub <- snps[sub_inds, ]
  pos_1 <- snp_sub$orig_pos[1:(nrow(snp_sub)-1)]
  pos_2 <- snp_sub$orig_pos[2:(nrow(snp_sub))]
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
  locus_names <- paste('locus', locus_pre, sc, 
    unlist(lapply(group_list, function(x) x[1])), sep = '_')
  # make dataframe for BED file
  sc_bed_df <- data.frame(chr = snp_sub[1,1], start = seq_start_pos, 
    end = seq_end_pos, name = locus_names, score = 1, strand = '+', 
    stringsAsFactors = F)
  bed_list[[sc]] <- sc_bed_df
  # add snp info to datafram
  for(i in seq(length(group_list))){
    add_new_inds <- snps$orig_pos %in% group_list[[i]]
    snps[[new_loc_col]][add_new_inds] <- locus_names[i]
    tmp_v4_start <- group_list[[i]][1] - seq_edge
    snps[[new_loc_start_col]][add_new_inds] <- tmp_v4_start
    snps[[new_loc_start_dist_col]][add_new_inds] <- (group_list[[i]] - 
      tmp_v4_start)
    tmp_v4_end <- group_list[[i]][length(group_list[[i]])] + seq_edge
    snps[[new_loc_end_col]][add_new_inds] <- tmp_v4_end
    snps[[new_loc_end_dist_col]][add_new_inds] <- (tmp_v4_end - group_list[[i]])
  }
}

if(length(bed_list) == 0){
  full_bed_df <- c()
}

if(length(bed_list) >= 1){
  full_bed_df <- bed_list[[1]]
}

if(length(bed_list) > 1){
  for(scl in c(2:length(bed_list))){
    full_bed_df <- rbind(full_bed_df, bed_list[[scl]])
  }
}

write.table(full_bed_df, file = bed_out, quote = F, sep = '\t', row.names = F,
  col.names = F)
saveRDS(snps, file = snp_df_rds_out)
write.table(snps, file = snp_df_txt_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)

quit(save = 'no')


