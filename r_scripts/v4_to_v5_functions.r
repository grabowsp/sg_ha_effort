# Functions used to transfer of v4 to v5 positions for exome-capture data

assign_v5_adj <- function(locus_name, snp_pos_df, locus_pre = ''){
  # Function to calculate the position adjustment for the v5 position. Some
  #   loci have insertions or deletions when looking at v4 vs v5, so need
  #   to calculate the adjustments necessary to get correct v5 positions
  #  This uses the CIGAR score outputed from bwa to get the positions of
  #   insertions or deletions in v4 seq compared to v5 seq
  #  'M' = nucleotide present in both ref and read, though can be mismatch
  #  'D' = nucleotide is present in ref but NOT in the read;  need
  #     to add the size of the D to the v4 position to get the corrected 
  #     relative v5 position
  #  'I' = nucleotide presnet in read but NOT in ref; need to subtract size of
  #     I from v4 position to get corrected relative v5 position
  #
  ######
  sep_char <- ''
  if(length(locus_pre) > 0){
    sep_char <- '_'
  }
  locus_col <- paste(locus_pre, 'locus', sep = sep_char)
  dist_start_col <- paste(locus_pre, 'dist_to_v4_start', sep = sep_char)
  dist_end_col <- paste(locus_pre, 'dist_to_v4_end', sep = sep_char)
  #
  snp_test_inds <- which(snp_pos_df[[locus_col]] == locus_name)
  dist_to_v4_start_vec <- snp_pos_df[[dist_start_col]][snp_test_inds]
  dist_to_v4_end_vec <- snp_pos_df[[dist_end_col]][snp_test_inds]
  snp_adj_vec <- rep(0, times = length(snp_test_inds))
  mask_adj_vec <- rep(0, times = length(snp_test_inds))
  sv_adj_vec <- rep(0, times = length(snp_test_inds))
  #
  snp_cigar_col <- paste(locus_pre, 'v5_cigar', sep = sep_char)
  snp_strand_col <- paste(locus_pre, 'v5_strand_flag', sep = sep_char)
  #
  test_cigar <- unique(snp_pos_df[[snp_cigar_col]][snp_test_inds])
  hit_strand_flag <- snp_pos_df[[snp_strand_col]][snp_test_inds[1]]
  #
  if(test_cigar == '*'){
    snp_adj_vec <- rep(NA, time = length(snp_test_inds))
    return(snp_adj_vec)
    stop()
  }
  #
  c_length_vec <- c(0, as.numeric(
    unlist(strsplit(test_cigar, split = '[A-Z]'))))
  c_code_vec <- unlist(strsplit(test_cigar, split = '[0-9]'))
  c_code_vec <- c('M', c_code_vec[-which(c_code_vec == '')])
  #
  test_m_inds <- which(c_code_vec == 'M')
  adj_mask_vec <- c('S', 'H')
  adj_m_dir_vec <- c(-1,-1)
  adj_code_vec <- c('D', 'I') 
  adj_dir_vec <- c(1, -1)
  # If beginning of seq is masked, need to subract that from rest of the 
  #   position info
  if(c_code_vec[2] == 'S'){
    mask_adjustment <- c_length_vec[2]
    mask_adj_vec <- mask_adj_vec - mask_adjustment
  }
  #
  adjusted_dist_to_start_vec <- dist_to_v4_start_vec + mask_adj_vec 
  adjusted_dist_to_end_vec <- dist_to_v4_end_vec + mask_adj_vec
  #
  if(c_code_vec[2] == 'S'){
    if(hit_strand_flag == 0){
      snp_adj_vec[which(dist_to_v4_start_vec <= c_length_vec[2])] <- NA
    }
    if(hit_strand_flag == 16){
      snp_adj_vec[which(dist_to_v4_end_vec <= c_length_vec[2])] <- NA
    }
  }
  c_last <- length(c_code_vec)
  if(c_code_vec[c_last] == 'S'){
    if(hit_strand_flag == 0){
      snp_adj_vec[which(dist_to_v4_end_vec <= c_length_vec[c_last])] <- NA
    }
    if(hit_strand_flag == 16){
      snp_adj_vec[which(dist_to_v4_start_vec <= c_length_vec[c_last])] <- NA
    } 
  }
  #
  for(acv in seq(length(adj_code_vec))){
    tmp_adj_code <- adj_code_vec[acv]
    tmp_adj <- adj_dir_vec[acv]
    test_indel_inds <- which(c_code_vec == tmp_adj_code)
    for(id_ind in test_indel_inds){
      tmp_indel_length <- c_length_vec[id_ind]*tmp_adj
      tmp_loc_pos <- sum(c_length_vec[test_m_inds[which(test_m_inds < id_ind)]])
      tmp_adj_inds <- which(adjusted_dist_to_start_vec > tmp_loc_pos)
      if(hit_strand_flag == 16){
        tmp_adj_inds <- which(adjusted_dist_to_end_vec > tmp_loc_pos)
      }
      sv_adj_vec[tmp_adj_inds] <- sv_adj_vec[tmp_adj_inds] + tmp_indel_length
    }
  }
  snp_adj_vec <- snp_adj_vec + mask_adj_vec + sv_adj_vec
  #
  return(snp_adj_vec)
}

