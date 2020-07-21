# Combine Results from 'geosamps' STRUCTURE results into single table

#bash
#source activate R_analysis

### LOAD INPUTS ###
# all paths are for HA cluster
c_res_file_1 <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo.clumpp.processed'
c_res_file_2 <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo.clumpp.processed'
c_res_file_3 <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo.clumpp.processed'

file_vec <- c(c_res_file_1, c_res_file_2, c_res_file_3)

res_list <- list()
for(i in seq(length(file_vec))){
  res_list[[i]] <- read.table(file_vec[i], header = F, sep = '\t', 
    stringsAsFactors = F)
}

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t',
  quote = "", comment.char = '$')

### SET OUTPUTS ####

out_file_pre <- 'upgeo'
out_file_short <- paste(out_file_pre, 'structure_combined_results.txt', 
  sep = '_')
out_dir <- dirname(c_res_file_1)

out_file <- file.path(out_dir, out_file_short)

##########

tot_df <- data.frame(res_list[[1]][,1], stringsAsFactors = F)
names(tot_df) <- 'LIB'

ploidy_ord <- c()
meta_ord <- c()
for(i in res_list[[1]][,1]){
  tmp_ploid <- which(ploidy$lib == i)
  tmp_meta <- which(meta$LIBRARY == i)
  ploidy_ord <- c(ploidy_ord, tmp_ploid)
  meta_ord <- c(meta_ord, tmp_meta)
}

tot_df <- data.frame(LIB = res_list[[1]][,1], 
  PLANT_ID = meta$PLANT_ID[meta_ord], 
  LATITUDE = meta$LATITUDE[meta_ord], LONGITUDE = meta$LONGITUDE[meta_ord],
  STATE = meta$STATE[meta_ord], PLOIDY_v1 = ploidy$total_ploidy[ploidy_ord],
  PLOIDY_v2 = ploidy$total_ploidy_2[ploidy_ord],
  OLD_SUBPOP = meta$SUBPOP_SNP[meta_ord], OLD_PLOIDY = meta$PLOIDY[meta_ord],
  stringsAsFactors = F)


for(i in seq(length(res_list))){
  tmp_df <- res_list[[i]]
  rownames(tmp_df) <- tmp_df[,1]
  colnames(tmp_df)[1] <- 'LIB'
  n_k <- ncol(tmp_df)-1
  tmp_colnames <- paste('K', n_k, '_', seq(n_k), sep = '')
  colnames(tmp_df)[-1] <- tmp_colnames
  tot_df[, tmp_colnames] <- tmp_df[tot_df$LIB, tmp_colnames]
}

write.table(tot_df, file = out_file, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')



