# Analysis of K=3 STRUCTURE results for `expandsamps`

## Calculate corelations and make plots
```
bash
source activate R_analysis

data_dir <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis/'

N_K <- 3

perm_files <- system(paste('ls ', data_dir, 'K', N_K, '_permute_clumpp*', 
  sep = ''), intern=T)

perm_list <- list()

for(pf in seq(length(perm_files))){
  perm_list[[pf]] <- read.table(perm_files[pf], header = F, 
    stringsAsFactors = F)
}

cor(perm_list[[1]][,6], perm_list[[3]][,6])

cor_mat_1 <- matrix(NA, nrow = 4, ncol = 4)
cor_mat_2 <- matrix(NA, nrow = 4, ncol = 4)
cor_mat_3 <- matrix(NA, nrow = 4, ncol = 4)

for(i in seq(4)){
  for(j in seq(4)){
    cor_mat_1[i,j] <- cor(perm_list[[i]][,6], perm_list[[j]][,6])
    cor_mat_2[i,j] <- cor(perm_list[[i]][,7], perm_list[[j]][,7])
    cor_mat_3[i,j] <- cor(perm_list[[i]][,8], perm_list[[j]][,8])
  }
}

cor_mat_1
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.9998298 0.9997375 0.9990131
[2,] 0.9998298 1.0000000 0.9998924 0.9990244
[3,] 0.9997375 0.9998924 1.0000000 0.9989213
[4,] 0.9990131 0.9990244 0.9989213 1.0000000

cor_mat_2
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.9999104 0.9998574 0.9994272
[2,] 0.9999104 1.0000000 0.9999577 0.9994413
[3,] 0.9998574 0.9999577 1.0000000 0.9994000
[4,] 0.9994272 0.9994413 0.9994000 1.0000000

cor_mat_3
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.9999524 0.9998440 0.9997395
[2,] 0.9999524 1.0000000 0.9998741 0.9997489
[3,] 0.9998440 0.9998741 1.0000000 0.9996479
[4,] 0.9997395 0.9997489 0.9996479 1.0000000


library(ggplot2)
library(patchwork)
library(reshape)

# Get library info

struc_info_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_3.1_q'
struc_info <- read.table(struc_info_file, header = F, stringsAsFactors = F)

samp_info <- data.frame(LIB = struc_info[,1], stringsAsFactors = F)

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

ploidy_inds <- c()

for(i in samp_info$LIB){
  tmp_ind <- which(ploidy$lib == i)
  ploidy_inds <- c(ploidy_inds, tmp_ind)
}

samp_info$subpop <- as.factor(ploidy$SUBPOP_SNP[ploidy_inds])
samp_info$ploidy <- as.factor(ploidy$total_ploidy_2[ploidy_inds])

# POP1

combo_list_1 <- list()

combo_list_1[[1]] <- data.frame(samp_info, p1 = perm_list[[1]][,6], 
  p2 = perm_list[[2]][,6], comp = 'tet_v_dip', stringsAsFactors = F)
combo_list_1[[2]] <- data.frame(samp_info, p1 = perm_list[[1]][,6], 
  p2 = perm_list[[3]][,6], comp = 'tet_v_NA', stringsAsFactors = F)
combo_list_1[[3]] <- data.frame(samp_info, p1 = perm_list[[1]][,6], 
  p2 = perm_list[[4]][,6], comp = 'tet_v_hap', stringsAsFactors = F)

combo_list_1[[4]] <- data.frame(samp_info, p1 = perm_list[[2]][,6], 
  p2 = perm_list[[3]][,6], comp = 'dip_v_NA', stringsAsFactors = F)
combo_list_1[[5]] <- data.frame(samp_info, p1 = perm_list[[2]][,6], 
  p2 = perm_list[[4]][,6], comp = 'dip_v_hap', stringsAsFactors = F)

combo_list_1[[6]] <- data.frame(samp_info, p1 = perm_list[[3]][,6], 
  p2 = perm_list[[4]][,6], comp = 'NA_v_hap', stringsAsFactors = F)

combo_df_1 <- data.frame(p1 = unlist(lapply(combo_list_1, function(x) x$p1)),
  p2 = unlist(lapply(combo_list_1, function(x) x$p2)),
  comp = unlist(lapply(combo_list_1, function(x) x$comp)),
  stringsAsFactors = T
)

# POP2

combo_list_2 <- list()

combo_list_2[[1]] <- data.frame(samp_info, p1 = perm_list[[1]][,7], 
  p2 = perm_list[[2]][,7], comp = 'tet_v_dip', stringsAsFactors = F)
combo_list_2[[2]] <- data.frame(samp_info, p1 = perm_list[[1]][,7], 
  p2 = perm_list[[3]][,7], comp = 'tet_v_NA', stringsAsFactors = F)
combo_list_2[[3]] <- data.frame(samp_info, p1 = perm_list[[1]][,7], 
  p2 = perm_list[[4]][,7], comp = 'tet_v_hap', stringsAsFactors = F)

combo_list_2[[4]] <- data.frame(samp_info, p1 = perm_list[[2]][,7], 
  p2 = perm_list[[3]][,7], comp = 'dip_v_NA', stringsAsFactors = F)
combo_list_2[[5]] <- data.frame(samp_info, p1 = perm_list[[2]][,7], 
  p2 = perm_list[[4]][,7], comp = 'dip_v_hap', stringsAsFactors = F)

combo_list_2[[6]] <- data.frame(samp_info, p1 = perm_list[[3]][,7], 
  p2 = perm_list[[4]][,7], comp = 'NA_v_hap', stringsAsFactors = F)

combo_df_2 <- data.frame(p1 = unlist(lapply(combo_list_2, function(x) x$p1)),
  p2 = unlist(lapply(combo_list_2, function(x) x$p2)),
  comp = unlist(lapply(combo_list_2, function(x) x$comp)),
  stringsAsFactors = T
)

# POP3

combo_list_3 <- list()

combo_list_3[[1]] <- data.frame(samp_info, p1 = perm_list[[1]][,8],
  p2 = perm_list[[2]][,8], comp = 'tet_v_dip', stringsAsFactors = F)
combo_list_3[[2]] <- data.frame(samp_info, p1 = perm_list[[1]][,8],
  p2 = perm_list[[3]][,8], comp = 'tet_v_NA', stringsAsFactors = F)
combo_list_3[[3]] <- data.frame(samp_info, p1 = perm_list[[1]][,8],
  p2 = perm_list[[4]][,8], comp = 'tet_v_hap', stringsAsFactors = F)

combo_list_3[[4]] <- data.frame(samp_info, p1 = perm_list[[2]][,8],
  p2 = perm_list[[3]][,8], comp = 'dip_v_NA', stringsAsFactors = F)
combo_list_3[[5]] <- data.frame(samp_info, p1 = perm_list[[2]][,8],
  p2 = perm_list[[4]][,8], comp = 'dip_v_hap', stringsAsFactors = F)

combo_list_3[[6]] <- data.frame(samp_info, p1 = perm_list[[3]][,8],
  p2 = perm_list[[4]][,8], comp = 'NA_v_hap', stringsAsFactors = F)

combo_df_3 <- data.frame(p1 = unlist(lapply(combo_list_3, function(x) x$p1)),
  p2 = unlist(lapply(combo_list_3, function(x) x$p2)),
  comp = unlist(lapply(combo_list_3, function(x) x$comp)),
  stringsAsFactors = T
)

# Plot correlations for all genotypes

gg_combo_list <- list()

gg_combo_list[[1]] <- ggplot(combo_df_1, aes(x = p1, y = p2, col = comp)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  ggtitle('K=3 Pop1 cor')

gg_combo_list[[2]] <- ggplot(combo_df_2, aes(x = p1, y = p2, col = comp)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  ggtitle('K=3 Pop2 cor')

gg_combo_list[[3]] <- ggplot(combo_df_3, aes(x = p1, y = p2, col = comp)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  ggtitle('K=3 Pop3 cor')

out_file <- paste(data_dir, 'structure_K3_allpops_cor.pdf', sep = '')

pdf(out_file, height = 5, width = 6*3)
wrap_plots(gg_combo_list, ncol = 3)
dev.off()

gg_list <- list()

for(i in seq(length(combo_list_1))){
  gg_list[[((6*i)-5)]] <- ggplot(combo_list_1[[i]], 
    aes_string(x = 'p1', y = 'p2', col = 'subpop')) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    ggtitle(paste(combo_list_1[[i]][1,'comp'], 'Pop1', sep = ' '))
  gg_list[[((6*i)-4)]] <- ggplot(combo_list_2[[i]], 
    aes_string(x = 'p1', y = 'p2', col = 'subpop')) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    ggtitle(paste(combo_list_2[[i]][1,'comp'], 'Pop2', sep = ' '))
  gg_list[[((6*i)-3)]] <- ggplot(combo_list_3[[i]], 
    aes_string(x = 'p1', y = 'p2', col = 'subpop')) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    ggtitle(paste(combo_list_3[[i]][1,'comp'], 'Pop3', sep = ' '))
#
  gg_list[[((6*i)-2)]] <- ggplot(combo_list_1[[i]], 
    aes_string(x = 'p1', y = 'p2', col = 'ploidy')) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    ggtitle(paste(combo_list_1[[i]][1,'comp'], 'Pop1', sep = ' '))
  gg_list[[((6*i)-1)]] <- ggplot(combo_list_2[[i]], 
    aes_string(x = 'p1', y = 'p2', col = 'ploidy')) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    ggtitle(paste(combo_list_2[[i]][1,'comp'], 'Pop2', sep = ' '))
  gg_list[[(6*i)]] <- ggplot(combo_list_3[[i]],
    aes_string(x = 'p1', y = 'p2', col = 'ploidy')) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    ggtitle(paste(combo_list_3[[i]][1,'comp'], 'Pop3', sep = ' '))
}

out_file_2 <- paste(data_dir, 'structure_K3_allPops_singleComps_cor.pdf', 
  sep = '')

pdf(out_file_2, height = 5*12, width = 6*3)
wrap_plots(gg_list, ncol = 3)
dev.off()





```

