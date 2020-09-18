# Analysis of K=2 STRUCTURE results for `expandsamps`

## Calculate corelations and make plots
```
bash
source activate R_analysis

data_dir <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis/'

N_K <- 2

perm_files <- system(paste('ls ', data_dir, 'permute_clumpp*K', N_K, sep = ''),
  intern=T)

perm_list <- list()

for(pf in seq(length(perm_files))){
  perm_list[[pf]] <- read.table(perm_files[pf], header = F, 
    stringsAsFactors = F)
}

cor(perm_list[[1]][,6], perm_list[[3]][,6])

cor_mat <- matrix(NA, nrow = 4, ncol = 4)

for(i in seq(4)){
  for(j in seq(4)){
    cor_mat[i,j] <- cor(perm_list[[i]][6], perm_list[[j]][6])
  }
}

> cor_mat
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.9999659 0.9998732 0.9997005
[2,] 0.9999659 1.0000000 0.9998743 0.9997033
[3,] 0.9998732 0.9998743 1.0000000 0.9995763
[4,] 0.9997005 0.9997033 0.9995763 1.0000000

library(ggplot2)
library(patchwork)
library(reshape)

combo_list <- list()

combo_list[[1]] <- data.frame(p1 = perm_list[[1]][6], p2 = perm_list[[2]][6], 
  comp = 'tet_v_dip', stringsAsFactors = F)
combo_list[[2]] <- data.frame(p1 = perm_list[[1]][6], p2 = perm_list[[3]][6], 
  comp = 'tet_v_NA', stringsAsFactors = F)
combo_list[[3]] <- data.frame(p1 = perm_list[[1]][6], p2 = perm_list[[4]][6],
  comp = 'tet_v_hap', stringsAsFactors = F)

combo_list[[4]] <- data.frame(p1 = perm_list[[2]][6], p2 = perm_list[[3]][6],
  comp = 'dip_v_NA', stringsAsFactors = F)
combo_list[[5]] <- data.frame(p1 = perm_list[[2]][6], p2 = perm_list[[4]][6],
  comp = 'dip_v_hap', stringsAsFactors = F)

combo_list[[6]] <- data.frame(p1 = perm_list[[3]][6], p2 = perm_list[[4]][6],
  comp = 'NA_v_hap', stringsAsFactors = F)

tot_c1 <- unlist(lapply(combo_list, function(x) x[,1]))
tot_c2 <- unlist(lapply(combo_list, function(x) x[,2]))
tot_c3 <- unlist(lapply(combo_list, function(x) x$comp))

combo_df <- data.frame(p1 = tot_c1, p2 = tot_c2, comp = tot_c3)

gg_1 <- ggplot(combo_df, aes(x = p1, y = p2, col = comp)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')

out_file <- paste(data_dir, 'structure_K2_P1_cor.pdf', sep = '')

pdf(out_file)
gg_1
dev.off()

gg_list <- list()

for(i in seq(length(combo_list))){
  gg_list[[i]] <- ggplot(combo_list[[i]], aes_string(x = 'V6', y = 'V6.1')) +
    geom_point() + geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    ggtitle(combo_list[[i]][1,'comp'])
}

out_file_2 <- paste(data_dir, 'structure_K2_P1_singleComps_cor.pdf', sep = '')

pdf(out_file_2, height = 5*3, width = 6*2)
wrap_plots(gg_list, ncol = 2)
dev.off()





```

