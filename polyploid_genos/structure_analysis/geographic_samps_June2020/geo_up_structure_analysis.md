# Steps and notes for analyzing STRUCTURE results for `up_geo`

## Overview
* This round of STRUCTURE run using tetrasomic genotypes for all samples
* delta-K results show K=2 is best and K=3 and K=4 are decent
* CLUMPP results used here to make graph

## K =2
### Make bargraph of CLUMPP results
```
bash
source activate R_analysis

library(ggplot2)

struc_func_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/struc_functions.r'
source(struc_func_file)

c_res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo.clumpp.processed'

c_res <- read.table(c_res_file, header = F, sep = '\t', stringsAsFactors = F)

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t',
  quote = "", comment.char = '$')

### evaluate what the groups are

pop_info_list <- list()

test_order_mat <- apply(c_res[, -1], 1, order, decreasing = T)

for(i in seq(ncol(c_res)-1)){
  test_pop <- i
  pop_name <- paste('Population', i)
  test_libs <- c_res[which(apply(test_order_mat, 2, 
    function(x) x[1] == test_pop)), 1]
  ploidy_inds <- which(ploidy$lib %in% test_libs)
  meta_inds <- which(meta$LIBRARY %in% test_libs)
  tmp_list <- list()
  tmp_list[['Chloroplast Ecotype']] <- table(ploidy$ECOTYPE_SNP_CHLR
    [ploidy_inds])
  tmp_list[['Old Subpop']] <- table(ploidy$SUBPOP_SNP[ploidy_inds])
  tmp_list[['Ploidy']] <- table(ploidy$total_ploidy_2[ploidy_inds])
  tmp_list[['State']] <- table(meta$STATE[meta_inds])
  pop_info_list[[pop_name]] <- tmp_list
}

pop_info_list

# pop 1 = mostly 8X: 84 8X, 13 4X
# pop 2 = all have Midwest designation, 178 4X, 2 8X

pop_order <- c(2,1)

######

test_samp_orders <- order_struc_samps(clumpp_result_df = c_res, 
  pop_order = pop_order, zero_cut = 0.1)

#####

res_ord <- c_res[test_samp_orders, ]

samp_vec <- rep(res_ord[,1], times = ncol(c_res)-1)
group_vec <- rep(seq(ncol(c_res)-1), each = nrow(res_ord))

res_vec <- c()
for(j in c(2:ncol(res_ord))){
  res_vec <- c(res_vec, res_ord[,j])
}

plot_df <- data.frame(LIB = samp_vec, GROUP = group_vec, MEMB = res_vec,
  stringsAsFactors = F)

plot_df$LIB <- factor(plot_df$LIB, levels = res_ord[,1])

plot_df$GROUP <- factor(plot_df$GROUP, levels = pop_order)

group_col_vec <- c()
group_col_vec['1'] <- 'darkorchid2'
group_col_vec['2'] <- 'blue2'

group_lab_vec <- c()
group_lab_vec['1'] <- 'Mainly 8X'
group_lab_vec['2'] <- 'Midwest 4X'
group_lab_vec <- group_lab_vec[pop_order]

#group_palette <- scale_colour_manual(name = 'K_Group', values = group_col_vec)

gg_bar <- ggplot(plot_df, aes(y = MEMB, x = LIB, fill = GROUP)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual('Group', values = group_col_vec, 
    labels = group_lab_vec) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = 'Group Membership')

out_file <- gsub('clumpp.processed', 'STRUCTURE.memb.pdf', c_res_file)

pdf(out_file, width = 12*(nrow(c_res)/125), height = 4)
gg_bar
dev.off()
```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo.STRUCTURE.memb.pdf`

########################

## K =3
### Make bargraph of CLUMPP results
```
bash
source activate R_analysis

library(ggplot2)

struc_func_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/struc_functions.r'
source(struc_func_file)

c_res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo.clumpp.processed'

c_res <- read.table(c_res_file, header = F, sep = '\t', stringsAsFactors = F)

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t',
  quote = "", comment.char = '$')

### evaluate what the groups are

pop_info_list <- list()

test_order_mat <- apply(c_res[, -1], 1, order, decreasing = T)

for(i in seq(ncol(c_res)-1)){
  test_pop <- i
  pop_name <- paste('Population', i)
  test_libs <- c_res[which(apply(test_order_mat, 2, 
    function(x) x[1] == test_pop)), 1]
  ploidy_inds <- which(ploidy$lib %in% test_libs)
  meta_inds <- which(meta$LIBRARY %in% test_libs)
  tmp_list <- list()
  tmp_list[['Chloroplast Ecotype']] <- table(ploidy$ECOTYPE_SNP_CHLR
    [ploidy_inds])
  tmp_list[['Old Subpop']] <- table(ploidy$SUBPOP_SNP[ploidy_inds])
  tmp_list[['Ploidy']] <- table(ploidy$total_ploidy_2[ploidy_inds])
  tmp_list[['State']] <- table(meta$STATE[meta_inds])
  pop_info_list[[pop_name]] <- tmp_list
}

pop_info_list

# pop 1 = 4X Midwest, 61 of 77 from North Dakota
# pop 2 = 4X Midwest, Great Lakes
# pop 3 = Mainly 8X - 84 8X, 10 4X

pop_order <- c(1,2,3)

######

test_samp_orders <- order_struc_samps(clumpp_result_df = c_res, 
  pop_order = pop_order, zero_cut = 0.1)

#####

res_ord <- c_res[test_samp_orders, ]

samp_vec <- rep(res_ord[,1], times = ncol(c_res)-1)
group_vec <- rep(seq(ncol(c_res)-1), each = nrow(res_ord))

res_vec <- c()
for(j in c(2:ncol(res_ord))){
  res_vec <- c(res_vec, res_ord[,j])
}

plot_df <- data.frame(LIB = samp_vec, GROUP = group_vec, MEMB = res_vec,
  stringsAsFactors = F)

plot_df$LIB <- factor(plot_df$LIB, levels = res_ord[,1])

plot_df$GROUP <- factor(plot_df$GROUP, levels = pop_order)

group_col_vec <- c()
group_col_vec['1'] <- 'blue2'
group_col_vec['2'] <- 'steelblue2'
group_col_vec['3'] <- 'darkorchid2'

group_lab_vec <- c()
group_lab_vec['1'] <- '4X Dakota'
group_lab_vec['2'] <- '4X Great Lakes'
group_lab_vec['3'] <- 'Mainly 8X'
group_lab_vec <- group_lab_vec[pop_order]

#group_palette <- scale_colour_manual(name = 'K_Group', values = group_col_vec)

gg_bar <- ggplot(plot_df, aes(y = MEMB, x = LIB, fill = GROUP)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual('Group', values = group_col_vec, 
    labels = group_lab_vec) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = 'Group Membership')

out_file <- gsub('clumpp.processed', 'STRUCTURE.memb.pdf', c_res_file)

pdf(out_file, width = 12*(nrow(c_res)/125), height = 4)
gg_bar
dev.off()

```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo.STRUCTURE.memb.pdf`

########################

## K =4
### Make bargraph of CLUMPP results
```
bash
source activate R_analysis

library(ggplot2)

struc_func_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/struc_functions.r'
source(struc_func_file)

c_res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo.clumpp.processed'

c_res <- read.table(c_res_file, header = F, sep = '\t', stringsAsFactors = F)

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t',
  quote = "", comment.char = '$')

### evaluate what the groups are

pop_info_list <- list()

test_order_mat <- apply(c_res[, -1], 1, order, decreasing = T)

for(i in seq(ncol(c_res)-1)){
  test_pop <- i
  pop_name <- paste('Population', i)
  test_libs <- c_res[which(apply(test_order_mat, 2, 
    function(x) x[1] == test_pop)), 1]
  ploidy_inds <- which(ploidy$lib %in% test_libs)
  meta_inds <- which(meta$LIBRARY %in% test_libs)
  tmp_list <- list()
  tmp_list[['Chloroplast Ecotype']] <- table(ploidy$ECOTYPE_SNP_CHLR
    [ploidy_inds])
  tmp_list[['Old Subpop']] <- table(ploidy$SUBPOP_SNP[ploidy_inds])
  tmp_list[['Ploidy']] <- table(ploidy$total_ploidy_2[ploidy_inds])
  tmp_list[['State']] <- table(meta$STATE[meta_inds])
  pop_info_list[[pop_name]] <- tmp_list
}

pop_info_list

# pop 1 = 32 8X, 15 4X, cosmopolitan
# pop 2 = 4X Midwest, Dakota
# pop 3 = 4X Midwest, Great Lakes
# pop 4 = 8X, MidEast

pop_order <- c(2,3,1,4)

######

test_samp_orders <- order_struc_samps(clumpp_result_df = c_res, 
  pop_order = pop_order, zero_cut = 0.1)

#####

res_ord <- c_res[test_samp_orders, ]

samp_vec <- rep(res_ord[,1], times = ncol(c_res)-1)
group_vec <- rep(seq(ncol(c_res)-1), each = nrow(res_ord))

res_vec <- c()
for(j in c(2:ncol(res_ord))){
  res_vec <- c(res_vec, res_ord[,j])
}

plot_df <- data.frame(LIB = samp_vec, GROUP = group_vec, MEMB = res_vec,
  stringsAsFactors = F)

plot_df$LIB <- factor(plot_df$LIB, levels = res_ord[,1])

plot_df$GROUP <- factor(plot_df$GROUP, levels = pop_order)

group_col_vec <- c()
group_col_vec['1'] <- 'darkorchid2'
group_col_vec['2'] <- 'blue2'
group_col_vec['3'] <- 'steelblue2'
group_col_vec['4'] <- 'turquoise4'

group_lab_vec <- c()
group_lab_vec['1'] <- 'Mainly 8X - Cosmo'
group_lab_vec['2'] <- '4X Dakotas'
group_lab_vec['3'] <- '4X Great Lakes'
group_lab_vec['4'] <- '8X MidEast'
group_lab_vec <- group_lab_vec[pop_order]

#group_palette <- scale_colour_manual(name = 'K_Group', values = group_col_vec)

gg_bar <- ggplot(plot_df, aes(y = MEMB, x = LIB, fill = GROUP)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual('Group', values = group_col_vec, 
    labels = group_lab_vec) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = 'Group Membership')

out_file <- gsub('clumpp.processed', 'STRUCTURE.memb.pdf', c_res_file)

pdf(out_file, width = 12*(nrow(c_res)/125), height = 4)
gg_bar
dev.off()

```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo.STRUCTURE.memb.pdf`






