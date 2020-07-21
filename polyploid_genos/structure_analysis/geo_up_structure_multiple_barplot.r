# Generate bargraphs for STRUCTURE results
##  Includes steps for evaluating the groups and deciding in what order
##    to plot groups in the plot

# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)
library(patchwork)

struc_func_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/struc_functions.r'
# * change to appropriate path
source(struc_func_file)

### LOAD INPUTS ###
# all paths are for HA cluster
res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_structure_combined_results.txt'
# * change to appropriate file
res <- read.table(res_file, header = T, sep = '\t', stringsAsFactors = F)

### SET INITIAL VARIABLES ###
k_vec <- c(2,3,4)
# * set the K's to include

# the K that will be used for ordering samples
k_order <- 3
# * set the K used for ording samples

### SET OUTPUTS ###
out_suf <- 'k_2_3_4_v2'
# * set the unique suffix for the figure file
out_suf_full <- paste(out_suf, '_barplots.pdf', sep = '')
out_file <- gsub('combined_results.txt', out_suf_full, res_file)

base_width <- 12
base_height <- 4

tot_height <- base_height * length(k_vec)

# the amount that the number of samples are divided by to scale the width of
#   the plot
width_factor <- 125

tot_width <- base_width * nrow(res)/width_factor

### SET REMAINING VARIABLES ###
group_col_list <- list()
group_col_list[[1]] <- c('darkorchid2', 'blue2')
group_col_list[[2]] <- c('blue2', 'steelblue2', 'darkorchid2')
group_col_list[[3]] <- c('darkorchid2', 'blue2', 'steelblue2', 'turquoise4')
names(group_col_list) <- as.character(k_vec)
# * input the group colors, as determined previously

group_name_list <- list()
group_name_list[[1]] <- c('Mainly 8X','Midwest 4X')
group_name_list[[2]] <- c('4X Dakotas', '4X Great\nLakes', 'Mainly 8X')
group_name_list[[3]] <- c('Mainly 8X:\nCosmo', '4X Dakotas', '4X Great\nLakes',
  '8X MidEast')
names(group_name_list) <- as.character(k_vec)
# * input the group names, as determined previously

pop_order_list <- list()
pop_order_list[[1]] <- c(2,1)
pop_order_list[[2]] <- c(1,2,3)
pop_order_list[[3]] <- c(2,3,1,4)
names(pop_order_list) <- as.character(k_vec)
# * input the pop orders, as determined previously

##############
###############
# make list of potential column names
col_name_list <- list()
for(i in seq(10)){
  col_name_list[[i]] <- paste('K', i, '_', seq(i), sep = '')
}
names(col_name_list) <- as.character(seq(10))

# make df of data with K results that will be used for ordering the samples
test_res <- res[, c('LIB', col_name_list[[as.character(k_order)]])]

# set order of samples
test_samp_orders <- order_struc_samps(clumpp_result_df = test_res, 
  pop_order = pop_order_list[[as.character(k_order)]], zero_cut = 0.1)

res_ord <- res[test_samp_orders, ]
########

barplot_list <- list()

for(i in k_vec){
  tmp_colnames <- col_name_list[[as.character(i)]]
  tmp_samp_vec <- rep(res_ord$LIB, times = length(tmp_colnames))
  tmp_group_vec <- rep(seq(length(tmp_colnames)), each = nrow(res_ord))
  tmp_res_vec <- c()
  for(j in tmp_colnames){
    tmp_res_vec <- c(tmp_res_vec, res_ord[, j])
  }
  tmp_plot_df <- data.frame(LIB = tmp_samp_vec, GROUP = tmp_group_vec, 
    MEMB = tmp_res_vec, stringsAsFactors = F)
  tmp_plot_df$LIB <- factor(tmp_plot_df$LIB, levels = res_ord$LIB)
  tmp_plot_df$GROUP <- factor(tmp_plot_df$GROUP, 
    levels = pop_order_list[[as.character(i)]])
  #
  tmp_group_col_vec <- group_col_list[[as.character(i)]]
  names(tmp_group_col_vec) <- as.character(seq(i))
  tmp_group_name_vec <- group_name_list[[as.character(i)]]
  names(tmp_group_name_vec) <- as.character(seq(i))
  tmp_group_name_vec <- tmp_group_name_vec[pop_order_list[[as.character(i)]]]
  #
  barplot_list[[as.character(i)]] <- ggplot(tmp_plot_df, 
      aes_string(x = 'LIB', y = 'MEMB', fill = 'GROUP')) +
    geom_bar(position = 'fill', stat = 'identity') +
    scale_fill_manual(name = 'Group', 
      values = tmp_group_col_vec, 
      labels = tmp_group_name_vec) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = 'Group Membership')
}

pdf(out_file, height = tot_height, width = tot_width)
wrap_plots(barplot_list, nrow = length(k_vec))
dev.off()

quit(save = 'no')


