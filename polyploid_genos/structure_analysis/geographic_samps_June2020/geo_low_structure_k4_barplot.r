# Generate bargraphs for STRUCTURE results on the `geosamps`
##  Includes steps for evaluating the groups and deciding in what order
##    to plot groups in the plot

#bash
#source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)

struc_func_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/struc_functions.r'
# This is the path on HA cluster
source(struc_func_file)

### LOAD INPUTS ###
# all paths are for HA cluster
c_res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo/lowgeo_k4_combo.clumpp.processed'
c_res <- read.table(c_res_file, header = F, sep = '\t', stringsAsFactors = F)

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t',
  quote = "", comment.char = '$')

### SET OUTPUTS ###

out_file <- gsub('clumpp.processed', 'STRUCTURE.memb.pdf', c_res_file)

### SET VARIABLES ###

base_width <- 12
base_height <- 4

# the amount that the number of samples are divided by to scale the width of
#   the plot
width_factor <- 125

tot_width <- base_width * nrow(c_res)/width_factor

##############
### evaluate what the groups are and set group plotting variables

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

# pop 1 = North Eastcoast NY/NJ/RI: 122 4X, 1 8X
# pop 2 = Texas; 131 4X, 3 6X, 4 8X
# pop 3 = Gulfcoast: 68 4X, 3 8X
# pop 4 = South Eastcoast: 147 4X, 1 6X, 3 8X

pop_order <- c(2,3,4,1)

group_col_vec <- c()
group_col_vec['1'] <- 'orange2'
group_col_vec['2'] <- 'red2'
group_col_vec['3'] <- 'salmon2'
group_col_vec['4'] <- 'yellow3'

group_lab_vec <- c()
group_lab_vec['1'] <- 'North East Coast'
group_lab_vec['2'] <- 'Texas'
group_lab_vec['3'] <- 'Gulf Coast'
group_lab_vec['4'] <- 'South East Coast'
group_lab_vec <- group_lab_vec[pop_order]

##########

# re-order samples
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

gg_bar <- ggplot(plot_df, aes(y = MEMB, x = LIB, fill = GROUP)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual('Group', values = group_col_vec, 
    labels = group_lab_vec) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = 'Group Membership')

pdf(out_file, width = tot_width, height = base_height)
gg_bar
dev.off()

quit(save = 'no')


