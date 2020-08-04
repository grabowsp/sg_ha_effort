# Code for generating r^2 figures for sim8X samples

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)
library(patchwork)
library(reshape)

general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)

### LOAD DATA ###

sim_res_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_r2_results'

sim_res_dir <- add_slash(sim_res_dir)

sim_combo_short <- 'combined.sim8X_subgroups.hi_r2.10bpwindow.rds'
sim_combo_file <- paste(sim_res_dir, sim_combo_short, sep = '')

sim_combo_r2 <- readRDS(sim_combo_file)

####

geo_res_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results'

geo_res_dir <- add_slash(geo_res_dir)

geo_combo_short <- 'combined.geo_subgroup.hi_r2.10bpwindow.rds'
geo_combo_file <- paste(geo_res_dir, geo_combo_short, sep = '')

geo_combo_r2 <- readRDS(geo_combo_file)

###############
# Separate panel for each sim8X group

for(i in names(sim_combo_r2)){
  colnames(sim_combo_r2[[i]])[2] <- paste(i, '_per_high_ld', sep = '')
}

combo_df <- sim_combo_r2[[1]][, c(1,2)]
for(j in c(2:length(sim_combo_r2))){
  combo_df[, colnames(sim_combo_r2[[j]][2])] <- sim_combo_r2[[j]][,2]
}

combo_sep_ls <- list()
for(i in seq(length(sim_combo_r2))){
  combo_sep_ls[[i]] <- ggplot(combo_df, aes_string(x = 'window_pos', 
    y = colnames(combo_df)[i+1])) +
  geom_line() + 
  xlab('distance between SNPs') +
  ylab('% r^2 above 0.9') +
  ylim(0,0.5) +
   ggtitle(paste('% r^2 above 0.9 in ', names(sim_combo_r2)[i],
    '\ngenome-wide 10bp windows'))
}

combo_fig_file <- paste(sim_res_dir, 'combo_sim8X_subgroup_hi.r2_fig_1.pdf',
  sep = '')

pdf(combo_fig_file, height = 5*(length(sim_combo_r2)/2), width = 6*2)
wrap_plots(combo_sep_ls, nrow = (length(sim_combo_r2)/2))
dev.off()

##################
# All groups in same panel

sim_combo_r2 <- readRDS(sim_combo_file)

wind_r2_df <- data.frame(window_pos = sim_combo_r2[[1]]$window_pos,
  stringsAsFactors = F)

for(i in seq(length(sim_combo_r2))){
  wind_r2_df[, names(sim_combo_r2)[i]] <- sim_combo_r2[[i]]$per_high_ld
}

wind_melt <- melt(wind_r2_df[, c(2:ncol(wind_r2_df))])

wind_melt <- melt(wind_r2_df, 
 measure.vars = colnames(wind_r2_df)[c(2:ncol(wind_r2_df))])

sim_group_names <- names(sim_combo_r2)

group_col_vec <- c()
group_col_vec['up_tet_1_8X'] <- 'blue1'
group_col_vec['up_tet_2_8X'] <- 'blue2'
#group_col_vec['up_oct_1'] <- 'black'
#group_col_vec['up_oct_2'] <- 'gray30'
group_col_vec['low_tx_1_8X'] <- 'magenta1'
group_col_vec['low_tx_2_8X'] <- 'magenta2'
group_col_vec['low_gc_1_8X'] <- 'green3'
group_col_vec['low_gc_2_8X'] <- 'green4'
group_col_vec['low_ec_1_8X'] <- 'darkorange1'
group_col_vec['low_ec_2_8X'] <- 'darkorange3'

group_col_vec['up_tet_1_up_tet_2_8X'] <- 'blue3'
group_col_vec['low_tx_1_low_tx_2_8X'] <- 'magenta3'
group_col_vec['low_gc_1_low_gc_2_8X'] <- 'green2'
group_col_vec['low_ec_1_low_ec_2_8X'] <- 'darkorange4'

group_col_vec[sim_group_names[intersect(grep('up', sim_group_names), 
  grep('low', sim_group_names))]] <- 'gray30'

group_col_vec[setdiff(sim_group_names, names(group_col_vec))] <- 'red3'

group_palette <- scale_colour_manual(name = 'Group', values = group_col_vec)

gg_tot_groups <- ggplot(wind_melt, aes(x = window_pos, y = value)) +
  geom_line(aes(color = variable)) + group_palette +
  xlab('distance between SNPs') +
  ylab('% r^2 above 0.9') +
  ggtitle('% r^2 above 0.9 in 10bp windows in sim8X subgroups')
  
combo_tot_file <- paste(sim_res_dir, 
  'combo_sim8X_subgroup_total_hi.r2_fig_1.pdf', sep = '')

pdf(combo_tot_file, width = 10, height = 5)
gg_tot_groups
dev.off()

###############
# Panels comparing 4X and 8X populations

sim_combo_r2 <- readRDS(sim_combo_file)

geo_combo_r2 <- readRDS(geo_combo_file)

sim_group_names <- names(sim_combo_r2)

geo_group_names <- names(geo_combo_r2)

tet_groups <- setdiff(geo_group_names, grep('oct', geo_group_names))

group_col_vec <- c()
group_col_vec['up_tet_1'] <- 'blue1'
group_col_vec['up_tet_2'] <- 'blue2'
group_col_vec['up_oct_1'] <- 'black'
group_col_vec['up_oct_2'] <- 'black'
group_col_vec['low_tx_1'] <- 'magenta1'
group_col_vec['low_tx_2'] <- 'magenta2'
group_col_vec['low_gc_1'] <- 'green3'
group_col_vec['low_gc_2'] <- 'green4'
group_col_vec['low_ec_1'] <- 'darkorange1'
group_col_vec['low_ec_2'] <- 'darkorange3'

group_col_vec['up_tet_1_8X'] <- 'blue1'
group_col_vec['up_tet_2_8X'] <- 'blue2'
#group_col_vec['up_oct_1'] <- 'black'
#group_col_vec['up_oct_2'] <- 'gray30'
group_col_vec['low_tx_1_8X'] <- 'magenta1'
group_col_vec['low_tx_2_8X'] <- 'magenta2'
group_col_vec['low_gc_1_8X'] <- 'green3'
group_col_vec['low_gc_2_8X'] <- 'green4'
group_col_vec['low_ec_1_8X'] <- 'darkorange1'
group_col_vec['low_ec_2_8X'] <- 'darkorange3'

group_col_vec['up_tet_1_up_tet_2_8X'] <- 'blue3'
group_col_vec['low_tx_1_low_tx_2_8X'] <- 'magenta3'
group_col_vec['low_gc_1_low_gc_2_8X'] <- 'green2'
group_col_vec['low_ec_1_low_ec_2_8X'] <- 'darkorange4'

group_col_vec[sim_group_names[intersect(grep('up', sim_group_names),
  grep('low', sim_group_names))]] <- 'gray30'

group_col_vec[setdiff(sim_group_names, names(group_col_vec))] <- 'red3'

group_palette <- scale_colour_manual(name = 'Group', values = group_col_vec)

####
group_lin_vec <- c()

group_lin_vec[geo_group_names] <- 'solid'

group_lin_vec[paste(tet_groups, '_8X', sep = '')] <- 'longdash'

group_lin_vec['up_tet_1_up_tet_2_8X'] <- 'dashed'
group_lin_vec['low_tx_1_low_tx_2_8X'] <- 'dashed'
group_lin_vec['low_gc_1_low_gc_2_8X'] <- 'dashed'
group_lin_vec['low_ec_1_low_ec_2_8X'] <- 'dashed'

group_lin_vec[sim_group_names[intersect(grep('up', sim_group_names),
  grep('low', sim_group_names))]] <- 'dotted'

group_lin_vec[setdiff(sim_group_names, names(group_lin_vec))] <- 'dashed'

group_line_set <- scale_linetype_manual(name = 'Group', values = group_lin_vec)
####


sim_tet_gg_list <- list()

for(tg in geo_group_names){
  test_group <- tg
  auto_8X <- paste(tg, '_8X', sep = '')
  tmp_eco <- unlist(strsplit(tg, split = '_'))[1]  
  begin_string <- paste('^', tmp_eco, sep = '')
  begin_inds <- grep(begin_string, sim_group_names)
  mid_string <- paste('_', tmp_eco, '_', sep = '')
  mid_inds <- grep(mid_string, sim_group_names)
  same_eco_8X <- sim_group_names[intersect(begin_inds, mid_inds)]
  up_inds <- grep('up', sim_group_names)
  low_inds <- grep('low', sim_group_names)
  inter_eco_8X <- sim_group_names[intersect(up_inds, low_inds)]
  #
  oct_groups <- c(auto_8X, same_eco_8X, inter_eco_8X)
  #
  wind_r2_df <- data.frame(window_pos = sim_combo_r2[[1]]$window_pos,
    stringsAsFactors = F)
  wind_r2_df[, tg] <- geo_combo_r2[[tg]]$per_high_ld

  for(i in oct_groups){
    wind_r2_df[, i] <- sim_combo_r2[[i]]$per_high_ld
  }
  wind_melt <- melt(wind_r2_df,
    measure.vars = colnames(wind_r2_df)[c(2:ncol(wind_r2_df))])
  sim_tet_gg_list[[tg]] <- ggplot(wind_melt, aes_string(x = 'window_pos',
      y = 'value')) +
    geom_line(aes(color = variable, linetype = variable)) + group_palette +
    group_line_set +
    xlab('distance between SNPs') +
    ylab('% r^2 above 0.9') +
    ylim(0, 0.4) +
    ggtitle('% r^2 above 0.9 in 10bp windows in 4X and associated sim8X subgroups')
}

sim_tet_file <- paste(sim_res_dir,
  'combo_4X_subgroup_withsim8X_subgroup_total_hi.r2_fig_1.pdf', sep = '')

pdf(sim_tet_file, height = 5*5, width = 6*3)
wrap_plots(sim_tet_gg_list, nrow = 5)
dev.off()

quit(save = 'no')
