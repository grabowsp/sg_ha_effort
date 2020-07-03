# Code for generating r^2 figures for geo sampls

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)
library(patchwork)
library(reshape)

general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)

### LOAD DATA ###

res_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results'

res_dir <- add_slash(res_dir)

combo_short <- 'combined.geo_subgroup.mean_r2.10bpwindow.rds'
combo_file <- paste(res_dir, combo_short, sep = '')

combo_r2 <- readRDS(combo_file)

chrom_short <- 'Chromosome.geo_subgroup.mean_r2.10bpwindow.rds'
chrom_file <- paste(res_dir, chrom_short, sep = '')

# I'm going to hard-code this for now
#   a good possibility is putting all values in same data frame, putting
#   individual values as ggplots using 'aes_string', then using
#   wrap_plots(...) to make the plot

###############
# Separate panel for each group

for(i in names(combo_r2)){
  colnames(combo_r2[[i]])[2] <- paste(i, '_window_avg', sep = '')
}

combo_df <- combo_r2[[1]][, c(1,2)]
for(j in c(2:length(combo_r2))){
  combo_df[, colnames(combo_r2[[j]][2])] <- combo_r2[[j]][,2]
}

combo_sep_ls <- list()
for(i in seq(length(combo_r2))){
  combo_sep_ls[[i]] <- ggplot(combo_df, aes_string(x = 'window_pos', 
    y = colnames(combo_df)[i+1])) +
  geom_line() + 
  xlab('distance between SNPs') +
  ylab('mean r^2') +
  ylim(0,0.5) +
   ggtitle(paste('r^2 in ', names(combo_r2)[i],
    '\ngenome-wide 10bp windows'))
}

combo_fig_file <- paste(res_dir, 'combo_geo_subgroup_mean.r2_fig_1.pdf',
  sep = '')

pdf(combo_fig_file, height = 5*5, width = 6*2)
wrap_plots(combo_sep_ls, nrow = 5)
dev.off()

##################
# All groups in same panel

wind_r2_df <- data.frame(window_pos = combo_r2[[1]]$window_pos,
  stringsAsFactors = F)

for(i in seq(length(combo_r2))){
  wind_r2_df[, names(combo_r2)[i]] <- combo_r2[[i]]$window_avg
}

wind_melt <- melt(wind_r2_df[, c(2:ncol(wind_r2_df))])

wind_melt <- melt(wind_r2_df, 
 measure.vars = colnames(wind_r2_df)[c(2:ncol(wind_r2_df))])

group_col_vec <- rep(NA, times = length(unique(wind_melt$variable)))
group_col_vec['up_tet_1'] <- 'blue1'
group_col_vec['up_tet_2'] <- 'blue2'
group_col_vec['up_oct_1'] <- 'black'
group_col_vec['up_oct_2'] <- 'gray30'
group_col_vec['low_tx_1'] <- 'magenta1'
group_col_vec['low_tx_2'] <- 'magenta2'
group_col_vec['low_gc_1'] <- 'green3'
group_col_vec['low_gc_2'] <- 'green4'
group_col_vec['low_ec_1'] <- 'darkorange1'
group_col_vec['low_ec_2'] <- 'darkorange3'

group_palette <- scale_colour_manual(name = 'Group', values = group_col_vec)

gg_tot_groups <- ggplot(wind_melt, aes(x = window_pos, y = value)) +
  geom_line(aes(color = variable)) + group_palette +
  xlab('distance between SNPs') +
  ylab('mean r^2') +
  ggtitle('Mean r^2 in 10bp windows in geo_samps subgroups')
  
combo_tot_file <- paste(res_dir, 'combo_geo_subgroup_total_mean.r2_fig_1.pdf',
  sep = '')

pdf(combo_tot_file, width = 6, height = 5)
gg_tot_groups
dev.off()

#####################
# Chromosomes separately, all groups in same chromosome panel

chrom_r2 <- readRDS(chrom_file)

group_col_vec <- rep(NA, times = length(unique(wind_melt$variable)))
group_col_vec['up_tet_1'] <- 'blue1'
group_col_vec['up_tet_2'] <- 'blue2'
group_col_vec['up_oct_1'] <- 'black'
group_col_vec['up_oct_2'] <- 'gray30'
group_col_vec['low_tx_1'] <- 'magenta1'
group_col_vec['low_tx_2'] <- 'magenta2'
group_col_vec['low_gc_1'] <- 'green3'
group_col_vec['low_gc_2'] <- 'green4'
group_col_vec['low_ec_1'] <- 'darkorange1'
group_col_vec['low_ec_2'] <- 'darkorange3'

group_palette <- scale_colour_manual(name = 'Group', values = group_col_vec)

chrom_names <- names(chrom_r2)

for(chr_nm in names(chrom_r2)){
  for(pop_name in names(chrom_r2[[chr_nm]])){
    colnames(chrom_r2[[chr_nm]][[pop_name]])[2] <- paste(chr_nm, 
      '_window_avg', sep = '')
    chrom_r2[[chr_nm]][[pop_name]]$pop <- pop_name
  }
}

chrom_unified <- lapply(chrom_r2, gen_unified_df)

chr_uni_df <- chrom_unified[[1]][, c(1,4,2)]
for(cn in c(2:length(chrom_names))){
  chr_uni_df[, paste(chrom_names[cn], '_window_avg', sep = '')] <- 
    chrom_unified[[chrom_names[cn]]][,2]
}

chr_plot_ls <- list()

for(i in seq(length(chrom_names))){
  chr_plot_ls[[i]] <- ggplot(chr_uni_df, aes_string(x = 'window_pos', 
    y = colnames(chr_uni_df)[i + 2])) +
  geom_line(aes(color = pop)) + group_palette +
  xlab('distance between SNPs') +
  ylab('mean r^2') + 
  ggtitle(paste('Mean r^2 in 10bp windows in ', chrom_names[i], 
    '\ngeo_samps subgroups'))
#  print(chr_plot_ls[[i]])
}

chr_plot_out <- paste(res_dir, 'chromosome_geo_subgroup_mean.r2_fig_1.pdf', 
  sep = '')

pdf(chr_plot_out, height = 5*9, width = 6*2)
wrap_plots(chr_plot_ls, nrow = 9)
dev.off()


