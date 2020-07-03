# scratch analysis space

#module load python/3.7-anaconda-2019.07
#source activate R_analysis

#source activate r_adegenet_env

general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)

chr_r2_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/Chromosome.geo_subgroup.r2.rds'

combo_r2_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/combined.geo_subgroup.r2.rds'

pop_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt'

pop_libs <- read.table(pop_lib_file, header = T, sep = '\t', 
  stringsAsFactors = F)

tx_2_libs <- pop_libs$LIB[which(pop_libs$POP=='low_tx_2')]
tx_1_libs <- pop_libs$LIB[which(pop_libs$POP=='low_tx_1')]

up4_2_libs <- pop_libs$LIB[which(pop_libs$POP == 'up_tet_2')]
up4_1_libs <- pop_libs$LIB[which(pop_libs$POP == 'up_tet_1')]

up8_1_libs <- pop_libs$LIB[which(pop_libs$POP == 'up_oct_1')]
up8_2_libs <- pop_libs$LIB[which(pop_libs$POP == 'up_oct_2')]

gc_1_libs <- pop_libs$LIB[which(pop_libs$POP=='low_gc_1')]
gc_2_libs <- pop_libs$LIB[which(pop_libs$POP=='low_gc_2')]

ec_1_libs <- pop_libs$LIB[which(pop_libs$POP=='low_ec_1')]
ec_2_libs <- pop_libs$LIB[which(pop_libs$POP=='low_ec_2')]

library(adegenet)
library(ggplot2)
library(reshape)

chr02K_genlight_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr02K.polyploid.CDS.geosamps.genlight.rds'

chr02K_gl <- readRDS(chr02K_genlight_file)

snp_ind <- grep('Chr02K_15271495', locNames(chr02K_gl))

snp_ind <- grep('Chr02K_16092256', locNames(chr02K_gl))

snp_ind_1 <- grep('Chr02K_17954779', locNames(chr02K_gl))
snp_ind_2 <- grep('Chr02K_17955025', locNames(chr02K_gl))
snp_ind_3 <- grep('Chr02K_17956516', locNames(chr02K_gl))
snp_ind_4 <- grep('Chr02K_17965022', locNames(chr02K_gl))

snp_vec <- c(snp_ind_1, snp_ind_2, snp_ind_3, snp_ind_4)

as.matrix(chr02K_gl[, snp_vec])[tx_2_libs, ]

as.matrix(chr02K_gl[, snp_vec])[up4_2_libs, ]

as.matrix(chr02K_gl[, snp_vec])[up4_1_libs, ]

as.matrix(chr02K_gl[, snp_vec])[tx_1_libs, ]

as.matrix(chr02K_gl[, snp_vec])[up8_1_libs, ]

as.matrix(chr02K_gl[, snp_vec])[up8_2_libs, ]

as.matrix(chr02K_gl[, snp_vec])[gc_1_libs, ]

as.matrix(chr02K_gl[, snp_vec])[gc_2_libs, ]

as.matrix(chr02K_gl[, snp_vec])[ec_1_libs, ]
as.matrix(chr02K_gl[, snp_vec])[ec_2_libs, ]

chr_r2 <- readRDS(chr_r2_file)

chr02K_r2 <- chr_r2[[3]]

tx2_r2 <- chr_r2[[3]][['low_tx_2']]

tx2_r2$test_pos <- as.numeric(unlist(lapply(strsplit(tx2_r2$test_snp, 
  split = '_'), function(x) x[2])))

tx2_r2 <- tx2_r2[order(tx2_r2$test_pos), ]

weird_inds <- intersect(which(tx2_r2$comp_dist > 9000), which(tx2_r2$r2 > 0.8))

tx2_weird <- tx2_r2[weird_inds, ]

# want to make a density plot of "weird_inds" along the chromosome to see if
#  there are enriched regions...

gg_r2 <- ggplot(tx2_weird, aes(x = test_pos)) +
  geom_density() +
  xlim(min(tx2_r2$test_pos), max(tx2_r2$test_pos)) +
  xlab('position on chromosome') +
  ggtitle('Density of SNPs > 9kb apart with r^2 > 0.8')

density_plot_out <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/',
  'tx_2_Chr01K_weird_density.pdf', sep = '')

pdf(density_plot_out, width = 12, height = 4)
gg_r2
dev.off()

max_dist <- max(tx2_r2$test_pos)
window_size <- 100000  

dist_mult_vec <- seq(ceiling(max_dist/window_size))

wind_dists <- lapply(dist_mult_vec, function(x) 
 ((x-1)*window_size) + c(1:window_size))

window_inds <- lapply(wind_dists, function(x)
    which(tx2_weird$test_pos %in% x))

num_window_inds <- unlist(lapply(window_inds, length))

summary(num_window_inds)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   0.000    0.000    0.000    7.785    1.000 3619.000 
which(num_window_inds > 100)
#[1] 121 180
which(num_window_inds > 1000)
#[1] 180
179*window_size
#[1] 17900000
sum((tx2_weird$test_pos > 17900000) & (tx2_weird$test_pos < 18000000))
#[1] 3619
weird_interval <- which((tx2_weird$test_pos > 17900000) & 
  (tx2_weird$test_pos < 18000000))

tot_weird_int <- which((tx2_r2$test_pos > 17954777) & 
  (tx2_r2$test_pos < 17966732))

summary(tx2_r2$r2[tot_weird_int])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.7639  1.0000  0.8342  1.0000  1.0000 
length(unique(tx2_r2$test_pos[tot_weird_int]))
#[1] 234


#######

pop_names <- names(chr02K_r2)

for(pn in pop_names){
  chr02K_r2[[pn]]$pop <- pn
}

tot_chr02K_df <- gen_unified_df(chr02K_r2)

tot_chr02K_df$test_pos <- as.numeric(unlist(lapply(strsplit(
  tot_chr02K_df$test_snp, split = '_'), function(x) x[2])))

tot_chr02K_df$pop <- as.factor(tot_chr02K_df$pop)

weird_df_inds <- intersect(which(tot_chr02K_df$comp_dist > 9000), 
  which(tot_chr02K_df$r2 > 0.8))

tot_weird_df <- tot_chr02K_df[weird_df_inds,]

group_col_vec <- rep(NA, times = length(unique(tot_weird_df$pop)))
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

group_col_vec <- group_col_vec[-which(is.na(group_col_vec))]

group_palette <- scale_colour_manual(name = 'Group', values = group_col_vec)

gg_r2_tot <- ggplot(tot_weird_df, aes(x = test_pos)) +
  geom_density(aes(color = pop)) + group_palette +
  xlim(min(tot_chr02K_df$test_pos), max(tot_chr02K_df$test_pos)) +
  xlab('position on chromosome') +
  ggtitle('Density of SNPs > 9kb apart with r^2 > 0.8')

tot_density_plot_out <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/',
  'allpops_Chr01K_weird_density.pdf', sep = '')

pdf(tot_density_plot_out, width = 12, height = 4)
gg_r2_tot
dev.off()

##########

# next: make figure showing relative patterns for each population on their
#  own scale

pop_names <- names(chr02K_r2)

weird_list <- list()

for(pn in pop_names){
  tmp_data <- chr02K_r2[[as.character(pn)]]
  tmp_data$test_pos <- as.numeric(unlist(lapply(strsplit(
    tmp_data$test_snp, split = '_'), function(x) x[2])))
  tmp_weird_inds <- intersect(which(tmp_data$comp_dist > 9000),
    which(tmp_data$r2 > 0.8))
  tmp_weird <- tmp_data[tmp_weird_inds,]
  gg_tmp <- ggplot(tmp_weird, aes(x = test_pos)) +
    geom_density() +
    xlim(min(tmp_data$test_pos), max(tmp_data$test_pos)) + 
    xlab('position on chromosome') + 
    ggtitle(paste(pn, ' Density of SNPs > 9kb apart with r^2 > 0.8', sep = ''))
  weird_list[[pn]] <- gg_tmp
}

pop_density_plot_out <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/',
  'pop_split_Chr01K_weird_density.pdf', sep = '')

pdf(pop_density_plot_out, height = 3*10, width = 12)
wrap_plots(weird_list, nrow = 10)
dev.off()

########

# Figure with up_tet_2, low_tx_2, and up_tet_1 and low_tx_1 as controls

sub_weird_list <- list

sub_pops <- c('low_tx_2', 'up_tet_2', 'low_tx_1', 'up_tet_1')

group_col_vec <- rep(NA, times = length(pop_names))
names(group_col_vec) <- pop_names
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

sub_list <- list()

for(pn in sub_pops){
  tmp_data <- chr02K_r2[[as.character(pn)]]
  tmp_data$test_pos <- as.numeric(unlist(lapply(strsplit(
    tmp_data$test_snp, split = '_'), function(x) x[2])))
  tmp_weird_inds <- intersect(which(tmp_data$comp_dist > 9000),
    which(tmp_data$r2 > 0.8))
  tmp_weird <- tmp_data[tmp_weird_inds,]
  n_comps <- nrow(tmp_weird)
  gg_tmp <- ggplot(tmp_weird, aes(x = test_pos)) +
    geom_density(color = group_col_vec[as.character(pn)]) +
#    scale_color_manual(values = group_col_vec[pn]) + 
    xlim(min(tmp_data$test_pos), max(tmp_data$test_pos)) +
    xlab('position on chromosome') +
    ggtitle(paste(pn, ' Density of SNPs > 9kb apart with r^2 > 0.8, ',
      n_comps, ' comps fit criteria', sep = ''))
  sub_list[[pn]] <- gg_tmp
}

subpop_density_plot_out <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/',
  '4_pop_split_Chr01K_weird_density.pdf', sep = '')

pdf(subpop_density_plot_out, height = 3*4, width = 12)
wrap_plots(sub_list, nrow = 4)
dev.off()




