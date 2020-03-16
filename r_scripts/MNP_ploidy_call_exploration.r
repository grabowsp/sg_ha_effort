# Call ploidy from MNP results

### LOAD PACKAGES ###
library(ggplot2)
library(reshape2)

### INPUT DATA ###
mnp_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v1.0.txt'
mnp_res <- read.table(mnp_res_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUTS ###
fig_out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_03/', sep = '')
plot_height <- 5
plot_width <- 5.5

### SET VARIABLES ###
cds_cut <- 100
genic_cut <- 400
thresh_per <- 0.1
cov_cut <- 1e10

##################
# Figures to visualize the MNP genotype calling results results

gg_temp <- ggplot(mnp_res) +
  geom_point(aes(x = seq_cov, y = cds_mnp_stand, color = nquire_c20_ploidy,
    shape = PLD_SNP)) +
  geom_hline(yintercept = cds_cut) +
  geom_hline(yintercept = (cds_cut + 10), linetype = 'dashed') +
  geom_hline(yintercept = (cds_cut - 10), linetype = 'dashed')

tmp_fig_name <- paste(fig_out_dir, 'cds_standMNP_vs_cov_nquireploidy.pdf', 
  sep = '')

pdf(tmp_fig_name, height = plot_height, width = plot_width)
gg_temp
dev.off()

gg_temp <- ggplot(mnp_res) +
  geom_point(aes(x = seq_cov, y = cds_mnp_stand, 
    color = nquire_c20cluster_ploidy, shape = PLD_SNP)) +
  geom_hline(yintercept = cds_cut) +
  geom_hline(yintercept = (cds_cut + 10), linetype = 'dashed') +
  geom_hline(yintercept = (cds_cut - 10), linetype = 'dashed')

tmp_fig_name <- paste(fig_out_dir, 'cds_standMNP_vs_cov_nquireploidy_2.pdf', 
  sep = '')

pdf(tmp_fig_name, height = plot_height, width = plot_width)
gg_temp
dev.off()

gg_temp <- ggplot(mnp_res) +
  geom_point(aes(x = seq_cov, y = genic_mnp_stand, color = nquire_c20_ploidy,
    shape = PLD_SNP)) +
  geom_hline(yintercept = genic_cut) +
  geom_hline(yintercept = (1+thresh_per)*genic_cut, linetype = 'dashed') +
  geom_hline(yintercept = (1-thresh_per)*genic_cut, linetype = 'dashed')

tmp_fig_name <- paste(fig_out_dir, 'genic_standMNP_vs_cov_nquireploidy.pdf', 
  sep = '')

pdf(tmp_fig_name, height = plot_height, width = plot_width)
gg_temp
dev.off()

gg_temp <- ggplot(mnp_res) +
  geom_point(aes(x = seq_cov, y = genic_mnp_stand, 
    color = nquire_c20cluster_ploidy, shape = PLD_SNP)) +
  geom_hline(yintercept = genic_cut) +
  geom_hline(yintercept = (1+thresh_per)*genic_cut, linetype = 'dashed') +
  geom_hline(yintercept = (1-thresh_per)*genic_cut, linetype = 'dashed')

tmp_fig_name <- paste(fig_out_dir, 'genic_standMNP_vs_cov_nquireploidy_2.pdf',
  sep = '')

pdf(tmp_fig_name, height = plot_height, width = plot_width)
gg_temp
dev.off()

####

gg_temp <- ggplot(mnp_res) +
  geom_point(aes(x = seq_cov, y = cds_mnp_stand, color = PLD_SNP, 
    shape = cds_mnp_ploidy)) +
  geom_hline(yintercept = cds_cut) +
  geom_hline(yintercept = (cds_cut + 10), linetype = 'dashed') + 
  geom_hline(yintercept = (cds_cut - 10), linetype = 'dashed')

tmp_fig_name <- paste(fig_out_dir, 'cds_standMNP_vs_cov_ploidy.pdf', sep = '')

pdf(tmp_fig_name, height = plot_height, width = plot_width)
gg_temp
dev.off()

###

gg_temp <- ggplot(mnp_res) +
  geom_point(aes(x = seq_cov, y = genic_mnp_stand, color = PLD_SNP, 
    shape = genic_mnp_ploidy)) +
  geom_hline(yintercept = genic_cut) +
  geom_hline(yintercept = (genic_cut + 40), linetype = 'dashed') + 
  geom_hline(yintercept = (genic_cut - 40), linetype = 'dashed')

tmp_fig_name <- paste(fig_out_dir, 'genic_standMNP_vs_cov_ploidy.pdf', sep = '')

pdf(tmp_fig_name, height = plot_height, width = plot_width)
gg_temp
dev.off()

###




quit(save = 'no')
