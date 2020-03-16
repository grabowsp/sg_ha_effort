# Script for analysis and plots for explaining calling ploidy

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)

### LOAD DATA ###
ploidy_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v2.0.txt'
ploidy_res <- read.table(ploidy_res_file, header = T, stringsAsFactors = F,
  sep = '\t')

c30_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire/CDS_nquire_c30_results_total.txt'
c30_res <- read.table(c30_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

c40_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire/CDS_nquire_c40_results_total.txt'
c40_res <- read.table(c40_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

### SET OUTPUTS ###
out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/', 
  'ploidy_calling_figs/', sep = '')

### SET VARIABLES ###

fig_width = 4.5
fig_height = 3

#####################################

### Tables of Values

table(ploidy_res$total_ploidy)
#  4X  8X 
# 760 141

table(ploidy_res$tot_ploid_confidence)
#   1   2   3   4 
#  15  25  86 775

table(ploidy_res$tot_ploid_confidence[ploidy_res$total_ploidy == '4X'])
table(ploidy_res$tot_ploid_confidence[ploidy_res$total_ploidy == '4X'])/sum(ploidy_res$total_ploidy == '4X') * 100
#   1   2     3     4 
#   5   7    72   676
# 0.66 0.92  9.5  88.9

table(ploidy_res$tot_ploid_confidence[ploidy_res$total_ploidy == '8X'])
table(ploidy_res$tot_ploid_confidence[ploidy_res$total_ploidy == '8X'])/sum(ploidy_res$total_ploidy == '8X') * 100
#  1    2      3      4 
# 10   18     14     99
# 7.1  12.8  9.93   70.2

table(ploidy_res$PLOIDY[intersect(ploidy_res$total_ploidy == '4X', ploidy_res$tot_ploid_confidence == '4')])

##########
# plots of nquire_d_dip_portion_20 vs c20 nSNPs

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_dip_portion_20) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 20', 
    y = 'nQuire diploid\nvalue at coverage = 20')

tmp_fig_name <- paste(out_dir, 'nQuire_c20dip_vs_c20nSNPS.pdf', sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_dip_portion_20) +
  geom_point(aes(color = PLOIDY)) +
  labs(x = 'number of SNPS at coverage = 20', 
    y = 'nQuire diploid\nvalue at coverage = 20', color = 'Previous Ploidy')

tmp_fig_name <- paste(out_dir, 'nQuire_c20dip_vs_c20nSNPS_bySujan_ploidy.pdf', 
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

# Plots of nquire_d_tet_portion_20 vs c20 nSNPs

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_tet_portion_20) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 20', 
    y = 'nQuire tetraploid\nvalue at coverage = 20')

tmp_fig_name <- paste(out_dir, 'nQuire_c20tet_vs_c20nSNPS.pdf', sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_tet_portion_20) +
  geom_point(aes(color = PLOIDY)) +
  labs(x = 'number of SNPS at coverage = 20',
    y = 'nQuire tetraploid\nvalue at coverage = 20', color = 'Previous Ploidy')

tmp_fig_name <- paste(out_dir, 'nQuire_c20tet_vs_c20nSNPS_bySujan_ploidy.pdf', 
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

##############################

# Add info about c30 and c40 results

c30_res$lib <- gsub('_c30_q30-bedcc.bin', '', c30_res$file, fixed = T)

ploidy_res$nSNPS_30 <- NA
ploidy_res$d_dip_portion_30 <- NA
ploidy_res$d_tet_portion_30 <- NA

for(LIB in seq(nrow(ploidy_res))){
  tmp_ind <- which(c30_res$lib == ploidy_res$lib[LIB])
  ploidy_res$nSNPS_30[LIB] <- c30_res$nSNPs[tmp_ind]
  ploidy_res$d_dip_portion_30[LIB] <- c30_res$d_dip_portion[tmp_ind]
  ploidy_res$d_tet_portion_30[LIB] <- c30_res$d_tet_portion[tmp_ind]
}

c40_res$lib <- gsub('_c40_q30-bedcc.bin', '', c40_res$file, fixed = T)

ploidy_res$nSNPS_40 <- NA
ploidy_res$d_dip_portion_40 <- NA
ploidy_res$d_tet_portion_40 <- NA

for(LIB in seq(nrow(ploidy_res))){
  tmp_ind <- which(c40_res$lib == ploidy_res$lib[LIB])
  ploidy_res$nSNPS_40[LIB] <- c40_res$nSNPs[tmp_ind]
  ploidy_res$d_dip_portion_40[LIB] <- c40_res$d_dip_portion[tmp_ind]
  ploidy_res$d_tet_portion_40[LIB] <- c40_res$d_tet_portion[tmp_ind]
}

ploidy_res$c30_dip_dif <- (ploidy_res$d_dip_portion_30 - 
  ploidy_res$nquire_d_dip_portion_20)

ploidy_res$c30_tet_dif <- (ploidy_res$d_tet_portion_30 - 
  ploidy_res$nquire_d_tet_portion_20)

ploidy_res$c40_dip_dif <- (ploidy_res$d_dip_portion_40 - 
  ploidy_res$nquire_d_dip_portion_20)

ploidy_res$c40_tet_dif <- (ploidy_res$d_tet_portion_40 -
  ploidy_res$nquire_d_tet_portion_20)

ploidy_res$c50_dip_dif <- (ploidy_res$nquire_d_dip_portion_50 -
  ploidy_res$nquire_d_dip_portion_20)

ploidy_res$c50_tet_dif <- (ploidy_res$nquire_d_tet_portion_50 -
  ploidy_res$nquire_d_tet_portion_20)


# Plot of c20-to-c30 difference vs c30 nSNPs

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nSNPS_30, y = c30_dip_dif) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 30',
    y = 'Change in diploid value\nfrom 20 to 30 coverage')

tmp_fig_name <- paste(out_dir, 'nQuire_changeDip_c20c30_vs_c30nSNPS.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nSNPS_30, y = c30_tet_dif) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 30',
    y = 'Change in tetraploid value\nfrom 20 to 30 coverage')

tmp_fig_name <- paste(out_dir, 'nQuire_changeTet_c20c30_vs_c30nSNPS.pdf',
  sep = '')
  
pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

# Plots of c20-to-c40 difference vs c40 nSNPs

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nSNPS_40, y = c40_dip_dif) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 40',
    y = 'Change in diploid value\nfrom 20 to 40 coverage')

tmp_fig_name <- paste(out_dir, 'nQuire_changeDip_c20c40_vs_c40nSNPS.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nSNPS_40, y = c40_tet_dif) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 40',
    y = 'Change in tetraploid value\nfrom 20 to 40 coverage')

tmp_fig_name <- paste(out_dir, 'nQuire_changeTet_c20c40_vs_c40nSNPS.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

# Plots of c20-to-c50 difference vs c50 nSNPs

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_50, y = c50_dip_dif) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 50',
    y = 'Change in diploid value\nfrom 20 to 50 coverage')

tmp_fig_name <- paste(out_dir, 'nQuire_changeDip_c20c50_vs_c50nSNPS.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_50, y = c50_tet_dif) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'number of SNPS at coverage = 50',
    y = 'Change in tetraploid value\nfrom 20 to 50 coverage')

tmp_fig_name <- paste(out_dir, 'nQuire_changeTet_c20c50_vs_c50nSNPS.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

# Plots of diploid and tetraploid vs c20 nSNPs color coded by subpop

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_dip_portion_20) +
  geom_point(aes(color = SUBPOP_SNP)) +
  labs(x = 'number of SNPS at coverage = 20',
    y = 'nQuire diploid\nvalue at coverage = 20')

tmp_fig_name <- paste(out_dir, 'nQuire_c20dip_vs_c20nSNPS_bySubpop.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_tet_portion_20) +
  geom_point(aes(color = SUBPOP_SNP)) +
  labs(x = 'number of SNPS at coverage = 20',
    y = 'nQuire tetraploid\nvalue at coverage = 20')

tmp_fig_name <- paste(out_dir, 'nQuire_c20tet_vs_c20nSNPS_bySubpop.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

##########

### MNP Plots ###

# CDS vs seq output

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = cds_mnp_stand) +
  geom_point(aes(color = PLOIDY)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nCDS MNP Count', color = 'Metadata\nPloidy')

tmp_fig_name <- paste(out_dir, 'CDS_MNP_vs_SeqOutput_by_MetaPloidy.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = cds_mnp_stand) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nCDS MNP Count', color = 'nQuire\nPloidy')

tmp_fig_name <- paste(out_dir, 'CDS_MNP_vs_SeqOutput_by_nQuirePloidy.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

# Genic vs seq outpu

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = genic_mnp_stand) +
  geom_point(aes(color = PLOIDY)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nGenic MNP Count', color = 'Metadata\nPloidy')

tmp_fig_name <- paste(out_dir, 'genic_MNP_vs_SeqOutput_by_MetaPloidy.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = genic_mnp_stand) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nGenic MNP Count', color = 'nQuire\nPloidy')

tmp_fig_name <- paste(out_dir, 'genic_MNP_vs_SeqOutput_by_nQuirePloidy.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

### MNP figs to demonstrate ploidy calling decisions

# MNP variables
cds_cut <- 100
genic_cut <- 400
thresh_per <- 0.1
cov_cut <- 1e10

tmp_gg <- ggplot(ploidy_res) + 
  aes(x = seq_cov, y = cds_mnp_stand) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'Total Sequencing Output', 
    y = 'Standardized\nCDS MNP Count', color = 'nQuire\nPloidy') +
  geom_hline(yintercept = cds_cut) +
  geom_hline(yintercept = (cds_cut + 10), linetype = 'dashed') +
  geom_hline(yintercept = (cds_cut - 10), linetype = 'dashed') +
  geom_vline(xintercept = cov_cut, linetype = 'dashed', color = 'red')

tmp_fig_name <- paste(out_dir, 'CDS_MNP_vs_SeqOutput_for_PloidyCalling.pdf',
  sep = '')

pdf(tmp_fig_name, height = 5, width = 5.5)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = genic_mnp_stand) +
  geom_point(aes(color = nquire_c20_ploidy)) +
  labs(x = 'Total Sequencing Output', 
    y = 'Standardized\nGenic MNP Count', color = 'nQuire\nPloidy') +
  geom_hline(yintercept = genic_cut) +
  geom_hline(yintercept = (genic_cut + 40), linetype = 'dashed') +
  geom_hline(yintercept = (genic_cut - 40), linetype = 'dashed') +
  geom_vline(xintercept = cov_cut, linetype = 'dashed', color = 'red')

tmp_fig_name <- paste(out_dir, 'genic_MNP_vs_SeqOutput_for_PloidyCalling.pdf',
  sep = '')

pdf(tmp_fig_name, height = 5, width = 5.5)
tmp_gg
dev.off()

### Figures showing nQuire results with final ploidy calls and confidence levels

ploidy_res$tot_ploid_confidence <- as.character(ploidy_res$tot_ploid_confidence)

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_dip_portion_20) +
  geom_point(aes(color = total_ploidy, shape = tot_ploid_confidence)) +
  labs(x = 'number of SNPS at coverage = 20',
    y = 'nQuire diploid\nvalue at coverage = 20', 
    color = 'Combined\nPloidy Call', shape = 'Ploidy\nConfidence')

tmp_fig_name <- paste(out_dir, 'nQuire_c20dip_vs_c20nSNPS_byTotPloidy_Confidence.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_dip_portion_20) +
  geom_point(aes(shape = total_ploidy, color = tot_ploid_confidence)) +
  labs(x = 'number of SNPS at coverage = 20',
    y = 'nQuire diploid\nvalue at coverage = 20', 
    shape = 'Combined\nPloidy Call', color = 'Ploidy\nConfidence')

tmp_fig_name <- paste(out_dir, 'nQuire_c20dip_vs_c20nSNPS_byTotPloidy_Confidence_2.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_tet_portion_20) +
  geom_point(aes(color = total_ploidy, shape = tot_ploid_confidence)) +
  labs(x = 'number of SNPS at coverage = 20',
    y = 'nQuire tetraploid\nvalue at coverage = 20', 
    color = 'Combined\nPloidy Call', shape = 'Ploidy\nConfidence')

tmp_fig_name <- paste(out_dir, 'nQuire_c20tet_vs_c20nSNPS_byTotPloidy_Confidence.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = nquire_nSNPS_20, y = nquire_d_tet_portion_20) +
  geom_point(aes(shape = total_ploidy, color = tot_ploid_confidence)) +
  labs(x = 'number of SNPS at coverage = 20',
    y = 'nQuire tetraploid\nvalue at coverage = 20',
    shape = 'Combined\nPloidy Call', color = 'Ploidy\nConfidence')

tmp_fig_name <- paste(out_dir, 'nQuire_c20tet_vs_c20nSNPS_byTotPloidy_Confidence_2.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

### MNP showing ploidy and confidence

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = cds_mnp_stand) +
  geom_point(aes(color = total_ploidy, shape = tot_ploid_confidence)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nCDS MNP Count', color = 'Combined\nPloidy Call',
    shape = 'Ploidy\nConfidence') 
#  geom_hline(yintercept = genic_cut) +
#  geom_hline(yintercept = (genic_cut + 40), linetype = 'dashed') +
#  geom_hline(yintercept = (genic_cut - 40), linetype = 'dashed') +
#  geom_vline(xintercept = cov_cut, linetype = 'dashed', color = 'red')

tmp_fig_name <- paste(out_dir, 'CDS_MNP_vs_SeqOutput_byTotPloidy_Confidence.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = cds_mnp_stand) +
  geom_point(aes(shape = total_ploidy, color = tot_ploid_confidence)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nCDS MNP Count', shape = 'Combined\nPloidy Call',
    color = 'Ploidy\nConfidence') 

tmp_fig_name <- paste(out_dir, 'CDS_MNP_vs_SeqOutput_byTotPloidy_Confidence_2.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = genic_mnp_stand) +
  geom_point(aes(color = total_ploidy, shape = tot_ploid_confidence)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nGenic MNP Count', color = 'Combined\nPloidy Call',
    shape = 'Ploidy\nConfidence') 

tmp_fig_name <- paste(out_dir, 'genic_MNP_vs_SeqOutput_byTotPloidy_Confidence.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()

tmp_gg <- ggplot(ploidy_res) +
  aes(x = seq_cov, y = genic_mnp_stand) +
  geom_point(aes(shape = total_ploidy, color = tot_ploid_confidence)) +
  labs(x = 'Total Sequencing Output',
    y = 'Standardized\nGenic MNP Count', shape = 'Combined\nPloidy Call',
    color = 'Ploidy\nConfidence')

tmp_fig_name <- paste(out_dir, 'genic_MNP_vs_SeqOutput_byTotPloidy_Confidence_2.pdf',
  sep = '')

pdf(tmp_fig_name, height = fig_height, width = fig_width)
tmp_gg
dev.off()





quit(save = 'no')
