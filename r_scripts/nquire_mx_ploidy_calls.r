# Script for exploring the nQuire results for the switchgrass samples

# Ploidy designation
# c20_ploidy = ploidy with lowest proportion
# c20_clust_ploidy = ploidy calls based on (visually) clustering the
#   samples from the c20 proportions and for slope_vs_c50 proportions; these
#   clusters end up being identicle

# module load python/3.7-anaconda-2019.07
# source activate R_analysis

### LOAD PACKAGES ###
library(ggplot2)
library(reshape2)
library(patchwork)

### LOAD DATA ###
nquire_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/', 
  'mx_CDS_nquire/mx_nQuire_results_summary_total.txt', sep = '')
nquire_res_0 <- read.table(nquire_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

samtools_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/',
  'mx_CDS_nquire/Mexican_samples_samtools_stats', sep = '')
samtools_res <- read.table(samtools_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)
for(i in c(3:8)){
  samtools_res[, i] <- as.numeric(gsub(',', '', samtools_res[,i]))
}

struc_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/',
  'mx_CDS_nquire/Mexican_Diversity_structure_combined.txt', sep = '')
struc_res <- read.table(struc_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

struc_res <- struc_res[, seq(which(colnames(struc_res) == 'C3_2.74'))]

### SET OUTPUTS ###

out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/',
  'mx_CDS_nquire/', sep = '')

tot_res_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/',
  'mx_CDS_nquire/mx_nquire_and_Sujan_tot_v1.0.txt', sep = '')

### SET VARIABLES ###

fig_height <- 4
fig_width <- 4.5

####################

# Make combined dataframe

nquire_res_ord <- order(nquire_res_0$samp)
samtools_ord <- order(samtools_res$LIB)
struc_ord <- order(struc_res$LIB)

sum(nquire_res_0$samp[nquire_res_ord] == struc_res$LIB[struc_ord])

tot_mx_res <- data.frame(nquire_res_0[nquire_res_ord, ], 
  samtools_res[samtools_ord, c(2:ncol(samtools_res))], 
  struc_res[struc_ord, c(2:ncol(struc_res))], stringsAsFactors = F)

c10_mins <- apply(tot_mx_res[, c(2:4)], 1, function(x) which(x == min(x)))
c20_mins <- apply(tot_mx_res[, c(5:7)], 1, function(x) which(x == min(x)))
c40_mins <- apply(tot_mx_res[, c(8:10)], 1, function(x) which(x == min(x)))

tot_mx_res$c20_ploidy <- paste((c20_mins + 1)*2, 'X', sep = '')
tot_mx_res$c10_ploidy <- paste((c10_mins + 1)*2, 'X', sep = '')
tot_mx_res$c40_ploidy <- paste((c40_mins + 1)*2, 'X', sep = '')

write.table(tot_mx_res, file = tot_res_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

#######

# make plots of nSNPs_20 vs Pct_bases_mapped colored by different variables

tot_mx_res$subpop2 <- as.factor(tot_mx_res$subpop2)

gg_1 <- ggplot(tot_mx_res) +
  aes(y = nSNPS_20, x = Pct_bases_mapped) +
  geom_point()

gg_2 <- ggplot(tot_mx_res) +
  aes(y = nSNPS_20, x = Pct_bases_mapped) +
  geom_point(aes(color = subpop2))

gg_3 <- ggplot(tot_mx_res) +
  aes(y = nSNPS_20, x = Pct_bases_mapped) +
  geom_point(aes(color = Mex_notes))

gg_4 <- ggplot(tot_mx_res) +
  aes(y = nSNPS_20, x = Pct_bases_mapped) +
  geom_point(aes(color = c20_ploidy))

tmp_fig_name <- paste(out_dir, 'Mx_nSNPs20_vs_PerBasesMapped_combo.pdf', 
  sep = '')

pdf(tmp_fig_name, height = 8, width = 9)
(gg_1 | gg_2) / (gg_3 | gg_4)
dev.off()


