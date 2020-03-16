# Script for exploring the nQuire results for the switchgrass samples

### LOAD PACKAGES ###
library(ggplot2)
library(reshape2)

### LOAD DATA ###
nquire_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'sg_nquire/sg_CDS_nquire/sg_reseq_nQuire_results_summary_total.txt', 
  sep = '')
nquire_res_0 <- read.table(nquire_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

mnp_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/', 
  'pos_03/pos_03_MNP_results_with_meta_info.txt', sep = '')
mnp_res <- read.table(mnp_res_file, header = T, sep = '\t', 
  stringsAsFactors = F)

### SET OUTPUTS ###

out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/',
  'pos_03/', sep = '')

######

nquire_good_libs <- which(nquire_res_0$samp %in% mnp_res$lib)

nquire_res <- nquire_res_0[nquire_good_libs, ]

#######

nq_melt_1 <- melt(nquire_res[, c(1:4)], id.vars = 'samp')
nq_melt_1$nSNPS_20 <- rep(nquire_res$nSNPS_20, times = 3)

gg_c20_portion <- ggplot(nq_melt_1, aes(x = nSNPS_20, y = value)) + 
  geom_point(aes(color = variable))

tmp_fig_name <- paste(out_dir, 'nquire_c20portionVnSNPs.pdf', sep = '')

pdf(tmp_fig_name)
gg_c20_portion
dev.off()

####

nq_melt_2 <- melt(nquire_res[, c(1, 5:7)], id.vars = 'samp')
nq_melt_2$nSNPS_50 <- rep(nquire_res$nSNPS_50, times = 3)

gg_c50_portion <- ggplot(nq_melt_2, aes(x = nSNPS_50, y = value)) +
  geom_point(aes(color = variable))

tmp_fig_name <- paste(out_dir, 'nquire_c50portionVnSNPs.pdf', sep = '')

pdf(tmp_fig_name)
gg_c50_portion
dev.off()

####

nq_melt_3 <- data.frame(samp = nq_melt_1$samp, variable = nq_melt_1$variable,
  c20_value = nq_melt_1$value, c50_value = nq_melt_2$value, 
  stringsAsFactors = F)

gg_c20vc50_portion <- ggplot(nq_melt_3, aes(x = c20_value, y = c50_value)) + 
  geom_point(aes(color = variable))

tmp_fig_name <- paste(out_dir, 'nquire_c20portVc50port.pdf', sep = '')

pdf(tmp_fig_name)
gg_c20vc50_portion
dev.off()

####

nq_melt_4 <- melt(nquire_res[, 
  c('samp', 'slope_dip', 'slope_tri', 'slope_tet')], id.vars = 'samp')
nq_melt_4$nSNPS_20 <- rep(nquire_res$nSNPS_20, times = 3)
nq_melt_4$nSNPS_50 <- rep(nquire_res$nSNPS_50, times = 3)

gg_slope_c20 <- ggplot(nq_melt_4, aes(x = nSNPS_20, y = value)) + 
  geom_point(aes(color = variable))

tmp_fig_name <- paste(out_dir, 'nquire_slopeVc20nSNPs.pdf', sep = '')

pdf(tmp_fig_name)
gg_slope_c20
dev.off()


####

nq_melt_5 <- data.frame(samp = nq_melt_1$samp, 
  ploidy_variable_1 = nq_melt_1$variable,
  c20_portion = nq_melt_1$value, c50_portion = nq_melt_2$value, 
  slope = nq_melt_4$value, ploidy_variable_2 = nq_melt_4$variable,
  stringsAsFactors = F)

gg_slopeVc20portion <- ggplot(nq_melt_5, aes(x = c20_portion, y = slope)) + 
  geom_point(aes(color = ploidy_variable_2))

tmp_fig_name <- paste(out_dir, 'nquire_slopeVc20portion.pdf', sep = '')

pdf(tmp_fig_name)
gg_slopeVc20portion
dev.off()

####

c50_sep_slope <- (0.0075-(-0.01225))/0.6


gg_slopeVc50port_dip <- ggplot(nquire_res, aes(x = d_dip_portion_50, 
  y = slope_dip)) + 
  geom_point() + 
  geom_abline(intercept = -0.01225, slope = c50_sep_slope)

tmp_fig_name <- paste(out_dir, 'nquire_slopeVc50portion_diploid.pdf', sep = '')

pdf(tmp_fig_name)
gg_slopeVc50port_dip
dev.off()

gg_slopeVc50port_tet <- ggplot(nquire_res, aes(x = d_tet_portion_50, 
  y = slope_tet)) +
  geom_point() + 
  geom_abline(intercept = -0.01, slope = c50_sep_slope)

tmp_fig_name <- paste(out_dir, 'nquire_slopeVc50portion_tetraploid.pdf', 
  sep = '')

pdf(tmp_fig_name)
gg_slopeVc50port_tet
dev.off()

gg_slopeVc50port_tri <- ggplot(nquire_res, aes(x = d_tri_portion_50,
  y = slope_tri)) +
  geom_point() +  
  geom_abline(intercept = -0.011, slope = c50_sep_slope)

tmp_fig_name <- paste(out_dir, 'nquire_slopeVc50portion_triploid.pdf', 
  sep = '')

pdf(tmp_fig_name)
gg_slopeVc50port_tri
dev.off()

gg_slopeVc50portion <- ggplot(nq_melt_5, aes(x = c50_portion, y = slope)) +
  geom_point(aes(color = ploidy_variable_2)) + 
  geom_abline(intercept = -0.01225, slope = c50_sep_slope, color = 'red') + 
  geom_abline(intercept = -0.01, slope = c50_sep_slope, color = 'blue2') + 
  geom_abline(intercept = -0.011, slope = c50_sep_slope, color = 'green3')

tmp_fig_name <- paste(out_dir, 'nquire_slopeVc50portion.pdf', sep = '')

pdf(tmp_fig_name)
gg_slopeVc50portion
dev.off()

####

# assign ploidy using slope and c50 value

nquire_res$slope_dip_call <- NA

non_dip_inds <- which(nquire_res$slope_dip < 
  (nquire_res$d_dip_portion_50 * c50_sep_slope - 0.01225))
dip_inds <- which(nquire_res$slope_dip > 
  (nquire_res$d_dip_portion_50 * c50_sep_slope - 0.01225))

nquire_res$slope_dip_call[dip_inds] <- 'dip'
nquire_res$slope_dip_call[non_dip_inds] <- 'not_dip'
 
nquire_res$slope_tet_call <- NA
non_tet_inds <- which(nquire_res$slope_tet < 
  (nquire_res$d_tet_portion_50 * c50_sep_slope - 0.01))
tet_inds <- which(nquire_res$slope_tet >
  (nquire_res$d_tet_portion_50 * c50_sep_slope - 0.01))

nquire_res$slope_tet_call[tet_inds] <- 'tet'
nquire_res$slope_tet_call[non_tet_inds] <- 'not_tet'


nquire_res$slope_tri_call <- NA
non_tri_inds <- which(nquire_res$slope_tri <
  (nquire_res$d_tri_portion_50 * c50_sep_slope - 0.011))
tri_inds <- which(nquire_res$slope_tri >
  (nquire_res$d_tri_portion_50 * c50_sep_slope - 0.01))

nquire_res$slope_tri_call[tri_inds] <- 'tri'
nquire_res$slope_tri_call[non_tri_inds] <- 'not_tri'

length(intersect(dip_inds, tet_inds))
# [1] 0

length(intersect(non_dip_inds, non_tet_inds))
# [1] 0

length(intersect(tet_inds, tri_inds))
# [1] 140

length(intersect(dip_inds, tri_inds))
# [1] 4

nquire_res$samp[intersect(non_dip_inds, non_tet_inds)]
# [1] "IENM" "IIBU" "IICE"
## these are 3 of the 4 weird samples with all the private alleles in the
##   the Gulf Coast group that had to be removed
## [this is leftover from analysis with all 1034 libraries]

nquire_res$samp[intersect(dip_inds, tri_inds)]
# there's nothing notable about these samples

####
# assign ploidy based on c20 call

c20_dip_inds <- which(apply(nquire_res[, 
  c('d_dip_portion_20', 'd_tri_portion_20', 'd_tet_portion_20')], 1, 
  function(x) min(x) == x[1]))

c20_tri_inds <- which(apply(nquire_res[,
  c('d_dip_portion_20', 'd_tri_portion_20', 'd_tet_portion_20')], 1, 
  function(x) min(x) == x[2]))

c20_tet_inds <- which(apply(nquire_res[,
  c('d_dip_portion_20', 'd_tri_portion_20', 'd_tet_portion_20')], 1,
  function(x) min(x) == x[3]))

#

c50_dip_inds <- which(apply(nquire_res[,
  c('d_dip_portion_50', 'd_tri_portion_50', 'd_tet_portion_50')], 1, 
  function(x) min(x) == x[1]))

c50_tri_inds <- which(apply(nquire_res[,
  c('d_dip_portion_50', 'd_tri_portion_50', 'd_tet_portion_50')], 1,
  function(x) min(x) == x[2]))

c50_tet_inds <- which(apply(nquire_res[,
  c('d_dip_portion_50', 'd_tri_portion_50', 'd_tet_portion_50')], 1,
  function(x) min(x) == x[3]))

nquire_res <- cbind(nquire_res, mnp_res[, -1])

# Note: I looked at nSNPS_20 and there didn't seem to be any pattern on 
#  matching/mismatching calss at c20 and c50

length(intersect(c20_dip_inds, c50_dip_inds))
# 552
## a brief look at metadata looks like these are called 4X
table(nquire_res$PLOIDY[intersect(c20_dip_inds, c50_dip_inds)])
#  4X  8X 
# 547   5
summary(nquire_res$cds_mnp_stand[intersect(c20_dip_inds, c50_dip_inds)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  12.38   58.19   65.74   62.95   71.59   86.42
### the SNP ploidy shows 5 8X in these samples, but the MNP samples
###   show only 4x (all under 100)

length(intersect(c20_tri_inds, c50_tri_inds))
# 2
# 1 4X, 1 8X
nquire_res$PLOIDY[intersect(c20_tri_inds, c50_tri_inds)]
# [1] "8X" "4X"
nquire_res$cds_mnp_stand[intersect(c20_tri_inds, c50_tri_inds)]
# [1] 89.16524 72.32954
nquire_res$seq_cov[intersect(c20_tri_inds, c50_tri_inds)]
# [1] 27,658,276,275 34,568,578,412
### The SNP ploidy shows 1 4X and 1 8X, but MNP shows both are 4X

length(intersect(c20_tet_inds, c50_tet_inds))
# 126
## MOST are 8X, but at least a few are called 4X (eg: IIJR, IILC)
table(nquire_res$PLOIDY[intersect(c20_tet_inds, c50_tet_inds)])
#  4X  8X 
#  11 115
summary(nquire_res$cds_mnp_stand[intersect(c20_tet_inds, c50_tet_inds)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  56.36  116.75  135.29  134.18  152.07  340.67

intersect(which(nquire_res$cds_mnp_stand < 100), intersect(c20_tet_inds, c50_tet_inds))
### 13 samples in this group have MNP counts below 100: 7 have low coverage,
###   5 are close to 100, and 2 are 8X Kanlow? (WZGG, XXAY) - if very 
###   recent 8X, then there wouldn't be time for MNPs to arise

## MISMATCHES

length(intersect(c20_dip_inds, c50_tri_inds))
# 171
## most are 4X
table(nquire_res$PLOIDY[intersect(c20_dip_inds, c50_tri_inds)])
#  4X 
# 171 
summary(nquire_res$cds_mnp_stand[intersect(c20_dip_inds, c50_tri_inds)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  31.23   56.38   68.87   66.51   76.23   94.50 
### All these have MNP values below 100
### They all seem to be 4X, so c20 value is accurate

length(intersect(c20_dip_inds, c50_tet_inds))
# 34
table(nquire_res$PLOIDY[intersect(c20_dip_inds, c50_tet_inds)])
# 4X 
# 34
summary(nquire_res$cds_mnp_stand[intersect(c20_dip_inds, c50_tet_inds)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  29.13   67.09   74.93   73.28   83.85  110.43

intersect(which(nquire_res$cds_mnp_stand > 100), intersect(c20_dip_inds, c50_tet_inds))
# 510 [IIMA]
# in metadata, this sample is flagged as weird
### All of these samples except for 1 seem to be 4X, and the one exception
###   seems to be a strange sample


length(intersect(c20_tri_inds, c50_dip_inds))
# 1
# is 4X
nquire_res$cds_mnp_stand[intersect(c20_tri_inds, c50_dip_inds)]
# [1] 33.89731
### This sample seems to be 4X, so c50 designation is right...

length(intersect(c20_tri_inds, c50_tet_inds))
# 2
nquire_res$PLOIDY[intersect(c20_tri_inds, c50_tet_inds)]
# [1] "4X" "4X"
nquire_res$cds_mnp_stand[intersect(c20_tri_inds, c50_tet_inds)]
# [1] 74.07652 85.52358
### both samples are 4X, so both designations are wrong...

length(intersect(c20_tet_inds, c50_dip_inds))
# 4
nquire_res$PLOIDY[intersect(c20_tet_inds, c50_dip_inds)]
# [1] "4X" "8X" "8X" "8X"
nquire_res$cds_mnp_stand[intersect(c20_tet_inds, c50_dip_inds)]
# [1]  61.85693 118.08364 130.59325 155.91916
### one sample is 4X and the rest at 8X

length(intersect(c20_tet_inds, c50_tri_inds))
# 9
table(nquire_res$PLOIDY[intersect(c20_tet_inds, c50_tri_inds)])
# 4X 8X 
#  1  8
summary(nquire_res$cds_mnp_stand[intersect(c20_tet_inds, c50_tri_inds)])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   72.85  110.94  117.89  120.73  133.60  167.07 

intersect(which(nquire_res$cds_mnp_stand < 100), intersect(c20_tet_inds, c50_tri_inds))
# [1] 128 874

## TAKE HOME ##
# Using the c20 designations seems to be most accurate
## there are only a handful of c20 designations that don't make sense with
##   the MNP data


# some exclued, 1 4X, rest 8X
summary(nquire_res$nSNPS_20[intersect(c20_tet_inds, c50_tri_inds)])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  59912  700852  919760  904330 1212384 1511376

