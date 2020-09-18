# Goals
# Evaluate if the number of fixed HOM-ALT and HET SNPs is consistent
#   with 8X switchgrass being auto-allopolyploid or allo-allopolyploid

# Final thought after a lot of time looking at results: The data are too 
#  messy for this analysis. 
#  There are a bunch of SNPs that show weird patterns, like being 
#  fixed-het in some/many subsets of samples
#  These are most likely assembly/mapping artifacts, like CNVs
#  Unfortunately, they add too much noise to try to model what's going on
#  with fixed-het SNPs in 8X populations, because are looking for new 
#  fixed-het SNPs, but there are too many artifacts.


# Goals
# Calculate number of Upland-fixed HOM-ALT
# Calculate number of Midwest-4X fixed HOM-ALT
# Calculate number of 8X fixed HOM-ALT
# Calculate estimated/predicted number of 8X fixed HOM-ALT
#   Not sure how to do that yet
# Calculate number of Midwest-4X fixed HET as error estimate
# Calculate number of 8X fixed HET
# Calculate prediced number of 8X fixed HET
#  based on number of Midwest-4X fixed HOM-ALT

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/R_analysis

### INPUT DATA ###

#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/'
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo/'

#vcf_in_short <- 'Chr01K.polyploid.CDS.expandgeosamps.vcf_00'
vcf_in_short <- 'Chr01K_geoSim8Xmerge_tet.vcf_00'

vcf_in <- paste(data_dir, vcf_in_short, sep = '')

vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

# Expand VCF to include more SNPs
subfile_vec <- c('01', '02', '03', '04')
for(sfv in subfile_vec){
  tmp_short <- paste('Chr01K.polyploid.CDS.upexpand.vcf_', sfv, sep = '')
  tmp_file_in <- paste(data_dir, tmp_short, sep = '')
  tmp_vcf <- read.table(tmp_file_in, header = F, stringsAsFactors = F, 
    sep = '\t')
  vcf <- rbind(vcf, tmp_vcf)
}

# Add header
#head_in_short <- 'CDS.expandgeosamps.vcf.header.txt'
vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo/geoSim8Xmerg_samps.txt'

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf) <- vcf_header

# Load library names of 4X and 8X samples
tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs <- intersect(tet_libs_0, colnames(vcf))

oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs <- intersect(oct_libs_0, colnames(vcf))

# Load file with upland/Midwest samples based on PCA
upland_samp_name_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/expandupland_libs_July2020.txt'
upland_samps <- read.table(upland_samp_name_file, header = F, 
  stringsAsFactors = F)[,1]

# Load file with lowland samples based on PCA
lowland_samp_name_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/lowland_libs_June2020.txt'
lowland_samps <- read.table(lowland_samp_name_file, header = F,
  stringsAsFactors = F)[,1]

# Load info about the simulated 8X samples
sim8X_samp_info_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_samp_sim8X_lib_combos.txt'
sim8X_samp_info <- read.table(sim8X_samp_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

# Load info about subpopulations used for r2 and simulated 8X construction
subpop_samp_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt', sep = '')
subpop_samp_info <- read.table(subpop_samp_file, header = T, sep = '\t', 
  stringsAsFactors = F)

######################

# 1) Set up matrices with different samples to be used for calculating
#     fixed genotype rates
#  - convert missing genotypes to NA's
#  - remove SNPs with high missing data across MW

## all MW, all MW-4X and all MW-8X
up_tet <- intersect(upland_samps, tet_libs)
up_oct <- intersect(upland_samps, oct_libs)

up_tot <- c(up_tet, up_oct)

up_df <- vcf[, up_tot]
up_df[up_df == './.'] <- NA

# find snps with hi missing data in uplands
up_na_bysnp <- apply(up_df, 1, function(x) sum(is.na(x)))
up_hi_miss <- which((up_na_bysnp/length(up_tot)) > 0.2)

# find snps with no alt alleles in uplands
up_n_altsamps <- apply(up_df, 1, function(x) sum((x == '3/1'| x == '2/2' | 
  x == '1/3' | x == '0/4'), na.rm = T))
up_no_alts <- which(up_n_altsamps == 0)

# find snps where the ALT allele is actually the major allele (won't work
#   for this analysis)
# need to use both upland and lowland samples for this
low_df <- vcf[, lowland_samps]
low_df[low_df == './.'] <- NA

up_tmp_a <- apply(up_df, 1, function(x) sum(x == '4/0', na.rm = T))
up_tmp_b <- apply(up_df, 1, function(x) sum(x == '0/4', na.rm = T))
up_alt_maj_snps <- which(up_tmp_b > up_tmp_a)

low_tmp_a <- apply(low_df, 1, function(x) sum(x == '4/0', na.rm = T))
low_tmp_b <- apply(low_df, 1, function(x) sum(x == '0/4', na.rm = T))
low_alt_maj_snps <- which(low_tmp_b > low_tmp_a)

alt_maj_snps <- intersect(up_alt_maj_snps, low_alt_maj_snps)

# decided to use all SNPs for the simulated and subpop analyses
#  will deflate the expected number of fixed HomALT and HET though
# bad_inds <- sort(union(up_hi_miss, up_no_alts))
# bad_inds <- up_hi_miss
bad_inds <- sort(union(up_hi_miss, alt_maj_snps))

up_tet_df <- vcf[-bad_inds, up_tet]
up_oct_df <- vcf[-bad_inds, up_oct]

up_tet_df[up_tet_df == './.'] <- NA
up_oct_df[up_oct_df == './.'] <- NA

## Simulated 8X groups
sim8X_groups <- unique(sim8X_samp_info$sim_pop)
sim8X_up_groups <- sim8X_groups[grep('up', sim8X_groups)]
sub_up_groups <- c('up_tet_1_8X', 'up_tet_2_8X')

sim8X_vcf_list <- list()
for(sug in sim8X_up_groups){
  tmp_samps <- sim8X_samp_info$sim_samp_name[
    which(sim8X_samp_info$sim_pop == sug)]
  tmp_vcf <- vcf[-bad_inds, tmp_samps]
  tmp_vcf[tmp_vcf == './.'] <- NA
  sim8X_vcf_list[[sug]] <- tmp_vcf
}

# 4X samples used to construct simulated 8X populations
up_sub_vcf_list <- list()
for(subp in sub_up_groups){
  tmp_name <- sub('8X', '4X', subp)
  tmp_inds <- which(sim8X_samp_info$sim_pop == subp)
  tmp_samps <- unique(c(sim8X_samp_info$keep_1[tmp_inds], 
    sim8X_samp_info$keep_2[tmp_inds]))
  tmp_vcf <- vcf[-bad_inds, tmp_samps]
  tmp_vcf[tmp_vcf == './.'] <- NA
  up_sub_vcf_list[[tmp_name]] <- tmp_vcf
}

## Subpopulations
subpop_groups <- unique(subpop_samp_info$POP)
sub_oct_groups <- c('up_oct_1', 'up_oct_2')

subpop_vcf_list <- list()
for(subpop in subpop_groups){
  tmp_samps <- subpop_samp_info$LIB[which(subpop_samp_info$POP == subpop)]
  tmp_vcf <- vcf[-bad_inds, tmp_samps]
  tmp_vcf[tmp_vcf == './.'] <- NA
  subpop_vcf_list[[subpop]] <- tmp_vcf
}

# 2) Calculate NA's and "normalization vectors"
# Note: need to add subpop stuff for remaining steps

up_tet_na_bysnp <- apply(up_tet_df, 1, function(x) sum(is.na(x)))
up_tet_norm_snp_vec <- 1/(1-(up_tet_na_bysnp/ncol(up_tet_df)))

up_oct_na_bysnp <- apply(up_oct_df, 1, function(x) sum(is.na(x)))
up_oct_norm_snp_vec <- 1/(1-(up_oct_na_bysnp/ncol(up_oct_df)))

sim8X_na_bysnp_list <- lapply(sim8X_vcf_list, function(x) 
  apply(x, 1, function(y) sum(is.na(y))))
sim8X_norm_snp_list <- lapply(sim8X_na_bysnp_list, function(x)
  1/(1-(x/16)))

sub_na_bysnp_list <- lapply(up_sub_vcf_list, function(x)
  apply(x, 1, function(y) sum(is.na(y))))
sub_norm_snp_list <- list()
for(subp in names(up_sub_vcf_list)){
  n_samps <- ncol(up_sub_vcf_list[[subp]])
  sub_norm_snp_list[[subp]] <- 1/(1-(sub_na_bysnp_list[[subp]]/n_samps))
}

subpop_na_bysnp_list <- lapply(subpop_vcf_list, function(x)
  apply(x, 1, function(y) sum(is.na(y))))
subpop_norm_snp_list <- list()
for(subpop in names(subpop_na_bysnp_list)){
  n_samps <- ncol(subpop_vcf_list[[subpop]])
  subpop_norm_snp_list[[subpop]] <- 1/(
    1-(subpop_na_bysnp_list[[subpop]]/n_samps))
}

# 3) Find/Quantify MW-4X, MW-8X, and MW-tot fixed Homozygous-ALT
# - includes look at levels of missing data to evaluate the effect of
#     missing data on quantities
up_tet_homA_bysnp <- apply(up_tet_df, 1, function(x) sum(x == '0/4', na.rm = T))
up_tet_homA_bysnp_1 <- up_tet_homA_bysnp * up_tet_norm_snp_vec
mw4_homA_fixed_snps <- which(up_tet_homA_bysnp_1 > (ncol(up_tet_df)*0.95))
length(mw4_homA_fixed_snps)
# 118 SNPs 
sum(up_tet_norm_snp_vec[mw4_homA_fixed_snps] > (0.95/0.8))
# 0

up_oct_homA_bysnp <- apply(up_oct_df, 1, function(x) sum(x == '0/4', na.rm = T))
up_oct_homA_bysnp_1 <- up_oct_homA_bysnp * up_oct_norm_snp_vec
up_oct_homA_fixed_snps <- which(up_oct_homA_bysnp_1 > (ncol(up_oct_df)*0.95))
length(up_oct_homA_fixed_snps)
# 53 SNPs
sum(up_oct_norm_snp_vec[up_oct_homA_fixed_snps] > (0.95/0.8))
# 0

up_homA_fixed_snps <- intersect(mw4_homA_fixed_snps, up_oct_homA_fixed_snps)
length(up_homA_fixed_snps)
# 53

mw4_only_homA_fixed_snps <- setdiff(mw4_homA_fixed_snps, up_oct_homA_fixed_snps)
length(mw4_only_homA_fixed_snps)
# 65 SNPs
sum(up_tet_norm_snp_vec[mw4_only_homA_fixed_snps]  > (0.95/0.8))
# 0

up_oct_only_homA_fixed_snps <- setdiff(up_oct_homA_fixed_snps, 
  mw4_homA_fixed_snps)
length(up_oct_only_homA_fixed_snps)
# 0 SNPs
sum(up_oct_norm_snp_vec[up_oct_only_homA_fixed_snps] > (0.95/0.8))
# 0

p_mw4_homA_fixed <- length(mw4_only_homA_fixed_snps)/nrow(up_tet_df)
#[1] 0.0006674676

expected_oct_homA_fixed <- (p_mw4_homA_fixed^2) * nrow(up_oct_df)
# 0.0433854 : the expected number of 8X homA-fixed SNPs if allo-allopolyploid
#   based on the Midwest-4X number of homA-fixed SNPs

# 4) Find/Quantify HOM-ALT in simulated 8X, samples used for simulated 8X, 
#    and subpopulations from r2 analysis
# - includes look at levels of missing data to evaluate the effect of
#     missing data on quantities

# Function to look at missing data
show_missdata_info <- function(fix_ind_list, norm_vec_list){
  for(tn in names(fix_ind_list)){
    print(tn)
    tmp_l <- length(fix_ind_list[[tn]])
    print(tmp_l)
    tmp_hi <- sum(norm_vec_list[[tn]][fix_ind_list[[tn]]] > (0.95/0.8))
    print(tmp_hi/tmp_l)
    tmp_vhi <- sum(norm_vec_list[[tn]][fix_ind_list[[tn]]] > 2)
    print(tmp_vhi/tmp_l)
    print('           ')
  }
}

sim8_homA_list <- lapply(sim8X_vcf_list, function(x) 
  apply(x, 1, function(y) sum(y == '0/4', na.rm = T)))

sim8_homA_1_list <- list()
for(sug in sim8X_up_groups){
  sim8_homA_1_list[[sug]] <- sim8_homA_list[[sug]] * sim8X_norm_snp_list[[sug]]
}

sim8X_homA_fixed_list <- lapply(sim8_homA_1_list, function(x) 
  which(x > (16*0.95)))
sim8X_pop_fixed_homA_list <- lapply(sim8X_homA_fixed_list, function(x) 
  setdiff(x, up_homA_fixed_snps))
# lapply(sim8X_pop_fixed_homA_list, length)
show_missdata_info(sim8X_pop_fixed_homA_list, sim8X_norm_snp_list)
# up_tet_1_8X: 3%/0% of 93 fixed homA with hi/very hi missing data
# up_tet_2_8X: 3.6%/2.7% of 111 with hi/very hi missing data
# up_tet_1_up_tet_2_8X: 0%/0% of 60 with hi/very hi missing data

sub_homA_list <- lapply(up_sub_vcf_list, function(x)
  apply(x, 1, function(y) sum(y == '0/4', na.rm = T)))
sub_homA_1_list <- list()
for(subp in names(sub_homA_list)){
  sub_homA_1_list[[subp]] <- sub_homA_list[[subp]] * sub_norm_snp_list[[subp]]
}

sub_homA_fixed_list <- list()
for(subp in names(sub_homA_1_list)){
  sub_homA_fixed_list[[subp]] <- which(sub_homA_1_list[[subp]] > 
    (ncol(up_sub_vcf_list[[subp]]) * 0.95))
}
sub_pop_fixed_homA_list <- lapply(sub_homA_fixed_list, function(x)
  setdiff(x, up_homA_fixed_snps))
# lapply(sub_pop_fixed_homA_list, length)
# 113 and 122
show_missdata_info(sub_pop_fixed_homA_list, sub_norm_snp_list)
# 2.5%/0% of 113 and 122 fixed homA in sublists show hi/very hi missing data

# length(intersect(sub_pop_fixed_homA_list[[1]], sub_pop_fixed_homA_list[[2]]))
# 82 shared by both subpops


subpop_homA_list <- lapply(subpop_vcf_list, function(x)
  apply(x, 1, function(y) sum(y == '0/4', na.rm = T)))
subpop_homA_1_list <- list()
for(subpop in names(subpop_homA_list)){
  subpop_homA_1_list[[subpop]] <- subpop_homA_list[[subpop]] * 
    subpop_norm_snp_list[[subpop]]
}

subpop_homA_fixed_list <- list()
for(subpop in names(subpop_homA_1_list)){
  subpop_homA_fixed_list[[subpop]] <- which(subpop_homA_1_list[[subpop]] >
    (ncol(subpop_vcf_list[[subpop]]) * 0.95))
}
subpop_pop_homA_list <- lapply(subpop_homA_fixed_list, function(x)
  setdiff(x, up_homA_fixed_snps))
# lapply(subpop_pop_homA_list, length)
show_missdata_info(subpop_pop_homA_list, subpop_norm_snp_list)
# 2%/0% of 98 fixed homALT in up_tet_1 hi/very hi missing data
# 2.5%/0% of 121 in up_tet_2
# 3.2%/0% of 62 in up_oct_1
# 0 in up_oct_2 - kind of odd...

# 5) Find/Quantify fixed HET sites
up_tet_het_bysnp <- apply(up_tet_df, 1, 
  function(x) sum(x == '3/1' | x == '2/2' | x == '1/3', na.rm = T))
up_tet_het_bysnp_1 <- up_tet_het_bysnp * up_tet_norm_snp_vec

up_oct_het_bysnp <- apply(up_oct_df, 1, 
  function(x) sum(x == '3/1' | x == '2/2' | x == '1/3', na.rm = T))
up_oct_het_bysnp_1 <- up_oct_het_bysnp * up_oct_norm_snp_vec

up_tet_hi_het_snps <- which(up_tet_het_bysnp_1 > (ncol(up_tet_df)*0.8))
length(up_tet_hi_het_snps)
# 299
# should be 0 unless there is balancing selection, so these are probably errors

up_oct_hi_het_snps <- which(up_oct_het_bysnp_1 > (ncol(up_oct_df)*0.95))
length(up_oct_hi_het_snps)
# 149
sum(up_oct_norm_snp_vec[up_oct_hi_het_snps] > (0.95/0.8))
# 1

up_fixed_hi_het_snps <- intersect(up_tet_hi_het_snps, up_oct_hi_het_snps)
length(up_fixed_hi_het_snps)
# 144

up_octOnly_hi_het_snps <- setdiff(up_oct_hi_het_snps, up_fixed_hi_het_snps)
length(up_octOnly_hi_het_snps)
# 5
up_tet_df[up_octOnly_hi_het_snps,]
# very high heterozygosity in the 4X samples, so these probably aren't 
#   specific to 8X samples

up_oct_22_bysnp <- apply(up_oct_df, 1,
  function(x) sum(x == '2/2', na.rm = T))
up_oct_22_bysnp_1 <- up_oct_22_bysnp * up_oct_norm_snp_vec
up_oct_22_fixed <- which(up_oct_22_bysnp_1 > (ncol(up_oct_df)*0.95))
length(up_oct_22_fixed)
# 0 show up as 2/2 in all samples, which is not a surprise

sim8X_het_bysnp <- lapply(sim8X_vcf_list, function(x)
  apply(x, 1, function(y) sum(y == '3/1' | y == '2/2' | y == '1/3', na.rm = T)))
sim8X_het_bysnp_1 <- list()
for(sug in sim8X_up_groups){
  sim8X_het_bysnp_1[[sug]] <- (sim8X_het_bysnp[[sug]] * 
    sim8X_norm_snp_list[[sug]])
}
sim8X_hi_het_snps <- lapply(sim8X_het_bysnp_1, function(x) which(x > (16*0.95)))
# lapply(sim8X_hi_het_snps, length)
show_missdata_info(sim8X_hi_het_snps, sim8X_norm_snp_list)
# 15%, 10%, and 12% of up_tet_1_8X, up_tet_2_8X, and up_tet_1_up_tet_2_8X 
#    hi_het_snps show hi missing data

sim8X_22_bysnp <- lapply(sim8X_vcf_list, function(x)
  apply(x, 1, function(y) sum(y == '2/2', na.rm = T)))
sim8X_22_bysnp_1 <- list()
for(sug in names(sim8X_22_bysnp)){
  sim8X_22_bysnp_1[[sug]] <- (sim8X_22_bysnp[[sug]] * 
    sim8X_norm_snp_list[[sug]])
}
sim8X_22_fixed <- lapply(sim8X_22_bysnp_1, function(x) which(x > (16*0.95)))
# lapply(sim8X_22_fixed, length)
show_missdata_info(sim8X_22_fixed, sim8X_norm_snp_list)
# 5.5%, 4.7%, and 6.6% of up_tet_1_8X, up_tet_2_8X, and up_tet_1_up_tet_2_8X 
#  fixed 2/2 SNPs show hi missing data

for(s8n in names(sim8X_22_fixed)){
  print(s8n)
  print(length(sim8X_22_fixed[[s8n]])/length(sim8X_hi_het_snps[[s8n]]))
}
# 33%, 17%, and 18% of hi_het_SNPs are 2/2 in up_tet_1_8X, up_tet_2_8X, and
#  up_tet_1_up_tet_2_8X simulated populations

sub_het_bysnp <- lapply(up_sub_vcf_list, function(x)
  apply(x, 1, function(y) sum(y == '3/1' | y == '2/2' | y == '1/3', na.rm = T)))
sub_het_bysnp_1 <- list()
for(subp in names(sub_het_bysnp)){
  sub_het_bysnp_1[[subp]] <- sub_het_bysnp[[subp]] * sub_norm_snp_list[[subp]]
}
sub_hi_het_snps <- list()
for(subp in names(sub_het_bysnp_1)){
  sub_hi_het_snps[[subp]] <- which(sub_het_bysnp_1[[subp]] > 
  (ncol(up_sub_vcf_list[[subp]]) * 0.95))
}
# lapply(sub_hi_het_snps, length)
# 186 and 151
# length(intersect(sub_hi_het_snps[[1]], sub_hi_het_snps[[2]]))
# 116 overlap
show_missdata_info(sub_hi_het_snps, sub_norm_snp_list)
# 10% and 13% show hi missing data

subpop_het_bysnp <- lapply(subpop_vcf_list, function(x)
  apply(x, 1, function(y) sum(y == '3/1' | y == '2/2' | y == '1/3', na.rm = T)))
subpop_het_bysnp_1 <- list()
for(subpop in names(subpop_het_bysnp)){
  subpop_het_bysnp_1[[subpop]] <- subpop_het_bysnp[[subpop]] * 
    subpop_norm_snp_list[[subpop]]
}
subpop_hi_het_snps <- list()
for(subpop in names(subpop_het_bysnp_1)){
  subpop_hi_het_snps[[subpop]] <- which(subpop_het_bysnp_1[[subpop]] >
  (ncol(subpop_vcf_list[[subpop]]) * 0.95))
}
# lapply(subpop_hi_het_snps, length)

# sim8X fixed HOM-ALTs
sim8X_pop_fixed_homA_list
# subpop_subsamples fixed HOM-ALTS
sub_pop_fixed_homA_list
# subpop fixed HOM-ALTS
subpop_pop_homA_list

# sim8X fixed 2/2s
sim8X_22_fixed
# subpop fixed Hets
subpop_hi_het_snps

length(intersect(subpop_hi_het_snps[['up_tet_1']], 
  subpop_hi_het_snps[['up_tet_2']]))
# 101

length(sim8X_22_fixed[['up_tet_1_up_tet_2_8X']])
# 90
length(intersect(up_fixed_hi_het_snps, 
  sim8X_22_fixed[['up_tet_1_up_tet_2_8X']]))
# 56; 34 are not fixed across all uplands

length(intersect(sim8X_22_fixed[['up_tet_1_up_tet_2_8X']], 
  intersect(subpop_hi_het_snps[['up_tet_1']],
  subpop_hi_het_snps[['up_tet_2']])))
# 78
# 12 are not fixed HET in both progenitors

length(intersect(sim8X_22_fixed[['up_tet_1_up_tet_2_8X']], unique(
  c(subpop_hi_het_snps[['up_tet_1']], subpop_hi_het_snps[['up_tet_2']]))))
# 83
# 7 are not fixed HET in either progenitor

length(intersect(sim8X_22_fixed[['up_tet_1_up_tet_2_8X']], unique(
  c(subpop_hi_het_snps[['up_tet_1']], subpop_hi_het_snps[['up_tet_2']],
  up_fixed_hi_het_snps))))
# 85
# 5 are new fixed HET SNPs (well, not really)

subpop_vcf_list[['up_tet_2']][
  setdiff(sim8X_22_fixed[['up_tet_1_up_tet_2_8X']], unique(
  c(subpop_hi_het_snps[['up_tet_1']], subpop_hi_het_snps[['up_tet_2']],
  up_fixed_hi_het_snps)))
, ]
# The remaining 5 are essentially fixed-HET in the two subpops, but below
#  the 95% cutoff

######

length(sim8X_22_fixed[['up_tet_1_low_tx_1_8X']])
# 33

length(intersect(subpop_hi_het_snps[['up_tet_1']],
  subpop_hi_het_snps[['low_tx_1']]))
# 7

subpop_vcf_list[['low_tx_1']][
  setdiff(sim8X_22_fixed[['up_tet_1_low_tx_1_8X']], 
    intersect(subpop_hi_het_snps[['up_tet_1']],
    subpop_hi_het_snps[['low_tx_1']]))
, ]
# 13 are homALT in up_tet_1 and homREF in low_tx_1; 15 fixed-HET in both
#  subpops but don't reach the 95% threshold

####
length(sim8X_22_fixed[['up_tet_1_low_ec_1_8X']])
# 63

length(intersect(subpop_hi_het_snps[['up_tet_1']],
  subpop_hi_het_snps[['low_ec_1']]))
# 59

length(setdiff(sim8X_22_fixed[['up_tet_1_low_ec_1_8X']],
intersect(subpop_hi_het_snps[['up_tet_1']],
  subpop_hi_het_snps[['low_ec_1']])
))
# 49

subpop_vcf_list[['low_ec_1']][
  setdiff(sim8X_22_fixed[['up_tet_1_low_ec_1_8X']],
    intersect(subpop_hi_het_snps[['up_tet_1']],
    subpop_hi_het_snps[['low_ec_1']])
  )
, ]

# 6 are fixed-HETs in the subpops, remaining 43 are "real" 




length(intersect(subpop_hi_het_snps[['up_tet_1']], 
  subpop_hi_het_snps[['up_oct_1']]))
# 106



length(intersect(subpop_hi_het_snps[['up_tet_1']], 
  subpop_hi_het_snps[['up_oct_2']]))
# 63

length(intersect(subpop_hi_het_snps[['up_tet_2']],
  subpop_hi_het_snps[['up_oct_1']]))
# 64

length(intersect(subpop_hi_het_snps[['up_tet_2']],
  subpop_hi_het_snps[['up_oct_2']]))
# 58

length(intersect(subpop_hi_het_snps[['up_oct_1']],
  subpop_hi_het_snps[['up_oct_2']]))
# 73

####
# KEEP GOING FROM HERE
####

length(intersect(sim8X_22_fixed[[1]], sub_pop_fixed_homA_list[[1]]))
# 0
length(intersect(sim8X_22_fixed[[1]], sub_hi_het_snps[[1]]))
# all of the fixed-het 8X loci  are fixed-het in the 4X population

# up_tet_1 and low_gc_2
length(intersect(sim8X_22_fixed[['up_tet_1_low_gc_2_8X']], 
  subpop_pop_homA_list[['up_tet_1']]))
# 3
length(intersect(sim8X_22_fixed[['up_tet_1_low_gc_2_8X']],
  subpop_pop_homA_list[['low_gc_2']]))
# 0

subpop_vcf_list[['up_tet_1']][sim8X_22_fixed[['up_tet_1_low_gc_2_8X']],]
subpop_vcf_list[['low_gc_2']][sim8X_22_fixed[['up_tet_1_low_gc_2_8X']],]

sim8X_vcf_list[['up_tet_1_low_gc_2_8X']][subpop_pop_homA_list[['up_tet_1']][1:20],]


# load annotation
exon_info_file <- '/global/homes/g/grabowsp/data/switchgrass/ref_stuff/Pvirgatum_516_v5.1/Pvirgatum_516_v5.1.gene_exons.gff3'

exon_info <- read.table(exon_info_file, stringsAsFactors = F)

exon_chr01k <- exon_info[which(exon_info$V1 == 'Chr01K'),]
cds_chr01k <- exon_chr01k[which(exon_chr01k$V3 == 'CDS'),]

####

up_tet_hi_df <- data.frame(POS = sort(vcf$POS[up_tet_hi_het_snps]),
  stringsAsFactors = F)

up_tet_hi_df$gff_ind <- unlist(sapply(up_tet_hi_df$POS, 
  function(x) intersect(which(cds_chr01k$V4 < x), which(cds_chr01k$V5 > x))[1]))

up_tet_hi_df$gene <- sapply(up_tet_hi_df$gff_ind, function(x) sub('ID=', '',
  paste(unlist(strsplit(cds_chr01k$V9[x], split = '.', fixed = T))[c(1,2)],
  collapse = '.')))

up_oct_hi_df <- data.frame(POS = sort(vcf$POS[up_oct_hi_het_snps]),
  stringsAsFactors = F)

up_oct_hi_df$gff_ind <- unlist(sapply(up_oct_hi_df$POS, 
  function(x) intersect(which(cds_chr01k$V4 < x), which(cds_chr01k$V5 > x))[1]))

up_oct_hi_df$gene <- sapply(up_oct_hi_df$gff_ind, function(x) sub('ID=', '',
  paste(unlist(strsplit(cds_chr01k$V9[x], split = '.', fixed = T))[c(1,2)],
  collapse = '.')))

up_tet_hi_genes <- unique(up_tet_hi_df$gene)
up_oct_hi_genes <- unique(up_oct_hi_df$gene)

length(up_tet_hi_genes)
#[1] 60

length(up_oct_hi_genes)
#[1] 94

length(intersect(up_tet_hi_genes, up_oct_hi_genes))
# [1] 58

length(setdiff(up_oct_hi_genes, up_tet_hi_genes))
# [1] 36

up_tet_overlap_inds <- which(up_tet_hi_df$gene %in% up_oct_hi_df$gene)

up_tet_unique_df <- up_tet_hi_df[-up_tet_overlap_inds, ]
# 2 SNPs, 2 genes

up_oct_overlap_inds <- which(up_oct_hi_df$gene %in% up_tet_hi_df$gene)

up_oct_hi_unique_df <- up_oct_hi_df[-up_oct_overlap_inds, ]
# 106 SNPs, 36 genes; 17 singletons (genes with just 1 SNP)

estimated_up_oct_fixHET <- (p_mw4_homA_fixed * (1-p_mw4_homA_fixed) * 
  nrow(up_oct_df))*2
# 264 estimated fixed HET SNPs for allo-autopolyploid 8X if use 
#   prob of HOM-ALT from Midwest-4X
# see 15.5X (17) fewer singleton 8X Fixed-HET SNPs





