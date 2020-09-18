# Goals
# Calculate homozygosity - % homozygous genotypes
# Calculate het-dosage - sum het genotypes, with RA = RRAA = 0.5;
#   RAAA = RRRA = 0.25

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/R_analysis


#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/'
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/'

#vcf_in_short <- 'Chr01K.polyploid.CDS.expandgeosamps.vcf_00'
vcf_in_short <- 'Chr01K.polyploid.CDS.upexpand.vcf_00'

vcf_in <- paste(data_dir, vcf_in_short, sep = '')

vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

#head_in_short <- 'CDS.expandgeosamps.vcf.header.txt'
head_in_short <- 'CDS.upexpand.vcf.header.txt'

vcf_header_file <- paste(data_dir, head_in_short, sep = '')

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf) <- vcf_header

tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs <- intersect(tet_libs_0, colnames(vcf))

oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs <- intersect(oct_libs_0, colnames(vcf))



oct_df <- vcf[, oct_libs]
tet_df <- vcf[, tet_libs]

oct_df[oct_df == './.'] <- NA
tet_df[tet_df == './.'] <- NA

geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
# note: the dosage vectors are wonky and strictly for this script
oct_dosage_vec <- c('0', '0.5', '1', '0.5', '0')
tet_dosage_vec <- c('0', '1', '1', '1', '0')

for(i in seq(length(geno_vec))){
  oct_df[oct_df == geno_vec[i]] <- oct_dosage_vec[i]
}

for(i in seq(length(geno_vec))){
  tet_df[tet_df == geno_vec[i]] <- tet_dosage_vec[i]
}

for(i in seq(ncol(oct_df))){
  oct_df[, i] <- as.numeric(oct_df[, i])
}

for(i in seq(ncol(tet_df))){
  tet_df[, i] <- as.numeric(tet_df[, i])
}

geno_df <- cbind(oct_df, tet_df)

geno_mat <- matrix(unlist(geno_df), ncol = nrow(geno_df),
  byrow = T)

rownames(geno_mat) <- colnames(geno_df)

na_vec <- apply(geno_mat, 1, function(x) sum(is.na(x)))

norm_vec <- 1/(1-(na_vec/ncol(geno_mat)))

hom_vec <- apply(geno_mat, 1, function(x) sum(x == 0, na.rm = T))
hom_vec_1 <- hom_vec * norm_vec
hom_vec_2 <- hom_vec_1/ncol(geno_mat)

het_vec <- apply(geno_mat, 1, function(x) sum(x, na.rm = T))
het_vec_1 <- het_vec * norm_vec
het_vec_2 <- het_vec_1/ncol(geno_mat)

#####################

# allo vs analysis

oct_na_bysnp <- apply(oct_df, 1, function(x) sum(is.na(x)))
oct_hom_bysnp <- apply(oct_df, 1, function(x) sum(x == 0, na.rm = T))

oct_norm_snp_vec <- 1/(1-(oct_na_bysnp/ncol(oct_df)))

oct_hom_bysnp_1 <- oct_hom_bysnp * oct_norm_snp_vec

sum(oct_hom_bysnp_1 < (ncol(oct_df)*.5))/nrow(oct_df)
#[1] 0.02750138

tet_na_bysnp <- apply(tet_df, 1, function(x) sum(is.na(x)))
tet_hom_bysnp <- apply(tet_df, 1, function(x) sum(x == 0, na.rm = T))

tet_norm_snp_vec <- 1/(1-(tet_na_bysnp/ncol(tet_df)))

tet_hom_bysnp_1 <- tet_hom_bysnp * tet_norm_snp_vec

sum(tet_hom_bysnp_1 < (ncol(tet_df)*.5))/nrow(tet_df)
# 0.01112056

tet_hi_het <- which(tet_hom_bysnp_1 < (ncol(tet_df)*.2))
oct_hi_het <- which(oct_hom_bysnp_1 < (ncol(oct_df)*.2))

exon_info_file <- '/global/homes/g/grabowsp/data/switchgrass/ref_stuff/Pvirgatum_516_v5.1/Pvirgatum_516_v5.1.gene_exons.gff3'

exon_info <- read.table(exon_info_file, stringsAsFactors = F)

exon_chr01k <- exon_info[which(exon_info$V1 == 'Chr01K'),]
cds_chr01k <- exon_chr01k[which(exon_chr01k$V3 == 'CDS'),]

tet_hi_df <- data.frame(POS = sort(vcf$POS[tet_hi_het]), 
  stringsAsFactors = F)

tet_hi_df$gff_ind <- unlist(sapply(tet_hi_df$POS, function(x) intersect(which(cds_chr01k$V4 < x), which(cds_chr01k$V5 > x))[1]))

tet_hi_df$gene <- sapply(tet_hi_df$gff_ind, function(x) sub('ID=', '', 
  paste(unlist(strsplit(cds_chr01k$V9[x], split = '.', fixed = T))[c(1,2)], 
  collapse = '.')))

oct_hi_df <- data.frame(POS = sort(vcf$POS[oct_hi_het]),
  stringsAsFactors = F)

oct_hi_df$gff_ind <- unlist(sapply(oct_hi_df$POS, function(x) intersect(which(cds_chr01k$V4 < x), which(cds_chr01k$V5 > x))[1]))

oct_hi_df$gene <- sapply(oct_hi_df$gff_ind, function(x) sub('ID=', '',
  paste(unlist(strsplit(cds_chr01k$V9[x], split = '.', fixed = T))[c(1,2)],
  collapse = '.')))

tet_hi_genes <- unique(tet_hi_df$gene)
oct_hi_genes <- unique(oct_hi_df$gene)

length(tet_hi_genes)
#[1] 65
length(oct_hi_genes)
#[1] 141
length(intersect(tet_hi_genes, oct_hi_genes))
# [1] 63
length(setdiff(oct_hi_genes, tet_hi_genes))
# [1] 78
oct_overlap_inds <- which(oct_hi_df$gene %in% tet_hi_df$gene)

oct_hi_unique_df <- oct_hi_df[-oct_overlap_inds, ]

table(table(oct_hi_unique_df$gene))
#  1  2  3  4  6  7  8 12 37 
# 40 20  6  3  4  2  1  1  1

table(table(tet_hi_df$gene))
#  1  2  3  4  5  6  7 10 11 13 16 17 22 
# 25 10  6  7  2  5  2  1  3  1  1  1  1

oct_df_2 <- vcf[, oct_libs]
tet_df_2 <- vcf[, tet_libs]

oct_df_2[oct_df_2 == './.'] <- NA
tet_df_2[tet_df_2 == './.'] <- NA

oct_4A_counts <- apply(oct_df_2, 1, function(x) sum(x == '0/4', na.rm = T))

oct_4A_counts_2 <- apply(oct_df_2, 2, function(x) sum(x == '0/4', na.rm = T))

oct_na_bysamp <- apply(oct_df_2, 2, function(x) sum(is.na(x)))
oct_norm_samp_vec <- 1 / (1-(oct_na_bysamp/nrow(oct_df_2)))

oct_4A_counts_norm <- oct_4A_counts_2*oct_norm_samp_vec

tet_4A_counts_2 <- apply(tet_df_2, 2, function(x) sum(x == '0/4', na.rm = T))
tet_na_bysamp <- apply(tet_df_2, 2, function(x) sum(is.na(x)))
tet_norm_samp_vec <- 1 / (1-(tet_na_bysamp/nrow(tet_df_2)))

tet_4A_counts_norm <- tet_4A_counts_2*tet_norm_samp_vec

p_tet_2A <- median(tet_4A_counts_norm/nrow(tet_df_2))


min_4A_count <- min(oct_4A_counts_norm)
max_4A_count <- max(oct_4A_counts_norm)

per_4A <- oct_4A_counts_norm / nrow(oct_df_2)

p_2A <- sqrt(per_4A)

expected_2A_0B <- median(p_2A * (1-p_2A)*nrow(oct_df_2))

summary(p_2A)

oct_norm_snp_vec <- 1/(1-(oct_na_bysnp/nrow(oct_df)))

oct_hom_bysnp_1 <- oct_hom_bysnp * oct_norm_snp_vec



geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
# note: the dosage vectors are wonky and strictly for this script
oct_dosage_vec <- c('0', '0.5', '1', '0.5', '0')
tet_dosage_vec <- c('0', '1', '1', '1', '0')

for(i in seq(length(geno_vec))){
  oct_df[oct_df == geno_vec[i]] <- oct_dosage_vec[i]
}

for(i in seq(length(geno_vec))){
  tet_df[tet_df == geno_vec[i]] <- tet_dosage_vec[i]
}

for(i in seq(ncol(oct_df))){
  oct_df[, i] <- as.numeric(oct_df[, i])
}


