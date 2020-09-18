# Goals
# Evaluate if the number of fixed HOM-ALT and HET SNPs is consistent
#   with 8X switchgrass being auto-allopolyploid or allo-allopolyploid

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


#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/'
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/'

#vcf_in_short <- 'Chr01K.polyploid.CDS.expandgeosamps.vcf_00'
vcf_in_short <- 'Chr01K.polyploid.CDS.upexpand.vcf_00'

vcf_in <- paste(data_dir, vcf_in_short, sep = '')

vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

### Try using more SNPs
subfile_vec <- c('01', '02', '03', '04')
for(sfv in subfile_vec){
  tmp_short <- paste('Chr01K.polyploid.CDS.upexpand.vcf_', sfv, sep = '')
  tmp_file_in <- paste(data_dir, tmp_short, sep = '')
  tmp_vcf <- read.table(tmp_file_in, header = F, stringsAsFactors = F, 
    sep = '\t')
  vcf <- rbind(vcf, tmp_vcf)
}

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

# Calculate n_NA's and normalization vectors

tet_na_bysnp <- apply(tet_df, 1, function(x) sum(is.na(x)))
tet_norm_snp_vec <- 1/(1-(tet_na_bysnp/ncol(tet_df)))

oct_na_bysnp <- apply(oct_df, 1, function(x) sum(is.na(x)))
oct_norm_snp_vec <- 1/(1-(oct_na_bysnp/ncol(oct_df)))

# Find HOM-ALT SNPs
tet_homA_bysnp <- apply(tet_df, 1, function(x) sum(x == '0/4', na.rm = T))
tet_homA_bysnp_1 <- tet_homA_bysnp * tet_norm_snp_vec

mw4_homA_fixed_snps <- which(tet_homA_bysnp_1 > (ncol(tet_df)*0.8))
# 1077 SNPs
# 5084 with 5 subfiles

oct_homA_bysnp <- apply(oct_df, 1, function(x) sum(x == '0/4', na.rm = T))
oct_homA_bysnp_1 <- oct_homA_bysnp * oct_norm_snp_vec

oct_homA_fixed_snps <- which(oct_homA_bysnp_1 > (ncol(oct_df)*0.8))
# 859 SNPs
# 4045 with 5

up_homA_fixed_snps <- intersect(mw4_homA_fixed_snps, oct_homA_fixed_snps)
# 855 SNPs
# 4006

mw4_only_homA_fixed_snps <- setdiff(mw4_homA_fixed_snps, oct_homA_fixed_snps)
# 222 SNPs
# 1078

oct_only_homA_fixed_snps <- setdiff(oct_homA_fixed_snps, mw4_homA_fixed_snps)
# 4 SNPs
# 39

p_mw4_homA_fixed <- length(mw4_only_homA_fixed_snps)/nrow(tet_df)
# 0.002220111
# 0.002156022

expected_oct_homA_fixed <- (p_mw4_homA_fixed^2) * nrow(oct_df)
# 0.4928646 : the expected number of 8X homA-fixed SNPs if allo-allopolyploid
#   based on the Midwest-4X number of homA-fixed SNPs
# 2.324191

### Calculate fixed HET sites

tet_het_bysnp <- apply(tet_df, 1, function(x) sum(x == '3/1' | x == '2/2' | x == '1/3', na.rm = T))
tet_het_bysnp_1 <- tet_het_bysnp * tet_norm_snp_vec

oct_het_bysnp <- apply(oct_df, 1, function(x) sum(x == '3/1' | x == '2/2' | x == '1/3', na.rm = T))
oct_het_bysnp_1 <- oct_het_bysnp * oct_norm_snp_vec

tet_hi_het_snps <- which(tet_het_bysnp_1 > (ncol(tet_df)*0.8))
# 256
# should be 0 unless there is balancing selection, so these are probably errors
# 1801

oct_hi_het_snps <- which(oct_het_bysnp_1 > (ncol(oct_df)*0.8))
# 603
# 3073

length(intersect(tet_hi_het_snps, oct_hi_het_snps))
# 236
# 1647

# It is possible that cetain genes have assembly issues (or CNVs), so let's
#  see what's shared and unique to each clade

# load annotation
exon_info_file <- '/global/homes/g/grabowsp/data/switchgrass/ref_stuff/Pvirgatum_516_v5.1/Pvirgatum_516_v5.1.gene_exons.gff3'

exon_info <- read.table(exon_info_file, stringsAsFactors = F)

exon_chr01k <- exon_info[which(exon_info$V1 == 'Chr01K'),]
cds_chr01k <- exon_chr01k[which(exon_chr01k$V3 == 'CDS'),]

####

tet_hi_df <- data.frame(POS = sort(vcf$POS[tet_hi_het_snps]),
  stringsAsFactors = F)

tet_hi_df$gff_ind <- unlist(sapply(tet_hi_df$POS, 
  function(x) intersect(which(cds_chr01k$V4 < x), which(cds_chr01k$V5 > x))[1]))

tet_hi_df$gene <- sapply(tet_hi_df$gff_ind, function(x) sub('ID=', '',
  paste(unlist(strsplit(cds_chr01k$V9[x], split = '.', fixed = T))[c(1,2)],
  collapse = '.')))

oct_hi_df <- data.frame(POS = sort(vcf$POS[oct_hi_het_snps]),
  stringsAsFactors = F)

oct_hi_df$gff_ind <- unlist(sapply(oct_hi_df$POS, 
  function(x) intersect(which(cds_chr01k$V4 < x), which(cds_chr01k$V5 > x))[1]))

oct_hi_df$gene <- sapply(oct_hi_df$gff_ind, function(x) sub('ID=', '',
  paste(unlist(strsplit(cds_chr01k$V9[x], split = '.', fixed = T))[c(1,2)],
  collapse = '.')))

tet_hi_genes <- unique(tet_hi_df$gene)
oct_hi_genes <- unique(oct_hi_df$gene)

length(tet_hi_genes)
#[1] 65
# 295
length(oct_hi_genes)
#[1] 139
# 591
length(intersect(tet_hi_genes, oct_hi_genes))
# [1] 63
# 286
length(setdiff(oct_hi_genes, tet_hi_genes))
# [1] 76
# 305
tet_overlap_inds <- which(tet_hi_df$gene %in% oct_hi_df$gene)

tet_unique_df <- tet_hi_df[-tet_overlap_inds, ]
# 2 SNPs, 2 genes
# 12 SNPs, 9 genes, 7 singletons (genes with just 1 SNP)

oct_overlap_inds <- which(oct_hi_df$gene %in% tet_hi_df$gene)

oct_hi_unique_df <- oct_hi_df[-oct_overlap_inds, ]
# 194 SNPs, 76 genes; 43 singletons (genes with just 1 SNP)
# 803 SNPs, 305 genes; 170 singletons

estimated_oct_fixHET <- (p_mw4_homA_fixed * (1-p_mw4_homA_fixed) * 
  nrow(oct_df))*2
# 443 estimated fixed HET SNPs for allo-autopolyploid 8X if use 
#   prob of HOM-ALT from Midwest-4X
# see 10X (43) fewer singleton 8X Fixed-HET SNPs

# Estimated = 2151.352; see 170 (using 5 subfiles)



