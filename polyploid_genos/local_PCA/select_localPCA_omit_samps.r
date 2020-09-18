# Script for selecting samples to omit based on excess missing data

### LOAD DATA #####
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/'

#vcf_in_short <- 'Chr01K.polyploid.CDS.expandgeosamps.vcf_00'
vcf_in_short <- 'Chr01K.polyploid.CDS.expandgeosamps.vcf_00'

vcf_in <- paste(data_dir, vcf_in_short, sep = '')

vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

#head_in_short <- 'CDS.expandgeosamps.vcf.header.txt'
head_in_short <- 'CDS.expandgeosamps.vcf.header.txt'

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

### SET OUTPUT ###
bad_samp_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/local_pca_remove_samps_09152020.txt'

###############

oct_df <- vcf[, oct_libs]
tet_df <- vcf[, tet_libs]

oct_df[oct_df == './.'] <- NA
tet_df[tet_df == './.'] <- NA

geno_df <- cbind(oct_df, tet_df)

geno_mat <- matrix(unlist(geno_df), ncol = ncol(geno_df),
  byrow = F)

rownames(geno_mat) <- seq(nrow(geno_mat))
colnames(geno_mat) <- colnames(geno_df)

per_miss_samp <- apply(geno_mat, 2, function(x) sum(is.na(x)) / nrow(geno_mat))

bad_samps <- names(which(per_miss_samp > 0.05))

write.table(bad_samps, bad_samp_out_file, quote = F, sep = '\t', row.names = F,
  col.names = F)


