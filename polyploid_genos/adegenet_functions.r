# Custom functions to use during analysis with adegenet library

gen_gl_preobj <- function(vcf, oct_libs = c(), tet_libs = c(),
  maf_cut = 0.001){
  # Function for filtering vcf by MAF freq and generating necessary
  #  objects to be used for making a "genlight" object
  # INPUTS
  # vcf = inputted vcf with header; needs following columns:
  #        ID with snp names
  #        CHROM with chromosome name
  #        POS with SNP position
  #        REF with reference allele
  #        ALT with alternate allele
  #        Sample library names as column names
  # oct_libs = names of 8X libraries
  # tet_libs = names of 4X libraries
  # maf_cut = minor allele frequency cutoff
  # OUTPUT
  # list with following elements, elements of list are based on vcf with 
  #   SNPs with MAF below the cutoff removed
  # 'geno_mat' = matrix of genotypes showing the number of ALT alleles;
  #                maximum number is 2 for 4X and 4 for 8X libraries
  # 'loc.names' =  vector of the names of the retained SNPs
  # 'chromosome' = vector of the chromosome name(s)
  # 'position' = vector of the position of the SNPs
  # 'ploidy' = vector of the ploidy for each sample
  #              2 = 4X switchgrass (because of disomic inheritance)
  #              4 = 8X switchgrass (because of tetrasomic inheritance)
  # 'loc.all' = vector of the REF and ALT alleles, ex: 'g/a'
  #############
  geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
  oct_alt_dose_vec <- c('0', '1', '2', '3', '4')
  tet_alt_dose_vec <- c('0', '1', '1', '1', '2')
  #
  if(length(oct_libs) > 0){
    oct_df <- vcf[, oct_libs]
    oct_df[oct_df == './.'] <- NA
    for(i in seq(length(geno_vec))){
      oct_df[oct_df == geno_vec[i]] <- oct_alt_dose_vec[i]
    }
    for(j8 in seq(ncol(oct_df))){
      oct_df[, j8] <- as.numeric(oct_df[, j8])
    }
    sum_alt_8 <- apply(oct_df, 1, function(x) sum(unlist(x), na.rm = T))
  } else{sum_alt_8 <- c()}
  #
  if(length(tet_libs) > 0){
    tet_df <- vcf[, tet_libs]
    tet_df[tet_df == './.'] <- NA
    for(i in seq(length(geno_vec))){
      tet_df[tet_df == geno_vec[i]] <- tet_alt_dose_vec[i]
    }
    for(j4 in seq(ncol(tet_df))){
      tet_df[, j4] <- as.numeric(tet_df[, j4])
    }
    sum_alt_4 <- apply(tet_df, 1, function(x) sum(unlist(x), na.rm = T))
  } else{sum_alt_4 <- c()}
  #
  if(length(oct_libs) == 0){
    tot_df <- tet_df
    sum_alt_total <- sum_alt_4
  } else if(length(tet_libs) == 0){
    tot_df <- oct_df
    sum_alt_total <- sum_alt_8/2
  } else{
    tot_df <- cbind(oct_df, tet_df)
    sum_alt_total <- (sum_alt_8/2) + sum_alt_4
  }
  #
  n_nas <- apply(tot_df, 1, function(x) sum(is.na(unlist(x))))
  n_gsamps <- ncol(tot_df) - n_nas
  allele_freq <- sum_alt_total / (2*n_gsamps)
  allele_freq_1 <- 1 - allele_freq
  minor_af <- apply(cbind(allele_freq, allele_freq_1), 1, function(x) min(x)[1])
  remove_inds <- which(minor_af < maf_cut)
  #
  keep_inds <- setdiff(seq(nrow(vcf)), remove_inds)
  #
  preobj_ls <- list()
  preobj_ls[['geno_mat']] <- tot_df[keep_inds, ]
  preobj_ls[['loc.names']] <- vcf$ID[keep_inds]
  preobj_ls[['chromosome']] <- vcf$CHROM[keep_inds]
  preobj_ls[['position']] <- vcf$POS[keep_inds]
  preobj_ls[['ploidy']] <- c(rep(4, times = length(oct_libs)),
    rep(2, times = length(tet_libs)))
  preobj_ls[['loc.all']] <- paste(tolower(vcf$REF[keep_inds]),
    tolower(vcf$ALT[keep_inds]), sep = '/')
  #
  return(preobj_ls)
}


