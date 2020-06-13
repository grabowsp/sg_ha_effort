# Functions used for analysis of polyploid genomes

generate_dosage_df <- dataframe(vcf, oct_libs = c(), tet_libs = c(), R1 = T){
  # Function for generating vcf of allele dosages ranging from 0 to 2
  #  based on 0-4 allele genotypes
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
  # R1 = if T, REF alleles = 1; if F, ALT alleles = 1
  # OUTPUT
  # data.frame with dosages ranging from 0 to 2 based on ploidy
  #########
  geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
  if(R1){
    oct_dosage_vec <- c('2', '1.5', '1', '0.5', '0')
    tet_dosage_vec <- c('2', '1', '1', '1', '0')
  } else {
    oct_dosage_vec <- c('0', '0.5', '1', '1.5', '2')
    tet_dosage_vec <- c('0', '1', '1', '1', '2')
  }
  vcf[vcf == './.'] <- NA
  if(length(oct_libs) > 0){
    oct_df <- vcf[, oct_libs]
    for(i in seq(length(geno_vec))){
      oct_df[oct_df == geno_vec[i]] <- oct_dosage_vec[i]
    }
    for(j8 in seq(ncol(oct_df))){
      oct_df[, j8] <- as.numeric(oct_df[, j8])
    }
  }
  if(length(tet_libs) > 0){
    tet_df <- vcf[, tet_libs]
    for(i in seq(length(geno_vec))){
      tet_df[tet_df == geno_vec[i]] <- tet_dosage_vec[i]
    }
    for(j4 in seq(ncol(tet_df))){
      tet_df[, j4] <- as.numeric(tet_df[, j4])
    }
  }
  if(length(oct_libs) == 0){
    tot_df <- tet_df
  } else if(length(tet_libs) == 0){
    tot_df <- oct_df
  } else{
    tot_df <- cbind(oct_df, tet_df)
  }
  return(tot_df)
}


