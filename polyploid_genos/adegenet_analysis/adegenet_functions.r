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

gen_gl_object <- function(preobj_list){
  # Function to generate 'genlight' object from the 'preobj_list' generated
  #   using the 'gen_gl_preobj' function
  # INPUTS
  # preobj_list = list generated by 'gen_gl_preobj' that contains elements
  #                for making the 'genlight' object
  # OUTPUT
  # 'genlight' object that includes genotypes, sample ploidy (2 = 4X 
  #    switchgrass, 4 = 8X switchgrass), SNP name, SNP chromosome, 
  #   SNP position, and alleles at each SNP
  ######
  gl <- new('genlight', gen = t(as.matrix(preobj_list[['geno_mat']])),
    ploidy = preobj_list[['ploidy']],
    loc.names = preobj_list[['loc.names']],
    chromosome = preobj_list[['chromosome']],
    position = preobj_list[['position']],
    loc.all = preobj_list[['loc.all']]
  )
  return(gl)
}

gen_struc_genos <- function(genotypes, lib_name, ploidy, 
  geno_type = 'all_tet'){
  # Function to generate genotypes for STRUCTURE in different formats
  #  Can generate:
  #    1) all_tet = tetrasomic genotypes for all samples: 4 lines for all 
  #        samples; 4X (disomic) switchgrass have genotypes replicated to 
  #        fill the 4 lines
  #    2) all_dip = disomic genotypes for all samples: 2 lines for all samples; 
  #        8X (tetrasomic) switchgrass have 2 alleles randomly selected from
  #        their genotypes
  #    3) part_NA = 4 lines for all samples, but 4X (disomic) switchgrass have
  #        NA in lines 3 and 4; 
  #        I suspect this type of genotype will cause an error in STRUCTURE
  #          but I'm not sure yet
  #    4) pseudohap = 1 line for each sample, allele is chosen randomly from
  #        alleles in the genotype
  # INPUTS
  # genotypes = vector of genotypes, typically from genlight object; 
  #   for ploidy = 2: 0 = HOM REF; 2 = HOM ALT
  #   for ploidy = 4: 0 = HOM REF; 4 = HOM ALT
  # lib_name = the name of the library/sample
  # ploidy = "inheritance" ploidy; 4X (disomic) switchgrass have ploidy = 2; 
  #            8x (tetrasomic) switchgrass have ploidy = 4 
  # geno_type = the type of genotypes to generate for STRUCTURE.
  #             Must be 'all_tet', 'all_dip', or 'part_NA'
  # OUTPUT
  # data.frame with first column being the library name and rest of columns 
  #   representing a SNP.
  # If geno_type = 'all_tet': sample has 4 lines, disomic samples have 2
  #   copies of each allele
  # If geno_type = 'all_dip': sample has 2 lines, tetrasomic samples have 2
  #   alleles randomly sampled from the overall genotype
  # If geno_type = 'part_NA': sample has 4 lines, if disomic, lines 3 and 4
  #   are all NA
  ##########
  if(geno_type == 'all_dip'){
    tmp_struc <- matrix(nrow = 2, ncol = length(genotypes), NA)
  } else if(geno_type == 'pseudohap'){
    tmp_struc <- rep(NA, times = length(genotypes))
  } else{
    tmp_struc <- matrix(nrow = 4, ncol = length(genotypes), NA)
  }
  tmp_na_inds <- which(is.na(genotypes))
  if(geno_type == 'pseudohap'){
    tmp_struc[tmp_na_inds] <- -9
  } else{
    tmp_struc[, tmp_na_inds] <- -9
  }
  #
  hom_ref <- which(genotypes == 0)
  A1_inds <- which(genotypes == 1)
  A2_inds <- which(genotypes == 2)
  A3_inds <- which(genotypes == 3)
  A4_inds <- which(genotypes == 4)
  #
  if(geno_type == 'pseudohap'){
    tmp_struc[hom_ref] <- 1
    if(ploidy == 2){
      for(i in A1_inds){
        tmp_struc[i] <- sample(c(1,2), size = 1)
      }
      tmp_struc[A2_inds] <- 2
    }
    if(ploidy == 4){
      for(i in A1_inds){
        tmp_struc[i] <- sample(c(1,1,1,2), size = 1)
      }
      for(j in A2_inds){
        tmp_struc[j] <- sample(c(1,1,2,2), size = 1)
      }
      for(k in A3_inds){
        tmp_struc[k] <- sample(c(1,2,2,2), size = 1)
      }
      tmp_struc[A4_inds] <- 2
    }
  } else {
  if(ploidy == 2){
    if(geno_type == 'all_tet'){
      tmp_struc[ , hom_ref] <- 1
      tmp_struc[c(1:2), A1_inds] <- 1
      tmp_struc[c(3:4), A1_inds] <- 2
      tmp_struc[ , A2_inds] <- 2
    }
    if(geno_type == 'all_dip'){
      tmp_struc[ , hom_ref] <- 1
      tmp_struc[1, A1_inds] <- 1
      tmp_struc[2, A1_inds] <- 2
      tmp_struc[ , A2_inds] <- 2
    }
    if(geno_type == 'part_NA'){
      tmp_struc[c(1:2), hom_ref] <- 1
      tmp_struc[1, A1_inds] <- 1
      tmp_struc[2, A1_inds] <- 2
      tmp_struc[c(1:2), A2_inds] <- 2
      tmp_struc[c(3:4), ] <- -9
    }
  }
  if(ploidy == 4){
    if(geno_type == 'all_dip'){
      tmp_struc[ , hom_ref] <- 1
      tmp_struc[ , A4_inds] <- 2
      for(i in A1_inds){
        tmp_sub <- sample(c(1,1,1,2), size = 2)
        tmp_struc[, i] <- tmp_sub
      }
      for(j in A2_inds){
        tmp_sub <- sample(c(1,1,2,2), size = 2)
        tmp_struc[, j] <- tmp_sub
      }
      for(k in A3_inds){
        tmp_sub <- sample(c(1,2,2,2), size = 2)
        tmp_struc[, k] <- tmp_sub
      }
    } else{
      tmp_struc[ , hom_ref] <- 1
      tmp_struc[c(1:3), A1_inds] <- 1
      tmp_struc[4, A1_inds] <- 2
      tmp_struc[c(1:2), A2_inds] <- 1
      tmp_struc[c(3:4), A2_inds] <- 2
      tmp_struc[1, A3_inds] <- 1
      tmp_struc[c(2:4), A3_inds] <- 2
      tmp_struc[ , A4_inds] <- 2
    }
  }
  }
  if(geno_type == 'pseudohap'){
    struc_df <- c(lib_name, as.character(tmp_struc))
  } else{
    struc_df <- data.frame(lib = lib_name, tmp_struc, stringsAsFactors = F)
  }
  return(struc_df)
}


