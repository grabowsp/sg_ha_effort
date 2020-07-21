gen_struc_genos <- function(genotypes, lib_name, ploidy,
  geno_type = 'all_tet'){
  if(geno_type == 'all_dip'){
    tmp_struc <- matrix(nrow = 2, ncol = length(genotypes), NA)
  } else if(geno_type == 'pseudohap'){
    tmp_struc <- matrix(nrow = 1, ncol = length(genotypes), NA)
  } else {
    tmp_struc <- matrix(nrow = 4, ncol = length(genotypes), NA)
  }
  tmp_na_inds <- which(is.na(genotypes))
  tmp_struc[, tmp_na_inds] <- -9
  #
  hom_ref <- which(genotypes == 0)
  A1_inds <- which(genotypes == 1)
  A2_inds <- which(genotypes == 2)
  A3_inds <- which(genotypes == 3)
  A4_inds <- which(genotypes == 4)
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
    if(geno_type == '')
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
  struc_df <- data.frame(lib = lib_name, tmp_struc, stringsAsFactors = F)
  return(struc_df)
}






