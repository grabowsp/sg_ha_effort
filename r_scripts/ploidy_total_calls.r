# Make ploidy calls using both nQuire and MNP data

### LOAD PACKAGES  ###


### LOAD DATA ###
ploidy_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v1.0.txt'
ploidy_res <- read.table(ploidy_res_file, header = T, stringsAsFactors = F, 
  sep = '\t')

### SET OUTPUT ###

ploidy_res_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v2.0.txt'


### SET VARIABLES ###


###########

great_4X <- intersect(
  intersect(which(ploidy_res$cds_mnp_ploidy == '4X'), 
    which(ploidy_res$genic_mnp_ploidy == '4X')),
  intersect(which(ploidy_res$nquire_c20_ploidy == '4X'), 
    which(ploidy_res$nquire_c20cluster_ploidy == '4X'))
  )

great_8X <- intersect(
  intersect(which(ploidy_res$cds_mnp_ploidy == '8X'), 
    which(ploidy_res$genic_mnp_ploidy == '8X')),
  intersect(which(ploidy_res$nquire_c20_ploidy == '8X'), 
    which(ploidy_res$nquire_c20cluster_ploidy == '8X'))
  )

### CONTINUE FROM HERE###
## assign ploidy_confidence:
## base = 1, good = 2, strong = 3, great = 4

# strong = c20cluster and cds_mnp or genic_mnp (and other NA)
# good = c20cluster and both cds_mnp and genic_mnp NA
# base = c20cluster (or outliers?)

## outliers:
# c20cluster above/below cds_mnp 100
# c20cluster with wrong genic_mnp genotype
# cds_mnp and genic_mnp don't match
 
strong_4X <- intersect(which(ploidy_res$nquire_c20cluster_ploidy == '4X'), 
  union(
    intersect(which(ploidy_res$cds_mnp_ploidy == '4X'), 
      which(ploidy_res$genic_mnp_ploidy == '?X')), 
    intersect(which(ploidy_res$cds_mnp_ploidy == '?X'), 
      which(ploidy_res$genic_mnp_ploidy == '4X'))
  )
)

strong_8X <- intersect(which(ploidy_res$nquire_c20cluster_ploidy == '8X'), 
  union(
    intersect(which(ploidy_res$cds_mnp_ploidy == '8X'),
      which(ploidy_res$genic_mnp_ploidy == '?X')),
    intersect(which(ploidy_res$cds_mnp_ploidy == '?X'),
      which(ploidy_res$genic_mnp_ploidy == '8X'))
  )
)

good_4X <- intersect(which(ploidy_res$nquire_c20cluster_ploidy == '4X'), 
  intersect(which(ploidy_res$cds_mnp_ploidy == '?X'), 
    which(ploidy_res$genic_mnp_ploidy == '?X')
  )
)

good_8X <- intersect(which(ploidy_res$nquire_c20cluster_ploidy == '8X'), 
  intersect(which(ploidy_res$cds_mnp_ploidy == '?X'),
    which(ploidy_res$genic_mnp_ploidy == '?X')
  )
)

base_4X <- setdiff(which(ploidy_res$nquire_c20cluster_ploidy == '4X'),
  c(great_4X, great_8X, strong_4X, strong_8X, good_4X, good_8X)
)

base_8X <- setdiff(which(ploidy_res$nquire_c20cluster_ploidy == '8X'),
  c(great_4X, great_8X, strong_4X, strong_8X, good_4X, good_8X)
)

###
overall_4X <- c(great_4X, strong_4X, good_4X, base_4X)
overall_8X <- c(great_8X, strong_8X, good_8X, base_8X)

great_inds <- c(great_4X, great_8X)
strong_inds <- c(strong_4X, strong_8X)
good_inds <- c(good_4X, good_8X)
base_inds <- c(base_4X, base_8X)

ploidy_res$total_ploidy <- '?X'
ploidy_res$total_ploidy[overall_4X] <- '4X'
ploidy_res$total_ploidy[overall_8X] <- '8X'

ploidy_res$tot_ploid_confidence <- 0
ploidy_res$tot_ploid_confidence[great_inds] <- 4
ploidy_res$tot_ploid_confidence[strong_inds] <- 3
ploidy_res$tot_ploid_confidence[good_inds] <- 2
ploidy_res$tot_ploid_confidence[base_inds] <- 1

write.table(ploidy_res, file = ploidy_res_out_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

##############
# lists of outliers
hex_inds <- which(ploidy_res$nquire_c20_ploidy == '6X')
# 5; all with 1 confidence
## IELY - "NOTE_PLOIDY" in metadata says "Investigate"; From Texas;
### no subpop because Sujan call is 8X; Old ploidy call is 6X; 
### MNP is 4X; nquire is 6X or 8X; 6X is lowest portion for c20 through c50
## IJBQ - Eastcoast; from North Carolina; 
### 4X in all other calls except the basic nQuire result;
### The other sample from this population is strongly 4X;
### my guess is this is a true 4X...; fairly low coverage
## INFT - from Florida; Gulf Coast Admixed
### MNP is 4X; nQuire cluster is 8X; nQuire is 6X for c20-c40
## PYSU - is Kanlow; MNP STRONGLY 4X; nquire cluster is 4X;
### at c20, 4X and 6X are essentially the same; at c30-c50, is 4X
### by guess is a true 4X
## XXAW - Texas_Admixed; F1 of AP13_vs_VS16; pretty strongly 4X in MNP;
### 4X in nQuire cluster; nQuire 6X at c20 but 4X at higher coverage
### is probably 4X BUT not important for natural pop analysis

c20_8X_wrongMNP_inds <- intersect(
  which(ploidy_res$nquire_c20cluster_ploidy == '8X'),
  union(which(ploidy_res$cds_mnp_ploidy == '4X'), 
    which(ploidy_res$genic_mnp_ploidy == '4X')
  )
)
# 10; all with 1 confidence
## ABHM - Texas; West Bee Caves; WBC3a; Sujan called 4X; Strongly 4X with MNP;
### nQuire 8X for c20-c40, 4X if used c50 results but very few c50 SNPs
## ASHOA - Texas; somewhat low coverage, but not terrible;
### Called as 8X by Sujan so no Subpop; 8X using CDS MNPs;
### 4X using genic MNPs; 8X using nQuire;
### Probably is 8X
## IELY - see 6X
## IEYM - Eastcoast; North Carolina; MNP indicated 4X (weakly); nquire says 8X
## IICH - From Pennsylvania; 4X via MNP; 8X via nQuire;
### there are two other samples from sampe pop: one can't be called by
### MNP and is 8X from nQuire, the other is 8X with all metrics.
### Is interesting that 2 samples from same pop have nQuire and MNP results
### that don't quite match
## INFT - See 6X
## ITAU - Gulfcoast Admixed; Texas; "Old" ploidy call is 8X; MNP is 4X;
### nQuire is 8X; other sample from pop is strongly 4X, but upland Chloroplast
## PZXO - North Carolina; strongly 4X from MNP; 8X with nQuire;
### Called 8X by sujan; other sample from pop is strongly 4X
### Interestingly, is nQuire 6X at c30-c50
## WZGG - Kanlow; strongly 4X with MNP, 8X with nQuire, called 8X by Sujan;
### Perhaps is a spontaneous 8X? There is plenty of sequencing coverage
## XXAY - "Kanlow Derived" in "NOTE_LATLONG"

c20_4X_wrongMNP_inds <- intersect(
  which(ploidy_res$nquire_c20cluster_ploidy == '4X'),
  union(which(ploidy_res$cds_mnp_ploidy == '8X'),
    which(ploidy_res$genic_mnp_ploidy == '8X')
  )
)
# 2; all with 1 confidence
## IIMA - Eastcoast admixed; New York; 4X with nQuire,
### Weakly 8X with MNP; metadata says this sample has high heterozygosity,
### so maybe that's why it has a lot of MNPs
## IRUY - Midwest; Illinois; CDS MNP = 4X; genic MNP = 8X; 
### nquire = 4X at all coverage levels 
### I'm guessing this is 4X and there's something weird with MNP results

MNP_mismatch <- union(
  intersect(which(ploidy_res$cds_mnp_ploidy == '8X'),
    which(ploidy_res$genic_mnp_ploidy == '4X')
  ),
  intersect(which(ploidy_res$cds_mnp_ploidy == '4X'),
    which(ploidy_res$genic_mnp_ploidy == '8X')
  )
)
# 2; all with 1 confidence
## ASHOA - see above
## IRUY - see above

tot_outliers <- unique(c(hex_inds, c20_8X_wrongMNP_inds, c20_4X_wrongMNP_inds, 
  MNP_mismatch))
# These include all the samples with confidence of 1

c20_4X_over100 <- intersect(
  which(ploidy_res$nquire_c20cluster_ploidy == '4X'),
  which(ploidy_res$cds_mnp_stand > 100)
)
# 1; with 1 confidence
## IIMA - see above

c20_8X_under100 <- intersect(
  which(ploidy_res$nquire_c20cluster_ploidy == '8X'),
  intersect(which(ploidy_res$cds_mnp_stand < 100),
    which(ploidy_res$seq_cov > 1e10)
  )
)
# 12; 1 with 3 confidence, 2 with 2 confidence, 9 with 1 confidence
# 3 Confidence
## INGC - Massachusetts; 8X call by Sujan so no subpop; 
### 8X with genic MNP and nQuire, no call with CDS MNP (91)
# 2 Confidence (399,403)
## IICA - Pennsylvania; is same pop as outlier IICH from above; 8X from nQuire;
### Weakly 4X using MNP (no call because within the no-call interval)
## IICG - Pennsylvania; 8X call by Sujan so no Subpop; 8X with nQuire;
### boarderline 4X with CDS MNP (97); boarderline 8X with geneic MNP (409);
### but is within the no-call intervals for MNP; other sample from pop is
### pretty strongly 8X 

quit(save = 'no')
