# Temporary code while Cori is down

## Identify 8X cultivars to remove

```
meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'

meta <- read.table(meta_file, sep = '\t',
  header = T, stringsAsFactors = F, quote = "", comment.char = '$')

# SOURCE_ID_1
c('Blackwell', 'BoMaster', 'Cave-in-Rock', 'EG1101', 'EG1102', 'Forestburg',
'Miami', 'Stuart', 'Sunburst', 'Trailblazer')
# additions from SOURCE_ID_2
c('Carthage', 'High Tide', 'KY 1625, 'KY1625', 'MIAMI', 'Pathfinder', 
  'Shawnee', 'SHAWNEE', 'Shelter', 'TRAILBLAZER', 'Timber', 'WS4U', 'WS8U')


# also grep "Dacotah", "AP13", "VS16", "NL94", "SL93"

# additional from SOURCE_ID_2

# grep "Flacon", "Summer", "Pathfinder", 'Alamo', 'ALAMO', 'Kanlow', 'KANLOW',

# COLLECTION_TYPE

c('Breeding Selection', 'Cultivar')


cultivar_names <- c('Blackwell', 'BoMaster', 'Cave-in-Rock', 'EG1101', 
  'EG1102', 'Forestburg', 'Miami', 'Stuart', 'Sunburst', 'Trailblazer',
  'Carthage', 'High Tide', 'KY 1625', 'KY1625', 'MIAMI', 'Pathfinder',
  'SHAWNEE', 'Shelter', 'TRAILBLAZER', 'Timber', 'WS4U', 'WS8U',
  'Dacotah', 'AP13', 'VS16', 'NL94', 'SL93', 'Falcon', 'Summer', 'Pathfinder',
  'Alamo', 'ALAMO', 'Kanlow', 'KANLOW')

for(cn in cultivar_names){
  tmp_inds <- union(grep(cn, meta$SOURCE_ID_1), grep(cn, meta$SOURCE_ID_2))
  print(cn)
  print(table(meta$COLLECTION_TYPE[tmp_inds]))
}

```

## Calculate r^2
```
p <- 0.9
q <- 1-p
geno_prob_0 <- c(p^4, 4*(p^3)*q, 6*(p^2)*(q^2), 4*p*(q^3), q^4)

na_prob <- 0.1

geno_prob_1 <- c(geno_prob_0 * (1-na_prob), na_prob)
geno_vec <- c(0,1,2,3,4,NA)

toy <- lapply(1:30, function(x) sample(geno_vec, 10000, prob = geno_prob_1,
  replace = TRUE))

# rows = SNP, samples = column
toy_mat <- matrix(data = unlist(toy), ncol = length(toy), byrow = F)

# coluns = SNP, rows = samples
toy_mat_1 <- matrix(data = unlist(toy), ncol = length(toy[[1]]), byrow = T)

#toy_cor_2 <- cor(toy_mat_1[, 1:1000], use = 'pairwise.complete.obs')

snp_pos <- c(1:10000)
snp_names <- paste('snp_', snp_pos, sep = '')

toy_df <- data.frame(snp_names = snp_names, snp_pos = snp_pos, toy_mat,
  stringsAsFactors = F)

COLNAME_1 <- 2
COLNAME_2 <- 1
FIRST_SAMP <- 3

#snp_dist <- dist(snp_pos, method = 'manhattan', diag = T, upper = T)

calc_snp_r2 <- function(vcf_df, max_dist){
  snp_pos <- vcf_df[, COLNAME_1]
  snp_names <- vcf_df[, COLNAME_2]
  #
  r2_list <- list()
  for(tm in seq(nrow(vcf_df))){
####
#for(tm in seq(1000)){
    tmp_dist <- snp_pos - snp_pos[tm]
    tmp_inds <- which((tmp_dist > 0) & (tmp_dist <= max_dist))
    if(length(tmp_inds) == 0){next}
    tmp_cors <- cor(t(vcf_df[tm,c(FIRST_SAMP:ncol(vcf_df))]), 
      t(vcf_df[tmp_inds, c(FIRST_SAMP:ncol(vcf_df))]),
      use = 'pairwise.complete.obs')
    tmp_r2 <- tmp_cors^2
    test_snp_vec <- rep(snp_names[tm], times = length(tmp_r2))
    comp_snp_vec <- snp_names[tmp_inds]
    comp_dist <- tmp_dist[tmp_inds]
    r2_list[[tm]] <- list()
    r2_list[[tm]][['test_snp_vec']] <- test_snp_vec
    r2_list[[tm]][['comp_snp_vec']] <- comp_snp_vec
    r2_list[[tm]][['comp_dist']] <- comp_dist
    r2_list[[tm]][['r2']] <- tmp_r2
####
  }
  tot_list <- list()
  tot_list[['test_snp']] <- unlist(lapply(r2_list, 
    function(x) x[['test_snp_vec']]))
  tot_list[['comp_snp']] <- unlist(lapply(r2_list, 
    function(x) x[['comp_snp_vec']]))
  tot_list[['comp_dist']] <- unlist(lapply(r2_list, 
    function(x) x[['comp_dist']]))
  tot_list[['r2']] <- unlist(lapply(r2_list, function(x) x[['r2']]))
  return(tot_list)
}

test <- calc_snp_r2(vcf_df = toy_df, max_dist = 500)



```


