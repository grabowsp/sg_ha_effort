# Script for identifying ploidy outliers from the pseudohaploid distances

# source activate r_phylo

# LOAD PACKAGES #

### LOAD DATA ###
data_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/',
  'filtered_vcfs/', 'pseudohapdists_v01_total.rds', sep = '')
data <- readRDS(data_file)

ploidy_info_file <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/', 
  'pseudohap/sg_ploidy_results_v2.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

### SET OUTPUTS ###
# manhattan-distance plots
#man_PCo1v2_fig <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/',
#  'pseudohap/filtered_vcfs/',
#  'pseudohapdists_v01_manhattan_PCo1vPCo2.pdf', sep = '')

# SET VARIABLES #

############

# Set Colors
table(ploidy_info$SUBPOP_SNP)
table(ploidy_info$ECOTYPE_SNP_CHLR)
table(ploidy_info$total_ploidy)

subpop_col_vec <- rep(NA, times = length(table(ploidy_info$SUBPOP_SNP)))
names(subpop_col_vec) <- names(table(ploidy_info$SUBPOP_SNP))
subpop_col_vec['?Subpop'] <- 'black'
subpop_col_vec['Eastcoast'] <- 'orange2'
subpop_col_vec['Eastcoast_Admixed'] <- 'orangered3'
subpop_col_vec['Gulfcoast'] <- 'springgreen3'
subpop_col_vec['Gulfcoast_Admixed'] <- 'darkgreen'
subpop_col_vec['Midwest'] <- 'royalblue2'
subpop_col_vec['Midwest_Admixed'] <- 'blue2'
subpop_col_vec['Texas'] <- 'violetred2'
subpop_col_vec['Texas_Admixed'] <- 'darkorchid2'

subpop_palette <- scale_colour_manual(name = 'Subpop', values = subpop_col_vec)

totPloid_col_vec <- rep(NA, times = length(table(ploidy_info$total_ploidy)))
names(totPloid_col_vec) <- names(table(ploidy_info$total_ploidy))
totPloid_col_vec['4X'] <- 'red2'
totPloid_col_vec['8X'] <- 'blue2'

ploidy_palette <- scale_colour_manual(name = 'Ploidy', values = totPloid_col_vec)

# Make Manhattan-distance PCoA plots 

man_dist_mat <- data[[2]][-902,-902]

man_cmd <- cmdscale(man_dist_mat, k = 200)
man_tot_var <- sum(apply(man_cmd, 2, var))
man_per_var <- (apply(man_cmd, 2, var)/man_tot_var)*100

man_dist_label_names <- gsub('^X', '', colnames(man_dist_mat))
man_dist_label_names[which(man_dist_label_names == '9001.3.BN389.69S')] <- '9001-3 BN389-69S'
man_dist_label_names[which(man_dist_label_names == '9020.5')] <- '9020-5'

meta_ord <- c()
for(DL in man_dist_label_names){
  meta_ord <- c(meta_ord, which(ploidy_info$PLANT_ID == DL))
}

subpop_plot_cols <- ploidy_info$SUBPOP_SNP[meta_ord]
for(SP in names(subpop_col_vec)){
  tmp_inds <- which(subpop_plot_cols == SP)
  subpop_plot_cols[tmp_inds] <- subpop_col_vec[SP]
}

totPloid_plot_cols <- ploidy_info$total_ploidy[meta_ord]
for(SP in names(totPloid_col_vec)){
  tmp_inds <- which(totPloid_plot_cols == SP)
  totPloid_plot_cols[tmp_inds] <- totPloid_col_vec[SP]
}

man_df <- data.frame(man_cmd, stringsAsFactors = F)
colnames(man_df) <- paste('PCo_', seq(ncol(man_df)), sep = '')

man_df$samp <- man_dist_label_names
man_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
man_df$totPloid <- ploidy_info$total_ploidy[meta_ord]


###################
# Euclidean distances

euc_dist_mat <- data[[3]][-902,-902]

euc_cmd <- cmdscale(euc_dist_mat, k = 200)
euc_tot_var <- sum(apply(euc_cmd, 2, var))
euc_per_var <- (apply(euc_cmd, 2, var)/euc_tot_var)*100

euc_dist_label_names <- gsub('^X', '', colnames(euc_dist_mat))
euc_dist_label_names[which(euc_dist_label_names == '9001.3.BN389.69S')] <- '9001-3 BN389-69S'
euc_dist_label_names[which(euc_dist_label_names == '9020.5')] <- '9020-5'

euc_df <- data.frame(euc_cmd, stringsAsFactors = F)
colnames(euc_df) <- paste('PCo_', seq(ncol(euc_df)), sep = '')

euc_df$samp <- euc_dist_label_names
euc_df$SUBPOP <- ploidy_info$SUBPOP_SNP[meta_ord]
euc_df$totPloid <- ploidy_info$total_ploidy[meta_ord]

######
# Find outliers

# PCo1-Lowland 8X
low_8X <- intersect(which(euc_df$PCo_1 < 0), which(euc_df$totPloid == '8X'))
maybe_low_8X <- setdiff(
  intersect(which(euc_df$PCo_1 < 1250), which(euc_df$totPloid == '8X')),
  low_8X)

odd_up_8X <- intersect(
  intersect(which(euc_df$PCo_4 < -250), which(euc_df$PCo_1 > 2500)),
  which(euc_df$totPloid == '8X'))


euc_df$samp[low_8X]
# 19 total

# [1] "J192.L1"   
# IEYM; Outlier - nQuire 8X, MNP 4X
# NC; EastCoast GP
# PCo1 = -2790.183; PCo2 = -2033.917

# "J021.A"    
# IHYU
# TX; Own cluster with J021.B (below), but near GulfCoast and TX_admixed
# PCo1 = -862.049; PCo2 = 1984.251

# "J462.C"    
# IIAK
# Louisiana, clusters with some GulfCoast Admixed
# PCo1 = -329.4186; PCo2 = 1271.319

# "J021.B"    
# IIBT
# TX; clusters with J021.A (above) near GulfCoast and TX_admixed
# PCo1 = -884.2043; PCo2 = 1977.77

# "J462.B"    
# IIDN
# Louisiana, clusters with some GulfCoast Admixed
# PCo1 = -203.2655; PCo2 = 1269.107

# "NFGA16_02"
# PZXO - Outlier, nQuire 8X; MNP 4X
# NC, clusters close to TX and TX_admixed; maybe a sample switch?;
# PCo1 = -2507.592; PCo2 = 2032.792
# The other samples that are supposed to come from the same population do 
#   not cluster with this or together

# [7] "2302"      
# WZGG - Outlier, nQuire 8X; MNP 4X
# Kanlow, clusters pretty closely with Texas
# PCo1 = -2221.461; PCo2 = 4448.141

# "2732"      
# XXAY - Outlier: nQuire 8X, MNP 4X
# Kanlow, clusters with Texas but slightly outside the cluster
# PCo1 = -1913.338; PCo2 = 4252.269

# "J249.A"    
# IELY - Outlier: nQuire 6X
# TX, clusters pretty strongly with TX
# PCo1 = -1990.89; PCo2 = 3483.804

# "WBC3a"     
# ABHM - Outlire; nQuire 8X, MNP 4X (strongly)
# TX, clusters with Texas 
# PCo1 = -2094.213; PCo2 = 3665.348

# "J302.A"    
# ASHOA - Outlier; nQuire 8X, MNP 8X (CDS), MNP 4X (genic)
# TX, clusters close to TX_admix and GC_Admix
# PCo1 = -1432.592; PCo2 = 2389.572

# "J615.A"   
# INFT - Outlier: nQuire 6X
# FL, GC_Admix; Clusters with TX_admix
# PCo1 = -1865.616; PCo2 = 2571.523

#[13] "J179.L1"   
# IEXZ
# FL, clusers close to EC_admix and GC_admix
# PCo1 = -2256.296; PCo2 = -971.6153

# "J180.L2"   
# IEYD
# FL, EastCoast; clusters with a bunch of EastCoast including 4X from same pop
# PCo1 = -2699.502; PCo2 =  -1599.151

# "J188.L1"   
# IEYI
# FL, GulfCoast Subpop, Clusters with GC and GC_admix
# PCo1 = -1950.557; PCo2 =  861.2388 

# "J530.B"    
# IFBB
# NY; EC genepool; Clusters sort of all by itself near EC_admix sample
# PCo1 = -436.7661; PCo2 = -1288.419

# "J462.A"    
# ITAT
# Louisiana, clusters with some GC_admix
# PCo1 = -307.6137; PCo2 = 1224.844 

# "J464.A"   
# ITAU, Outlier - nQuire 8X, MNP 4X
# TX, GC_admix, clusters with some other GC_admix
# PCo1 = -548.7651; PCo2 =  1332.142

#[19] "J439.A"   
# IIKL
# Cave-in-Rock? ; Clusters right in the middle of EastCoast gene pool, maybe
#  a sample switch?
# PCo1 = -2747.183; PCo2 = -2391.417

euc_df$samp[maybe_low_8X]
# 5 samples

# [1] "J579.A" 
# IICP
# Georgia; in a cluster with other non-assigned samps; in NJ tree, is near
#   TX_admix, EC_admix, and MW_admix samples
# PCo1 = 790.2757; PCo2 = -233.2951

# "J579.C" 
# IICX
# Georgia, in a cluster with J589.A (above) and other non-assigned samps
# PCo1 = 706.2056; PCo2 = -286.9408

# "2785b"  
# WZGH
# Kanlow x Summer hybrid; Clusters with a bunch of TX_admix samples
# PCo1 = 741.9168; PCo2 = 2400.332

# "J199.A" 
# IELU
# Georgia, clusters by itself but kind of close to the other GA samples in this
#   group of 8X samples in PCoA and NJ tree
# PCo1 = 134.8305; PCo2 =  -417.1597

# "J579.B"
# IRUX
# Georgia; in cluster with other samples from same population
# PCo1 = 762.7978; PCo2 = -277.3008

euc_df$samp[odd_up_8X]
# 4 samples

# [1] "J045.A" 
# IIDF
# Summer (?), is kind of by itself in the plot
# PCo1 =3997.249; PCo4 = -702.0549

# "2719"   
# XXAZ
# Summer cross; clusters with 4X samples
# PCo1 = 4541.763; PCo4 = -726.8848

# "J400.A" 
# IGRS
# NY (Buffalo?); Clusters close but not exactly within 4X samples
# PCo1 = 4377.015; PCo4 = -452.3318

# "J148.A"
# IJAF
# North Dacotah; Clusters with 4X samples, also in NJ tree
# PCo1 = 4600.994; PCo4 = -910.3999

quit(save = 'no')

