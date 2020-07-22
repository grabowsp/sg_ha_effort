# Steps and notes for running STRUCTURE on the `geosamps`

## Test running using different types of genotypes
* all samples with 4 lines
* all samples with 2 lines
* all samples with 4 lines but 4X (disomic) SG have NA (-9) in lines 3 and 4
### Generate test genotypes
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.rds

OUT_NAME=test_100_alltet.strucgenos.txt
N_SNPS=100
GENO_TYPE=all_tet

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE

OUT_NAME=test_100_alldip.strucgenos.txt
GENO_TYPE=all_dip

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE

OUT_NAME=test_100_partNA.strucgenos.txt
GENO_TYPE=part_NA

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```
### Make Generic param files
* Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/structure_tests`
#### mainparam files
* Tetrasomic file
  * `tet_generic.mainparams`
* Disomic file
  * `dip_generic.mainparams`
#### extraparams
* `generic.extraparams`
### Run tests
```
module load python/3.7-anaconda-2019.07
source activate structure_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/structure_tests

structure -m tet_generic.mainparams -e generic.extraparams -K 3 -L 100 \
-N 826 -i alltet_test -o at_test_2

structure -m tet_generic.mainparams -e generic.extraparams -K 3 -L 100 \
-N 826 -i partNA_test -o pNA_test_1

structure -m dip_generic.mainparams -e generic.extraparams -K 3 -L 100 \
-N 826 -i alldip_test -o ap_test_1
```

## Full run for geo_samps
* 20k SNPs, 10k burnin, 30k run
### Generate Genotypes
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

#IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.rds

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/combo.sub.polyploid.CDS.geosamps.genlight.rds

#OUT_NAME=geosamps_20k_alltet.strucgenos.txt
OUT_NAME=geosamps_v2_20k_alltet.strucgenos.txt
N_SNPS=20000
GENO_TYPE=all_tet

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE

#OUT_NAME=geosamps_20k_alldip.strucgenos.txt
OUT_NAME=geosamps_v2_20k_alldip.strucgenos.txt
GENO_TYPE=all_dip

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```
#### Move genotypes to new directory
```
#cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_20k_struc_files
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_v2_struc_files

mv ../geosamps_v2_20k*strucgenos.txt .
cp ../structure_tests/*params .

tar -czvf geosamp_v2_files.tar.gz *params *strucgenos.txt
```
### Transfer to HA
```
#cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_20k_struc_files

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_v2_struc_files

scp geosamp_v2_files.tar.gz grabowsk@pants.hagsc.org:/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/

# on HA
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2
tar -xzvf geosamp_v2_files.tar.gz

#ln -s geosamps_20k_alldip.strucgenos.txt geo_20k_alldip.txt
#ln -s geosamps_20k_alltet.strucgenos.txt geo_20k_alltet.txt

```
### Run K=1 to K=10
#### Generate submission scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps
#cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

bash
for SUB_FILE in *sh;
  do
  sed 's/struc\/geo_samps/struc\/geo_v2/g; s/geosamps_20k/geosamps_v2_20k/g; s/OUT_NAME=geo20k_/OUT_NAME=geo_v2_/g; s/N_SAMPS=826/N_SAMPS=772/g' \
  $SUB_FILE > ../geo_v2/$SUB_FILE;
  done
n {1..10};
  do
  for KR in {1..3};
    do
    qsub upgeo_structure_k$KT'.'$KR'.sh';
  done; ``
#### Submit jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

bash
for KT in {1..10};
  do
  qsub geo_structure_k$KT'.1.sh';
  done;

bash
for KT in {1..10};
  do
  qsub geo_structure_k$KT'.2.sh';
  done;

for KT in {1..10};
  do
  qsub geo_structure_k$KT'.3.sh';
  done;
```

### Get Estimated Ln Prob of Data
```
#cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' geo_v2_$KT'.'$KR'_f' >> \
      geo_v2_LnProbData.txt;
    done;
  done;
```
### Estimate delta-K and other Evanno et al metrics
```
bash
source activate R_analysis

nKs <- 10
nreps <- 3

ln_prob_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_LnProbData.txt'

ln_prob <- read.table(ln_prob_in, header = F, sep = '\t', stringsAsFactors = F)

ln_prob_vals <- unlist(lapply(strsplit(ln_prob[,1], split = ' '), 
  function(x) as.numeric(rev(x)[1])))

ln_prob_df <- data.frame(K = rep(seq(nKs), times = 1, each = nreps), 
  rep = rep(seq(nreps), times = nKs), ln_prob = ln_prob_vals, 
  stringsAsFactors = F)

ln_prob_df$ord_1 <- NA
ln_prob_df$ord_2 <- NA
for(i in seq(nreps)){
  tmp_inds <- which(ln_prob_df$rep == i)
  tmp_vals <- ln_prob_df$ln_prob[tmp_inds]
  tmp_change_1 <- tmp_vals[2:nKs] - tmp_vals[1:(nKs-1)]
  change_1 <- c(NA, tmp_change_1)
  ln_prob_df$ord_1[tmp_inds] <- change_1
  tmp_change_2 <- abs(tmp_change_1[2:(nKs-1)] - tmp_change_1[1:(nKs-2)])
  change_2 <- c(NA, tmp_change_2, NA)
  ln_prob_df$ord_2[tmp_inds] <- change_2
}

mean_ln_prob <- tapply(ln_prob_df$ln_prob, ln_prob_df$K, mean)
sd_ln_prob <- tapply(ln_prob_df$ln_prob, ln_prob_df$K, sd)

mean_change1_K <- tapply(ln_prob_df$ord_1, ln_prob_df$K, mean)
sd_change1_K <- tapply(ln_prob_df$ord_1, ln_prob_df$K, sd)

mean_change2_K <- tapply(ln_prob_df$ord_2, ln_prob_df$K, mean)
sd_change2_K <- tapply(ln_prob_df$ord_2, ln_prob_df$K, sd)

delta_K <- mean_change2_K / sd_ln_prob

delta_df <- data.frame(K = seq(nKs), mean_ln_prob, sd_ln_prob, mean_change1_K,
  sd_change1_K, mean_change2_K, sd_change2_K, delta_K, stringsAsFactors = F)

library(ggplot2)
library(patchwork)

gg_lnprob <- ggplot(delta_df, aes(x = K, y = mean_ln_prob)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_ln_prob - sd_ln_prob, 
    ymax = mean_ln_prob + sd_ln_prob)) +
  xlab('K') +
  ylab('ln Prob Data') +
  ggtitle('Ln Prob Data for Geographic Samples (v2)')

gg_change1 <- ggplot(delta_df, aes(x = K, y = mean_change1_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change1_K - sd_change1_K, 
    ymax = mean_change1_K + sd_change1_K)) +
  xlab('K') +
  ylab("L'(K)") +
  ggtitle("L'(K) for Geographic Samples (v2)")

gg_change2 <- ggplot(delta_df, aes(x = K, y = mean_change2_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change2_K - sd_change2_K,
    ymax = mean_change2_K + sd_change2_K)) +
  xlab('K') +
  ylab("|L''(K)|") +
  ggtitle("|L''(K)| for Geographic Samples (v2)")

gg_deltaK <- ggplot(delta_df, aes(x = K, y = delta_K)) +
  geom_point() +
  geom_line() +
  xlab('K') +
  ylab("delta K") +
  ggtitle("delta K for Geographic Samples (v2)")

out_fig_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_structure_evanno.pdf'

pdf(out_fig_file, width = 12, height = 10)
(gg_lnprob + gg_change1) / (gg_change2 + gg_deltaK)
dev.off()

out_tab_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_structure_lnprob.txt'

write.table(ln_prob_df, file = out_tab_file, quote = F, sep= '\t', 
  row.names = F, col.names = T)

```
## Run CLUMPP on results
### K=3
#### Concatenate results from the 3 replicates
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

cat geo_v2_3.*_q > geo_v2_k3_combo_q_tmp
```
#### Format input file for CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo_q_tmp'

in_res <- read.table(in_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(in_res)/3

samp_int_vec <- rep(seq(nsamps), times = 3)

third_vec <- rep('(x)', times = nrow(in_res))
#third_vec <- rep(paste('(', in_res[,1], ')', sep = ''), times = 3)

pop_vec <- rep(NA, times = nrow(in_res))
fifth_vec <- rep(':', times = nrow(in_res))

out_res_tmp <- data.frame( C1=samp_int_vec, C2 = samp_int_vec, 
  C3 = third_vec, C4 = samp_int_vec, C5 = fifth_vec,
  stringsAsFactors = F)

out_res <- cbind(out_res_tmp, in_res[, c(2:ncol(in_res))])

out_file <- gsub('_q_tmp', '.clumpp.indfile', in_file)

write.table(out_res, out_file, quote=F, sep = '\t', row.names = F, 
  col.names = F)
```
#### Run CLUMPP
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2
qsub geo_v2_k3_clumpp.sh

```
#### Process output from CLUMPP
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo_q_tmp'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo.clumpp.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('clumpp.out', 'clumpp.processed', res_file)
write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
#### Make bargraph of CLUMPP results
```
bash
source activate R_analysis

library(ggplot2)

struc_func_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/struc_functions.r'
source(struc_func_file)

c_res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo.clumpp.processed'

c_res <- read.table(c_res_file, header = F, sep = '\t', stringsAsFactors = F)

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t',
  quote = "", comment.char = '$')

### evaluate what the groups are

pop_info_list <- list()

test_order_mat <- apply(c_res[, -1], 1, order, decreasing = T)

for(i in seq(ncol(c_res)-1)){
  test_pop <- i
  pop_name <- paste('Population', i)
  test_libs <- c_res[which(apply(test_order_mat, 2, 
    function(x) x[1] == test_pop)), 1]
  ploidy_inds <- which(ploidy$lib %in% test_libs)
  meta_inds <- which(meta$LIBRARY %in% test_libs)
  tmp_list <- list()
  tmp_list[['Cloroplast Ecotype']] <- table(ploidy$ECOTYPE_SNP_CHLR
    [ploidy_inds])
  tmp_list[['Old Subpop']] <- table(ploidy$SUBPOP_SNP[ploidy_inds])
  tmp_list[['Ploidy']] <- table(ploidy$total_ploidy_2[ploidy_inds])
  tmp_list[['State']] <- table(meta$STATE[meta_inds])
  pop_info_list[[pop_name]] <- tmp_list
}

# pop 1 = EC lowland
# pop 2 = upland
# pop 3 = TX lowland

pop_order <- c(2,1,3)

######

test_samp_orders <- order_struc_samps(clumpp_result_df = c_res, 
  pop_order = pop_order, zero_cut = 0.1)

#####

res_ord <- c_res[test_samp_orders, ]

samp_vec <- rep(res_ord[,1], times = ncol(c_res)-1)
group_vec <- rep(seq(ncol(c_res)-1), each = nrow(res_ord))

res_vec <- c()
for(j in c(2:ncol(res_ord))){
  res_vec <- c(res_vec, res_ord[,j])
}

plot_df <- data.frame(LIB = samp_vec, GROUP = group_vec, MEMB = res_vec,
  stringsAsFactors = F)

plot_df$LIB <- factor(plot_df$LIB, levels = res_ord[,1])

plot_df$GROUP <- factor(plot_df$GROUP, levels = pop_order)

group_col_vec <- c()
group_col_vec['1'] <- 'yellow3'
group_col_vec['2'] <- 'blue2'
group_col_vec['3'] <- 'red2'

group_lab_vec <- c()
group_lab_vec['1'] <- 'EC'
group_lab_vec['2'] <- 'MW'
group_lab_vec['3'] <- 'TX'
group_lab_vec <- group_lab_vec[pop_order]

#group_palette <- scale_colour_manual(name = 'K_Group', values = group_col_vec)

gg_bar <- ggplot(plot_df, aes(y = MEMB, x = LIB, fill = GROUP)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual('Group', values = group_col_vec, 
    labels = group_lab_vec) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = 'Group Membership')

out_file <- gsub('clumpp.processed', 'STRUCTURE.memb.pdf', c_res_file)

pdf(out_file, width = 12*6, height = 4)
gg_bar
dev.off()


```
