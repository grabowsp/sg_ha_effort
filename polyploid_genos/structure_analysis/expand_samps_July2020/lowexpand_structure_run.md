# Steps and notes for running STRUCTURE on 'lowexpand', the upland samples
#   in the expanded geographic 4X and 8X samples

## Test robustness of results by  using different types of genotypes
* all samples with 4 lines
* all samples with 2 lines
* all samples with 4 lines but 4X (disomic) SG have NA (-9) in lines 3 and 4
* all samples with pseudhaploid genotypes

### Location of Generic param files
* NERSC Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/structure_tests`
* HA Directory
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc`
#### mainparam files
* Tetrasomic file
  * `tet_generic.mainparams`
* Disomic file
  * `dip_generic.mainparams`
* Pseudohaploid file
  * `hap_generic.mainparams`
#### extraparams
* `generic.extraparams`

## Generate Genotypes
* 25k SNPs, 10k burnin, 30k run
### Make Genotypes with R
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/lowexpand_samps

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/lowexpand_samps/combo.sub.polyploid.CDS.lowexpand.genlight.rds

OUT_NAME=lowexpand_25k.strucgenos.txt
N_SNPS=25000
GENO_TYPE=everything

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```
### Transfer genotype files to HA
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/lowexpand_samps

scp lowexpand_25k.strucgenos.txt* grabowsk@pants.hagsc.org:/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/lowexpand
```

## Run `all_tet` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/lowexpand/all_tet

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/tet_generic.mainparams ./tet_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```

### Generate submission scripts
# CONTINUE FROM HERE - EDIT FROM N_SAMPS and beyond
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet 

bash
for SUB_FILE in upexpand_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $3}'`
  sed 's/struc\/upexpand/struc\/lowexpand/g; s/upexpand_25k.strucgenos.txt_all_tet/lowexpand_25k.strucgenos.txt_all_tet/g; s/OUT_NAME=upexpand_alltet_/OUT_NAME=lowexpand_alltet_/g; s/N_SAMPS=772/N_SAMPS=785/g; s/N_SNPS=20000/N_SNPS=25000/g; s/-N struc_/-N alltet_struc_/g' \
  $SUB_FILE > ../expand_geo/all_tet/expandgeo_structure_alltet_$FILE_SUF;
  done
```
### Add random number step to script
* had issue with same seed used for multiple reps, resulting in identical results
* had to add a random number generating step to avoid that
```
bash

for TK in 2 3;
  do
  sed 's/alltet_struc_k1.1/'alltet_struc_k1."$TK"'/g; s/TEST_REP=1/'TEST_REP="$TK"'/g ' \
expandgeo_structure_alltet_k1.1.sh > expandgeo_structure_alltet_k1.$TK'.sh';
  done

for TK in {2..10};
  do
  for TR in {1..3};
    do
    sed 's/alltet_struc_k1/'alltet_struc_k"$TK"'/g; s/TEST_K=1/'TEST_K="$TK"'/g' \
expandgeo_structure_alltet_k1.$TR'.sh' > \
expandgeo_structure_alltet_k$TK'.'$TR'.sh';
    done;
  done
```

### Submit Jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

bash
for KT in {1..10};
  do
# for KR in {1..3}; I had to re-run reps 2 and 3 because of SEED issues
  for KR in {2..3};
    do
    qsub expandgeo_structure_alltet_k$KT'.'$KR'.sh';
  done;
done
```

## Run `all_dip` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/dip_generic.mainparams ./dip_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate Submission Scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

bash
for SUB_FILE in expandgeo_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`
  sed 's/struc\/expand_geo\/all_tet/struc\/expand_geo\/all_dip/g; s/expandgeo_25k.strucgenos.txt_all_tet/expandgeo_25k.strucgenos.txt_all_dip/g; s/OUT_NAME=expandgeo_alltet_/OUT_NAME=expandgeo_alldip_/g; s/-N alltet_struc_/-N alldip_struc_/g; s/tet_generic.mainparams/dip_generic.mainparams/g' \
  $SUB_FILE > ../all_dip/expandgeo_structure_alldip_$FILE_SUF;
  done
```
### Submit Jobs
* NEED TO DO THIS
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub expandgeo_structure_alldip_k$KT'.'$KR'.sh';
  done;
done
```

## Run `part_NA` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/tet_generic.mainparams ./tet_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate Submission Scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

bash
for SUB_FILE in expandgeo_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`
  sed 's/struc\/expand_geo\/all_tet/struc\/expand_geo\/part_NA/g; s/expandgeo_25k.strucgenos.txt_all_tet/expandgeo_25k.strucgenos.txt_part_NA/g; s/OUT_NAME=expandgeo_alltet_/OUT_NAME=expandgeo_partNA_/g; s/-N alltet_struc_/-N partNA_struc_/g' \
  $SUB_FILE > ../part_NA/expandgeo_structure_partNA_$FILE_SUF;
  done
```
### Submit Jobs
* NEED TO DO THIS

## Run `pseudohap` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/hap_generic.mainparams ./hap_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate Submission Scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

bash
for SUB_FILE in expandgeo_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`
  sed 's/struc\/expand_geo\/all_tet/struc\/expand_geo\/pseudohap/g; s/expandgeo_25k.strucgenos.txt_all_tet/expandgeo_25k.strucgenos.txt_pseudohap/g; s/OUT_NAME=expandgeo_alltet_/OUT_NAME=expandgeo_pseudohap_/g; s/-N alltet_struc_/-N pseudohap_struc_/g; s/tet_generic.mainparams/hap_generic.mainparams/g' \
  $SUB_FILE > ../pseudohap/expandgeo_structure_pseudohap_$FILE_SUF;
  done
```
### Submit Jobs
* NEED TO DO THI


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
