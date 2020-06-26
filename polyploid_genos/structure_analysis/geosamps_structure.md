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
### Run K=1 to K=8
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
```
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
### CONTINUE RE-RUN FROM HERE
### Get Estimated Ln Prob of Data
```
#cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' geo20k_$KT'.'$KR'_f' >> \
      geo20k_LnProbData.txt;
    done;
  done;
```
* in R
```
bash
source activate R_analysis

ln_prob_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps/geo20k_LnProbData.txt'

ln_prob <- read.table(ln_prob_in, header = F, sep = '\t', stringsAsFactors = F)

ln_prob_vals <- unlist(lapply(strsplit(ln_prob[,1], split = ' '), 
  function(x) as.numeric(rev(x)[1])))

ln_prob_df <- data.frame(K = rep(seq(8), times = 1, each = 3), 
  rep = rep(seq(3), times = 8), ln_prob = ln_prob_vals, stringsAsFactors = F)

ln_prob_df$ord_1 <- NA
ln_prob_df$ord_2 <- NA
for(i in seq(1:3)){
  tmp_inds <- which(ln_prob_df$rep == i)
  tmp_vals <- ln_prob_df$ln_prob[tmp_inds]
  tmp_change_1 <- tmp_vals[2:8] - tmp_vals[1:7]
  change_1 <- c(NA, tmp_change_1)
  ln_prob_df$ord_1[tmp_inds] <- change_1
  tmp_change_2 <- abs(tmp_change_1[2:7] - tmp_change_1[1:6])
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

delta_df <- data.frame(K = seq(8), mean_ln_prob, sd_ln_prob, mean_change1_K,
  sd_change1_K, mean_change2_K, sd_change2_K, delta_K, stringsAsFactors = F)

library(ggplot2)
library(patchwork)

gg_lnprob <- ggplot(delta_df, aes(x = K, y = mean_ln_prob)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_ln_prob - sd_ln_prob, 
    ymax = mean_ln_prob + sd_ln_prob)) +
  xlab('K') +
  ylab('ln Prob Data') +
  ggtitle('Ln Prob Data for Geographic Samples (v1)')

gg_change1 <- ggplot(delta_df, aes(x = K, y = mean_change1_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change1_K - sd_change1_K, 
    ymax = mean_change1_K + sd_change1_K)) +
  xlab('K') +
  ylab("L'(K)") +
  ggtitle("L'(K) for Geographic Samples (v1)")

gg_change2 <- ggplot(delta_df, aes(x = K, y = mean_change2_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change2_K - sd_change2_K,
    ymax = mean_change2_K + sd_change2_K)) +
  xlab('K') +
  ylab("|L''(K)|") +
  ggtitle("|L''(K)| for Geographic Samples (v1)")

gg_deltaK <- ggplot(delta_df, aes(x = K, y = delta_K)) +
  geom_point() +
  geom_line() +
  xlab('K') +
  ylab("delta K") +
  ggtitle("delta K for Geographic Samples (v1)")

out_fig_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps/geo_v1_structure_evanno.pdf'

pdf(out_fig_file, width = 12, height = 10)
(gg_lnprob + gg_change1) / (gg_change2 + gg_deltaK)
dev.off()

out_tab_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps/geo_v1_structure_lnprob.txt'

write.table(ln_prob_df, file = out_tab_file, quote = F, sep= '\t', 
  row.names = F, col.names = T)

```
