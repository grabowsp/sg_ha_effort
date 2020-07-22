# Steps and notes for running STRUCTURE on the `geosamps`

## Test running using different types of genotypes
* all samples with 4 lines
* all samples with 2 lines
* all samples with 4 lines but 4X (disomic) SG have NA (-9) in lines 3 and 4

## Generic param files
* Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/structure_tests`
#### mainparam files
* Tetrasomic file
  * `tet_generic.mainparams`
* Disomic file
  * `dip_generic.mainparams`
#### extraparams
* `generic.extraparams`

## Full run for `low_geo`
* 20k SNPs, 10k burnin, 30k run
  * 482 samples
### Generate Genotypes
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/combo.sub.polyploid.CDS.lowgeosamps.genlight.rds

OUT_NAME=lowgeo_20k_alltet.strucgenos.txt
N_SNPS=20000
GENO_TYPE=all_tet

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE

OUT_NAME=lowgeo_20k_alldip.strucgenos.txt
GENO_TYPE=all_dip

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```
#### Move genotypes to new directory
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/low_geo_struc_files

mv ../lowgeo_20k*strucgenos.txt .
cp ../../geo_samps/structure_tests/*params .

tar -czvf lowgeo_struc_files.tar.gz *params *strucgenos.txt
```
### Transfer to HA
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/low_geo_struc_files

scp lowgeo_struc_files.tar.gz grabowsk@pants.hagsc.org:/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo/

# on HA
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo
tar -xzvf lowgeo_struc_files.tar.gz
```
### Run K=1 to K=10
#### Generate submission scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

bash
for SUB_FILE in *sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $3}'`
  sed 's/struc\/geo_v2/struc\/low_geo/g; s/geosamps_v2_20k/lowgeo_20k/g; s/OUT_NAME=geo_v2_/OUT_NAME=lowgeo_/g; s/N_SAMPS=772/N_SAMPS=482/g' \
  $SUB_FILE > ../low_geo/lowgeo_structure_$FILE_SUF;
  done
```
#### Submit jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub lowgeo_structure_k$KT'.'$KR'.sh';
  done;
done

```

### Get Estimated Ln Prob of Data
```
#cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' lowgeo_$KT'.'$KR'_f' >> \
      lowgeo_LnProbData.txt;
    done;
  done;
```
* in R
```
bash
source activate R_analysis

nKs <- 10
nreps <- 3

ln_prob_in <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo/lowgeo_LnProbData.txt'

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
  ggtitle('Ln Prob Data for Lowland Geographic Samples')

gg_change1 <- ggplot(delta_df, aes(x = K, y = mean_change1_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change1_K - sd_change1_K, 
    ymax = mean_change1_K + sd_change1_K)) +
  xlab('K') +
  ylab("L'(K)") +
  ggtitle("L'(K) for LowlandGeographic Samples")

gg_change2 <- ggplot(delta_df, aes(x = K, y = mean_change2_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change2_K - sd_change2_K,
    ymax = mean_change2_K + sd_change2_K)) +
  xlab('K') +
  ylab("|L''(K)|") +
  ggtitle("|L''(K)| for Lowland Geographic Samples")

gg_deltaK <- ggplot(delta_df, aes(x = K, y = delta_K)) +
  geom_point() +
  geom_line() +
  xlab('K') +
  ylab("delta K") +
  ggtitle("delta K for Lowland Geographic Samples")

out_fig_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo/low_geo_structure_evanno.pdf'

pdf(out_fig_file, width = 12, height = 10)
(gg_lnprob + gg_change1) / (gg_change2 + gg_deltaK)
dev.off()

out_tab_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/low_geo/low_geo_structure_lnprob.txt'

write.table(ln_prob_df, file = out_tab_file, quote = F, sep= '\t', 
  row.names = F, col.names = T)

```

