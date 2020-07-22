# Steps and notes for running CLUMPP on STRUCTURE results for `geosamps`

## Test running using different types of genotypes
* all samples with 4 lines
* all samples with 2 lines
* all samples with 4 lines but 4X (disomic) SG have NA (-9) in lines 3 and 4

## Overview
* This round of STRUCTURE run using tetrasomic genotypes for all samples 
* Based on the delta-K results, k=2 is best results, but k=3 is decent

## CLUMPP Process for K=2
### Concatenate K=2 STRUCTURE results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

cat geo_v2_2.*_q > geo_v2_k2_combo_q_tmp
```
### Format input file for K=2 CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k2_combo_q_tmp'

in_res <- read.table(in_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(in_res)/3

samp_int_vec <- rep(seq(nsamps), times = 3)

third_vec <- rep('(x)', times = nrow(in_res))

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k2_combo.clumpp.indfile`

### Run CLUMPP program for K=2
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2
qsub geo_v2_k2_clumpp.sh
```
### Process output from CLUMPP for K=2
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k2_combo_q_tmp'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k2_combo.clumpp.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('clumpp.out', 'clumpp.processed', res_file)
write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k2_combo.clumpp.processed`


## CLUMPP Process for K=3 
### Concatenate STRUCTURE results for K=3
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

cat geo_v2_3.*_q > geo_v2_k3_combo_q_tmp
```
### Format input file for CLUMPP for K = 3
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
### Run CLUMPP program for k=3
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2
qsub geo_v2_k3_clumpp.sh

```
### Process output from CLUMPP for K=3
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

