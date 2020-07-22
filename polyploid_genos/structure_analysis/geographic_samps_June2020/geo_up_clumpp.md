# Steps and notes for running CLUMPP on STRUCTURE results for `geo_up`

## Overview
* This round of STRUCTURE run using tetrasomic genotypes for all samples 
* Based on the delta-K results, k=2 is best results, but k=3 and k=4  is decent
* 277 samples

## CLUMPP Process for K=2
### Concatenate K=2 STRUCTURE results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo

cat upgeo_2.*_q > upgeo_k2_combo_q_tmp
```
### Format input file for K=2 CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo_q_tmp'

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo.clumpp.indfile`

### Run CLUMPP program for K=2
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo

bash
source activate structure_env

PARAMFILE=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile

ln -s $PARAMFILE ./paramfile

DATA_PRE=upgeo_k2_combo.clumpp
INDFILE=$DATA_PRE'.indfile'
OUTFILE=$DATA_PRE'.out'
MISCFILE=$DATA_PRE'.miscfile'

N_K=2
N_C=277
N_R=3

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $N_K -c $N_C -r $N_R
```
* `upgeo_k2_combo.clumpp.out`

### Process output from CLUMPP for K=2
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo_q_tmp'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo.clumpp.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('clumpp.out', 'clumpp.processed', res_file)
write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k2_combo.clumpp.processed`

################

## CLUMPP Process for K=3
### Concatenate K=3 STRUCTURE results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo

cat upgeo_3.*_q > upgeo_k3_combo_q_tmp
```
### Format input file for K=3 CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo_q_tmp'

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo.clumpp.indfile`

### Run CLUMPP program for K=3
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo

bash
source activate structure_env

PARAMFILE=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile

ln -s $PARAMFILE ./paramfile

DATA_PRE=upgeo_k3_combo.clumpp
INDFILE=$DATA_PRE'.indfile'
OUTFILE=$DATA_PRE'.out'
MISCFILE=$DATA_PRE'.miscfile'

N_K=3
N_C=277
N_R=3

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $N_K -c $N_C -r $N_R
```
* `upgeo_k3_combo.clumpp.out`

### Process output from CLUMPP for K=3
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo_q_tmp'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo.clumpp.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('clumpp.out', 'clumpp.processed', res_file)
write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k3_combo.clumpp.processed`

################

## CLUMPP Process for K=4
### Concatenate K=3 STRUCTURE results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo

cat upgeo_4.*_q > upgeo_k4_combo_q_tmp
```
### Format input file for K=3 CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo_q_tmp'

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo.clumpp.indfile`

### Run CLUMPP program for K=4
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo

bash
source activate structure_env

PARAMFILE=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile

ln -s $PARAMFILE ./paramfile

DATA_PRE=upgeo_k4_combo.clumpp
INDFILE=$DATA_PRE'.indfile'
OUTFILE=$DATA_PRE'.out'
MISCFILE=$DATA_PRE'.miscfile'

N_K=4
N_C=277
N_R=3

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $N_K -c $N_C -r $N_R
```
* `upgeo_k4_combo.clumpp.out`

### Process output from CLUMPP for K=4
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo_q_tmp'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo.clumpp.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('clumpp.out', 'clumpp.processed', res_file)
write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_k4_combo.clumpp.processed`


