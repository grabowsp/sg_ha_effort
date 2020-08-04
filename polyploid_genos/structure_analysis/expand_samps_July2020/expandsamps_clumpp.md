# Steps and notes for running CLUMPP on STRUCTURE results for `expand_samps`

## Overview
* Run CLUMPP on STRUCTURE results from 4 types of genotypes:
  * alltet = all samples get tetrasomic genotypes
  * alldip = all samples get disomic genotypes
  * partNA = 4X switchgrass get 2 lines of NA
  * pseudohap = all samples get haplotype genotype
* 785 samples

## CLUMPP for `alltet`
### K=2
#### Concatenate results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

cat expandgeo_alltet_2*_q > expandgeo_alltet_combo_q_tmp
```
#### Format input file for CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo_q_tmp'

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo.clumpp.indfile`

#### Run CLUMPP program for K=2
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet 

bash
source activate structure_env

PARAMFILE=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile

ln -s $PARAMFILE ./paramfile

DATA_PRE=expandgeo_alltet_combo.clumpp
INDFILE=$DATA_PRE'.indfile'
MISCFILE=$DATA_PRE'.miscfile'

N_K=2
N_C=785
N_R=3

OUTFILE=$DATA_PRE'_k'$N_K'.out'

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $N_K -c $N_C -r $N_R
```
* `expandgeo_alltet_combo.clumpp_k2.out`

#### Process output from CLUMPP for K=2
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo.clumpp.indfile'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo.clumpp_k2.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('.clumpp_', '_', gsub('out', 'clumpp.processed', res_file), 
  fixed = T)
write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo_k2.clumpp.processed`

### K=3
#### Concatenate results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

cat expandgeo_alltet_3*_q > expandgeo_alltet_combo_q_tmp
```
#### Format input file for CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo_q_tmp'

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo.clumpp.indfile`
#### Run CLUMPP program for K=3
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet 

bash
source activate structure_env

PARAMFILE=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile

ln -s $PARAMFILE ./paramfile

DATA_PRE=expandgeo_alltet_combo.clumpp
INDFILE=$DATA_PRE'.indfile'
MISCFILE=$DATA_PRE'.miscfile'

N_K=3
N_C=785
N_R=3

OUTFILE=$DATA_PRE'_k'$N_K'.out'

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $N_K -c $N_C -r $N_R
```
* `expandgeo_alltet_combo.clumpp_k3.out`

#### Process output from CLUMPP for K=3
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo.clumpp.indfile'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo.clumpp_k3.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('.clumpp_', '_', gsub('out', 'clumpp.processed', res_file),
  fixed = T)

write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo_k3.clumpp.processed`

