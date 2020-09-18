# CLUMPP for comparing the STRUCTURE results for `expandsamps`

## K=2 Results
### Concatenate results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis

bash

TEST_K=2

ALLTET_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo_k'$TEST_K'.clumpp.processed'
ALLDIP_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip/expandgeo_alldip_combo_k'$TEST_K'.clumpp.processed'
PARTNA_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA/expandgeo_partNA_combo_k'$TEST_K'.clumpp.processed'
HAP_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap/expandgeo_pseudohap_2_combo_k'$TEST_K'.clumpp.processed'

OUT_FILE='expandgeo_combo_res_k'$TEST_K'_q_tmp'

cat $ALLTET_RES $ALLDIP_RES $PARTNA_RES $HAP_RES > $OUT_FILE

```
### Format input for CLUMPP
```
bash
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis/expandgeo_combo_res_k2_q_tmp'

in_res <- read.table(in_file, header = F, stringsAsFactors = F)

nreps <- 4
nsamps <- nrow(in_res)/nreps

samp_int_vec <- rep(seq(nsamps), times = nreps)

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis/expandgeo_combo_res_k2.clumpp.indfile`

### Run CLUMPP program for K=2
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis

bash
source activate structure_env

PARAMFILE=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile_2

ln -s $PARAMFILE ./paramfile

DATA_PRE=expandgeo_combo_res_k2.clumpp
INDFILE=$DATA_PRE'.indfile'
MISCFILE=$DATA_PRE'.miscfile'

N_K=2
N_C=785
N_R=4

OUTFILE=$DATA_PRE'_k'$N_K'.out'

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $N_K -c $N_C -r $N_R

for CF in permute_clumpp*;
do
  mv $CF $CF'_K'$N_K;
done
```

###########

## K=3 Results
### Concatenate results
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis

bash

TEST_K=3

ALLTET_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_combo_k'$TEST_K'.clumpp.processed'
ALLDIP_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip/expandgeo_alldip_combo_k'$TEST_K'.clumpp.processed'
PARTNA_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA/expandgeo_partNA_combo_k'$TEST_K'.clumpp.processed'
HAP_RES='/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap/expandgeo_pseudohap_'$TEST_K'_combo_k'$TEST_K'.clumpp.processed'

OUT_FILE='expandgeo_combo_res_k'$TEST_K'_q_tmp'

cat $ALLTET_RES $ALLDIP_RES $PARTNA_RES $HAP_RES > $OUT_FILE

```
### Format input for CLUMPP
```
bash
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis/expandgeo_combo_res_k3_q_tmp'

in_res <- read.table(in_file, header = F, stringsAsFactors = F)

nreps <- 4
nsamps <- nrow(in_res)/nreps

samp_int_vec <- rep(seq(nsamps), times = nreps)

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
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis/expandgeo_combo_res_k3.clumpp.indfile`

### Run CLUMPP program for K=3
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/combo_analysis

bash
source activate structure_env

PARAMFILE=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile_2

ln -s $PARAMFILE ./paramfile

DATA_PRE=expandgeo_combo_res_k3.clumpp
INDFILE=$DATA_PRE'.indfile'
MISCFILE=$DATA_PRE'.miscfile'

N_K=3
N_C=785
N_R=4

OUTFILE=$DATA_PRE'_k'$N_K'.out'

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $N_K -c $N_C -r $N_R

for CF in permute_clumpp*;
do
  mv $CF 'K'$N_K'_'$CF;
done
```


