# Steps and notes for running CLUMPP on structure results for `up_expand`

## Overview 
* Run CLUMPP on STRUCTURE results for 4 types of genotypes
  * k=2 and k=5
* 287 samples

## CLUMPP for `alltet`
```
bash
source activate structure_env

cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet

#k=2
RES_PRE=upexpand_alltet_
K_VAL=2
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

#k=5

RES_PRE=upexpand_alltet_
K_VAL=5
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R
```
* output files
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet/upexpand_alltet_2_combo.clumpp_k2.out`
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet/upexpand_alltet_5_combo.clumpp_k5.out`

## CLUMPP for `alldip`
```
bash
source activate structure_env

cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_dip

#k=2
RES_PRE=upexpand_alldip_
K_VAL=2
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

#k=5
RES_PRE=upexpand_alldip_
K_VAL=5
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R
```
* Output files
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_dip/upexpand_alldip_2_combo.clumpp_k2.out`
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_dip/upexpand_alldip_5_combo.clumpp_k5.out`

## CLUMPP for `partNA`
```
bash
source activate structure_env

cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/part_NA

#k=2
RES_PRE=upexpand_partNA_
K_VAL=2
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

#k=5
RES_PRE=upexpand_partNA_
K_VAL=5
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R
```
* Output files
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/part_NA/upexpand_partNA_2_combo.clumpp_k2.out`
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/part_NA/upexpand_partNA_5_combo.clumpp_k5.out`

## CLUMPP for 'pseudohap'
```
bash
source activate structure_env

cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/pseudohap

#k=2
RES_PRE=upexpand_pseudohap_
K_VAL=2
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

#k=5
RES_PRE=upexpand_pseudohap_
K_VAL=5
N_C=287
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R
```
* Output
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/pseudohap/upexpand_pseudohap_2_combo.clumpp_k2.out`
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/pseudohap/upexpand_pseudohap_5_combo.clumpp_k5.out`

