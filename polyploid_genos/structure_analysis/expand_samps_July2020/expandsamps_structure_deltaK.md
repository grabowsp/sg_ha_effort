# Steps and looking at deltaK and other metrics used to evaluate K as 
#    demonstrated in the Evanno paper

## Test robustness of results by  using different types of genotypes
* all samples with 4 lines
* all samples with 2 lines
* all samples with 4 lines but 4X (disomic) SG have NA (-9) in lines 3 and 4
* all samples with pseudhaploid genotypes

## `expand_geo` `all_tet` Results 
### Extract values from STRUCTURE output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

bash

OUT_PRE=expandgeo_alltet

for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' $OUT_PRE'_'$KT'.'$KR'_f' >> \
      $OUT_PRE'_LnProbData.txt';
    done;
  done;
```
### Calculate delta-K and other Evanno et al metrics
```
bash
source activate R_analysis

# ADJUST THESE VARIABLES TO CORRECT NAMES AND VALUES
PARENT_DIR=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/

FILE_SHORT=expandgeo_alltet_LnProbData.txt

NKS=10
NREPS=3

SAMP_SET_NAME=expandgeo_alltet

####
LN_DATA_FILE=$PARENT_DIR$FILE_SHORT

Rscript /home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/gen_deltaK_plot.r $LN_DATA_FILE $NKS $NREPS $SAMP_SET_NAME

```
### Output files
* Interpretation:
  * K=2 is best; K=3 is okay, other K's don't show much support
* Figure
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_LnProbData_deltaK.pdf`
* Table
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_LnProbData_deltaK.txt`

## `expand_geo` `all_dip` Results
### Extract values from STRUCTURE output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip

bash

OUT_PRE=expandgeo_alldip

for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' $OUT_PRE'_'$KT'.'$KR'_f' >> \
      $OUT_PRE'_LnProbData.txt';
    done;
  done;
```
### Calculate delta-K and other Evanno et al metrics
```
bash
source activate R_analysis

# ADJUST THESE VARIABLES TO CORRECT NAMES AND VALUES
PARENT_DIR=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip/

FILE_SHORT=expandgeo_alldip_LnProbData.txt

NKS=10
NREPS=3

SAMP_SET_NAME=expandgeo_alldip

####
LN_DATA_FILE=$PARENT_DIR$FILE_SHORT

Rscript /home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/gen_deltaK_plot.r $LN_DATA_FILE $NKS $NREPS $SAMP_SET_NAME

```

## `expand_geo` `part_NA` Results
### Extract values from STRUCTURE output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA

bash

OUT_PRE=expandgeo_partNA

for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' $OUT_PRE'_'$KT'.'$KR'_f' >> \
      $OUT_PRE'_LnProbData.txt';
    done;
  done;
```
### Calculate delta-K and other Evanno et al metrics
```
bash
source activate R_analysis

# ADJUST THESE VARIABLES TO CORRECT NAMES AND VALUES
PARENT_DIR=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA/

FILE_SHORT=expandgeo_partNA_LnProbData.txt

NKS=10
NREPS=3

SAMP_SET_NAME=expandgeo_partNA

####
LN_DATA_FILE=$PARENT_DIR$FILE_SHORT

Rscript /home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/gen_deltaK_plot.r $LN_DATA_FILE $NKS $NREPS $SAMP_SET_NAME

```
### Output
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA/expandgeo_partNA_LnProbData_deltaK.pdf`
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA/expandgeo_partNA_LnProbData_deltaK.txt`


######################
## `expand_geo` `pseudohap` Results
### Extract values from STRUCTURE output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap

bash

OUT_PRE=expandgeo_pseudohap

for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' $OUT_PRE'_'$KT'.'$KR'_f' >> \
      $OUT_PRE'_LnProbData.txt';
    done;
  done;
```
### Calculate delta-K and other Evanno et al metrics
```
bash
source activate R_analysis

# ADJUST THESE VARIABLES TO CORRECT NAMES AND VALUES
PARENT_DIR=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap/

FILE_SHORT=expandgeo_pseudohap_LnProbData.txt

NKS=10
NREPS=3

SAMP_SET_NAME=expandgeo_partNA

####
LN_DATA_FILE=$PARENT_DIR$FILE_SHORT

Rscript /home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/gen_deltaK_plot.r $LN_DATA_FILE $NKS $NREPS $SAMP_SET_NAME

```
### Output
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap/expandgeo_pseudohap_LnProbData_deltaK.pdf`
* `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap/expandgeo_pseudohap_LnProbData_deltaK.txt`




