# Steps and notes for running STRUCTURE on 'expandsamps', the expanded 
#    4X and 8X samples

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

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/combo.sub.polyploid.CDS.expandgeosamps.genlight.rds

OUT_NAME=expandgeo_25k.strucgenos.txt
N_SNPS=25000
GENO_TYPE=everything

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```
### Transfer genotype files to HA
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps

scp expandgeo_25k.strucgenos.txt* grabowsk@pants.hagsc.org:/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo
```

## Run `all_tet` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/tet_generic.mainparams ./tet_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate submission scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2
q
bash
for SUB_FILE in geo_structure*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $3}'`
  sed 's/struc\/geo_v2/struc\/expand_geo\/all_tet/g; s/geosamps_v2_20k_alltet.strucgenos.txt/expandgeo_25k.strucgenos.txt_all_tet/g; s/OUT_NAME=geo_v2_/OUT_NAME=expandgeo_alltet_/g; s/N_SAMPS=772/N_SAMPS=785/g; s/N_SNPS=20000/N_SNPS=25000/g; s/-N struc_/-N alltet_struc_/g' \
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
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/part_NA

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub expandgeo_structure_partNA_k$KT'.'$KR'.sh';
  done;
done
```


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
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub expandgeo_structure_pseudohap_k$KT'.'$KR'.sh';
  done;
done
```

