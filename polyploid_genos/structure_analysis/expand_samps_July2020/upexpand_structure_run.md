# Steps and notes for running STRUCTURE on 'upexpand', the upland samples
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

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/combo.sub.polyploid.CDS.upexpand.genlight.rds

OUT_NAME=upexpand_25k.strucgenos.txt
N_SNPS=25000
GENO_TYPE=everything

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```
### Transfer genotype files to HA
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

scp upexpand_25k.strucgenos.txt* grabowsk@pants.hagsc.org:/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand
```

## Run `all_tet` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/tet_generic.mainparams ./tet_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate submission scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip

bash
for SUB_FILE in expandgeo_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`
  sed 's/expand_geo\/all_dip/upexpand\/all_tet/g; s/expandgeo_25k.strucgenos.txt_all_dip/upexpand_25k.strucgenos.txt_all_tet/g; s/OUT_NAME=expandgeo_alldip_/OUT_NAME=upexpand_alltet_/g; s/N_SAMPS=785/N_SAMPS=287/g; s/N_SNPS=25000/N_SNPS=25000/g; s/MAIN_PARAM=dip_generic.mainparams/MAIN_PARAM=tet_generic.mainparams/g; s/-N alldip_/-N alltet_/g' \
  $SUB_FILE > ../../upexpand/all_tet/upexpand_structure_alltet_$FILE_SUF;
  done
```
### Submit Jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet

bash
for KT in {1..10};
  do
for KR in {1..3};
    do
    qsub upexpand_structure_alltet_k$KT'.'$KR'.sh';
  done;
done
```

## Run `all_dip` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_dip

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/dip_generic.mainparams ./dip_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate Submission Scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet

bash

for SUB_FILE in upexpand_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`
  sed 's/upexpand\/all_tet/upexpand\/all_dip/g; s/upexpand_25k.strucgenos.txt_all_tet/upexpand_25k.strucgenos.txt_all_dip/g; s/OUT_NAME=upexpand_alltet_/OUT_NAME=upexpand_alldip_/g; s/tet_generic.mainparams/dip_generic.mainparams/g' \
  $SUB_FILE > ../all_dip/upexpand_structure_alldip_$FILE_SUF;
  done
```
### Submit Jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_dip

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub upexpand_structure_alldip_k$KT'.'$KR'.sh';
  done;
done
```

## Run `part_NA` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/part_NA 

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/tet_generic.mainparams ./tet_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate Submission Scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet

bash

for SUB_FILE in upexpand_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`
  sed 's/upexpand\/all_tet/upexpand\/part_NA/g; s/upexpand_25k.strucgenos.txt_all_tet/upexpand_25k.strucgenos.txt_part_NA/g; s/OUT_NAME=upexpand_alltet_/OUT_NAME=upexpand_partNA_/g; s/tet_generic.mainparams/tet_generic.mainparams/g; s/-N alldip_/-N partNA_/g' \
  $SUB_FILE > ../part_NA/upexpand_structure_partNA_$FILE_SUF;
  done
```
### Submit Jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/part_NA

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub upexpand_structure_partNA_k$KT'.'$KR'.sh';
  done;
done
```

## Run `pseudohap` K=1 to K=10
### Soft link the parameter files
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/pseudohap 

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/hap_generic.mainparams ./hap_generic.mainparams
ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```
### Generate Submission Scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet

bash

for SUB_FILE in upexpand_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`
  sed 's/upexpand\/all_tet/upexpand\/pseudohap/g; s/upexpand_25k.strucgenos.txt_all_tet/upexpand_25k.strucgenos.txt_pseudohap/g; s/OUT_NAME=upexpand_alltet_/OUT_NAME=upexpand_pseudohap_/g; s/tet_generic.mainparams/hap_generic.mainparams/g; s/-N alldip_/-N pseudohap_/g' \
  $SUB_FILE > ../pseudohap/upexpand_structure_pseudohap_$FILE_SUF;
  done
```
### Submit Jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/pseudohap

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub upexpand_structure_pseudohap_k$KT'.'$KR'.sh';
  done;
done
```


