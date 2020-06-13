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

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.rds

OUT_NAME=geosamps_20k_alltet.strucgenos.txt
N_SNPS=20000
GENO_TYPE=all_tet

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE

OUT_NAME=geosamps_20k_alldip.strucgenos.txt
GENO_TYPE=all_dip

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```
#### Move genotypes to new directory
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_20k_struc_files

mv ../geosamps_20k*strucgenos.txt .
cp ../structure_tests/*params .

tar -czvf geosamp_20k_files.tar.gz *params *strucgenos.txt
```
### Transfer to HA
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_20k_struc_files

scp geosamp_20k_files.tar.gz grabowsk@pants.hagsc.org:/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps/

# on HA
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps
tar -xzvf geosamp_20k_files.tar.gz

ln -s geosamps_20k_alldip.strucgenos.txt geo_20k_alldip.txt
ln -s geosamps_20k_alltet.strucgenos.txt geo_20k_alltet.txt

```
### Run K=1 to K=8
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_samps
qsub geo_structure_k1.1.sh

bash
for KT in {2..8};
  do
  qsub geo_structure_k$KT'.1.sh';
  done;

bash
for KT in {1..8};
  do
  qsub geo_structure_k$KT'.2.sh';
  done;

for KT in {1..8};
  do
  qsub geo_structure_k$KT'.3.sh';
  done;
```
