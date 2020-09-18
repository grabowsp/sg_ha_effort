# Instructions and Explanations for pipeline for STRUCTURE Analysis

## Steps
1) Generate genotype file
2) Make or copy parameter files
3) Copy genotype and parameter files to HA Cluster
4) Generate submit scripts
5) Submit STRUCTURE jobs
6) Evaluate best K(s) via deltaK
7) CLUMPP to aggregate results from replicates
8) Generate barplot figure
9) Generate map pieplot figure

## Generate Genotype File
### Overview
* Start with genlight object that has SNPs already filtered for MAF
  * Look at "Manipulate Genomtype Objects" Section for details on generating
the input file
* Choose the number of SNPs to include in the STRUCTURE genotype file
  * 20k-50k should be more than adequate
* Decide on the "style" of genotypes
  * See below
### Required Resources
* `r_adegenet_env` conda environment
  * includes R and adegenet and ggplot2 (and other) packages
* `generate_structure_input.r` script from github
  * `/PATHTOREPOSITORY/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r`
### Genotype "styles"
* STRUCTURE assumes all samples have the same ploidy, so need to choose one
of 4 solutions:
  * all_tet = all samples get 4 lines; 4X (disomic) samples have their 
genotypes repeated to make their 4 lines
  * all_dip = all samples get 2 lines; 8X (tetrasomic) samples have 2 alleles
randomly chosen from their genotypes to make their 2 lines
  * part_NA = all samples get 4 lines; 4X (disomic) samples have NA's for 
their 3rd and 4th lines
  * pseudohap = all samples get 1 line; all samples have 1 allele randomly 
chosen from their genotype
* if chose "everything", then make all 4 styles of genotype files
### Run Script
* Example for running script on Cori to generate genotypes
* Is fast enough to run in interactive session
```
# load conda environment
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# go to directory where file will be saved
cd /DATA/DIRECTORY

# the genlight genotype file
IN_FILE=/PATH/TO/GENLIGHT_FILE.genlight.rds

# the name of the outputted genotype file
OUT_NAME=sampleset_nSNPs.strucgenos.txt

# the number of SNPs to have in outputted genotype file
N_SNPS=25000

# choose the style of the outputted genotype file; 
#   options are: all_tet, all_dip, part_NA, pseudohap, everything
GENO_TYPE=all_tet

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/generate_structure_input.r \
$IN_FILE $OUT_NAME $N_SNPS $GENO_TYPE
```

## Make or Copy Parameter Files
### Overview
* STRUCTURE needs a 'mainparams' and 'extraparams' file to run
* There are generic 'mainparams' and 'extraparams' files that can be used 
for most uses, and command line flags can be used to update the necessary
parameters
* 
### Selecting 'mainparams' file
* `tet_generic.mainparams` is used for all_tet and part_NA genotype styles
* `dip_generic.mainparams` is used for all_dip genotype style
* `hap_generic.mainparams` is used for pseudohap genotype style
### Location of files
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

## Copy genotype and parameter files to HA Cluster
### Overview
* Some STRUCTURE runs take longer than the allowed time on NERSC, so run
STRUCTURE on the HudsonAlpha GSC cluster instead
* Generally generate genotype file at NERSC and transfer genotype file to 
HA GSC cluster
* Need to make soft-link of parameter files into the directory that will 
contain the genotype file and STRUCTURE outputs
  * Genotype file and links to param files need to be in same directory for
scripts to work properly
  * I usually make a separate directory for each sample set to make it easier  
for troubleshooting and keeping track of results
### Example of transfering genotype file
```
NERSC_GENO_FILE=/PATH/TO/GENOFILE.genlight.rds
HA_DIR=/PATH/ON/HAGSC/SAMPSET_DIR

scp $NERSC_GENO_FILE grabowspgrabowsk@pants.hagsc.org:$HA_DIR
```
### Example of making soft links
```
HA_DIR=/PATH/ON/HAGSC/SAMPSET_DIR

cd $HA_DIR

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/tet_generic.mainparams ./tet_generic.mainparams

ln -s /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/generic.extraparams ./generic.extraparams
```

## Generate Submit Script
### Description
* Submit script includes:
  * the location of input files (and where the output should be saved)
  * name of genotype file
  * name of structure parameter files
  * what to name the generated files
  * the number of samples
  * the number of SNPs
* The script should include a step for generating the starting seed
  * otherwise all jobs that start at same time use same seed
* Easiest way to generate new submit scripts is to copy-and-edit existing files
  * example of how to copy-and-edit is below
### Example of script
```
#!/bin/bash

#$ -q "mecat.q"
#$ -pe smp 1
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae
#$ -V
#$ -cwd
#$ -N pseudohap_struc_k2.1

TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
OUT_DIR=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap/

cd $TMP_DIR

source /home/raid2/LINUXOPT/miniconda2/bin/activate /home/grabowsky/.conda/envs/structure_env

IN_FILE=expandgeo_25k.strucgenos.txt_pseudohap
MAIN_PARAM=hap_generic.mainparams
EX_PARAM=generic.extraparams

TEST_K=2
TEST_REP=1

OUT_NAME=expandgeo_pseudohap_$TEST_K.$TEST_REP

N_SAMPS=785
N_SNPS=25000

R_SEED=$(($RANDOM + $TEST_REP))

cp $OUT_DIR$IN_FILE $TMP_DIR'/struc_in'
cp $OUT_DIR$MAIN_PARAM $TMP_DIR'/struc_mp'
cp $OUT_DIR$EX_PARAM $TMP_DIR'/struc_ep'

structure -m struc_mp -e struc_ep -K $TEST_K -L $N_SNPS -N $N_SAMPS \
-i struc_in -o $OUT_NAME -D $R_SEED

rm ./struc_in
rm ./struc_mp
rm ./struc_ep

mv ./seed.txt ./seed_$TEST_K'_'$TEST_REP'.txt'

/usr/bin/rsync -avuP $TMP_DIR/* $OUT_DIR

/bin/rm -rf $TMP_DIR
```
### How to copy-and-edit previous scripts
* Steps
  * Go to directory with sh submit scripts for K's and reps
  * Cycle through all scripts and use sed to edit the files
  * Save edited files to appropriate directory
#### Example
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

## Submit Jobs
```
cd /DIRECTORY/WITH/FILES

bash
for KT in {1..10};
  do
  for KR in {1..3};
    do
    qsub expandgeo_structure_alltet_k$KT'.'$KR'.sh';
  done;
done
```

## Calculate delta-K
* Calculate delta-K statistics from Evanno et al. for determing best K
* STEPS
  * Concatenate the Ln Prob Data info from all structure results
  * Use Rscript to generate Evanno delta-K statistics, output figure and table
* Steps are put into wrapper script
  * `~/sg_ha_effort/polyploid_genos/structure_analysis/run_deltaK_steps.sh`
### Using wrapper script
* Requires 3 inputs
  1. prefix of the structure results; should include the underscore at end
    * ex: expandgeo_pseudohap_
  2. The maximum K run for structure
    * ex: 10
  3. The number of replicate structure runs for each K
    * ex: 3
#### Example
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap

RES_PRE=expandgeo_pseudohap_
MAX_K=10
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_deltaK_steps.sh $RES_PRE $MAX_K $N_R
```

## Run CLUMPP
* Use the program CLUMPP to average the runs that use the same K
* Steps:
  1. Concatenate structure output for K value
  2. Generate CLUMPP input
  3. Run CLUMPP
  4. Process CLUMPP output
* Steps run as wrapper script
* Requires 4 inputs
  1. Prefix of the structure results; best to include underscore at end
    * ex: expandgeo_pseudohap_
  2. K value
    * ex: 2
  3. Number of samples
    * ex: 785
  4. Number of structure replicate runs for K
    * ex: 3
### Example of running script
```
bash

cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/pseudohap

#k=2
RES_PRE=expandgeo_pseudohap_
K_VAL=2
N_C=785
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

```

