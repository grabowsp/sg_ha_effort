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

