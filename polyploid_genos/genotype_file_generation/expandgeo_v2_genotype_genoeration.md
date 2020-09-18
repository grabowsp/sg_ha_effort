# Template to use for the generation and manipulation of genotypes

* Change variables or lines that begin with #$&#

## Repeated Variables
* Data directory
#$&#  * ``
* Name for sample set when used in file names (SAMPSET variable)
#$&#  * ``
  
## Overview
* Generate CDS VCFs for
  * Filter by MAF
* Split CDS VCFs into 100k-line subfiles
* Generate header file to use with subfiles
* Generate `genlight` objects for each chromsome using adegenet
* Generate genome-wide, subsampled `genlight` object

## Location of Output Files
* Directory with `expand_v2` genotype files
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2`

## location of input files
* chromosome cds vcfs for all libraries - USE SORTED VERSIONS
  *  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
#$&# * Library names of the 785 "expand" libraries I consider "natural" samples
#$&#  * `/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_expand_785_lib_names.txt`

## Generate CDS VCFs and split into 100k lines
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

sbatch generate_Chr01_Chr05_CDS_expand_v2_vcf.sh
sbatch generate_Chr06_Chr09_CDS_expand_v2_vcf.sh
```
#### Example of submit file
* From `generate_Chr01_Chr05_CDS_expand_v2_vcf.sh`
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

# Directory where job info (and new VCFs) will be saved
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

# Directory with base CDS VCFs
DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/

# File of Library Names to keep
SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_expand_785_lib_names.txt

# Directory where new VCFs will be saved
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/

# sample set name to be used in fame files
SAMPSET=expandv2

for CHROM in 01K 01N 02K 02N 03K 03N 04K 04N 05K 05N;
  do
  vcftools --gzvcf $DATA_DIR'Chr'$CHROM'.polyploid.CDS.sorted.vcf.gz' \
  --max-missing 0.8 \
  --stdout --keep $SAMP_FILE --recode --recode-INFO-all | \
  gzip -c > \
  $OUT_DIR'Chr'$CHROM'.polyploid.CDS.'$SAMPSET'.vcf.gz';
  gunzip -kc Chr$CHROM'.polyploid.CDS.'$SAMPSET'.vcf.gz' | \
  split -l 100000 -d - Chr$CHROM'.polyploid.CDS.'$SAMPSET'.vcf_';
  done
```

## Generate VCF header file for subfiles
### Output
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/CDS.expandv2.vcf.header.txt`
### Commands to generate Output
```
# go to directory with genotype files
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/ 

SAMPSET=expandv2

head -24 Chr01K.polyploid.CDS.$SAMPSET'.vcf_00' | tail -1 > \
CDS.$SAMPSET'.vcf.header.txt'
```

## Generate `genlight` Object
* use 0.002 as MAF cutoff
  * between 3 and 6 copies of minor allele at MAF
### Submit job
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

# NEED TO RUN THIS
#$&# sbatch gen_expandv2_genlight.sh
```
### Example of submit script
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/r_adegenet_env

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/

cd $DATA_DIR

HEADER_FILE=$DATA_DIR'CDS.expandv2.vcf.header.txt'

SAMPSET=expandv2

MAF_CUT=0.002

for CHR_N in {01..09};
  do
  for CHR_T in K N;
    do
    TEST_CHR=$CHR_N$CHR_T
    SEARCH_STRING=Chr$TEST_CHR'.polyploid.CDS.'$SAMPSET'.vcf_'
    Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/make_Chr_genlight_objs.r \
    $DATA_DIR $SEARCH_STRING'*' $HEADER_FILE $MAF_CUT;
    done;
  done;
```

## Generate genome-wide subsampled `genlight` object
### Output
* on Cori
#$&#  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/combo.sub.polyploid.CDS.geosamps.genlight.rds`
#$&#    * 772 samples
#$&#    * 592,150 SNPs
### Select Subsampling level
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# R
library(adegenet)

# adjust these variables
###
#$&# data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/'

#$&# goal_nsnps <- 5e5
###

file_vec <- system(paste('ls ', data_dir, 'Chr*genlight.rds', sep = ''), 
  intern = T)

n_snp_vec <- c()
for(i in file_vec){
  print(i)
  tmp_in <- readRDS(i)
  tmp_nsnps <- nLoc(tmp_in)
  n_snp_vec <- c(n_snp_vec, tmp_nsnps)
}

tot_snps <- sum(n_snp_vec)

goal_sub <- goal_nsnps/tot_snps

print(goal_sub)
#$&# [1] 0.05910562
# I deciede to use 0.07
```
### Generate object
```
# load conda environment
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# go to directory where want to run script
#$&# cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

# directory with chromosome genlight objects
#$&# DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

# how sample set is used in file names
#$&# SAMPSET=geosamps

# percentage of total SNPs to select
#$&# PER_SUBSAMP=0.07

# FILE_SUB may need to be hardcoded (without the $SAMPSET variable) - try that
#   if there is an issue
FILE_SUB=$SAMPSET'.genlight.rds'
OUT_NAME=combo.sub.polyploid.CDS.$SAMPSET'.genlight.rds'

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/subsample_genlight.r \
$DATA_DIR '*'$FILE_SUB $OUT_NAME $PER_SUBSAMP
```


