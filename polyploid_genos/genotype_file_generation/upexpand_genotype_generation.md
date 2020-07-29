# Steps for generation and manipulation of `upexpand` genotypes

* Change variables or lines that begin with #$&#

## Repeated Variables
* Data directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps`
* Name for sample set when used in file names (SAMPSET variable)
  * `upexpand`
  
## Overview
* Generate CDS VCFs for
  * Filter by MAF
* Split CDS VCFs into 100k-line subfiles
* Generate header file to use with subfiles
* Generate `genlight` objects for each chromsome using adegenet
* Generate genome-wide, subsampled `genlight` object

## Location of Output Files
* Directory with `upexpand` genotype files
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps`

## location of input files
* chromosome cds vcfs for all libraries
  *  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
 * Library names of the 287 "expandgeo" libraries I consider upland
   * `/global/cscratch1/sd/grabowsp/sg_ploidy/expandupland_libs_July2020.txt`

## Generate CDS VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

sbatch generate_Chr01_Chr05_CDS_upexpand_vcf.sh
sbatch generate_Chr06_Chr09_CDS_upexpand_vcf.sh
```
#### Example of submit file
* From `generate_Chr01_Chr05_CDS_upexpand_vcf.sh`
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

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/expandupland_libs_July2020.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/

SAMPSET=upexpand

for CHROM in 01K 01N 02K 02N 03K 03N 04K 04N 05K 05N;
  do
  vcftools --gzvcf $DATA_DIR'Chr'$CHROM'.polyploid.CDS.vcf.gz' \
  --max-missing 0.8 \
  --stdout --keep $SAMP_FILE --recode --recode-INFO-all | \
  gzip -c > \
  $OUT_DIR'Chr'$CHROM'.polyploid.CDS.'$SAMPSET'.vcf.gz';
  done
```

## Generate 100k-line subfiles
```
# Go to directory with VCFs to split
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

SAMPSET=upexpand

for CN in {01..09};
do
gunzip -kc Chr$CN'K.polyploid.CDS.'$SAMPSET'.vcf.gz' | \
split -l 100000 -d - Chr$CN'K.polyploid.CDS.'$SAMPSET'.vcf_';
done

for CN in {01..09};
do
gunzip -kc Chr$CN'N.polyploid.CDS.'$SAMPSET'.vcf.gz' | \
split -l 100000 -d - Chr$CN'N.polyploid.CDS.'$SAMPSET'.vcf_';
done
```

## Generate VCF header file for subfiles
### Output
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/CDS.upexpand.vcf.header.txt`
### Commands to generate Output
```
# go to directory with genotype files
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

SAMPSET=upexpand

head -5 Chr01K.polyploid.CDS.$SAMPSET'.vcf_00' | tail -1 > \
CDS.$SAMPSET'.vcf.header.txt'
```

## Generate `genlight` Object
#$&# * use 0.004 as MAF cutoff
  * between 3 and 5 copies of minor allele at MAF
### Submit job
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

sbatch gen_upexpand_genlight.sh
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

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/

cd $DATA_DIR

HEADER_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/CDS.upexpand.vcf.header.txt

SAMPSET=upexpand

MAF_CUT=0.004

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
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/combo.sub.polyploid.CDS.upexpand.genlight.rds`
    * 287 samples
    * 225,937 SNPs
### Select Subsampling level
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# R
library(adegenet)

# adjust these variables
###
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/'

goal_nsnps <- 287/772*5e5
# 186000
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
# [1] 0.04113422
# I deciede to use 0.05
```
### Generate object
```
# load conda environment
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# go to directory where want to run script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/ 

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/

# how sample set is used in file names
SAMPSET=upexpand

# percentage of total SNPs to select
PER_SUBSAMP=0.05

# FILE_SUB may need to be hardcoded (without the $SAMPSET variable) - try that
#   if there is an issue
FILE_SUB=$SAMPSET'.genlight.rds'
OUT_NAME=combo.sub.polyploid.CDS.$SAMPSET'.genlight.rds'

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/subsample_genlight.r \
$DATA_DIR '*'$FILE_SUB $OUT_NAME $PER_SUBSAMP
```


