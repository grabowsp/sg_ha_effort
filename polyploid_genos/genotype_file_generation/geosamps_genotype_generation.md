# Steps for the generation and manipulation of genotypes for `geosamps`

## Overview
* Generate CDS VCFs for `geosamps`
  * Filter by MAF
* Split CDS VCFs into 100k-line subfiles
* Generate header file to use with subfiles
* Generate `genlight` objects for each chromsome using adegenet
* Generate genome-wide, subsampled `genlight` object

## Location of Output Files
* Directory with `geosamps` genotype files
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps`

## location of input files
* chromosome cds vcfs for all libraries
  *  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
* Library names of the 772 "geographic" libraries I consider "natural" samples
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_772_lib_names.txt`

## Generate CDS VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

sbatch generate_Chr01_Chr05_CDS_geosamps_vcf.sh
sbatch generate_Chr06_Chr09_CDS_geosamps_vcf.sh
```
#### Example of submit file
* From `generate_Chr01_Chr05_CDS_geosamps_vcf.sh`
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

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/
SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_772_lib_names.txt
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/
SAMPSET=geosamps

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
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

gunzip -kc Chr01K.polyploid.CDS.geosamps.vcf.gz | \
split -l 100000 -d - Chr01K.polyploid.CDS.geosamps.vcf_

gunzip -kc Chr01N.polyploid.CDS.geosamps.vcf.gz | \
split -l 100000 -d - Chr01N.polyploid.CDS.geosamps.vcf_

for CN in {02..09};
do
gunzip -kc Chr$CN'K.polyploid.CDS.geosamps.vcf.gz' | \
split -l 100000 -d - Chr$CN'K.polyploid.CDS.geosamps.vcf_';
done

for CN in {02..09};
do
gunzip -kc Chr$CN'N.polyploid.CDS.geosamps.vcf.gz' | \
split -l 100000 -d - Chr$CN'N.polyploid.CDS.geosamps.vcf_';
done

```

## Generate VCF header file for subfiles
### Output
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samp/CDS.geosamps.vcf.header.txt`
### Commands to generate Output
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps
head -5 Chr01K.polyploid.CDS.geosamps.vcf_00 | tail -1 > \
CDS.geosamps.vcf.header.txt
```

## Generate `genlight` Object
* use 0.002 as MAF cutoff
  * between 3 and 6 copies of minor allele at MAF
### Submit job
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

sbatch gen_geo_genlight.sh
```

## Generate genome-wide subsampled `genlight` object
### Output
* on Cori
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/combo.sub.polyploid.CDS.geosamps.genlight.rds`
    * 772 samples
    * 592,150 SNPs
### Select Subsampling level
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# R
library(adegenet)

# adjust these variables
###
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/'

goal_nsnps <- 5e5
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
[1] 0.05910562
# I deciede to use 0.07
```
### Generate object
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

FILE_SUB=Chr*geosamps.genlight.rds

OUT_NAME=combo.sub.polyploid.CDS.geosamps.genlight.rds

PER_SUBSAMP=0.07

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/subsample_genlight.r \
$DATA_DIR '*'$FILE_SUB $OUT_NAME $PER_SUBSAMP
```


