# Template to use for the generation and manipulation of genotypes

* Change variables or lines that begin with #$&#

## Repeated Variables
* Data directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps`
* Name for sample set when used in file names (SAMPSET variable)
  * `lowgeosamps`
  
## Overview
* Generate CDS VCFs for
  * Filter by MAF
* Split CDS VCFs into 100k-line subfiles
* Generate header file to use with subfiles
* Generate `genlight` objects for each chromsome using adegenet
* Generate genome-wide, subsampled `genlight` object

## Location of Output Files
* Directory with `lowgeosamps` genotype files
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps`
 
## location of input files
* chromosome cds vcfs for all libraries
  *  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
* Library names of the lowland "geographic" libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/lowland_libs_June2020.txt`

## Generate CDS VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps

sbatch generate_Chr01K_01N_CDS_lowgeosamps_vcf.sh
sbatch generate_Chr02_Chr05_CDS_lowgeosamps_vcf.sh
sbatch generate_Chr06_Chr09_CDS_lowgeosamps_vcf.sh
```
#### Example of submit file
* From `sbatch generate_Chr02_Chr05_CDS_lowgeosamps_vcf.sh`
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

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/
SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/lowland_libs_June2020.txt
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/

for CHROM in 02K 02N 03K 03N 04K 04N 05K 05N;
  do
  vcftools --gzvcf $DATA_DIR'Chr'$CHROM'.polyploid.CDS.vcf.gz' \
  --max-missing 0.8 \
  --stdout --keep $SAMP_FILE --recode --recode-INFO-all | \
  gzip -c > \
  $OUT_DIR'Chr'$CHROM'.polyploid.CDS.lowgeosamps.vcf.gz';
  done
```

## Generate 100k-line subfiles
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps

for CN in {01..09};
do
gunzip -kc Chr$CN'K.polyploid.CDS.lowgeosamps.vcf.gz' | \
split -l 100000 -d - Chr$CN'K.polyploid.CDS.lowgeosamps.vcf_';
done

for CN in {01..09};
do
gunzip -kc Chr$CN'N.polyploid.CDS.lowgeosamps.vcf.gz' | \
split -l 100000 -d - Chr$CN'N.polyploid.CDS.lowgeosamps.vcf_';
done
```

## Generate VCF header file for subfiles
### Output
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/CDS.lowgeosamps.vcf.header.txt` 
### Commands to generate Output
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps
head -5 Chr01K.polyploid.CDS.lowgeosamps.vcf_00 | tail -1 > \
CDS.lowgeosamps.vcf.header.txt
```

## Generate `genlight` Object
* use 0.005 as MAF cutoff
  * between 3 and 6 copies of minor allele at MAF
### Submit job
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps

sbatch gen_lowgeo_genlight.sh
```

## Generate genome-wide subsampled `genlight` object
### Overview
* Goal is ~370k SNPs
  * keeps roughly the same SNP/sample ratio as with the geo samp results
  * To get this:
    * get the SNP/sample number from the geo_samps subsampled genlight object
    * multiply the SNP/sample number by the number of upland samples
* Subsample 7.1% of the SNPs
  * to get this:
    * estimate total SNPs in vcf's by counting up number of sub_vcfs
    * for 1 chromosome, estimate the % of SNPs that made it from the sub_vcfs
to the genlight.rds
    * divide goal amount of SNPs by (total SNPs in sub_vcfs * % that made it
into the genlight objects)
### Output
* on Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/combo.sub.polyploid.CDS.lowgeosamps.genlight.rds`
  * 482 genotypes
  * 385,251 SNPs
### Generate object
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/

FILE_SUB=lowgeosamps.genlight.rds

#OUT_NAME=Combo.sub.polyploid.CDS.lowgeosamps.genlight.rds
OUT_NAME=combo.sub.polyploid.CDS.lowgeosamps.genlight.rds

PER_SUBSAMP=0.071

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/subsample_genlight.r \
$DATA_DIR '*'$FILE_SUB $OUT_NAME $PER_SUBSAMP
```


