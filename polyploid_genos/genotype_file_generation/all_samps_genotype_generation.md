# Steps for the generation and manipulation of genotypes for `all_samps`

## Overview
* Generate CDS VCFs for `all_samps` 
  * Filter by MAF
* Split CDS VCFs into 100k-line subfiles
* Generate header file to use with subfiles
* Did not go further (yet) for this sample set

## Location of Output Files
* Directory with `all_samps` genotype files
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/all_samps`

## Location of Input Files
* Chromosome CDS VCFs for all samples
  *  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
* Library names of the 901 "good" libraries used in Genome Paper
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/samp901_lib_names.txt`

## Generate CDS VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/all_samps

sbatch generate_Chr01K_01N_CDS_allsamp_vcf.sh
sbatch generate_Chr02_05_CDS_allsamp_vcf.sh
sbatch generate_Chr06_09_CDS_allsamp_vcf.sh
```
### Example of submit script
* From `generate_Chr01K_01N_CDS_allsamp_vcf.sh`
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

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/all_samps/

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/
SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/samp901_lib_names.txt
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/all_samps/

for CHROM in 01K 01N;
  do
  vcftools --gzvcf $DATA_DIR'Chr'$CHROM'.polyploid.CDS.vcf.gz' \
  --stdout --keep $SAMP_FILE --recode --recode-INFO-all | \
  gzip -c > \
  $OUT_DIR'Chr'$CHROM'.polyploid.CDS.allsamps.vcf.gz';
  done
``` 

## Generate VCF header file to be used with VCF subfiles
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/all_samps
head -5 Chr01K.polyploid.CDS.allsamps.vcf_00 | tail -1 > \
CDS.allsamps.vcf.header.txt
```

