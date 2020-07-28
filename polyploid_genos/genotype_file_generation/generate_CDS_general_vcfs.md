# Steps used for generating CDS VCFs for downstream analysis

## Overview
* Starting VCFs are tetrasomic genotypes generated for all samples
* Use VCFtools to extract CDS SNPs using a .bed file of CDS positions

## Location of Chromosome CDS VCFs
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
* Where the output files are stored

## Location of Input Files
### Raw VCF's (genotype calls, no read counts)
* Sujan's directory with VCF's
  * `/global/projectb/scratch/sujan/Pvirgatum_polyploid/`
* My directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs`
### Files to help process VCFs
* BED Files
  * CDS BED
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/sg_v5_CDS.bed`
  * Genes BED
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/sg_v5_genes.bed`

## Generate CDS VCFs for each chromosome
### Run jobs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs

sbatch generate_Chr01K_N_CDS_vcfs.sh

# check that works
#gunzip -kc Chr01K.polyploid.CDS.vcf.gz | head

sbatch generate_Chr02to05_CDS_vcfs.sh
sbatch generate_Chr06to09_CDS_vcfs.sh
```
* Moved to `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
### Example of submit script
* From `generate_Chr01K_N_CDS_vcfs.sh`
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

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs

DATA_DIR=/global/projectb/scratch/sujan/Pvirgatum_polyploid/
BED_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/sg_v5_CDS.bed
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/

for CHROM in 01K 01N;
  do
  bzip2 -dkc $DATA_DIR'Chr'$CHROM'.polyploid.vcf.bz2' | \
  vcftools --vcf - --bed $BED_FILE --recode --recode-INFO-all -c | \
  gzip -c > $OUT_DIR'Chr'$CHROM'.polyploid.CDS.vcf.gz';
  done
```
