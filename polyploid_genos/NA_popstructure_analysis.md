# Info about using NA's to analyze population structure

## Overview
* Convert genotypes to 1=genotype; 0=NA
* Generate distance matrices based on 1's and 0's
* Look at PCoA and NJ trees
* Look at --max-missing 0 (no filter), 0.2 (80% missing data), 0.5 (50% 
missing data) and 0.8 (20% missing data)

## Generate NA Distance Matrices
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/all_samps/

CHROM=01K
VCF_MEAT=.polyploid.CDS.allsamps.vcf_
SUB_NUM=00

VCF_IN=$DATA_DIR'Chr'$CHROM$VCF_MEAT$SUB_NUM

HEADER_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/full_vcf_header.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/na_popstructure/all_samps/no_filtering/

Rscript DIR/generate_NA_dist_mat.r $VCF_IN $HEADER_FILE $OUT_DIR


```


