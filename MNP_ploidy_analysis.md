# Ploidy Analysis Based on MNP Counts

## Overview
* Choose criteria for calling ploidy using MNP counts
* Select outliers based on previous ploidy desinations in metadata

## R script
* on HA
  * `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/MNP_ploidy_calls_and_outliers.r`

## Table of MNP results using depth of 3
* on NERSC
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_03/pos_03_MNP_results_with_meta_info.txt`

## Figures
* in this directory on NERSC with .pdf suffix
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_03/`

## Calculate stats from VCF for adjusting MNP results
* Turns out can't use these results
```
module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_ploidy/ld_prune_vcf
vcftools --vcf ldprune0.6_unimpute_4d_anc.vcf --out ldprune --het
```

