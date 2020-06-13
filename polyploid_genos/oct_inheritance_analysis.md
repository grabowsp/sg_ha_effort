# Analysis of inheritance patterns in Octoploid (8X) switchgrass

## Overview
### Analysis ideas
* measure LD via r^2 in different groups
* Look at levels of heterozygosity within individuals
  * compare levels of 1:3, 2:2, 3:1 genotypes
* Look for fixed HET SNPs and/or high numbers of HET loci across population

## LD analysis
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/Chr01K.polyploid.CDS.up_oct_1.vcf_00

HEAD_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/CDS.up_oct_1.vcf.header.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/

MAX_DIST=10000

MAF_CUT=0.05

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
$IN_FILE $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT
 
```

 
