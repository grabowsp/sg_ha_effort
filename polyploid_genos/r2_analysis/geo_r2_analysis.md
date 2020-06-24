# Process for evaluating r2 in the "Natural" Samples

## Generate Sample and Population dataframe
* R script for generating tables
  * `~/geo_r2_pop_table.r`
* Table of all subpopulation libraries
  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt`
  * 352 Libraries
  * 10 subpopulations
* Table of 4X subpopulation libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt`
  * 306 Libraries
  * 8 subpopulations
* Table of 4X subpopulations with 25+ libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large4X_subpop_libs_for_r2.txt`
  * 271 Libraries
  * 6 subpopulations

## Run r^2 analysis
### Test
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

VCF_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00

HEADER_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SAMP_POP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt

VCF_TYPE=allele_count

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

FILE_PRE=geo_subgroup

MAX_DIST=10000

MAF=0.1

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/r2_analysis/calc_snp_r2_v2.r

```


