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
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

```


