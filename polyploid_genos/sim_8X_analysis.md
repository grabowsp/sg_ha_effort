# (Preliminary) Info about analysis of Simulated 8X samples

## Generate combinations of samples for making simulated 8X samples
### Info about libraries and subpops used for combinations
* Files with Library and Sub-populations chosen for r2 analysis
  * All 4X sub-populations
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt`
  * 4X sub-populations with 25+ Libraries
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large4X_subpop_libs_for_r2.txt`
  * All (4X and 8X) sub-populations chosen for r2 analysis
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt`
  * List generated with this script:
    * `~/geo_r2_pop_table.r`
### Generate combinations
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

SUBPOP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt

OUT_PRE=geo_samp

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/make_sim_8X_combos.r $SUBPOP_FILE $OUT_PRE


```



