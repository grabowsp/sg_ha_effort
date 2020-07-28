# Explanation of choosing sample combinations for simulated 8X genotypes

## Overview
* Select pairwise library combinations within- and between 4X populations 
  for combining into simulated 8X genotypes
* Two sets of combinations
  * Using all 4X subpopulations
    * more simulated 8X populations, fewer simulated genotype per sim pop
  * Using 4X supopulations with 25+ samples
    * fewer simulated 8X pops, more simulated genotypes per sim pop 
* Uses 4X populations based on the `geo_samp` subpopulations chosen for r^2 
analysis
* Rules used for making sample combinations:
  * No library combinations are repeated in any of the population combinations
  * For within-population combinations, each library is included in at most 
2 combinations

## Sample combination explanation and output files
* Combinations using all 8 4X subpopulations
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_samp_sim8X_lib_combos.txt`
  * 16 combinations per simulated sub-population
  * 8 auto-8X sub-populations
  * 28 allo-8X sub-populations
    * samples chosen from 2 different sub-populations
* Combinations using 6 4X subpopulations with 25+ libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large_geo_samp_sim8X_lib_combos.txt`
  * 30 combinations per simulated sub-population
  * 6 auto-8X sub-populations
  * 15 allo-8X sub-populations

## Info about libraries and subpops used for combinations
* Libraries/samples come from sub-populations chosed for r2 analysis
* List of libraries generated with this script:
    * `~/geo_r2_pop_table.r`
* Files with Libraries and Sub-pops chosen for r2 analysis and used for sim 8X combinations
  * All 4X sub-populations
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt`
  * 4X sub-populations with 25+ Libraries
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large4X_subpop_libs_for_r2.txt`


## Make sample combinations
* Script used for generating combinations
  * `~/sg_ha_effort/polyploid_genos/make_sim_8X_combos.r`
### all 4X sub-populations
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

SUBPOP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt

OUT_PRE=geo_samp

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/make_sim_8X_combos.r $SUBPOP_FILE $OUT_PRE
```
### All large (25+ libraries) 4X sub-populations
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SUBPOP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large4X_subpop_libs_for_r2.txt

OUT_PRE=large_geo_samp

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/make_sim_8X_combos.r $SUBPOP_FILE $OUT_PRE
```


