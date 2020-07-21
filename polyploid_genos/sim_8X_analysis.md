# (Preliminary) Info about analysis of Simulated 8X samples

## Generate combinations of samples for making simulated 8X samples
### Output files
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
#### all 4X sub-populations
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

SUBPOP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt

OUT_PRE=geo_samp

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/make_sim_8X_combos.r $SUBPOP_FILE $OUT_PRE
```
#### All large (25+ libraries) 4X sub-populations
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

SUBPOP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large4X_subpop_libs_for_r2.txt

OUT_PRE=large_geo_samp

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/make_sim_8X_combos.r $SUBPOP_FILE $OUT_PRE
```

## Generate Simulated 8X VCFs
### Test
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

VCF_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00

VCF_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SIM_COMBO_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_samp_sim8X_lib_combos.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

OUT_PRE=geo_samp

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/make_sim_8X_genotypes.r $VCF_IN $VCF_HEAD $SIM_COMBO_IN $OUT_DIR $OUT_PRE

```
### Combos using all 4X sub-populations
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

sbatch make_geo_sim8X_geno_vcfs.sh
```
### Combos using 4X sub-pops with 25+ libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/largeSim8X_vcfs

sbatch make_geo_largeSim8X_geno_vcfs.sh
```

# NEXT STEPS
* r2 analysis
* combine vcfs with standard VCFs
* PCA including standard VCFs



