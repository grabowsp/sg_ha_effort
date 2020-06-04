# Notes about PCA-based population structure analysis using adegenet

## Overview
* using `adegenet` package in R to run PCA for 4X/8X population structure
analysis
* As of June 1 2020, the PCA function in `adegenet` was not working properly
on Cori (there was an error with the eigen() function), so current workflow is:
  * Generate genlight object on Cori (documented in other .md file)
  * Transfer genlight object to HA
  * Run PCA on HA
  * Transfer PCA results to Cori
  * PCA-based figures and analysis on Cori

## PCA on all geo samples
* uses the `adegenet_pca.r` script
### Location of genlight object
* on Cori
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.rds`
* on HA
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.595K.polyploid.CDS.geosamps.genlight.rds`
### Run PCA on HA
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/

qsub combo595K_PCA_submit.sh
```
### Location of PCA results
* on Cori
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds`
* on HA
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds`
### Generate PCA figures
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

DATA_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds

PLOT_PRE=Geographic_samples_595K_SNPs

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE

```
### Divide Samples into genetic ecotype groups
* upland = PC1 < -25
* lowland = PC1 > 0

