# Steps for PCA analysis of the 'expanded geographic samples' that include
#   8X cultivars

## Location of Files
### Genlight genotypes
* All samples
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/combo.sub.polyploid.CDS.expandgeosamps.genlight.rds`
  * 785 samples
  * 552,225 SNPs
  * SNPs in CDS

## Run PCA
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps

sbatch expandgeo_PCA.sh
```

## Generate PCA Figures
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps

DATA_FILE=combo.sub.polyploid.CDS.expandgeosamps.genlight.PCAresults.rds

PLOT_PRE=Expanded_geographic_samples

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE
```

## Set criteria for sub-divisions
* upland = PC1 < -20
* lowland = PC1 > -3

