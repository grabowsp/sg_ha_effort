# Steps for PCA analysis of the 'expanded geographic samples' that include
#   8X cultivars

## Location of Files
### Genlight genotypes
* All samples
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/combo.sub.polyploid.CDS.expandgeosamps.genlight.rds`
  * 785 samples
  * 552,225 SNPs
  * SNPs in CDS
* Upland samples
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/combo.sub.polyploid.CDS.upexpand.genlight.rds`
    * 287 samples    
    * 225,937 SNPs
* Lowland samples
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/lowexpand_samps/combo.sub.polyploid.CDS.lowexpand.genlight.rds`
    * 486 samples
    * 360,952 SNPs

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

## Run PCA on `upexpand` samples
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

sbatch upexpand_PCA.sh
```

## Generate `upexpand` PCA figures
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps/combo.sub.polyploid.CDS.upexpand.genlight.PCAresults.figs.pdf`
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/upexpand_samps

DATA_FILE=combo.sub.polyploid.CDS.upexpand.genlight.PCAresults.rds

PLOT_PRE=Upland_expanded_samples

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE
```

## Run PCA on `lowexpand` samples
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/lowexpand_samps

sbatch lowexpand_PCA.sh
```

## Generate `lowexpand` PCA figures
* `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/lowexpand_samps/combo.sub.polyploid.CDS.lowexpand.genlight.PCAresults.figs.pdf`
```
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/lowexpand_samps

DATA_FILE=combo.sub.polyploid.CDS.lowexpand.genlight.PCAresults.rds

PLOT_PRE=Lowland_expanded_samples

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE
```

