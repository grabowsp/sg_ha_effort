# Steps for analysis of 4X/8X data sets using adegenet

## Overview
* Generate 'genlight' object
  * Should be done ahead of time
* Run PCA
* Plot PCA results

## Starting Files
* genlight object
  * subsampled set of SNPs randomly selected from across genome
    * ex: 500k SNPs
  * Contains ploidy information
  * Steps for generating genlight objects and subsampling SNPs in 'genotype_manipulations' file/directory
  * Examples:
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps/combo.sub.polyploid.CDS.expandgeosamps.genlight.rds`

## Run PCA on NERSC
* Can take a while to run, so is best to submit job to cluster
### Requirements
* `r_adegenet_env`
  * conda environment with adegenet and other necessary packages
  * full path:
    * `/global/homes/g/grabowsp/.conda/envs/r_adegenet_env`
* `adegenet_pca.r`
  * R script for running PCA
### Example submit script
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 72:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/r_adegenet_env

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

DATA_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/combo.sub.polyploid.CDS.geosamps.genlight.rds

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/adegenet_pca.r $DATA_FILE

```

## Plot PCA Results
### Overview
* Use script to automatically generate plots
* Plots are PC1 vs PC2 through PC5 in same total figure
* Total figure include 2 versions of each plot color coded by A) gene pool or
B) ploidy
### Required resources
* `r_adegenet_env` conda environment
* `pca_figs_adegenet.r` script
  * `/PATH/TO/REPOSITORY/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r`
### Required files
* PCA results file
  * should end in '.PCAresults.rds'
* Ploidy info file
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v3.0.txt`
  * This is hard-coded into the script
### Variables and colors used in plotting
* The script has columns from the ploidy info file hard-coded into the process
* The current gene pool colors are based on old results
* May want to make a v2 of the results with new gene pools
### Run script
* Is fast - can just run in interactive session
* Figure name is same as PCA results file but with 'figs.pdf' replacing '.rda'
#### Example
* For running on Cori
```
# load environment
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# go to directory with files
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_geo_samps

# The name of the PCA results file
DATA_FILE=combo.sub.polyploid.CDS.expandgeosamps.genlight.PCAresults.rds

# Name of the sample set - used for making titles for each plot
PLOT_PRE=Expanded_geographic_samples

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/pca_figs_adegenet.r $DATA_FILE $PLOT_PRE
```
### Analysis of Results
* Can use plot to set criteria for dividing samples for hierarchical analysis
  * Example: use PC1 to divide upland and lowland samples

