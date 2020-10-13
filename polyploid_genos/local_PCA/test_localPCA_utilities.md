# Testing localPCA utilities

## Choosing bp windows for SNP windows
* Test script for choosing best bp-windows for desired number of SNPs in a window
* R script = `~/sg_ha_effort/polyploid_genos/local_PCA/pick_window_v2.r`
* test shell script
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/localPCA_window_test_1.sh`
### Submit test of just 1 chromosome
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2
sbatch localPCA_window_test_1.sh
```
### Example of submit script
```
#!/bin/bash

#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=16G
#SBATCH --qos=genepool_shared
#SBATCH -t 24:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/local_PCA

REPO_DIR=/global/homes/g/grabowsp/tools/
DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/
CHR_NAME=Chr01K
HEAD_SHORT=CDS.expandv2.vcf.header.txt
INPUT_FILES=$REPO_DIR'sg_ha_effort/polyploid_genos/local_PCA/standard_input_files.r'

R_SCRIPT_NAME=$REPO_DIR'sg_ha_effort/polyploid_genos/local_PCA/pick_window_v2.r '

cd $DATA_DIR

$R_SCRIPT_NAME $REPO_DIR $DATA_DIR $CHR_NAME $HEAD_SHORT $INPUT_FILES



```
