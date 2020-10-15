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
### Submit test for rest of chromosomes
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2
sbatch localPCA_window_test_2.sh

```
### Compile output
* This example is for old version of output that doesn't include Chr
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/'

sys_com <- paste('ls ', data_dir, '*windowtests.txt', sep = '')

in_files <- system(sys_com, intern = T)

for(infi in in_files){
  res <- read.table(infi, header = T, sep = '\t', stringsAsFactors = F)
  chr_name <- unlist(strsplit(sub(data_dir, '', infi), split = '.', 
    fixed = T))[1]
  res$chr <- chr_name
  if(infi == in_files[1]){
    tot_res <- res
  } else{
    tot_res <- rbind(tot_res, res)
  }
}

tot_med_windows <- tapply(tot_res$median_best_window, tot_res$snp_window, 
  median)

tot_n_windows <- tapply(tot_res$n_windows, tot_res$snp_window, sum)

tot_per_windows <- tapply(tot_res$median_percent_good_window, 
  tot_res$snp_window, mean)

tot_wind_df <- data.frame(SNP_window = names(tot_med_windows), 
  best_bp_window = tot_med_windows, stringsAsFactors = F)

tot_wind_df$n_good_windows <- tot_n_windows[rownames(tot_wind_df)]

tot_wind_df$per_good_windows <- tot_per_windows[rownames(tot_wind_df)]

tot_wind_df <- tot_wind_df[ order(as.numeric(
  sub('_SNPs', '', tot_wind_df$SNP_window))), ]

out_file <- paste(data_dir, 'tot_best_bp_window_test.txt', sep = '')


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

Rscript R_SCRIPT_NAME=$REPO_DIR'sg_ha_effort/polyploid_genos/local_PCA/pick_window_v2.r '

cd $DATA_DIR

$R_SCRIPT_NAME $REPO_DIR $DATA_DIR $CHR_NAME $HEAD_SHORT $INPUT_FILES
```

## Get info about window choice
* R script = `~/sg_ha_effort/polyploid_genos/local_PCA/get_window_info.r`
* test shell script
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/localPCA_windowinfo_test_1.sh`
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/
sbatch localPCA_windowinfo_test_1.sh
```


