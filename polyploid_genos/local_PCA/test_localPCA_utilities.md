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
### Test for larger windows
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2
sbatch localPCA_window_test_Chr01K_larger_windows.sh

```

### Submit test for rest of chromosomes
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2
sbatch localPCA_window_test_2.sh

```
### Compile output
* R script
  * `~/sg_ha_effort/polyploid_genos/local_PCA/combine_pick_window_res.r`
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

RES_SUF=polyploid.CDS.expandv2.windowtests.txt

OUT_PRE=CDS.expandv2

R_SCRIPT=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/local_PCA/combine_pick_window_res.r

Rscript $R_SCRIPT $DATA_DIR $RES_SUF $OUT_PRE

```

### Test jackknife script

```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

for tf in localPCA_jack*test.sh;
do
sbatch $tf;
done


```
### Look at jackknife output
```

module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2

RES_SUF=jackknifetest.txt

OUT_PRE=Chr01K.CDS.expandv2

R_SCRIPT=/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/local_PCA/combine_jackknife_res.r

Rscript $R_SCRIPT $DATA_DIR $RES_SUF $OUT_PRE

```

### Look at results
```
# module load python/3.7-anaconda-2019.07
# source activate local_PCA

library(data.table)
library(lostruct)
library(ggplot2)

data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/expand_v2/'

data_file_short <- 'Chr01K.polyploid.CDS.expandv2.300_SNP.20000_bp.localPCA.rds'

data_file <- paste(data_dir, data_file_short, sep = '')


tmp_res <- readRDS(data_file)

pcdist <- pc_dist(tmp_res, npc = 2)
tot_cmd <- cmdscale(pcdist, k = 50)
dist_tot_var <- sum(apply(tot_cmd, 2, var))
dist_per_var <- (apply(tot_cmd, 2 , var)/dist_tot_var)*100

tot_cmd_df <- data.frame(tot_cmd, stringsAsFactors = F)
colnames(tot_cmd_df) <- paste('PCo_', seq(ncol(tot_cmd_df)), sep = '')

pcs_to_keep <- c(1)

tot_res_cols <- c()
for(ptk in pcs_to_keep){
  tmp_cols <- colnames(tmp_res)[grep(paste('PC_', ptk, '_', sep = ''), 
    colnames(tmp_res))]
  tot_res_cols <- c(tot_res_cols, tmp_cols)
}

cols_to_keep <- c('total', paste('lam_', pcs_to_keep, sep = ''), tot_res_cols)

tmp_res_sub <- tmp_res[, cols_to_keep]

pcdist_sub <- pc_dist(tmp_res_sub, npc = 1)
sub_cmd <- cmdscale(pcdist_sub, k = 50)
sub_tot_var <- sum(apply(sub_cmd, 2, var))
sub_per_var <- (apply(sub_cmd, 2, var)/sub_tot_var)*100

sub_cmd_df <- data.frame(sub_cmd, stringsAsFactors = F)
colnames(sub_cmd_df) <- paste('PCo_', seq(ncol(sub_cmd_df)), sep = '')

gg_c1 <- ggplot(tot_cmd_df, aes(x = PCo_1, y = PCo_2)) + geom_point() +
  ggtitle('PCoA of Chr01K expand_samps local PCA\n300 SNPs, 20kb windows\nDistance using PC 1 and PC 2')

gg_sub_c1 <- ggplot(sub_cmd_df, aes(x = PCo_1, y = PCo_2)) + geom_point() +
  ggtitle('PCoA of Chr01K expand_samps local PCA\n300 SNPs, 20kb windows\nDistance using PC 1 only')

cmd_pdf_out_file_1 <- paste(data_dir, 
  'localPCA_Chr01K_300SNPs20kb_PC1_PC2_test.pdf', sep = '')

pdf(cmd_pdf_out_file_1, width = 7, height = 7)
gg_c1
dev.off()

cmd_pdf_out_file_2 <- paste(data_dir, 
  'localPCA_Chr01K_300SNPs20kb_PC1_only_test.pdf', sep = '')

pdf(cmd_pdf_out_file_2, width = 7, height = 7)
gg_sub_c1
dev.off()


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


