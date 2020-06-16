# Analysis of inheritance patterns in Octoploid (8X) switchgrass

## Overview
### Analysis ideas
* measure LD via r^2 in different groups
* Look at levels of heterozygosity within individuals
  * compare levels of 1:3, 2:2, 3:1 genotypes
* Look for fixed HET SNPs and/or high numbers of HET loci across population

## LD analysis
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups

IN_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/Chr01K.polyploid.CDS.up_oct_1.vcf_00

HEAD_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/CDS.up_oct_1.vcf.header.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/

MAX_DIST=10000

MAF_CUT=0.05

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
$IN_FILE $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT
 
```
### Chr01K analysis
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups

TOP_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/

HEAD_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/CDS.up_oct_1.vcf.header.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/

MAX_DIST=10000

MAF_CUT=0.05

for OCT_1_F in $TOP_DIR'Chr01K.polyploid.CDS.up_oct_1.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $OCT_1_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

for OCT_2_F in $TOP_DIR'Chr01K.polyploid.CDS.up_oct_2.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $OCT_2_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

for TET_1_F in $TOP_DIR'Chr01K.polyploid.CDS.up_tet_1.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $TET_1_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

for TET_2_F in $TOP_DIR'Chr01K.polyploid.CDS.up_tet_2.vcf_'*;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
calc_snp_r2.r \
  $TET_2_F $HEAD_FILE $OUT_DIR $MAX_DIST $MAF_CUT;
  done;

```
#### Consolidate results
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups

# in R

file_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/'

group_name_vec <- c('up_tet_1', 'up_tet_2', 'up_oct_1', 'up_oct_2')

for(j in seq(length(group_name_vec))){
  res_file_pre <- paste('Chr01K.polyploid.CDS.', group_name_vec[j], '_', 
    sep = '')
  r2_file_vec <- system(paste('ls ', file_dir, res_file_pre, '*', 'r2.rds', 
    sep = ''), intern = T)
  r2_tot <- readRDS(r2_file_vec[1])
  for(i in c(2:length(r2_file_vec))){
    tmp_res <- readRDS(r2_file_vec[i])
    r2_tot <- rbind(r2_tot, tmp_res)
  }
  out_file <- paste(file_dir, res_file_pre, 'tot_r2.rds', sep = '')
  saveRDS(r2_tot, out_file)
}
```
#### Generate Plot
* `~/r2_analysis_code.r`
* Plot locations
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/Chr01K.polyploid.CDS.upgroups_hiLD.pdf`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/up_subgroups/r2_results/Chr01K.polyploid.CDS.upgroups_r2.pdf`
