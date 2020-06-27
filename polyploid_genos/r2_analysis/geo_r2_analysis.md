# Process for evaluating r2 in the "Natural" Samples

## Generate Sample and Population dataframe
* R script for generating tables
  * `~/geo_r2_pop_table.r`
* Table of all subpopulation libraries
  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt`
  * 352 Libraries
  * 10 subpopulations
* Table of 4X subpopulation libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt`
  * 306 Libraries
  * 8 subpopulations
* Table of 4X subpopulations with 25+ libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large4X_subpop_libs_for_r2.txt`
  * 271 Libraries
  * 6 subpopulations

## Calculate Pairwise r^2 in CDS SNPs across genome
* Important Notes
  * max distance between SNPs = 10kbp
  * min MAF = 0.1
  * r^2 calculated separately within each of the 10 subgroups
* R script for generating sub-vcf objects
  * `~/calc_snp_r2_v2.r`
* Cori directory with results
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results`
  * There is a separate .rds file for each sub-vcf
* consolidated r^2 results separated by chromosomes
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/Chromosome.geo_subgroup.r2.rds`
* consolidated r^2 results combined across all chromosomes
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/combined.geo_subgroup.r2.rds` 
### Test script for 1 sub-vcf
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

VCF_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00

HEADER_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SAMP_POP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt

VCF_TYPE=allele_count

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

FILE_PRE=`basename $VCF_IN`'_geo_subgroup'

MAX_DIST=10000

MAF=0.1

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/r2_analysis/calc_snp_r2_v2.r $VCF_IN $HEADER_IN $SAMP_POP_FILE $VCF_TYPE $OUT_DIR \
$FILE_PRE $MAX_DIST $MAF

```
### Submit jobs for remaining sub-vcfs by chromosomes
* I tried using an interactive session to loop through the script for all
the sub-vcf results but was taking too long so stopped after Chr03K and Chr03N
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

sbatch calc_r2_Chr04_Chr06.sh
sbatch calc_r2_Chr07_Chr09.sh
```
### Consolidate Results
* Rscript for consolidating results from sub-vcf results
  * `~/consolidate_r2_results.r`
#### Can run interactively but takes a long time; next time I'd submit the job
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/

FILE_SUF=geo_subgroup

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
r2_analysis/consolidate_r2_results.r \
$DATA_DIR $FILE_SUF
```

## Generate r^ window value for windows

* Important Notes:
  * window size are 10bp, so 1-10bp, 11-20bp, etc distance between SNPs used
for pairwise comparison
  * For Hi-LD results, used r^2 > 0.9 as cutoff for high LD
  * generated window-based results for each of the 10 subgroups
* R script used for generating mean r^2 window results
  * `~/gen_mean_r2_windows.r`
* R script used for generating percentage of values in window with high r^2
  * `~/gen_hi_r2_windows.r`
### Submit jobs for generating window-based results
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

sbatch gen_chrom_mean_r2_windows.sh 
sbatch gen_chrom_hi_r2_windows.sh

sbatch gen_combo_mean_r2_windows.sh
sbatch gen_combo_hi_r2_windows.sh

```




