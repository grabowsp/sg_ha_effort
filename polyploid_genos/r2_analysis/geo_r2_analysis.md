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

## Run r^2 analysis
### Test
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
# Run on all subfiles
```

module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

# Variables for bash loop

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

LS_STRING=$DATA_DIR*vcf_*

# Variables for R script

HEADER_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SAMP_POP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt

VCF_TYPE=allele_count

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

MAX_DIST=10000

MAF=0.1

for VCF_IN in `ls $LS_STRING`;
  do
  echo $VCF_IN;
  FILE_PRE=`basename $VCF_IN`'_geo_subgroup';
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/r2_analysis/calc_snp_r2_v2.r \
  $VCF_IN $HEADER_IN $SAMP_POP_FILE $VCF_TYPE $OUT_DIR $FILE_PRE \
  $MAX_DIST $MAF;
  done;

```

## Errors

Error in if (sum(minor_af < maf_cut) > 0) { : 
  missing value where TRUE/FALSE needed
In addition: There were 50 or more warnings (use warnings() to see the first 50)
Execution halted


Files:
* /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_03
* /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01N.polyploid.CDS.geosamps.vcf_07
* /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr02K.polyploid.CDS.geosamps.vcf_01
* /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr02N.polyploid.CDS.geosamps.vcf_01
* /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr03K.polyploid.CDS.geosamps.vcf_04



