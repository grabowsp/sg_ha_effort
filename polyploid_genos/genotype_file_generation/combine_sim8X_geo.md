# Steps for combining simulated 8X and geosamps VCFs

## Overview
* To combine tetrasomic-style VCFs, need to convert them
* The main issue was needed additional ALT alleles so that genotypes with 
numbers above 1 are accepted
  * I decided to also convert the genotypes to a code using 0, 1, and 2 
so that would only need to add 1 additional character, N, to the ALT column
  
## Steps
* Convert subfiles to "coded" disomic genotypes
* Concatenate subfiles
* Sort, bgzip and tabix subfiles
* Merge sim8X and geosamps
* Generate new subfiles
* Convert new subfiles to tetrasomic genotypes

## Coded disomic genotypes
* TET genotype : DIP code
* 4/0 : 2/0
* 3/1 : 2/1
* 2/2 : 2/2
* 1/3 : 1/2
* 0/4 : 0/2




## Convert subfile to coded disomic genotypes 
* Steps:
  * convert tet genotypes to dipcode genotypes for subfiles
    * includes adding ALT alleles so that bcftools is happy
  * concatenate the subfiles from a chromosome
  * sort, bgzip, and tabix single chromosome
### `geo_samps` files 
* all steps included in single script
  * Ran Chr01K as step-by-step test
  * Ran Chr01N as sh script test
  * Ran remaining chromosomes in single script
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

sbatch tet_to_dipcode_Chr01N.sh

sbatch tet_to_dipcode_Chr02_to_Chr09_KN.sh
```
### Simulated 8X files
* Additional first step is to convert AltDosage genotypes to "standard" format
* Tried to as many steps in a single script as possible
  * Chr01K was run as step-by-step test
  * Chr01N was run as test of shell script
  * Rest of chromsomsomes run in single scripts
#### Convert AltDosage to standard file
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

sbatch convert_AltDosage_standard_Chr01N.sh
sbatch convert_AltDosage_standard_Chr02_Chr09_KN.sh
```
#### Convert to coded disomic file
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs

sbatch tet_to_dipcode_Chr01N.sh
sbatch tet_to_dipcode_Chr02_to_Chr09_KN.sh
```

## Merge and convert geo_samps and sim8X
* Steps:
  * merge dipcoded chromosome VCFs
  * split merged VCF into 100k line subfiles
  * convert dipcode subfiles to tetrasomic genotypes
  * delete dipcoded merged VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo

sbatch merge_convert_geo_sim8X_Chr01N.sh

sbatch  merge_convert_geo_sim8X_Chr02_09.sh
```

## Delete remaining dipcoded VCFs












### Test
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00

CONVERT_TYPE=1

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $CONVERT_TYPE $OUT_DIR

```
### Rest of Chr01K
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

bash

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

CONVERT_TYPE=1

for VN in {01..10};
do
echo $VN;
IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_$VN
Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $CONVERT_TYPE $OUT_DIR;
done

```
### combine converted subfiles
* includes VCF headers generated previously
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

FORMAT_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/vcf_format_header.txt

SAMP_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

OUT_FILE=Chr01K.polyploid.CDS.geosamps.vcf_dipcode.vcf

cat $FORMAT_HEAD $SAMP_HEAD \
Chr01K.polyploid.CDS.geosamps.vcf_dipcode_* > $OUT_FILE

```
### Sort, bgzip, and tabix the new vcf
```
module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

bcftools sort Chr01K.polyploid.CDS.geosamps.vcf_dipcode.vcf -Ov -o Chr01K_dipcode.vcf_sort

bgzip Chr01K_dipcode.vcf_sort

tabix -p vcf Chr01K_dipcode.vcf_sort.gz
```



### Test with sim8X vcf
### Change sim8X Chr01K vcfs to standard tet format
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

#IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_AltDosage.vcf_01

#Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_sim8X_VCF_to_standard.r $IN_VCF

for VN in {00..10};
do
echo $VN;
IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_AltDosage.vcf_$VN
Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_sim8X_VCF_to_standard.r $IN_VCF;
done

```
### Convert sim8X to coded
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs

CONVERT_TYPE=1
#IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_standard.vcf_00

#Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $OUT_DIR

for TN in {00..10};
do
echo $TN
IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_standard.vcf_$TN
Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $CONVERT_TYPE $OUT_DIR;
done

```
### Concatenate sim8X Chr01K subfiles

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs

FORMAT_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/vcf_format_header.txt

SIM8X_SAMP_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/CDS.sim8X.vcf.header.txt

OUT_FILE=Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_standard.vcf_dipcode.vcf

cat $FORMAT_HEAD $SIM8X_SAMP_HEAD Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_standard.vcf_dipcode_* > $OUT_FILE

### Sort, bgzip, and tabix the new vcf
```
module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs

bcftools sort Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_standard.vcf_dipcode.vcf -Ov -o Chr01K_sim8X_dipcode.vcf_sort

bgzip Chr01K_sim8X_dipcode.vcf_sort

tabix -p vcf Chr01K_sim8X_dipcode.vcf_sort.gz
```

### Merge sim8X and geo
```
module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo

SIM_TO_MERGE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs/Chr01K_sim8X_dipcode.vcf_sort.gz

GEO_TO_MERGE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs/Chr01K_dipcode.vcf_sort.gz

OUT_FILE=Chr01K_geoSim8Xmerge_dipcode.vcf.gz

bcftools merge $GEO_TO_MERGE $SIM_TO_MERGE -Oz -o $OUT_FILE

# save the header
bcftools merge $GEO_TO_MERGE $SIM_TO_MERGE --print-header > \
geoSim8Xmerg_header.txt

tail -n 1 geoSim8Xmerg_header.txt > geoSim8Xmerg_samps.txt


```

#### Make subfiles
```
module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo

bgzip -cd Chr01K_geoSim8Xmerge_dipcode.vcf.gz \
| split -l 100000 -d - Chr01K_geoSim8Xmerge_dipcode.vcf_
```

#### Convert files to tetrasomic genotypes
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

bash

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo

CONVERT_TYPE=2

for VN in {00..10};
do
echo $VN;
IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8Xgeo_combo/Chr01K_geoSim8Xmerge_dipcode.vcf_$VN
Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $CONVERT_TYPE $OUT_DIR;
done

```




