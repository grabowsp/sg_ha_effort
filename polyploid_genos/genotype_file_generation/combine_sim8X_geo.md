# Steps for combining simulated 8X and geosamps VCFs

## Proposed Steps
* Convert subfiles to "coded" disomic genotypes
* Concatenate subfiles
* Sort, bgzip and tabix subfiles
* Combine sim8X and geosamps
* Generate new subfiles
* Convert new subfiles to tetrasomic genotypes

## Coded disomic genotypes
* TET genotype : DIP code
* 4/0 : 2/0
* 3/1 : 2/1
* 2/2 : 2/2
* 1/3 : 1/2
* 0/4 : 0/2

## Test
* Convert geosamps Chr01K subfile(s) to coded disomic genotypes
* Concatenate subfiles
* sort, bgzip, and tabix
* Convert sim8X Chr01K subfiles to coded disomic genotypes
* concatenate subfiles
* sort, bgzip, and tabix
* merge the geosamps and sim8X Chr01K files
* generate merged subfiles
* convert merged subfiles back to tetrasomic genotypes

## Convert subfile to coded disomic genotypes
### Test
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $OUT_DIR

```
### Rest of Chr01K
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

bash

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs


for VN in {01..10};
do
echo $VN;
IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_$VN
Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $OUT_DIR;
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
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs
 
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/convert_vcfs

IN_VCF=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_standard.vcf_00

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/convert_tet_to_dipcoded_vcf.r $IN_VCF $OUT_DIR

```


