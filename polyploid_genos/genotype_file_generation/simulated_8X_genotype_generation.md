# Steps for generating simulated 8X genotypes

## Overview
* Two sets of VCFs
  * One set using all 4X subpopulations
    * more simulated 8X populations, fewer simulated genotype per sim pop
  * One set Using 4X supopulations with 25+ samples
    * fewer simulated 8X pops, more simulated genotypes per sim pop
* Sample come from 4X populations based on the `geo_samp` subpopulations chosen for r^2
analysis
* Rules used for making sample combinations used for genotypes:
  * No library combinations are repeated in any of the population combinations
  * For within-population combinations, each library is included in at most
2 combinations
* Genotypes are tetrasomic dosage (0-4) genotypes based on genotypes of the two
4X samples

## Location of Simulated 8X VCFs
* Combinations using all 4X subpopulations
  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs`
* Combinations of large 4X subpopulations
  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/largeSim8X_vcfs`

## Make simulated 8X Genotypes
* Script used for making genotypes
  * `~/sg_ha_effort/polyploid_genos/make_sim_8X_genotypes.r`
* Process takes a long time, so should run via submit script
### Submit jobs
#### All 4X sub-populations
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

sbatch make_geo_sim8X_geno_vcfs.sh
```
#### 4X sub-pops with 25+ libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/largeSim8X_vcfs

sbatch make_geo_largeSim8X_geno_vcfs.sh
```
### Example of submit script
* From `make_geo_sim8X_geno_vcfs.sh`
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

IN_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

VCF_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SIM_COMBO_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/geo_samp_sim8X_lib_combos.txt

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

OUT_PRE=geo_samp

for VCF_IN in `ls $IN_DIR*vcf_*`;
  do
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/make_sim_8X_genotypes.r \
  $VCF_IN $VCF_HEAD $SIM_COMBO_IN $OUT_DIR $OUT_PRE;
  done;

```

## Convert Allele-count VCFs to standard-style VCFs
### Steps
* Generate sim8X vcf
  * make header that includes Chr info
  * sort
  * tabix
  * convert to BCF
* Generate geo vcf
  * adjust header to include Chr info
  * sort
  * tabix
  * convert to BCF

### Generate sim8X sample header
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

head -1 Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_AltDosage.vcf_00 > CDS.sim8X.vcf.header.txt
# add '#' to beginning of line
```
### Generate vcf formatting header
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/
head -4 Chr01K.polyploid.CDS.geosamps.vcf_00 > vcf_format_header.txt 
```
* add Chr info to vcf_format_header.txt
```
##contig=<ID=Chr01K,length=56546490>
##contig=<ID=Chr01N,length=66475219>
##contig=<ID=Chr02K,length=68340493>
##contig=<ID=Chr02N,length=69773881>
##contig=<ID=Chr03K,length=63841035>
##contig=<ID=Chr03N,length=69486008>
##contig=<ID=Chr04K,length=47708993>
##contig=<ID=Chr04N,length=50282576>
##contig=<ID=Chr05K,length=61896298>
##contig=<ID=Chr05N,length=71719776>
##contig=<ID=Chr06K,length=48106579>
##contig=<ID=Chr06N,length=53044726>
##contig=<ID=Chr07K,length=52492986>
##contig=<ID=Chr07N,length=50739414>
##contig=<ID=Chr08K,length=55599782>
##contig=<ID=Chr08N,length=50668275>
##contig=<ID=Chr09K,length=70504668>
##contig=<ID=Chr09N,length=82442687>
```
:
### Generate adjusted VCF files and merge
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs

module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

SIM8X_VCF_IN=Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_standard.vcf_00

FORMAT_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/vcf_format_header.txt

SIM8X_SAMP_HEAD=CDS.sim8X.vcf.header.txt

cat $FORMAT_HEAD $SIM8X_SAMP_HEAD $SIM8X_VCF_IN > tmp_vcf

bcftools sort tmp_vcf -Ov -o tmp_vcf_sort
bgzip tmp_vcf_sort 
tabix -p vcf tmp_vcf_sort.gz

# make temporary geosamps file
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

GEO_VCF_IN=Chr01K.polyploid.CDS.geosamps.vcf_00

FORMAT_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/vcf_format_header.txt

GEO_SAMP_HEAD=CDS.geosamps.vcf.header.txt

tail -n+5 $GEO_VCF_IN > tmp_vcf_1

cat $FORMAT_HEAD tmp_vcf_1 > tmp_geo_vcf

bcftools sort tmp_geo_vcf -Ov -o tmp_geo_vcf_sort
bgzip tmp_geo_vcf_sort
tabix -p vcf tmp_geo_vcf_sort.gz

# try merging

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

SIM_TO_MERGE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/tmp_vcf_sort.gz

GEO_TO_MERGE=tmp_geo_vcf_sort.gz

bcftools merge $GEO_TO_MERGE $SIM_TO_MERGE -Oz -o tmp_merged.vcf.gz
```
