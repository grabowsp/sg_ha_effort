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

