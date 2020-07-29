# Steps for analysing r^2 in simulated 8X samples

## Generate sample lists
* lists have to be formatted correctly for the R-scripts
* made using `~/sim8X_r2_pop_table.r`
### Files
* All sim8X populations
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_libs_for_r2.txt`
    * 16 genotypes per populationl; 576 total
    * 8 auto-8X populations
    * 28 allo-8X populations
* Large (25+) sim8X populations
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large_sim8X_libs_for_r2.txt`
    * 30 genotypes per population; 630 total 
    * 6 auto-8X populations
    * 15 allo-8X populations

## Calculate r^2 for sim8X sub-vcf files
* Important notes
  * max distance between SNPs = 10kbp
  * min MAF = 0.1
  * r^2 calculated separately within each of the simulated populations
* sim8X VCFs are formatted differently than the other VCFs
  * they only show Alternate-allele counts, so aren't formatted like genotypes
in "regular" VCFs
* Rscript used for this
  * `~/calc_snp_r2_sim8X.r`
* Directory with results
  * All sim8X pops
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_r2_results`
  * Large sim8X pops
    * ``
### Submit jobs
#### All sim8X pops
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_r2_results

sbatch calc_r2_Chr01_Chr03.sh
sbatch calc_r2_Chr04_Chr06.sh
sbatch calc_r2_Chr07_Chr09.sh
```

### Consolidate Results
* Rscript
  * `~/consolidate_r2_results.r`
#### Submit job
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_r2_results

sbatch consolidate_sim8X_r2_results.sh
```
#### Submit file
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
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
source activate /global/homes/g/grabowsp/.conda/envs/R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_r2_results/

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_r2_results/

FILE_SUF=sim8X_subgroups

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/\
r2_analysis/consolidate_r2_results.r \
$DATA_DIR $FILE_SUF

```

* CONTINUE FROM HERE

