# Explanation of Steps for r^2 analysis

## Overview
* The steps going from VCF to r^2 results and figures

## VCF Inputs
* Pipeline is designed to be run on Chromosome-level VCF files divided into 
100k line subfiles, with all sub-vcf files in the same directory
  * the sub-vcf files end with '.vcf_' and a number
* Example of sub-files in this Cori directory:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps`
* Documentation of genotype file 'manipulation', including using 'split' to 
make the sub-vcf files
  * `~/../genotype_file_manipulations.md`

## Choose Samples for Sub-groups
* Want to select genetically-similar samples for examining r^2
* Visually select samples based on position of PCA graph
  * It seems easiest to select groups based on results from upland-only or
lowland-only PCA
* Example of process:
  * `~/../adegenet_analysis/select_geosamp_subgroups.md`

## Generate Sub-group file used for r^2 analysis
* Need to have properly formatted table to process sub-groups separately
* Example of code for generating table
  * `~/geo_r2_pop_table.r`
  * Generated 3 tables:
    * 1 with all subgroups for r^2 analysis
    * One only with 4X samples, for generating simulated 8X genotypes for
as many combindations of subgroups as possible
    * One with large 4X sub-groups for making larger populations of simulated
8X genotypes

## Calculate r^2 for sub-vcf files
* R script used for this
  * `~/calc_snp_r2_v2.r`
* Example on Cori of submitting jobs using this script
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/calc_r2_Chr04_Chr06.sh`
* Example of running code on inital 'geo_samps'
  * `~/geo_r2_analysis.md`

## Consolidate r^2 across sub-vcf results
* Generates Chromosome-level and genome-wide results
* R script used for his:
  * `~/consolidate_r2_results.r`
* Example of running code on inital 'geo_samps'
  * `~/geo_r2_analysis.md`

## Calculate window-based values
* Script for getting mean r^2 in distance-based windows
  * `~/gen_mean_r2_windows.r`
* Script for getting percentage of r^2 in window with high r^2
  * `~/gen_hi_r2_windows.r`
* Should do this separatelyfor chromosome-level and genome-wide consolidated
results
* Example of running code on inital 'geo_samps'
  * `~/geo_r2_analysis.md`
* Example of submitting job using script
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/gen_chrom_mean_r2_windows.sh`


