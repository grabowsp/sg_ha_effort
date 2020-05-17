# Information about Files Used/Needed for Analysis

## Switchgrass v5 files
### v5 reference genome sequence
* Official NERSC location
  * `/global/dna/projectdirs/plant/assembly_releases/Panicum_virgatum/20180630`
My NERSC location
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/ref_stuff/assembly`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/ref_stuff/assembly/Panicum_virgatum_var_AP13/sequences`
### v5 annotation
* Official NERSC location
  * `/global/dna/projectdirs/plant/annotation/www/Pvirgatum`
* My NERSC location
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/ref_stuff/annotation`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/ref_stuff/annotation/Pvirgatum_516_v5.1`

## Switchgrass Resequencing Metadata
* NERSC Directory
  * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata`
### Sample Metadata
* Metadata as of 9/3/2019 - Full
  * `PVDIV_Master_Metadata_File_9-3-2019.xlsx`
* Metadata as of 9/3/2019 - Edited for loading into R
  * `Reseq_Metadata_Sept_2019_Edited_for_R.tsv`
  * deleted the `LIB_NOTES` column because of missing values
  * use `read.table('Reseq_Metadata_Sept_2019_Edited_for_R.tsv', sep = '\t', \
header = T, stringsAsFactors = F, quote = "", comment.char = '$')` \
to load in R
* Metadata from John - May 2020
  * `genotype.metadata.May2020.rds`
### Switchgrass Cultivars in Reseq Data
* Table of cultivars I found in the reseq data
  * Original file with lots of typos:
    * `reseq_cultivars.tsv`
  * Updated file
    * `switchgrass_reseq_cultivars_corrected.tsv`

## Genotype VCFs
### LD-prunded VCF
* `/global/cscratch1/sd/grabowsp/sg_ploidy/ld_prune_vcf/ldprune0.6_unimpute_4d_anc.vcf`
  * This turned out to only have the punitive 4X (not totally correct) samples

## Ploidy info
* `/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v1.0.txt`

 
