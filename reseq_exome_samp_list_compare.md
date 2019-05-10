# Compare the sample lists from the resequencing and exome-capture analyses

## Sample info file locations
### Resequencing Samples
* Metadata directory:
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/`
  * last updated 4/18/2019
* Metadata files:
  * Metadata with Sujan's ploidy info: `/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/Pvirgatum_921samples_ploidy_FINAL.txt`
  * Original sample info from Sujan: `/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/sample_list_1005_to_Paul.txt`
  * Metadata from Jason in the Juenger lab (4/29/2019): `/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/PVDIV_Master_Metadata_File_5-1-2019.txt`
### Exome-capture
* Metadata directory: 
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/exome/metadata/`
* Metadata files:
  * Sample: `/home/t4c1/WORK/grabowsk/data/switchgrass/exome/metadata/combo_samp_metadata_v1.9.txt`
  * Population: `/home/t4c1/WORK/grabowsk/data/switchgrass/exome/metadata/combo_pop_metadata_v1.4.txt`
  * Last updated 4/18/2019
### Reseq metadata with exome-capture match info
* `/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/sg_reseq_metadata_with_exome_info.txt`
* Explanation of columns with info about exome-capture matches:
  * `exome_pop`: the name(s) of exome-capture population(s) that match the \
name and/or geographic info of the sample. 
  * `exome_ecotype_info`: A tally of the ecotype designation of all the \
samples in the exome-capture date that are part of the population(s) shown \
in `exome_pop`. 
    * For example: 'L.12.__U.0.__A.1' indicates 13 samples belong to that \
population: 12 Lowland, 0 Upland, 1 Admixed.  
  * `exome_ecotype`: Ecotype designation of population in `exome_pop`. This \
column will indicate a unanimous designation. 'MIXED' indicates that some \
samples from `exome_pop` have different ecotype designations. 
  * `exome_gp_info`: A tally of genepool designations for all the \
exome-capture samples from `exome_pop`. Similar format as \
`exome_ecotype_info`. 
  * `exome_gp`: Genepool designation for population in `exome_pop`. Similar \
format to `exome_ecotype`.
  * `exome_ploidy_info`: A tally of the ploidy designation for all the \
exome-capture samples from `exome_pop`. Similar format as `exome_ecotype_info`.
  * `exome_ploidy`: Ploidy designation for accession in 'exome_pop'. Similar \
format to `exome_ecotype`.

## Comparison
### Overview
* Looked for matches in reseq sample names and Jason's (Juenger lab) metadata \
to population names, accession names, and aliases in the exome-capture metadata
* Looked for geographic matches (within a certain lat/long rounding error)
  * reseq sample can match multiple exome-catpure populations if the \
exome-capture populations are so close that they can't be separated based \
on lat/long
* Prioritize name matches, then add geographic matches
### R script
* `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/samp_list_comparisons.r`
  * includes portion at the end getting info about ecotype, genepool, and \
ploidy for the matches from the exome-capture data




