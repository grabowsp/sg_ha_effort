# Overview of analysis of MNP analysis for estimating ploidy

## Look at Seq Depth and generate chromosome-depth file
### Output files
* Matrix containing the total coverage for each sample by chromsome
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth/sample_seq_depth_mat.rds`
### Script used for analysis
* On HA
  * `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/process_MNP_seq_depth_files.r`

## Analyse the Genic MNPs with 20+ depth per allele
### Figures
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/genic_MNP_20depth_*.pdf`
### Tally genic MNPs
#### Output
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/genic_MNP_count_list_20depth.rds`
#### Script
* R script on Cori to Count up the number of MNPs that are in genic regions \
for each chromsome in each sample
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/MNP_20depth_genic_counts.r`
* Submit script
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/count_genic_MNPs_20depth.sh`
  * Submit job because takes too long to run interactively
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20
sbatch count_genic_MNPs_20depth.sh
```  
### Analyse genic MNP numbers
* `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/analyze_genic_MNPs_20depth.r`

## Analyse the Genic MNPs with 10+ depth per allele
### Tally genic MNPs

###
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/MNP_20depth_genic_counts.r ./MNP_10depth_genic_counts.r
# adjust in vim

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/count_genic_MNPs_20depth.sh ./count_genic_MNPs_10depth.sh
# adjust in vim

cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10
sbatch count_genic_MNPs_10depth.sh
```

## Analyze All MNPs with 20+ depth per allele
* One big peak and then long tail - hard to separate 4X from 8X
### Figures
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth20/mnp_hist_20depth_*.pdf`
### Script
* On HA
  * `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/analyze_MNPs_20depth.r`

## Analyze All MNPs with 10+ depth per allele
* Don't get any resolution to separate ploidy levels
### Figures
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10/mnp_hist_10depth_*.pdf`
### Script
* On HA
  * `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/analyze_MNPs_10depth.r`

