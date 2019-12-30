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
### Analyse 20+ depth genic MNP numbers
* `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/analyze_genic_MNPs_20depth.r`


## CDS MNPs with 20+ depth per allele
### Tally CDS MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20
cp ./MNP_20depth_genic_counts.r MNP_20depth_CDS_counts.r
# adjust in vim

# copy and adjust submit script
cp ./count_genic_MNPs_20depth.sh count_CDS_MNPs_20depth.sh
# adjust in vim

sbatch count_CDS_MNPs_20depth.sh
```
### Analyse 20+ depth CDS MNPs
* on HA
* `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/analyze_CDS_MNPs_20depth.r`
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_genic_MNPs_20depth.r analyze_CDS_MNPs_20depth.r
```

## Analyse the Genic MNPs with 10+ depth per allele
### Tally genic MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/MNP_20depth_genic_counts.r ./MNP_10depth_genic_counts.r
# adjust in vim

# adjust submit script
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/count_genic_MNPs_20depth.sh ./count_genic_MNPs_10depth.sh
# adjust in vim

cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10
sbatch count_genic_MNPs_10depth.sh
```
### Analyse genic MNPs with 10+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_genic_MNPs_20depth.r analyze_genic_MNPs_10depth.r
```

## CDS MNPs with 10+ depth per allele
### Figures
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/CDS_MNP_10depth*pdf`
### Tally CDS MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10
cp ../pos_20/MNP_20depth_CDS_counts.r MNP_10depth_CDS_counts.r
# adjust in vim

# copy and adjust submit script
cp ../pos_20/count_CDS_MNPs_20depth.sh count_CDS_MNPs_10depth.sh
# adjust in vim

sbatch count_CDS_MNPs_10depth.sh
```
### Analyse CDS MNPs with 10+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_genic_MNPs_10depth.r analyze_CDS_MNPs_10depth.r
```

## Genic MNPs with 15+ depth per allele
### Tally genic MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_15
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/MNP_20depth_genic_counts.r ./MNP_15depth_genic_counts.r
# adjust in vim

# adjust submit script
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10/count_genic_MNPs_10depth.sh ./count_genic_MNPs_15depth.sh
# adjust in vim

cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_15
sbatch count_genic_MNPs_15depth.sh
```
### Analyse genic MNPs with 15+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_genic_MNPs_10depth.r analyze_genic_MNPs_15depth.r
```

## CDS MNPs with 15+ depth per allele
### Figures
### Tally CDS MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_15
cp ../pos_20/MNP_20depth_CDS_counts.r MNP_15depth_CDS_counts.r
# adjust in vim

# copy and adjust submit script
cp ../pos_10/count_CDS_MNPs_10depth.sh count_CDS_MNPs_15depth.sh
# adjust in vim

sbatch count_CDS_MNPs_15depth.sh
```
### Analyse CDS MNPs with 15+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_CDS_MNPs_10depth.r analyze_CDS_MNPs_15depth.r
```

## Genic MNPs with 25+ depth per allele
### Tally genic MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_25
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/MNP_20depth_genic_counts.r ./MNP_25depth_genic_counts.r
# adjust in vim

# adjust submit script
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10/count_genic_MNPs_10depth.sh ./count_genic_MNPs_25depth.sh
# adjust in vim

cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_25
sbatch count_genic_MNPs_25depth.sh
```
### Analyse genic MNPs with 25+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_genic_MNPs_10depth.r analyze_genic_MNPs_25depth.r
```

## CDS MNPs with 25+ depth per allele
### Figures
### Tally CDS MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_25
cp ../pos_20/MNP_20depth_CDS_counts.r MNP_25depth_CDS_counts.r
# adjust in vim

# copy and adjust submit script
cp ../pos_10/count_CDS_MNPs_10depth.sh count_CDS_MNPs_25depth.sh
# adjust in vim

sbatch count_CDS_MNPs_25depth.sh
```
### Analyse CDS MNPs with 25+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_CDS_MNPs_10depth.r analyze_CDS_MNPs_25depth.r
```

## Genic MNPs with 5+ depth per allele
### Tally genic MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_05
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/MNP_20depth_genic_counts.r ./MNP_05depth_genic_counts.r
# adjust in vim

# adjust submit script
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10/count_genic_MNPs_10depth.sh ./count_genic_MNPs_05depth.sh
# adjust in vim

cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_05
sbatch count_genic_MNPs_05depth.sh
```
### Analyse genic MNPs with 5+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_genic_MNPs_10depth.r analyze_genic_MNPs_05depth.r
```

## CDS MNPs with 5+ depth per allele
### Figures
### Tally CDS MNPs
#### Scripts
```
# copy and adjust R script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_05
cp ../pos_20/MNP_20depth_CDS_counts.r MNP_05depth_CDS_counts.r
# adjust in vim

# copy and adjust submit script
cp ../pos_10/count_CDS_MNPs_10depth.sh count_CDS_MNPs_05depth.sh
# adjust in vim

sbatch count_CDS_MNPs_05depth.sh
```
### Analyse CDS MNPs with 5+ depth
* on HA
```
cd /home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/
cp analyze_CDS_MNPs_10depth.r analyze_CDS_MNPs_05depth.r
# NEED TO RUN
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

