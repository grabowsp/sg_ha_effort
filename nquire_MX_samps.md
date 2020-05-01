# Running nQuire to determine ploidy in Mexican samples

## Generate soft link to the switchgrass bams
### Location of bams in Sujan's directory
* `/global/cscratch1/sd/sujan/Pvirg_mexican`
  * `/global/cscratch1/sd/sujan/Pvirg_mexican/*.bam`
### Make list of library IDs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/mx_bam_links

ls /global/cscratch1/sd/sujan/Pvirg_mexican/*.bam | cut -d '/' -f 7 | \
cut -d '.' -f 1 > Pvig_mex_libs.txt
```
### Make soft links
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/mx_bam_links

for LIB in `cat Pvig_mex_libs.txt`;
  do
  ln -s /global/cscratch1/sd/sujan/Pvirg_mexican/$LIB'.merged.sort.bam' \
  /global/cscratch1/sd/grabowsp/sg_ploidy/mx_bam_links/$LIB'.bam';
  done
```

## Generate Indexes for each bam
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/mx_bam_links

cp /global/cscratch1/sd/grabowsp/sg_ploidy/sg_bam_links/make_bam_index_00.sh \
make_bam_index_full.sh

sbatch make_bam_index_full.sh
```

## Run nQuire on all libraries
### Copy old submit script
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire

cp /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire/CDS_nquire_000.sh ./mx_CDS_nquire_000.sh
```
### Divide Libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/mx_bam_links

split -l 4 -d Pvig_mex_libs.txt Pv_mex_libs_sub4_
```
### Test with first 4 libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire

sbatch mx_CDS_nquire_000.sh
```
### Run on remaining libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire

for SUBSET in {01..23};
  do
  sed 's/sub4_00/'sub4_"$SUBSET"'/g' mx_CDS_nquire_000.sh > \
mx_CDS_nquire_0$SUBSET'.sh';
  done

for SUBSET in {01..23};
  do
  sbatch mx_CDS_nquire_0$SUBSET'.sh';
  done
```
## Combine results
* Combined result file on NERSC:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire/mx_nQuire_results_summary_total.txt`
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire

head -1 IYYS.nQuire_res_summary.txt | cat - *nQuire_res_summary.no_head.txt \
> mx_nQuire_results_summary_total.txt
```
## Get mapping and structure results from Sujan
* samtools mapping stats
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire/Mexican_samples_samtools_stats`
* Sujan's fastStructure results
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire/Mexican_Diversity_structure_combined.txt`

## Combine nQuire and Sujan results and make figures
* Directory with figures and data
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire/`
* Combined table
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/mx_CDS_nquire/mx_nquire_and_Sujan_tot_v1.0.txt`
* R Script
  * `~/sg_ha_effort/r_scripts/nquire_mx_ploidy_calls.r`
