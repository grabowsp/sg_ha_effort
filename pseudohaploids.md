# Steps for Pseudohaploid Analysis

## Overview


## Location of Files
### Sujan's Unimputed VCFs
  * `/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_V5_1035g_*.unfiltered_reheaded.sorted.vcf.gz`
    * Separate file for each chromosome
    * These filess are BIG
### Main Directory for Analysis
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap`

## Make Filtered VCFs to use to make pseudohaploids
### Directory
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs`
### Submit jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs
qsub gen_filtered_vcf_Chr01K.sh
qsub gen_filtered_vcf_Chr01N.sh

bash

for CHR in {02..09};
  do
  sed 's/01N/'"$CHR"K'/g' gen_filtered_vcf_Chr01N.sh > \
    gen_filtered_vcf_Chr$CHR'K.sh';
  done

for CHR in {02..09};
  do
  sed 's/01N/'"$CHR"N'/g' gen_filtered_vcf_Chr01N.sh > \
    gen_filtered_vcf_Chr$CHR'N.sh';
  done

####

for CHR in {02..09};
  do
  qsub gen_filtered_vcf_Chr$CHR'K.sh';
  qsub gen_filtered_vcf_Chr$CHR'N.sh';
  done
```

## Split VCFs into 100k line files
### Split Chr01K
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

split -l 100000 -d Chr01K_filt.recode.vcf Chr01K_filt_split_
```
### Generate vector of sample names
* `/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs/filtered_vcf_samp_names.txt`
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

head -5 Chr01K_filt.recode.vcf | tail -n 1 | \
cut --complement -f 1-9 > filtered_vcf_samp_names.txt
```
### Split remaining Chromosomes
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

bash
for NUM in {02..09};
  do
  split -l 100000 -d Chr$NUM'K_filt.recode.vcf' Chr$NUM'K_filt_split_';
  done

for NUM in {01..09};
  do
  split -l 100000 -d Chr$NUM'N_filt.recode.vcf' Chr$NUM'N_filt_split_';
  done

```

## Generate Allele Ratio Files
### Test for Chr01K
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs
qsub gen_ratios_Chr01K_00.sh

bash
for SF in {01..05};
do
sed 's/_00/'_"$SF"'/g' gen_ratios_Chr01K_00.sh > gen_ratios_Chr01K_$SF'.sh';
done

bash
for SF in {01..05};
do
qsub gen_ratios_Chr01K_$SF'.sh';
done
```
### Run on rest of chromosomes
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

for CHRNUM in {02..09};
  do
  for SUBFILE in Chr$CHRNUM'K_filt_split'*;
    do
    SUBNUM=`echo $SUBFILE | cut -d '_' -f 4`;
    sed 's/_00/'_"$SUBNUM"'/g; s/01K/'"$CHRNUM"K'/g' gen_ratios_Chr01K_00.sh > gen_ratios_Chr$CHRNUM'K_'$SUBNUM'.sh';
    done 
  done

for CHRNUM in {01..09};
  do
  for SUBFILE in Chr$CHRNUM'N_filt_split'*;
    do
    SUBNUM=`echo $SUBFILE | cut -d '_' -f 4`;
    sed 's/_00/'_"$SUBNUM"'/g; s/01K/'"$CHRNUM"N'/g' gen_ratios_Chr01K_00.sh > gen_ratios_Chr$CHRNUM'N_'$SUBNUM'.sh';
    done
  done

bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs
for CHRNUM in {02..09};
  do
  for SUBFILE in gen_ratios_Chr$CHRNUM'K'*;
    do
    qsub $SUBFILE
    done
  done

for CHRNUM in {01..09};
  do
  for SUBFILE in gen_ratios_Chr$CHRNUM'N'*;
    do
    qsub $SUBFILE
    done
  done
```
### Check that all files finished
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

ls *filt_split_* | wc -l
# 238
ls *alleleratios.rds | wc -l
# 119
```

## Generate pseudohaploid genotypes
### Test with Chr01K
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs
qsub gen_pseudohapgenos_Chr01K_00.sh

bash
for SF in {01..05};
do
sed 's/_00/'_"$SF"'/g' gen_pseudohapgenos_Chr01K_00.sh > \
gen_pseudohapgenos_Chr01K_$SF'.sh';
done

bash
for SF in {01..05};
do
qsub gen_pseudohapgenos_Chr01K_$SF'.sh';
done
```
### Run on rest of chromosomes
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

for CHRNUM in {02..09};
  do
  for SUBFILE in Chr$CHRNUM'K_filt_split'*;
    do
    SUBNUM=`echo $SUBFILE | cut -d '_' -f 4`;
    sed 's/_00/'_"$SUBNUM"'/g; s/01K/'"$CHRNUM"K'/g' \
gen_pseudohapgenos_Chr01K_00.sh > \
gen_pseudohapgenos_Chr$CHRNUM'K_'$SUBNUM'.sh';
    done
  done

for CHRNUM in {01..09};
  do
  for SUBFILE in Chr$CHRNUM'N_filt_split'*;
    do
    SUBNUM=`echo $SUBFILE | cut -d '_' -f 4`;
    sed 's/_00/'_"$SUBNUM"'/g; s/01K/'"$CHRNUM"N'/g' \
gen_pseudohapgenos_Chr01K_00.sh > \
gen_pseudohapgenos_Chr$CHRNUM'N_'$SUBNUM'.sh';
    done
  done

bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs
for CHRNUM in {02..09};
  do
  for SUBFILE in gen_pseudohapgenos_Chr$CHRNUM'K'*;
    do
    qsub $SUBFILE
    done
  done

for CHRNUM in {01..09};
  do
  for SUBFILE in gen_pseudohapgenos_Chr$CHRNUM'N'*;
    do
    qsub $SUBFILE
    done
  done
```

## Generate Distance Matrices
### Test with Chr01K
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

cp gen_pseudohapgenos_Chr01K_00.sh calc_dist_Chr01K_00.sh

qsub calc_dist_Chr01K_00.sh

bash
for SF in {01..05};
do
sed 's/_00/'_"$SF"'/g' calc_dist_Chr01K_00.sh > \
calc_dist_Chr01K_$SF'.sh';
done

bash
for SF in {01..05};
do
qsub calc_dist_Chr01K_$SF'.sh';
done
```
### Run on remaining chromosomes
#### Generate shell scripts
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

for CHRNUM in {02..09};
  do
  for SUBFILE in Chr$CHRNUM'K_filt_split'*'alleleratios.rds';
    do
    SUBNUM=`echo $SUBFILE | cut -d '_' -f 4`;
    sed 's/_00/'_"$SUBNUM"'/g; s/01K/'"$CHRNUM"K'/g' \
calc_dist_Chr01K_00.sh > \
calc_dist_Chr$CHRNUM'K_'$SUBNUM'.sh';
    done
  done

for CHRNUM in {01..09};
  do
  for SUBFILE in Chr$CHRNUM'N_filt_split'*'alleleratios.rds';
    do
    SUBNUM=`echo $SUBFILE | cut -d '_' -f 4`;
    sed 's/_00/'_"$SUBNUM"'/g; s/01K/'"$CHRNUM"N'/g' \
calc_dist_Chr01K_00.sh > \
calc_dist_Chr$CHRNUM'N_'$SUBNUM'.sh';
    done
  done
```
#### Run jobs
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs

bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs
for CHRNUM in {02..09};
  do
  for SUBFILE in calc_dist_Chr$CHRNUM'K'*;
    do
    qsub $SUBFILE
    done
  done

for CHRNUM in {01..09};
  do
  for SUBFILE in calc_dist_Chr$CHRNUM'N'*;
    do
    qsub $SUBFILE
    done
  done
```



