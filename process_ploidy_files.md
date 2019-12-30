# Process MNP Files that Sujan made

## Overview
* 2 Approaches
  1. Count up the total number of MNPs for different sequencing depths
  2. Save the positions of MNPs and filter based on annotation of that position
    * For example, count up MNPs in genic regions
    * Requires more steps for the analysis

## Location of Data
* Location of ATGC count files Sujan generated
  * `/global/cscratch1/sd/sujan/Pvirg_ATGC_counts/`
* Directory with ATCG Counts for chromosomes made for the ~200 libraries \
that have whole-genome files (rather than chromosome files)
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files`
  * I split the whole-genome files into chromsome files and saved the \
results in this directory
* Testing directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir`
* Parent directory for MNP counting results
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts`
* SHOULD make separate .md file for this: Location of sequencing depth files
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth`

## Example scripts for steps
### Divide whole-genome count files into chromosome-scale count files
* `/global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files/make_chr_files_sub_00.sh`
### Count up MNPs with 10+ reads
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10/chr_awk_depth10_00.sh`
### Get positions of MNPS with 10+ reads
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10/chr_MNP_pos_depth10_00.sh`
### Calculate Sequencing Depth
* `/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth/chr_file_get_seq_depth_00.sh`

## Generate Lists of libraries
* Some libraries have chromosome-specific files, others have full-genome files
  * need to make separate lists because these libraries are treated a little \
differently
    * The whole-genome files need to get divided into chromsome files
* in R
```
all_files <- system('ls /global/cscratch1/sd/sujan/Pvirg_ATGC_counts/*bz2',
  intern = T)

chr01k_files <- all_files[grep('Chr01K', all_files)] 

chr01_libs <- gsub('/global/cscratch1/sd/sujan/Pvirg_ATGC_counts/', '',
  gsub('.Chr01K.counts.bz2', '', chr01k_files))

full_files <- all_files[-grep('Chr', all_files)]
full_libs <- gsub('/global/cscratch1/sd/sujan/Pvirg_ATGC_counts/', '',
  gsub('.counts.bz2', '', full_files))
full_libs <- setdiff(full_libs, c('LIB', 'Pvirgatum_MNP.positions.bz2'))

full_lib_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/full_file_libs.txt'
write.table(full_libs, file = full_lib_out_file, quote = F, sep = '\t',
  row.names = F, col.names = F)

chr01_lib_out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/chr_file_libs.txt'
write.table(chr01_libs, file = chr01_lib_out_file, quote = F, sep = '\t',
  row.names = F, col.names = F)
```
### Divide the lists
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/test_dir

split -l 30 -d chr_file_libs.txt chr_file_libs_sub_
split -l 10 -d full_file_libs.txt full_file_libs_sub_
```
* some of the original files needed to be remade, so I need to adjust \
these subfiles
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/test_dir

head -9 full_file_libs_sub_18 > chr_file_libs_sub_28
echo INHR >> full_file_libs_sub_21
rm full_file_libs_sub_18
mv full_file_libs_sub_21 full_file_libs_sub_18
```

## Divide whole-genome count files
* Need to divide the ~200 whole-genome ATCG files into chromosome files \
otherwise the files are too big to process, and it's good to get chromosome \
level counts
* Results in this directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files/`
### Example of what the script does
```
bzip2 -dc $DATA_FILE | awk -v chrname="$CHR" '$0~chrname {print $0}' \
| bzip2 -c > $OUT_FILE
```
### Test with first 10 libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files
sbatch make_chr_files_sub_00.sh
```
### Generate rest of commands
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files

for LIB in {01..21};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' make_chr_files_sub_00.sh \
> make_chr_files_sub_$LIB'.sh';
  done 
```
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files

for LIB in {01..21};
  do
  sbatch make_chr_files_sub_$LIB'.sh';
  done
```

## Count up MNPs with 10+ reads
* Use `awk` to get line numbers with 10+ reads of A; repeat for other 3
  * output gets concatenated into single file
* Use `uniq` and `grep` to count the line numbers repeate 3 or 4 times
* Location of results
  * /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10
### Example of what the script does
```
bzip2 -dc $DATA_FILE | awk '$5 >= 10 {print FNR}' > $OUT_FILE;
bzip2 -dc $DATA_FILE | awk '$6 >= 10 {print FNR}' >> $OUT_FILE;
bzip2 -dc $DATA_FILE | awk '$7 >= 10 {print FNR}' >> $OUT_FILE;
bzip2 -dc $DATA_FILE | awk '$8 >= 10 {print FNR}' >> $OUT_FILE;
echo $CHR `sort $OUT_FILE | uniq -c | grep -c '^\s*[34]\s'` > $KEEP_FILE;
```
### Testing
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/test_dir
sbatch chr_libs_awk_count_TESTING.sh
# took 1hr 55min to go through 2 full libraries; about 1hr per library
```
### Run for Sujan-made Chr Libraries
* Libraries that Sujan made chromosome-scale files for
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10` 
#### First job
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10
sbatch chr_awk_depth10_00.sh
```
#### Generate and submit commands for remaining Sujan-made Chr libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10

for LIB in {01..27};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_awk_depth10_00.sh \
> chr_awk_depth10_$LIB'.sh';
  done

for LIB in {01..27};
  do
  sbatch chr_awk_depth10_$LIB'.sh';
  done
```
### Run for split-up whole-genome Libraries
* Libraries that I had to split up into chromosome-scale files otherwise \
the files were too big
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10

# cp chr_awk_depth10_00.sh full_split_depth10_00.sh
# edit full_split_depth10_00.sh in vim

for LIB in {01..21};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_split_depth10_00.sh \
> full_split_depth10_$LIB'.sh';
  done
```
#### Submit Commands
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10

for LIB in {00..21};
  do
  sbatch full_split_depth10_$LIB'.sh';
  done
```

## Count up MNPs with 20+ reads
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth20`
### Generate commands
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth20

cp ../depth10/chr_awk_depth10_00.sh ./chr_awk_depth20_00.sh
# edited with vim

for LIB in {01..28};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_awk_depth20_00.sh \
> chr_awk_depth20_$LIB'.sh';
  done

for LIB in {00..28};
  do
  sbatch chr_awk_depth20_$LIB'.sh';
  done

cp ../depth10/full_split_depth10_00.sh ./full_split_depth20_00.sh
# edited with vim

for LIB in {01..20};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_split_depth20_00.sh \
> full_split_depth20_$LIB'.sh';
  done

for LIB in {00..20};
  do
  sbatch full_split_depth20_$LIB'.sh';
  done
```

## Get Chromosome Positions of MNPs with 10+ reads
* Results Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10`
### Example of what the script does
```
bzip2 -dc $DATA_FILE | awk '$5 >= 10 {print $2}' > $OUT_FILE;
bzip2 -dc $DATA_FILE | awk '$6 >= 10 {print $2}' >> $OUT_FILE;
bzip2 -dc $DATA_FILE | awk '$7 >= 10 {print $2}' >> $OUT_FILE;
bzip2 -dc $DATA_FILE | awk '$8 >= 10 {print $2}' >> $OUT_FILE;
sort $OUT_FILE | uniq -c | grep '^\s*[34]\s' | \
            sed 's/^\s*[0-9]*\s//' | sed '$d' > $KEEP_FILE;
```
### Generate scripts
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10

cp /global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/chr_libs_MNP_position_TESTING.sh ./chr_MNP_pos_depth10_00.sh

sbatch chr_MNP_pos_depth10_00.sh

for LIB in {01..28};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_MNP_pos_depth10_00.sh \
> chr_MNP_pos_depth10_$LIB'.sh';
  done

for LIB in {01..28};
  do
  sbatch chr_MNP_pos_depth10_$LIB'.sh';
  done

cp chr_MNP_pos_depth10_00.sh full_MNP_pos_depth10_00.sh

for LIB in {01..20};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_MNP_pos_depth10_00.sh \
> full_MNP_pos_depth10_$LIB'.sh';
  done

for LIB in {00..20};
  do
  sbatch full_MNP_pos_depth10_$LIB'.sh';
  done
```
### Redo for libraries that were reprocessed
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10

cp chr_MNP_pos_depth10_00.sh chr_MNP_pos_depth10_redo01.sh
# adjust with vim

sbatch chr_MNP_pos_depth10_redo01.sh
```

## Get Chromsome Positions of MNPs with 20+ reads
* Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20`
### Generate Scripts
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10/chr_MNP_pos_depth10_00.sh ./chr_MNP_pos_depth20_00.sh

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_10/full_MNP_pos_depth10_00.sh ./full_MNP_pos_depth20_00.sh

for LIB in {01..28};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_MNP_pos_depth20_00.sh \
> chr_MNP_pos_depth20_$LIB'.sh';
  done

for LIB in {00..28};
  do
  sbatch chr_MNP_pos_depth20_$LIB'.sh';
  done

for LIB in {01..20};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_MNP_pos_depth20_00.sh \
> full_MNP_pos_depth20_$LIB'.sh';
  done

for LIB in {00..20};
  do
  sbatch full_MNP_pos_depth20_$LIB'.sh';
  done
```
### Rerun for the 9 libraries that did not get processed properly
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20
cp chr_MNP_pos_depth20_00.sh chr_MNP_pos_depth20_redo01.sh

# adjust in vim

sbatch chr_MNP_pos_depth20_redo01.sh
```

## Get Chromsome Positions of MNPs with 15+ reads
* Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_15`
### Generate Scripts
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_15

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/chr_MNP_pos_depth20_00.sh ./chr_MNP_pos_depth15_00.sh

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/full_MNP_pos_depth20_00.sh ./full_MNP_pos_depth15_00.sh

for LIB in {01..28};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_MNP_pos_depth15_00.sh \
> chr_MNP_pos_depth15_$LIB'.sh';
  done

for LIB in {00..28};
  do
  sbatch chr_MNP_pos_depth15_$LIB'.sh';
  done

for LIB in {01..20};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_MNP_pos_depth15_00.sh \
> full_MNP_pos_depth15_$LIB'.sh';
  done

for LIB in {00..20};
  do
  sbatch full_MNP_pos_depth15_$LIB'.sh';
  done
```

## Get Chromsome Positions of MNPs with 25+ reads
* Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_25`
### Generate Scripts
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_25
  
cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/chr_MNP_pos_depth20_00.sh ./chr_MNP_pos_depth25_00.sh

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/full_MNP_pos_depth20_00.sh ./full_MNP_pos_depth25_00.sh

for LIB in {01..28};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_MNP_pos_depth25_00.sh \
> chr_MNP_pos_depth25_$LIB'.sh';
  done

for LIB in {00..28};
  do
  sbatch chr_MNP_pos_depth25_$LIB'.sh';
  done

for LIB in {01..20};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_MNP_pos_depth25_00.sh \
> full_MNP_pos_depth25_$LIB'.sh';
  done

for LIB in {00..20};
  do
  sbatch full_MNP_pos_depth25_$LIB'.sh';
  done

# Generate two new script: one for a failed node and one for results I \
# deleted like a moron

cp chr_MNP_pos_depth25_00.sh chr_MNP_pos_depth25_finish20.sh

cp chr_MNP_pos_depth25_00.sh chr_MNP_pos_depth25_rerun_IGN.sh

sbatch chr_MNP_pos_depth25_finish20.sh
sbatch chr_MNP_pos_depth25_rerun_IGN.sh
```

## Get Chromsome Positions of MNPs with 5+ reads
* Directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_05`
### Generate Scripts
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_05

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/chr_MNP_pos_depth20_00.sh ./chr_MNP_pos_depth05_00.sh

cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/full_MNP_pos_depth20_00.sh ./full_MNP_pos_depth05_00.sh

for LIB in {01..28};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_MNP_pos_depth05_00.sh \
> chr_MNP_pos_depth05_$LIB'.sh';
  done

for LIB in {00..28};
  do
  sbatch chr_MNP_pos_depth05_$LIB'.sh';
  done

for LIB in {01..20};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_MNP_pos_depth05_00.sh \
> full_MNP_pos_depth05_$LIB'.sh';
  done

for LIB in {00..20};
  do
  sbatch full_MNP_pos_depth05_$LIB'.sh';
  done
```

###################

## Get Chromosome Positions for MNPs with 3+, 4+, 6+ and 7+ reads
* Will do as much of the prep in parallel as possible
* Directories
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_03`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_04`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_06`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_07`
### Generate Scripts
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts

# Copy scripts from previous depths
for DEPTH in 03 04 06 07;
  do
  cd ./pos_$DEPTH;
  cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/chr_MNP_pos_depth20_00.sh chr_MNP_pos_depth$DEPTH'_00.sh';
  cp /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/pos_20/full_MNP_pos_depth20_00.sh full_MNP_pos_depth$DEPTH'_00.sh';
  cd ..;
  done 
# adjust scripts with vim

# Make copies for each library name sub-file
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts

for DEPTH in 03 04 06 07;
  do
  cd ./pos_$DEPTH; 
  for LIB in {01..28};
    do 
    sed 's/sub_00/'sub_"$LIB"'/g' chr_MNP_pos_depth$DEPTH'_00.sh' \
      > 'chr_MNP_pos_depth'$DEPTH'_'$LIB'.sh';
    done;
  cd ..;
  done

for DEPTH in 03 04 06 07;
  do
  cd ./pos_$DEPTH;
  for LIB in {01..20};
    do
    sed 's/sub_00/'sub_"$LIB"'/g' full_MNP_pos_depth$DEPTH'_00.sh' \
      > 'full_MNP_pos_depth'$DEPTH'_'$LIB'.sh';
    done;
  cd ..;
  done

# submit jobs
cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts

for DEPTH in 03 04 06 07;
  do
  cd ./pos_$DEPTH;
  for LIB in {00..28};
    do
    sbatch 'chr_MNP_pos_depth'$DEPTH'_'$LIB'.sh'
    done;
  cd ..;
  done

cd /global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts

for DEPTH in 03 04 06 07;
  do
  cd ./pos_$DEPTH;
  for LIB in {00..20};
    do
    sbatch 'full_MNP_pos_depth'$DEPTH'_'$LIB'.sh'
    done;
  cd ..;
  done
```



#####################





## Tally total seq coverage from Sujans ATCG count files
* Location of results
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth`
### Example of what script does
```
echo $CHR `bzip2 -dc $DATA_FILE | awk '{sum+=$4} END {print sum}'` \
   > $OUT_FILE
```
### Test
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth

sbatch chr_file_get_seq_depth_00.sh
```
* works
### Generate and run scripts
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth

for LIB in {01..27};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_file_get_seq_depth_00.sh \
> chr_file_get_seq_depth_$LIB'.sh';
  done

for LIB in {01..27};
  do
  sbatch chr_file_get_seq_depth_$LIB'.sh';
  done

# cp chr_file_get_seq_depth_00.sh full_split_get_seq_depth_00.sh
# adjust full_split_get_seq_depth_00.sh with vim

for LIB in {01..21};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_split_get_seq_depth_00.sh \
> full_split_get_seq_depth_$LIB'.sh';
  done

# Submitted 12/2
for LIB in {00..21};
  do
  sbatch full_split_get_seq_depth_$LIB'.sh';
  done
```
### Rerun for the 9 files that did not get processed properly
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth

cp chr_file_get_seq_depth_00.sh chr_file_get_seq_depth_redo01.sh
# adjust commands with vim

sbatch chr_file_get_seq_depth_redo01.sh
```


