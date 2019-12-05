# Process Files that Sujan made

## Overview
* Look at Chromosome Chr01K and count up the number of positions with 3 or 4 \
alleles with a count greater than a set depth
  * starting with 10
* Need to normalize by sequencing depth
* Use `awk` to get the line numbers of that have the minimum depth for \
each allele, then use R to find the line numbers represented 3+ times

## Set up Directories
* Testing directory
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir`
* Directory with ATCG Counts for chromosomes made for the ~200 libraries \
that have whole-genome files (rather than chromosome files)
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files`
* Directory with MNP counts for different seq depth thresholds
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts`
### Copy some files for testing
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/test_dir

cp /global/cscratch1/sd/sujan/Pvirg_ATGC_counts/XXBC.Chr01K.counts.bz2 .

cp /global/cscratch1/sd/sujan/Pvirg_ATGC_counts/PZXN.counts.bz2 .
```

## Example scripts for steps
### Divide whole-genome count files into chromosome-scale count files
* `/global/cscratch1/sd/grabowsp/sg_ploidy/chr_count_files/make_chr_files_sub_00.sh`
### Count up MNPs with 10+ reads
* `/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10/chr_awk_depth10_00.sh`
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

## Divide whole-genome count files
* Need to divide the ~200 whole-genome ATCG files into chromosome files \
otherwise the files are too big to process, and it's good to get chromosome \
level counts
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

for LIB in {01..27};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' chr_awk_depth20_00.sh \
> chr_awk_depth20_$LIB'.sh';
  done

for LIB in {00..27};
  do
  sbatch chr_awk_depth20_$LIB'.sh';
  done

cp ../depth10/full_split_depth10_00.sh ./full_split_depth20_00.sh
# edited with vim

for LIB in {01..21};
  do
  sed 's/sub_00/'sub_"$LIB"'/g' full_split_depth20_00.sh \
> full_split_depth20_$LIB'.sh';
  done

for LIB in {00..21};
  do
  sbatch full_split_depth20_$LIB'.sh';
  done


```


## Tally total seq coverage from Sujans ATCG count files
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
### Combine results
* in R
```
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/seq_depth/'

depth_files <- system(paste('ls ', data_dir, '*tot_seq_depth.txt', sep = ''), 
  intern = T)

lib_names <- gsub(data_dir, '', gsub('_tot_seq_depth.txt', '', depth_files))

######
# check that all libraries are represented in the files
chr_libs <- read.table('/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/chr_file_libs.txt', header = F, stringsAsFactors = F)
full_libs <- read.table('/global/cscratch1/sd/grabowsp/sg_ploidy/test_dir/full_file_libs.txt', header = F, stringsAsFactors = F)

test_libs <- c(chr_libs[,1], full_libs[,1])
setdiff(test_libs, lib_names)
# character(0)
setdiff(lib_names, test_libs)
# character(0)
## all libraries are present

rm_libs <- c('INHH', 'INHI', 'INHJ', 'INHK', 'INHL', 'INHM', 'INHN', 
  'INHP', 'INHQ')

rm_file_inds <- c()
for(i in rm_libs){
  rm_file_inds <- c(rm_file_inds, grep(i, depth_files))
}

tmp_libs <- setdiff(lib_names, rm_libs)
tmp_files <- depth_files[-rm_file_inds]
#####

depth_list <- list()

# Temporary script until fix the few unfinished libraries

#for(i in seq(length(lib_names))){
#  depth_list[[lib_names[i]]] <- read.table(depth_files[i], header = F, 
#  stringsAsFactors = F, sep = ' ')
#  print(lib_names[i])
#}

for(i in seq(length(tmp_libs))){
  depth_list[[tmp_libs[i]]] <- read.table(tmp_files[i], header = F,
    stringsAsFactors = F, sep = ' ')
}

tot_coverage <- unlist(lapply(depth_list, function(x) sum(as.numeric(x[,2]))))

depth_mat <- matrix(
  data = unlist(lapply(depth_list, function(x) as.numeric(x[,2]))),
  nrow = length(depth_list), byrow = T
)
colnames(depth_mat) <- depth_list[[1]][,1]
rownames(depth_mat) <- names(depth_list)

# calculate the portion of coverage in each chromosome
port_mat <- t(apply(depth_mat, 1, function(x) x/sum(x)))

port_mean <- apply(port_mat, 2, mean)

port_sd <- apply(port_mat, 2, sd)

# look for outliers - chromosomes that are 4 SD away from the mean
#  of the coverage proportion for that chromosome across all samples
cov_outlier_list <- list()
for(i in seq(ncol(port_mat))){
  cov_outlier_list[[i]] <- which(
    abs(port_mat[,i] - port_mean[i]) - (port_sd[i]*4) > 0)
}

cov_out_tab <- table(unlist(cov_outlier_list))

weird_cov_libs <- names(depth_list)[
  as.numeric(names(cov_out_tab)[cov_out_tab > 2])]

### Input the MNP counts
count_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/MNP_counts/depth10/'

count_list <- list()

for(LIB in tmp_libs){
  tmp_count_file <- paste(count_dir, LIB, '_10_MNP_total_count.txt', sep = '')
  count_list[[LIB]] <- read.table(tmp_count_file, header = F, sep = ' ',
    stringsAsFactors = F)
}

tot_mnp <- unlist(lapply(count_list, function(x) sum(as.numeric(x[,2]))))

tot_mnp_corrected <- (tot_mnp / tot_coverage) * min(tot_coverage)

tot_mnp['IIAS'] / tot_coverage['IIAS'] * min(tot_coverage)
```


## Testing of saving the MNP positions
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/test_dir

# test_10k.txt
# chr_libs_MNP_position_TESTING.sh

awk '$5 >= 2 {print $2}' test_10k.txt > test_pos_output.txt

head -45 test_pos_output.txt > test_pos_short.txt

cat test_pos_output.txt test_pos_short.txt > test_pos_output_2.txt 

sort test_pos_output_2.txt | uniq -c | grep '^\s*2\s' | sed 's/^\s*[0-9]*\s//' | sed '$d'
```
