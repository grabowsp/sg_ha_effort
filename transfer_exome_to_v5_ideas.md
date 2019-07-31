# Plan for converting exome-capture SNPs to v5 genome

## Steps
1. Find/load exome-capture data
1. Find/load switchgrass genome
2. Get coordinates for a SNP
3. Extract surrounding 100, 150, 200, 500bp in both directions
* some options: samtools faidx; bedtools getfasta
4. Use blast (or something faster - blat?) to find hits in v5 genome
5. Find location of SNP in homologous sequence
6. Figure out how fast running on one SNP is.
7. Try writing steps in Python
8. Compare speeds
9.

## Data locations
### v4 exome snps
* `/home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/combined_all_1168_samples_filtered_matrix_renamed_chroms_v4_coord.txt`
  * 7.3G zipped, 38G unzipped
  * 5,596,351 SNPs
### v4 genome
* `/home/t4c1/WORK/grabowsk/data/switchgrass/ref_seqs/v4_seq/Pvirgatum_450_v4.0.fa`
  * 349M zipped, 1.6G unzipped
### v5 genome
* `/home/t4c1/WORK/grabowsk/data/switchgrass/ref_seqs/v5_seq/Pvirgatum_516_v5.0.fa`
  * 318M zipped, 1.1G unzipped

## Divide v4 exome snps into subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps
mkdir sub_files
split combined_all_1168_samples_filtered_matrix_renamed_chroms_v4_coord.txt \
-l 100000 -d v4_exome_snps_sub
mv v4_exome_snps_sub* ./sub_files
cd sub_files
tail -n 99999 v4_exome_snps_sub00 > tmp_sub00
mv tmp_sub00 v4_exome_snps_sub00
```

## Get header line from snp file
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps
head -1 combined_all_1168_samples_filtered_matrix_renamed_chroms_v4_coord.txt \
> v4_combo_snps_header.txt
```

## Generate v4 BED files for mapping and SNP info tables for subfiles
### Overview
* Load SNP file
* Use SNP positions to bunch SNPs into loci
* Generate BED file with start and end positions of loci to be used \
to extract v4 sequence using BEDtools and then mapping to v5 reference with \
bwa mem
* Make table with info about each SNP in relation to loci used for mapping \
to v5 - this is important for determining exact location from bwa results
### R script
* `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/make_v4_df_and_bed.r`
### Testing with subfile 1
* Test sh file
  * /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files/make_v4_BED_00.sh
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub make_v4_BED_00.sh
```
* works
### Make .sh for remaining subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./make_v4_BED_00.sh > ./make_v4_BED_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub make_v4_BED_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*v4BED.e*
rm sub*v4BED.o*
rm sub*v4BED.p*
```

## Get Loci v4 seq using BED files and bedtools
### Test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp make_v4_BED_00.sh bedtools_v4seq_00.sh
```
* adjusted commands in vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub bedtools_v4seq_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./bedtools_v4seq_00.sh > ./bedtools_v4seq_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub bedtools_v4seq_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*bedtools.e*
rm sub*bedtools.o*
rm sub*bedtools.p*
```

## Map v4 loci to v5 reference with bwa
### test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp bedtools_v4seq_00.sh bwa_v4v5_00.sh
```
* adjust with vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub bwa_v4v5_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./bwa_v4v5_00.sh > ./bwa_v4v5_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub bwa_v4v5_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*bwa_v4v5.e*
rm sub*bwa_v4v5.o*
rm sub*bwa_v4v5.p*
```

## Connect bwa mapping info to SNP info
### Test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp make_v4_BED_00.sh connect_v4_and_v5_00.sh
```
* adjusted script with vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub connect_v4_and_v5_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./connect_v4_and_v5_00.sh > \
./connect_v4_and_v5_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub connect_v4_and_v5_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*connect_v4v5.e*
rm sub*connect_v4v5.o*
rm sub*connect_v4v5.p*
```

# Test that there's no hard-masking
```
file_inds <- sprintf('%02d', c(0:55))
data_dir <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/exome/', 
  'v4_snps/sub_files/', sep = '')

for(i in file_inds){
  print(paste('importing', i))
  tmp_data_file <- paste(data_dir, 'snps_sub', i, '_v4_and_v5info.rds', 
  sep = '')
  tmp_data <- readRDS(tmp_data_file)
  h_inds <- grep('H', tmp_data$v5_cigar)
  if(length(h_inds) > 0){
    print(paste('H in', i))
  }
#  s_inds <- grep('S', tmp_data$v5_cigar)
#  if(length(s_inds) > 0){
#    print(paste('S in', i))
  }
}

for(i in file_inds){
  print(paste('importing', i))
  tmp_data_file <- paste(data_dir, 'snps_sub', i, '_v4_and_v5info.rds',
  sep = '')
  tmp_data <- readRDS(tmp_data_file)
  per_nas <- sum(is.na(tmp_data$v5_pos))/nrow(tmp_data)
  print(per_nas)
}
```
* no hard masking, so scripts works fine for now
* about 1-2% of SNPs are unmapped to v5 at this point
  * should run a second round with smaller loci

## Generate BED files and snp info for second round (try2) of mapping 
### R script
* `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/make_unmappped_bed.r`
### Testing with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp make_v4_BED_00.sh make_try2_BED_00.sh
```
* adjust submit script in vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub make_try2_BED_00.sh
```
* works
### Make .sh for remaining subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./make_try2_BED_00.sh > ./make_try2_BED_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub make_try2_BED_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*try2BED.e*
rm sub*try2BED.o*
rm sub*try2BED.p*
```

## Get v4 seq using try2 BED files and bedtools
### Test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp bedtools_v4seq_00.sh bedtools_try2_00.sh
```
* adjusted commands in vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub bedtools_try2_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./bedtools_try2_00.sh > ./bedtools_try2_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub bedtools_try2_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*bed_try2.e*
rm sub*bed_try2.o*
rm sub*bed_try2.p*
```

## Map try2 loci to v5 reference with bwa
### test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp bwa_v4v5_00.sh bwa_try2_00.sh
```
* adjust with vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub bwa_try2_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./bwa_try2_00.sh > ./bwa_try2_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub bwa_try2_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*bwa_try2.e*
rm sub*bwa_try2.o*
rm sub*bwa_try2.p*
```

## Connect try2 bwa mapping info to SNP info
### Test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp connect_v4_and_v5_00.sh connect_v4v5_try2_00.sh
```
* adjusted script with vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub connect_v4v5_try2_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./connect_v4v5_try2_00.sh > \
./connect_v4v5_try2_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub connect_v4v5_try2_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*connect_try2.e*
rm sub*connect_try2.o*
rm sub*connect_try2.p*
```
# Test that there's no hard-masking in try2 cigar scores
```
file_inds <- sprintf('%02d', c(0:55))
data_dir <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/exome/', 
  'v4_snps/sub_files/', sep = '')

for(i in file_inds){
  print(paste('importing', i))
  tmp_data_file <- paste(data_dir, 'sub', i, '_try2_v4_and_v5info.rds', 
  sep = '')
  tmp_data <- readRDS(tmp_data_file)
  h_inds <- grep('H', tmp_data$v5_cigar)
  if(length(h_inds) > 0){
    print(paste('H in', i))
  }
#  s_inds <- grep('S', tmp_data$v5_cigar)
#  if(length(s_inds) > 0){
#    print(paste('S in', i))
#  }
}

```
* no hard masking, so scripts works fine for now


## Generate BED files and snp info for second round (try3) of mapping 
### R script
* `/home/grabowsky/tools/workflows/sg_ha_effort/r_scripts/make_unmappped_bed.r`
### Testing with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp make_try2_BED_00.sh make_try3_BED_00.sh
```
* adjust submit script in vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub make_try3_BED_00.sh
```
* works
### Make .sh for remaining subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./make_try3_BED_00.sh > ./make_try3_BED_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub make_try3_BED_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*try3BED.e*
rm sub*try3BED.o*
rm sub*try3BED.p*
```

## Get v4 seq using try3 BED files and bedtools
### Test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp bedtools_try2_00.sh bedtools_try3_00.sh
```
* adjusted commands in vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub bedtools_try3_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./bedtools_try3_00.sh > ./bedtools_try3_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub bedtools_try3_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*bed_try3.e*
rm sub*bed_try3.o*
rm sub*bed_try3.p*
```

## Map try3 loci to v5 reference with bwa
### test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp bwa_try2_00.sh bwa_try3_00.sh
```
* adjust with vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub bwa_try3_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./bwa_try3_00.sh > ./bwa_try3_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub bwa_try3_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*bwa_try3.e*
rm sub*bwa_try3.o*
rm sub*bwa_try3.p*
```

## Connect try3 bwa mapping info to SNP info
### Test with subfile 1
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
cp connect_v4v5_try2_00.sh connect_v4v5_try3_00.sh
```
* adjusted script with vim
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
qsub connect_v4v5_try3_00.sh
```
### Make submit files for rest of subfiles
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do sed 's/sub00/sub'"$SF"'/g' ./connect_v4v5_try3_00.sh > \
./connect_v4v5_try3_$SF.sh;
done
```
### Submit jobs
```
bash
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
for SF in {01..55};
do qsub connect_v4v5_try3_$SF.sh;
done
```
### Clean up output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/exome/v4_snps/sub_files
rm sub*connect_try3.e*
rm sub*connect_try3.o*
rm sub*connect_try3.p*
```

## Test that there's no hard-masking in try3 cigar scores
```
file_inds <- sprintf('%02d', c(0:55))
data_dir <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/exome/', 
  'v4_snps/sub_files/', sep = '')

for(i in file_inds){
  print(paste('importing', i))
  tmp_data_file <- paste(data_dir, 'sub', i, '_try3_v4_and_v5info.rds', 
  sep = '')
  tmp_data <- readRDS(tmp_data_file)
  h_inds <- grep('H', tmp_data$v5_cigar)
  if(length(h_inds) > 0){
    print(paste('H in', i))
  }
#  s_inds <- grep('S', tmp_data$v5_cigar)
#  if(length(s_inds) > 0){
#    print(paste('S in', i))
#  }
}
```
* no hard masking

## Compare missing data between runs
```
file_inds <- sprintf('%02d', c(0:55))
data_dir <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/exome/',
  'v4_snps/sub_files/', sep = '')

file_add_vec_1 <- c('snps_', '', '')
file_add_vec_2 <- c('', '_try2', '_try3')

for(i in file_inds){
  print(paste('importing', i))
  for(j in c(1:3)){
    tmp_data_file <- paste(data_dir, file_add_vec_1[j], 'sub', i, 
      file_add_vec_2[j], '_v4_and_v5info.rds', sep = '')
    tmp_data <- readRDS(tmp_data_file)
    per_nas <- sum(is.na(tmp_data$v5_pos))/nrow(tmp_data)
    print(per_nas)
  }
}

```
* missing data decreases with each round

## Stats on mapping to v5
```
file_inds <- sprintf('%02d', c(0:55))
data_dir <- paste('/home/t4c1/WORK/grabowsk/data/switchgrass/exome/',
  'v4_snps/sub_files/', sep = '')

data_list <- list()

for(i in file_inds){
  print(paste('importing', i))
  tmp_data_file <- paste(data_dir, 'sub', i, '_try3_v4_and_v5info.rds',
    sep = '')
  tmp_data <- readRDS(tmp_data_file)
  # combine chromosome info from all runs
  t2_chr_inds <- which(is.na(tmp_data$try2_v5_chr) == F)
  tmp_data$v5_chr[t2_chr_inds] <- tmp_data$try2_v5_chr[t2_chr_inds]
  t3_chr_inds <- which(is.na(tmp_data$try3_v5_chr) == F)
  tmp_data$v5_chr[t3_chr_inds] <- tmp_data$try3_v5_chr[t3_chr_inds]
  tot_na_inds <- which(is.na(tmp_data$v5_pos))
  tmp_data$v5_chr[tot_na_inds] <- NA
  #
  data_list[[i]][['v4_chr']] <- tmp_data$v4_chr
  data_list[[i]][['v4_pos']] <- tmp_data$orig_pos
  data_list[[i]][['v5_chr']] <- tmp_data$v5_chr
  data_list[[i]][['v5_pos']] <- tmp_data$v5_pos
}

tot_v4_chr_vec <- unlist(lapply(data_list, function(x) x[['v4_chr']])) 
tot_v5_chr_vec <- unlist(lapply(data_list, function(x) x[['v5_chr']]))


sum(is.na(tot_v5_chr_vec))/length(tot_v5_chr_vec)
# [1] 0.01632957
# 1.6% of SNPs not mapped

sum(tot_v4_chr_vec == tot_v5_chr_vec, na.rm = T)/length(tot_v4_chr_vec)
# [1] 0.8766353
# 87.7% of SNPs mapped to same chromosome, so 10% of SNPs map to different
#  chromsome

chr_comp_vec <- paste(tot_v4_chr_vec, tot_v5_chr_vec, sep = '_')
table(chr_comp_vec)
# Long table - take home: see mapping to many other chromsomes, not just
#  homeologs, though majority of off-v4-chrom mapping is to homeologs

```



## Manually inspect quality of mapping
### v4
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/ref_seqs/v4_seq
samtools faidx Pvirgatum_450_v4.0.fa Chr01K:1657454-1657533
```
### v5
* forward strand
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/ref_seqs/v5_seq
samtools faidx Pvirgatum_516_v5.0.fa Chr01K:1798152-1798197
```
* reverse strand
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/ref_seqs/v5_seq
samtools faidx -i Pvirgatum_516_v5.0.fa Chr01N:1692297-1692376
```

