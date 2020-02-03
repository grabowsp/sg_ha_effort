# Test nQuire to Decide on Appropiate Parameters for Switchgrass

## Notes
* Need to use jamo to find the location of all the bams that Sujan made
* Tested on ABHM - the first library in the file
* Using all the SNPs seems to be too noisy - nQuire says ABHM is 8X while 
  all other methods indicate it is a 4X
* So, will make BED files for CDS and genic regions and re-run the tests that 
way
* I will also test a 4X and 8X
  * 4X = ABHM
  * 8X = IENS


## Test nQuire using CDS BED and 4X sample
* CDS BED location
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/sg_v5_CDS.bed`
### 4X
```

module load python/3.7-anaconda-2019.07
source activate poplar_align

ln -s /global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/ABHM.Pvirgatum.V5.merged.gatk.bam /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/ABHM.bam

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
samtools index ABHM.bam

## Default values: q = 1, c = 10
#/global/homes/g/grabowsp/tools/nQuire/nQuire create -b \
#/global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/ABHM.bam \
#-q 1 -c 10 \
#-o test_4X_CDS_interactive \
#-r /global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/sg_v5_CDS.bed -y


-o test_4X_CDS_interactive
#-o test_4X_CDS_c10_q01

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_4X_CDS_c10_q01.workflow.sh

sbatch test_4X_CDS_c20_q20.workflow.sh
sbatch test_4X_CDS_c20_q30.workflow.sh
sbatch test_4X_CDS_c30_q20.workflow.sh
sbatch test_4X_CDS_c30_q30.workflow.sh
sbatch test_4X_CDS_c40_q20.workflow.sh
sbatch test_4X_CDS_c40_q30.workflow.sh
sbatch test_4X_CDS_c50_q20.workflow.sh
sbatch test_4X_CDS_c50_q30.workflow.sh

#### Concatenate results
cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
for i in test_4X*lrd_results.txt; do tail -n 1 $i >> \
test_4X_combo_test_lrd_results.txt; done

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
for i in test_4X_CDS*bin; do echo $i; done

for i in test_4X_CDS*bin;
  do /global/homes/g/grabowsp/tools/nQuire/nQuire view $i | wc -l \
    >> test_4X_CDS_nSNPs.txt;
  done
 
```
### Process results in R
```
lrd_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_4X_combo_test_lrd_results.txt', sep = '')
lrd_res <- read.table(lrd_res_file, header = F, sep = '\t',
  stringsAsFactors = F)

lrd_head_0 <- system(paste('head -1 ',
  '/global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/',
  'test_full_1_lrd_results.txt', sep = ''), intern = T)

lrd_head <- unlist(strsplit(lrd_head_0, split = '\t'))

colnames(lrd_res) <- lrd_head

nquire_nsnps_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_4X_CDS_nSNPs.txt', sep = '')
nquire_nsnps <- read.table(nquire_nsnps_file, header = F, stringsAsFactors = F)

lrd_res$nSNPs <- nquire_nsnps[,1]

lrd_res$coverage <- c(10,10,20,20,20,20,30,30,30,30,40,40,40,40,50,50,50,50)
lrd_res$qual <- c(1,1,20,20,30,30,20,20,30,30,20,20,30,30,20,20,30,30)

lrd_res$d_dip_portion <- (lrd_res$d_dip /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tri_portion <- (lrd_res$d_tri /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tet_portion <- (lrd_res$d_tet /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

combo_processed_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_4X_lrd_results_processed.txt', sep = '')
write.table(lrd_res, file = combo_processed_file, quote = F, sep = '\t',
  row.names = F, col.names = T)


```
### 8X
```
module load python/3.7-anaconda-2019.07
source activate poplar_align

ln -s /global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/IENS.Pvirgatum.V5.merged.gatk.bam /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/IENS.bam

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
samtools index IENS.bam

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_8X_CDS_c10_q01.workflow.sh

sbatch test_8X_CDS_c20_q20.workflow.sh
sbatch test_8X_CDS_c20_q30.workflow.sh
sbatch test_8X_CDS_c30_q20.workflow.sh
sbatch test_8X_CDS_c30_q30.workflow.sh
sbatch test_8X_CDS_c40_q20.workflow.sh
sbatch test_8X_CDS_c40_q30.workflow.sh
sbatch test_8X_CDS_c50_q20.workflow.sh
sbatch test_8X_CDS_c50_q30.workflow.sh

#### Concatenate results
cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
for i in test_8X*lrd_results.txt; do tail -n 1 $i >> \
test_8X_combo_test_lrd_results.txt; done

for i in test_8X_CDS*bin;
  do /global/homes/g/grabowsp/tools/nQuire/nQuire view $i | wc -l \
    >> test_8X_CDS_nSNPs.txt;
  done

```
### Process results in R
```
lrd_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_8X_combo_test_lrd_results.txt', sep = '')
lrd_res <- read.table(lrd_res_file, header = F, sep = '\t',
  stringsAsFactors = F)

lrd_head_0 <- system(paste('head -1 ',
  '/global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/',
  'test_full_1_lrd_results.txt', sep = ''), intern = T)

lrd_head <- unlist(strsplit(lrd_head_0, split = '\t'))

colnames(lrd_res) <- lrd_head

nquire_nsnps_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_8X_CDS_nSNPs.txt', sep = '')
nquire_nsnps <- read.table(nquire_nsnps_file, header = F, stringsAsFactors = F)

lrd_res$nSNPs <- nquire_nsnps[,1]

lrd_res$coverage <- c(10,10,20,20,20,20,30,30,30,30,40,40,40,40,50,50,50,50)
lrd_res$qual <- c(1,1,20,20,30,30,20,20,30,30,20,20,30,30,20,20,30,30)

lrd_res$d_dip_portion <- (lrd_res$d_dip /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tri_portion <- (lrd_res$d_tri /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tet_portion <- (lrd_res$d_tet /
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

combo_processed_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_8X_lrd_results_processed.txt', sep = '')
write.table(lrd_res, file = combo_processed_file, quote = F, sep = '\t',
  row.names = F, col.names = T)


```
* Take-homes:
  * For 4X, the d_dip_portion increases with coverage while the d_tet_portion
decreases with coverage
  * For 8X, the d_dip_portion and d_tet_portion stay consistent with 
increasing coverage
  * Denoising does not improve the patterns much
  * using Quality of 30 (-q 30) is better than 20
  * The "true" pattern for the 4X sample comes out at depth of 50 but the change in d_dip_portion is apparent at lower coverage 

* Moving forward
  * Make plots showing the patterns with increasing coverage
  * Run nquire (without denoising) at 20, 30, 40, and 50 (can add intermediate 
points if need be
  * Check results and number of SNPs at depth of 50
  * Look at changes in d_PLOIDY_portion to confirm call

## Generate plots
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

library(ggplot2)

res_4X_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_4X_lrd_results_processed.txt', sep = '')

res_4X <- read.table(res_4X_file, header = T, sep = '\t', stringsAsFactors = F)


res_8X_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/test_8X_lrd_results_processed.txt', sep = '')

res_8X <- read.table(res_8X_file, header = T, sep = '\t', stringsAsFactors = F)

res_4X$ploidy <- '4X'
res_8X$ploidy <- '8X'

tot_res <- rbind(res_4X, res_8X)

raw_30_inds <- intersect(grep('bedcc.bin', tot_res$file), 
  which(tot_res$qual == 30))
denoise_30_inds <- intersect(grep('denoised', tot_res$file), 
  which(tot_res$qual == 30))

raw_30_df <- tot_res[raw_30_inds, ]
denoise_30_df <- tot_res[denoise_30_inds, ]

raw_30_dip <- data.frame(file = raw_30_df$file, coverage = raw_30_df$coverage,  
  d_loglik_portion = raw_30_df$d_dip_portion, ploidy_test = 'dip', 
  ploidy = raw_30_df$ploidy, stringsAsFactors = F)

raw_30_tri <- data.frame(file = raw_30_df$file, coverage = raw_30_df$coverage,
  d_loglik_portion = raw_30_df$d_tri_portion, ploidy_test = 'tri', 
  ploidy = raw_30_df$ploidy, stringsAsFactors = F)

raw_30_tet <- data.frame(file = raw_30_df$file, coverage = raw_30_df$coverage,
  d_loglik_portion = raw_30_df$d_tet_portion, ploidy_test = 'tet', 
  ploidy = raw_30_df$ploidy, stringsAsFactors = F)

raw_30_combo <- rbind(rbind(raw_30_dip, raw_30_tri), raw_30_tet)

raw_30_combo_4X <- raw_30_combo[raw_30_combo$ploidy == '4X', ]
raw_30_combo_8X <- raw_30_combo[raw_30_combo$ploidy == '8X', ]

gg_raw4X <- ggplot(raw_30_combo_4X, aes(x = coverage, y = d_loglik_portion)) + 
  geom_line(aes(color = factor(ploidy_test), group=ploidy_test)) + 
  labs(title = 'Portion of nQuire delta log-likelihood for 4X test sample', 
    x = 'Coverage cutoff', y = 'Portion of delta log-likelihood')

raw4X_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/nquire_test_4X_raw_q30_delta_loglik_results.pdf', sep = '')

pdf(file = raw4X_file)
gg_raw4X
dev.off()

gg_raw8X <- ggplot(raw_30_combo_8X, aes(x = coverage, y = d_loglik_portion)) +
  geom_line(aes(color = factor(ploidy_test), group=ploidy_test)) +
  labs(title = 'Portion of nQuire delta log-likelihood for 8X test sample', 
    x = 'Coverage cutoff', y = 'Portion of delta log-likelihood')

raw8X_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/nquire_test_8X_raw_q30_delta_loglik_results.pdf', sep = '')

pdf(file = raw8X_file)
gg_raw8X
dev.off()

raw_30_combo_2 <- raw_30_combo
raw_30_combo_2$comp_group <- paste(raw_30_combo_2$ploidy_test, 
  raw_30_combo_2$ploidy, sep = '_')

gg_raw_combo <- ggplot(raw_30_combo_2, aes(x = coverage, 
    y = d_loglik_portion)) +
  geom_line(aes(color = factor(ploidy_test),  lty = factor(ploidy), 
    group=comp_group)) +
  labs(title = paste('Portion of nQuire delta log-likelihood for test',
    ' samples\nUsing CDS SNPs', sep = ''),
    x = 'Coverage cutoff', y = 'Portion of delta log-likelihood')

raw_combo_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/nquire_test_combo_raw_q30_delta_loglik_results.pdf', sep = '')

pdf(file = raw_combo_file)
gg_raw_combo
dev.off()


gg_raw_coverage <- ggplot(raw_30_df, aes(x = coverage, y = nSNPs)) + 
  geom_line(aes(color = factor(ploidy), group = ploidy)) + 
  labs(title = paste('Number CDS SNPs used by nQuire at different',
    '\ncoverage cutoffs', sep = ''),
    x = 'Coverage cutoff', y = 'number CDS SNPs')

raw_coverage_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/nquire_test_combo_raw_q30_coverage.pdf', sep = '')

pdf(raw_coverage_file)
gg_raw_coverage
dev.off()


#############

den_30_dip <- data.frame(file = denoise_30_df$file, 
  coverage = denoise_30_df$coverage,
  d_loglik_portion = denoise_30_df$d_dip_portion, ploidy_test = 'dip',
  ploidy = denoise_30_df$ploidy, stringsAsFactors = F)

den_30_tri <- data.frame(file = denoise_30_df$file, 
  coverage = denoise_30_df$coverage,
  d_loglik_portion = denoise_30_df$d_tri_portion, ploidy_test = 'tri',
  ploidy = denoise_30_df$ploidy, stringsAsFactors = F)

den_30_tet <- data.frame(file = denoise_30_df$file, 
  coverage = denoise_30_df$coverage,
  d_loglik_portion = denoise_30_df$d_tet_portion, ploidy_test = 'tet',
  ploidy = denoise_30_df$ploidy, stringsAsFactors = F)

den_30_combo <- rbind(rbind(den_30_dip, den_30_tri), den_30_tet)

den_30_combo_2 <- den_30_combo
den_30_combo_2$comp_group <- paste(den_30_combo_2$ploidy_test, 
  den_30_combo_2$ploidy, sep = '_')

gg_denoise_combo <- ggplot(den_30_combo_2, aes(x = coverage, 
    y = d_loglik_portion)) +
  geom_line(aes(color = factor(ploidy_test),  lty = factor(ploidy),
    group=comp_group)) +
  labs(title = paste('Portion of nQuire delta log-likelihood for test', 
    ' samples\nUsing "Denoised" CDS SNPs', sep = ''),
    x = 'Coverage cutoff', y = 'Portion of delta log-likelihood')

den_combo_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/nquire_test_combo_denoise_q30_delta_loglik_results.pdf', 
   sep = '')

pdf(file = den_combo_file)
gg_denoise_combo
dev.off()

gg_denoise_coverage <- ggplot(denoise_30_df, aes(x = coverage, y = nSNPs)) +
  geom_line(aes(color = factor(ploidy), group = ploidy)) +
  labs(title = paste('Number "Denoise" CDS SNPs used by nQuire at different',
    '\ncoverage cutoffs', sep = ''),
    x = 'Coverage cutoff', y = 'number CDS SNPs')

den_coverage_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/nquire_test_combo_denoise_q30_coverage.pdf', sep = '')

pdf(den_coverage_file)
gg_denoise_coverage
dev.off()

######################

tot_res_2 <- tot_res

tot_res_2$lib_type <- 'raw'
tot_res_2$lib_type[grep('denoise', tot_res_2$file)] <- 'denoised'

tot_res_3 <- tot_res_2[which(tot_res_2$qual == 30), ]
tot_res_3$comp_group <- paste(tot_res_3$ploidy, tot_res_3$lib_type, sep = '_')

gg_tot_coverage <- ggplot(tot_res_3, aes(x = coverage, y = nSNPs)) + 
  geom_line(aes(color = factor(ploidy), lty = factor(lib_type), 
    group = comp_group)) + 
  labs(title = paste('Number CDS SNPs used by nQuire at different',
    '\ncoverage cutoffs', sep = ''),
    x = 'Coverage cutoff', y = 'number CDS SNPs')

combo_coverage_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/nquire_test_combo_combo_q30_coverage.pdf', sep = '')

pdf(combo_coverage_file)
gg_tot_coverage
dev.off()



```



## Test nQuire - using all SNPs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test

module load python/3.7-anaconda-2019.07
source activate poplar_align

# generate smaller bam file
samtools view -h /global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/ABHM.Pvirgatum.V5.merged.gatk.bam | head -n 100000 | samtools view -bS - > test_100k.bam

/global/homes/g/grabowsp/tools/nQuire/nQuire create -b test_100k.bam -o test_1

/global/homes/g/grabowsp/tools/nQuire/nQuire create -b test_100k.bam -o \
test_2 -x

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_1.bin
# produces a ASCII (printed to command line) histogram

# To look at base distribution for included SNPs
/global/homes/g/grabowsp/tools/nQuire/nQuire view test_1.bin
# for eXtended .bin, need to include original bam to get annotation info
/global/homes/g/grabowsp/tools/nQuire/nQuire view test_2.bin -a \
test_100k.bam | head

# Test the [denoising] option
/global/homes/g/grabowsp/tools/nQuire/nQuire denoise test_1.bin -o \
test_1_denoised
## Before: 6113
## After:  1806 (29.5%)

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_1_denoised.bin

# compare [denoising] to using more stringent mapping
/global/homes/g/grabowsp/tools/nQuire/nQuire create -b test_100k.bam \
-q 30 -o test_3

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_3.bin

/global/homes/g/grabowsp/tools/nQuire/nQuire denoise test_3.bin -o \
test_3_denoised

# Before: 4076
# After:  1326 (32.5%)
# Maintains higher percentage but fewer SNPs than when using min mapping 
#   score of 1

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_3_denoised.bin

# test using higher coverage requirement
/global/homes/g/grabowsp/tools/nQuire/nQuire create -b test_100k.bam \
-c 20 -o test_4

/global/homes/g/grabowsp/tools/nQuire/nQuire denoise test_4.bin -o \
test_4_denoised

# Before: 2576
# After:  860 (33.4%)

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_4_denoised.bin

# Take homes so far:
# 1) Should use `-q 30` to select for better mapping
# 2) Should us `-c 20` or something like that to select higher coverage 
# requirment

/global/homes/g/grabowsp/tools/nQuire/nQuire lrdmodel test_1_denoised.bin
/global/homes/g/grabowsp/tools/nQuire/nQuire modeltest test_1_denoised.bin
/global/homes/g/grabowsp/tools/nQuire/nQuire estmodel test_1_denoised.bin
/global/homes/g/grabowsp/tools/nQuire/nQuire histotest test_1_denoised.bin

# Full BAM, default parameters

/global/homes/g/grabowsp/tools/nQuire/nQuire create -b /global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/ABHM.Pvirgatum.V5.merged.gatk.bam -o test_full_1

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_1.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_1_lrd.sh

## This takes a long time; need to re-run again!
# /global/homes/g/grabowsp/tools/nQuire/nQuire denoise test_full_1.bin -o \
test_full_1_denoised

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_1_denoise.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_1_denoised.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_1_denoised_lrd.sh

### Full BAM, stringent mapping ###
/global/homes/g/grabowsp/tools/nQuire/nQuire create -b \
/global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/ABHM.Pvirgatum.V5.merged.gatk.bam \
-q 30 -c 40 -o test_full_2

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_2.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_2_lrd.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire denoise test_full_2.bin -o \
test_full_2_denoised

# Before: 961307
# After:  236133 (24.6%)

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_2_denoised.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_2_denoised_lrd.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire lrdmodel -t 10 \
test_full_2_denoised.bin

### FULL BAM, medium mapping ###
# need to test this
#/global/homes/g/grabowsp/tools/nQuire/nQuire create -b \
#/global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/ABHM.Pvirgatum.V5.merged.gatk.bam \
#-q 20 -c 20 -o test_full_3

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_3_create.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_3.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_3_lrd.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_3_denoise.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_3_denoised.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_3_denoised_lrd.sh

## FULL BAM, high coverage
#/global/homes/g/grabowsp/tools/nQuire/nQuire create -b \
#/global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/ABHM.Pvirgatum.V5.merged.gatk.bam \
#-q 30 -c 50 -o test_full_4

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_4_create.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_4.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_4_lrd.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_4_denoise.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_4_denoised.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_4_denoised_lrd.sh

# FULL BAM, 40 coverage, 20 quality
## -q 20 -c 40 -o test_full_5
sbatch test_full_5_create.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_5.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_5_lrd.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_5_denoise.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_5_denoised.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_5_denoised_lrd.sh

# FULL BAM, 20 coverage, 30 quality
## -q 30 -c 20 -o test_full_6
sbatch test_full_6_create.sh

/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_6.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_6_lrd.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_6_denoise.sh

## NEED TO RUN NEXT 2 PARTS
/global/homes/g/grabowsp/tools/nQuire/nQuire histo test_full_6_denoised.bin

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_full_6_denoised_lrd.sh

#### Concatenate results
cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
for i in *lrd_results.txt; do tail -n 1 $i >> combo_test_lrd_results.txt; done

cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
for i in {1..6}; do echo test_full_$i'_denoised.bin'; echo test_full_$i'.bin'; done

for i in {1..6};
do /global/homes/g/grabowsp/tools/nQuire/nQuire view \
test_full_$i'_denoised.bin' | wc -l >> combo_test_nSNPs.txt;
/global/homes/g/grabowsp/tools/nQuire/nQuire view \
test_full_$i'.bin' | wc -l >> combo_test_nSNPs.txt;
done

#### Analyze results

lrd_res_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'nquire_test/combo_test_lrd_results.txt', sep = '')
lrd_res <- read.table(lrd_res_file, header = F, sep = '\t',
  stringsAsFactors = F)

lrd_head_0 <- system(paste('head -1 ',
  '/global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/',
  'test_full_1_lrd_results.txt', sep = ''), intern = T)

lrd_head <- unlist(strsplit(lrd_head_0, split = '\t'))

colnames(lrd_res) <- lrd_head

nquire_nsnps_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'nquire_test/combo_test_nSNPs.txt', sep = '')
nquire_nsnps <- read.table(nquire_nsnps_file, header = F, stringsAsFactors = F)

lrd_res$nSNPs <- nquire_nsnps[,1]

lrd_res$coverage <- c(10,10,40,40,20,20,50,50,40,40,20,20)
lrd_res$qual <- c(1,1,30,30,20,20,30,30,20,20,30,30)

lrd_res$d_dip_portion <- (lrd_res$d_dip / 
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tri_portion <- (lrd_res$d_tri / 
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

lrd_res$d_tet_portion <- (lrd_res$d_tet / 
  (lrd_res$d_dip + lrd_res$d_tri + lrd_res$d_tet))

combo_processed_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'nquire_test/combo_test_lrd_results_processed.txt', sep = '')
write.table(lrd_res, file = combo_processed_file, quote = F, sep = '\t',
  row.names = F, col.names = T)

```

## Generate .bin files for each BAM
* will use -q 30 and -c 20 values 
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire
cut Pvirgatum_V5_bams_arcived.txt -d "." -f 1 > Pvirgatum_V5_libs.txt



for SGBAM in `cat Pvirgatum_V5_libs.txt`;
  do sed 's/TESTLIB/'"$SGBAM"'/g' TESTLIB_nquire_create.sh > \
$SGBAM'_nquire_create.sh';
done

$SGBAM'_test'; done 


cut 
```

