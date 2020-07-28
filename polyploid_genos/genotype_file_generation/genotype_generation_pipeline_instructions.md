# Instructions and descriptions of steps for generating genotype files
#  for different analyses of 4X/8X samples

## Background
* Current "raw" genotypes are tetrasomic genotypes for all samples generated
by Sujan
  * Good for 8X, not correct for 4X
  * 4X genotypes are adjusted when `genlight` objects are generated
* Use CDS SNPs because they seem to be highest quality based on previous 
experience
  * There seem to be fewer mapping issues with CDS SNPs
* 

## Steps
* Generate chromosome CDS VCFs for sample set
* Split chromsome CDS VCFs into 100k-line subfiles
* Generate header file to be used with subfiles
* Generate `genlight` object for each chromosome
* Make genome-wide, subsampled `genlight` object

## Files and Resource Used
### Needed Resources
#### conda environments
* `gen_bioinformatics` environment
  * includes VCFtools which is needed for filtering VCFs
  * `/global/homes/g/grabowsp/.conda/envs/gen_bioinformatics`
* `r_adegenet_env` environment
  * includeds R and adegenet package
  * `/global/homes/g/grabowsp/.conda/envs/r_adegenet_env`
#### Scripts from github
* `make_Chr_genlight_objs.r`
  * used for making chromosome 'genlight' genotype objects
  * `/PATH/TO/REPOSITORY/sg_ha_effort/polyploid_genos/adegenet_analysis/make_Chr_genlight_objs.r`
* `subsample_genlight.r`
  * used for making subsampled genome-wide 'genlight' object
  * `/PATH/TO/REPOSITORY/sg_ha_effort/polyploid_genos/adegenet_analysis/subsample_genlight.r`
#### R-Function files from github
* `general_functions.r`
  * contains functions used in many types of analyses
  * `/PATH/TO/REPOSITORY/sg_ha_effort/polyploid_genos/general_functions.r`
* `adegenet_functions.r`
  * contains functions specific to analyses using adegenet package
  * `/PATH/TO/REPOSITORY//home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/adegenet_analysis/adegenet_functions.r`

### Needed Files
#### Starting-point VCFs
* CDS SNP VCFs for all libraries for each chromosome
  * Found in this directory:
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs`
  * Generated using `generate_CDS_general_vcfs.md` script
#### Library names to include
* Make a file with the library names of the samples to be included
#### Library names for different ploidies
* Tetraploid libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt`
* Octoploid libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt`

## Generate chromosome CDS VCFs for sample set
* I generally use maximum 20% missing data (--max-missing 0.8)
* Requires file of library names to keep in a single column
* Easiest way to to submit jobs to cluster
### Submit job
```
cd /PATH/TO/DATA/BEING/USED

sbatch generate_Chrset_sampleset_vcf.sh
```
### Example of submit script
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

# Directory where job info (and new VCFs) will be saved
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

# Directory with base CDS VCFs
DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/

# File of Library Names to keep
SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_772_lib_names.txt

# Directory where new VCFs will be saved
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

# sample set name to be used in fame files
SAMPSET=geosamps

for CHROM in 01K 01N 02K 02N 03K 03N 04K 04N 05K 05N;
  do
  vcftools --gzvcf $DATA_DIR'Chr'$CHROM'.polyploid.CDS.vcf.gz' \
  --max-missing 0.8 \
  --stdout --keep $SAMP_FILE --recode --recode-INFO-all | \
  gzip -c > \
  $OUT_DIR'Chr'$CHROM'.polyploid.CDS.'$SAMPSET'.vcf.gz';
  done
```

## Split chromosome CDS VCFs into 100k-line subfiles
* For some analyses, is faster and practical to work on subsets of the genotypes
  * I chose 100k lines because not too many files are generated (~200 for 
first round of analysis) and files are processed reasonably fast in R
* Is reasonably fast to do this in an interactive session
  * can break it up into loops in separate sessions
### Example Code for splitting files
```
# Go to directory with VCFs to split
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

SAMPSET=geosamps

for CN in {01..09};
do
gunzip -kc Chr$CN'K.polyploid.CDS.'$SAMPSET'.vcf.gz' | \
split -l 100000 -d - Chr$CN'K.polyploid.CDS.'$SAMPSET'.vcf_';
done

for CN in {01..09};
do
gunzip -kc Chr$CN'N.polyploid.CDS.'$SAMPSET'.vcf.gz' | \
split -l 100000 -d - Chr$CN'N.polyploid.CDS.'$SAMPSET'.vcf_';
done
```

## Generate header file to be used with subfiles
* 100k-line VCF subfiles require header when loading them into R for analysis
* 1 header file can work for all the sub-files
  * VCF header information begins with # so isn't loaded into R, so can use
same approach for all subfiles, even those at beginning of chromsome that
have the standard VCF header
* Use the very first sub-file from any chromosome to make the generic header
file
  * Example: `Chr01K...vcf_00`
* The number of lines to be used with `head` depends on the amount of
information in the header section of the VCF
  * For current formatting, there are only 5 lines before the genotypes begin
### Example of generating file
```
# go to directory with genotype files
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

SAMPSET=geosamps

head -5 Chr01K.polyploid.CDS.$SAMPSET'.vcf_00' | tail -1 > \
CDS.$SAMPSET'.vcf.header.txt'
```

## Generate `genlight` object for each chromosome
* `Genlight` object is genotypes formated for `adegenet` R package
* Uses `make_Chr_genlight_objs.r` script to generate `genlight` object for
each chromosome
  * hardcoded to include 2 files that have the library names of 8X and 4X samples
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt`
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt`
* Have to decide what MAF cutoff to use
  * I typically aim for a MAF cutoff where there are ~5 copies of the allele
at the MAF, but using higher MAF is certainly defensible
### Example of submitting job
```
# go to directory with genotypes and submit shell script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

sbatch gen_geo_genlight.sh
```
### Example of submit script
* From `gen_geo_genlight.sh`
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/r_adegenet_env

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

cd $DATA_DIR

HEADER_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SAMPSET=geosamps

MAF_CUT=0.002

for CHR_N in {01..09};
  do
  for CHR_T in K N;
    do
    TEST_CHR=$CHR_N$CHR_T
    SEARCH_STRING=Chr$TEST_CHR'.polyploid.CDS.'$SAMPSET'.vcf_'
    Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/make_Chr_genlight_objs.r \
    $DATA_DIR $SEARCH_STRING'*' $HEADER_FILE $MAF_CUT;
    done;
  done;
```

## Generate genome-wide subsampled `genlight` object
* Subsample SNPs from all chromosomes to get desired number of genome-wide SNPs
### Chose subsampling level
* Decide on goal number of SNPs, calculate total number of SNPs in all 
chromosome 'genlight' objects, and calculate the percentage needed for goal
number of SNPs
```
# load conda environment
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# R
library(adegenet)

# adjust these variables
###
# directory where genolight objects are
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/'

# goal number of subsampled SNPs
goal_nsnps <- 5e5
###

file_vec <- system(paste('ls ', data_dir, 'Chr*genlight.rds', sep = ''),
  intern = T)

n_snp_vec <- c()
for(i in file_vec){
  print(i)
  tmp_in <- readRDS(i)
  tmp_nsnps <- nLoc(tmp_in)
  n_snp_vec <- c(n_snp_vec, tmp_nsnps)
}

tot_snps <- sum(n_snp_vec)
goal_sub <- goal_nsnps/tot_snps
print(goal_sub)
# [1] 0.05910562
```
### Generate object
* uses `subsample_genlight.r` script
* Takes a few minutes, but can be run in interactive session
#### Example
```
# load conda environment
module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

# go to directory where want to run script
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps

# directory with chromosome genlight objects
DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

# how sample set is used in file names
SAMPSET=geosamps

# percentage of total SNPs to select
PER_SUBSAMP=0.07

# FILE_SUB may need to be hardcoded (without the $SAMPSET variable) - try that
#   if there is an issue
FILE_SUB=Chr*$SAMPSET'.genlight.rds'
OUT_NAME=combo.sub.polyploid.CDS.$SAMPSET'.genlight.rds'

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/adegenet_analysis/subsample_genlight.r \
$DATA_DIR '*'$FILE_SUB $OUT_NAME $PER_SUBSAMP
```


