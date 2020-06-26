# Process for evaluating r2 in the "Natural" Samples

## Generate Sample and Population dataframe
* R script for generating tables
  * `~/geo_r2_pop_table.r`
* Table of all subpopulation libraries
  `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt`
  * 352 Libraries
  * 10 subpopulations
* Table of 4X subpopulation libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/4X_subpop_libs_for_r2.txt`
  * 306 Libraries
  * 8 subpopulations
* Table of 4X subpopulations with 25+ libraries
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/large4X_subpop_libs_for_r2.txt`
  * 271 Libraries
  * 6 subpopulations

## Run r^2 analysis
### Test for 1 sub-vcf
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

VCF_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00

HEADER_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SAMP_POP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt

VCF_TYPE=allele_count

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

FILE_PRE=`basename $VCF_IN`'_geo_subgroup'

MAX_DIST=10000

MAF=0.1

Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/r2_analysis/calc_snp_r2_v2.r $VCF_IN $HEADER_IN $SAMP_POP_FILE $VCF_TYPE $OUT_DIR \
$FILE_PRE $MAX_DIST $MAF

```
### Run interactively on all subfiles
* was taking too long so stopped after finishing Chr03K and Chr03N
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

# Variables for bash loop

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/

LS_STRING=$DATA_DIR*vcf_*

# Variables for R script

HEADER_IN=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt

SAMP_POP_FILE=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/subpop_libs_for_r2.txt

VCF_TYPE=allele_count

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

MAX_DIST=10000

MAF=0.1

for VCF_IN in `ls $LS_STRING`;
  do
  echo $VCF_IN;
  FILE_PRE=`basename $VCF_IN`'_geo_subgroup';
  Rscript /global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/r2_analysis/calc_snp_r2_v2.r \
  $VCF_IN $HEADER_IN $SAMP_POP_FILE $VCF_TYPE $OUT_DIR $FILE_PRE \
  $MAX_DIST $MAF;
  done;
```
### Submit jobs on chromosomes
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results

sbatch calc_r2_Chr04_Chr06.sh
sbatch calc_r2_Chr07_Chr09.sh
```
### Consolidate Results
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

chr_vec <- paste(rep('Chr0', times = 18), rep(seq(9), each = 2), 
  rep(c('K', 'N'), times = 9), sep = '')

data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/'

# add_slash

tot_r2_ls <- list()

for(CHR_NAME in chr_vec){
  print(CHR_NAME)
  ls_com <- paste('ls ', data_dir, CHR_NAME, '*geo_subgroup_r2.rds', sep = '')
  res_files <- system(ls_com, intern = T)
  tot_chr_r2 <- readRDS(res_files[1])
  if(length(res_files > 1)){
    for(RF in c(2:length(res_files))){
      tmp_r2 <- readRDS(res_files[RF])
      for(i in seq(length(tot_chr_r2))){
        tot_chr_r2[[i]] <- rbind(tot_chr_r2[[i]], tmp_r2[[i]])
      }
    }
  }
  tot_r2_ls[[CHR_NAME]] <- tot_chr_r2
}

tot_r2_out_file <- paste(data_dir, 'Chromosome.geo_subgroup.r2.rds', sep = '')

saveRDS(tot_r2_ls, file = tot_r2_out_file)

combo_r2 <- tot_r2_ls[[1]]

for(j in c(2:length(tot_r2_ls))){
  print(j)
  for(k in seq(length(combo_r2))){
    combo_r2[[k]] <- rbind(combo_r2[[k]], tot_r2_ls[[j]][[k]])
  }
}

combo_r2_out_file <- paste(data_dir, 'combined.geo_subgroup.r2.rds', sep = '')
saveRDS(combo_r2, file = combo_r2_out_file)
```
######
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

chr_r2_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results/Chromosome.geo_subgroup.r2.rds'

general_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(general_function_file)



r2_res <- readRDS(chr_r2_file)
chr_vec <- names(r2_res)

test <- r2_res[1:2]

r2_windows <- lapply(r2_res, function(y) 
  lapply(y, function(x) 
    generate_dist_window_df(dist_vec = x$comp_dist, value_vec = x$r2, 
      window_size = 10)
  ))

perHighr2_windows <- lapply(r2_res, function(y)
  lapply(y, function(x)
    generate_perHighVal_window_df(dist_vec = x$comp_dist, value_vec = x$r2,
      window_size = 10, max_val = 0.9)
))

```
### Generate r^ window value for windows
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/r2_results


```




# NEED to ADD the HIGH_VALUE/LD function to the General Funciton File, then
# run that

test <- lapply(tot_r2_ls[[1]], function(x) generate_dist_window_df(dist_vec = x$comp_dist, value_vec = x$r2, window_size = 10))




file_short <- 'Chr07N.polyploid.CDS.geosamps.vcf_06_geo_subgroup_r2.rds'

in_file <- paste(data_dir, file_short, sep = '')

r2_data <- readRDS(in_file)

tot_r2 <- r2_data

file_short_2 <- 'Chr07N.polyploid.CDS.geosamps.vcf_07_geo_subgroup_r2.rds'

in_file_2 <- paste(data_dir, file_short_2, sep = '')

r2_data <- readRDS(in_file_2)

for(i in seq(length(tot_r2))){
  tot_r2[[i]] <- rbind(tot_r2[[i]], r2_data[[i]])
}


``` 

