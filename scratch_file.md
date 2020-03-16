```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap

cp /home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_V5_1035g_Chr01K.unfiltered_reheaded.sorted.vcf.gz .

gunzip -c Pvirgatum_V5_1035g_Chr01K.unfiltered_reheaded.sorted.vcf.gz | \
head -100000 > test.vcf

# filter vcf based on missing data and maf

# missing data < 20%
# maf > 
# only select libraries in full analysis

vcftools --vcf test.vcf --out test_filt --maf 0.01 --max-missing 0.8 \
--keep samp901_PlantIDs.txt --indv '9001-3 BN389-69S' --bed sg_v5_genes.bed \
--recode --recode-INFO-all


vcftools --gzvcf Pvirgatum_V5_1035g_Chr01K.unfiltered_reheaded.sorted.vcf.gz \
--out Chr01K_filt --maf 0.01 --max-missing 0.8 \
--keep samp901_PlantIDs.txt --indv '9001-3 BN389-69S' --bed sg_v5_genes.bed \
--recode --recode-INFO-all




head -5 Chr01K_filt.recode.vcf | tail -n 1 | \
cut --complement -f 1-9 > 


/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs/filtered_vcf_samp_names.txt

#####

test_vcf_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/test_filt.recode.vcf'

test_vcf_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs/Chr01K_filt_split_00'

test_vcf <- read.table(test_vcf_file, header = F, stringsAsFactors = F)

test_counts <- apply(test_vcf[, c(10:ncol(test_vcf))], 1, function(x)
  unlist(lapply(strsplit(x, split = ':'), function(y) y[[2]]))
)


test_counts_2 <- apply(test_vcf[, c(10:ncol(test_vcf))], 2, function(x) 
  unlist(lapply(strsplit(x, split = ':'), function(y) y[[2]]))
)

test_counts_3 <- apply(test_vcf[c(1:50000), c(10:ncol(test_vcf))], 2, function(x)
  unlist(lapply(strsplit(x, split = ':'), function(y) y[[2]]))
)

test_ratios_3 <- apply(test_counts_3, 2, function(x)
  unlist(lapply(strsplit(x, split = ','), function(y)
  as.numeric(y[1])/sum(as.numeric(y))))
)

test_ratio_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/filtered_vcfs/Chr01K_filt_split_00_alleleratios.rds'

test_ratio <- readRDS(test_ratio_file)

pseudo_hap_genos <- apply(test_ratio[c(1:100), c(6:ncol(test_ratio))], 2, function(x)
  sapply(x, function(y) rbinom(n = 1, size = 1, prob = y))
)

snp_info <- test_ratio[c(1:100), c(1:5)]

pseudohap_tot <- data.frame(snp_info, pseudo_hap_genos, stringsAsFactors = F)

ratios_in <- test_ratio_file




test_ratios <- apply(test_counts, 2, function(x)
  unlist(lapply(strsplit(x, split = ','), function(y)
  as.numeric(y[1])/sum(as.numeric(y))))
)

test_ratios_2 <- apply(test_counts_2, 2, function(x)
  unlist(lapply(strsplit(x, split = ','), function(y)
  as.numeric(y[1])/sum(as.numeric(y))))
)


pseudohap_genos_3 <- apply(test_ratios_3, 2, function(x)
  sapply(x, function(y) rbinom(n = 1, size = 1, prob = y))
)

pseudohap_genos_1 <- apply(test_ratios, 2, function(x)
  sapply(x, function(y) rbinom(n = 1, size = 1, prob = y))
)

sapply(test_ratios[,1], function(x) rbinom(n = 1, size = 1, prob = x))

rbinom(n = 1, size = 1, prob = test_ratios[1,1])


unlist(lapply(strsplit(test_counts[,1], split = ','), function(y) 
  as.numeric(y[1])/sum(as.numeric(y))))



unlist(lapply(strsplit(test_vcf[, 10], split = ':'), function(x) x[[2]]))


## Generate list of samples to keep on Cori

# module load python/3.7-anaconda-2019.07
# source activate R_analysis


ploidy_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v2.0.txt'
ploidy_res <- read.table(ploidy_res_file, header = T, stringsAsFactors = F,
  sep = '\t')

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/samp901_lib_names.txt'

write.table(ploidy_res$lib, file = out_file, quote = F, sep = '\t', 
  row.names = F, col.names = F)

#######

scp /global/cscratch1/sd/grabowsp/sg_ploidy/samp901_lib_names.txt \
grabowsk@pants.hagsc.org:/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap

# Extract sample names

lib_name_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/samp901_lib_names.txt'
lib_keep_names <- read.table(lib_name_file, header = F, stringsAsFactors = F)

meta_samp_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/PVDIV_Master_Metadata_File_9-3-2019_Samps.txt'

meta_samps <- read.table(meta_samp_file, header = T, stringsAsFactors = F, 
  sep = '\t')

samps_to_keep <- c()

for(i in seq(nrow(lib_keep_names))){
  tmp_meta_ind <- which(meta_samps$LIBRARY == lib_keep_names[i, 1])
  samps_to_keep <- c(samps_to_keep, meta_samps$PLANT_ID[tmp_meta_ind])
}

out_df <- data.frame(plant_id = samps_to_keep, stringsAsFactors = F)

out_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/samp901_PlantIDs.txt'

write.table(out_df, file = out_file, quote = F, sep = '\t', row.names = F,
  col.names = F)

```
