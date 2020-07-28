# Explanation of 887 Hi Coverage Samples
* Start with 901 samples used in the Genome Paper
* Remove 14 samples with sequencing output below 1e10
* These samples are tricky to determine ploidy and call genotypes
* On Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/hi_cov_887_lib_names.txt`
### Script used for making this sample list
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

ploidy_info_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'ploidy_calling/sg_ploidy_results_v3.0.txt', sep = '')
ploidy_info <- read.table(ploidy_info_file, header = T, sep = '\t',
  stringsAsFactors = F)

low_cov_inds <- which(ploidy_info$seq_cov < 1e10)

out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'hi_cov_887_lib_names.txt', sep = '')

write.table(ploidy_info$lib[-low_cov_inds], file = out_file, quote = F,
  sep = '\t', row.names = F, col.names = F)
```

