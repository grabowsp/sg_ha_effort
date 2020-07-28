# Explanation of Lists of 4X and 8X Samples
* These library lists are based on best estimates of ploidy using nQuire
and MNP counts done in spring 2020.
* Potential 6X samples have been asigned to the 8X list
* On Cori
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt`
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

ploidy_meta_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v2.0.txt'

ploidy_meta <- read.table(ploidy_meta_file, header = T, stringsAsFactors = F,
  sep = '\t')

tet_libs <- ploidy_meta$lib[which(ploidy_meta$total_ploidy == '4X')]
oct_libs <- ploidy_meta$lib[which(ploidy_meta$total_ploidy == '8X')]

tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'

write.table(tet_libs, file = tet_lib_file, quote = F, row.names = F,
  col.names = F)
write.table(oct_libs, file = oct_lib_file, quote = F, row.names = F,
  col.names = F)
```


