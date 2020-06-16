# Explanations of the samples used for analyses

## 901 Reseq samples
* These are the resequencing samples used for the analyses in Genome Paper
  * Were resequenced
  * Duplicate libraries removed
  * "Bad Samples" removed
* On Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/samp901_lib_names.txt`

## 887 Hi Coverage Samples
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

## Samples to be used for Natural/Geographic Analysis
* 4X samples used for geographic analysis in genome paper and 8X samples \
that are "Natural Collections" (and "Cultivars" for expanded set)
* On Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_772_lib_names.txt`
    * included 8X "Natural Collections"
    * old file: `/global/cscratch1/sd/grabowsp/sg_ploidy/geo826_lib_names.txt`
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_expand_785_lib_names.txt`
    * included both 8X "Natural Collections" and 8X "Cultivars"
    * old file: `/global/cscratch1/sd/grabowsp/sg_ploidy/geo_expand_839_lib_names.txt`
### Code used for generating sample sets
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/genotype.metadata.May2020.rds'

meta <- readRDS(meta_file)
# use this to chose "Climate" samples for 4X samples

# "good" libraries
good_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/hi_cov_887_lib_names.txt'
good_libs <- read.table(good_lib_file, header = F, stringsAsFactors = F)

good_inds <- which(meta$LIBRARY %in% good_libs[,1])

# 4X samples chosen for climate analysis
tet_clim_libs <- which(meta$LIB_BIOCLIM == 'Y')

# presumptive 8X samples that come from natural collections
other_clim_libs <- intersect(
  which(meta$COLLECTION_TYPE == 'Natural Collection'),
  which(meta$LIB_CLIMATE == 'NA'))

good_clim_libs <- intersect(good_inds, union(tet_clim_libs, other_clim_libs))

#out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/geo826_lib_names.txt'

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_772_lib_names.txt'

write.table(meta$LIBRARY[good_clim_libs], file = out_file, quote = F,
  sep = '\t', row.names = F, col.names = F)

# presumptive 8X samples that come from Cultivars
other_cultivar_libs <- intersect(
  which(meta$COLLECTION_TYPE == 'Cultivar'),
  which(meta$LIB_CLIMATE == 'NA'))

expanded_clim_libs <- intersect(good_inds, union(tet_clim_libs,
  union(other_clim_libs, other_cultivar_libs)))

expand_out_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'geo_v2_expand_785_lib_names.txt', sep = '')

write.table(meta$LIBRARY[expanded_clim_libs], file = expand_out_file, quote = F,
  sep = '\t', row.names = F, col.names = F)
```

## Lists of 4X and 8X Samples
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
