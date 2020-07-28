## Samples from Original Sample Set used for Natural/Geographic Analysis
* 4X samples used for geographic analysis in genome paper and 8X samples \
that are "Natural Collections" 
  * 8X "Cultivars" for included in the expanded_geo set
* On Cori:
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_772_lib_names.txt`
    * included 8X "Natural Collections"
    * old file: `/global/cscratch1/sd/grabowsp/sg_ploidy/geo826_lib_names.txt`
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

out_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/geo_v2_772_lib_names.txt'

write.table(meta$LIBRARY[good_clim_libs], file = out_file, quote = F,
  sep = '\t', row.names = F, col.names = F)

```
