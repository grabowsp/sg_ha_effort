# Steps for combining simulated 8X and geosamps VCFs

## Proposed Steps
* Convert subfiles to "coded" disomic genotypes
* Concatenate subfiles
* Sort, bgzip and tabix subfiles
* Combine sim8X and geosamps
* Generate new subfiles
* Convert new subfiles to tetrasomic genotypes

## Coded disomic genotypes
* TET genotype : DIP code
* 4/0 : 2/0
* 3/1 : 2/1
* 2/2 : 2/2
* 1/3 : 1/2
* 0/4 : 0/2

## Test
* Convert geosamps Chr01K subfile(s) to coded disomic genotypes
* Concatenate subfiles
* sort, bgzip, and tabix
* Convert sim8X Chr01K subfiles to coded disomic genotypes
* concatenate subfiles
* sort, bgzip, and tabix
* merge the geosamps and sim8X Chr01K files
* generate merged subfiles
* convert merged subfiles back to tetrasomic genotypes

## Convert subfile to coded disomic genotypes
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(gen_function_file)


in_vcf_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00'

in_vcf <- read.table(in_vcf_file, header = F, stringsAsFactors = F, sep = '\t')

out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/convert_vcfs'

out_dir <- add_slash(out_dir)

out_file_tmp <- rev(unlist(strsplit(in_vcf_file, split = '/')))[1]

out_file_short <- gsub('vcf_', 'vcf_dipcode_', out_file_tmp)

out_file <- paste(out_dir, out_file_short, sep = '')

####

tet_genos <- c('4/0', '3/1', '1/3', '0/4')
dip_genos <- c('2/0', '2/1', '1/2', '0/2')

for(i in seq(length(tet_genos))){
  in_vcf[in_vcf == tet_genos[i]] <- dip_genos[i]
}





```
