# Trying to use adegenet

## Notes for using adegenet up until now
### * I am able to make and combine genlight objects with the genotypes using
###     a mix of disomic and tetrasomic genotypes
### * However, there is an error whenever I try to run glPca() with more than
###    around 210 samples -  i get the following error:
###    " Error in eigen(allProd, symmetric = TRUE, only.values = FALSE) : 
###       error code 1 from Lapack routine 'dsyevr' "
###   * I get the error regardless of the number of SNPs I use - I tested 10
###       through 10K SNPs
###   * It is not because of ploidy differences - it happens regardless of 
###       which samples I select
###   * It is not because of identity between samples - again, it happens
###       regardless of the samples I select
### * The error is an issue withthe package that underlies eigen(), and the 
###     manual says that to solve the problem requires detailed analysis
###     of the FORTRAN code, which I can't do

module load python/3.7-anaconda-2019.07
source activate r_adegenet_env

library(adegenet)
library(parallel)

# on HA: gen_function_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/general_functions.r'
# on Cori
gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(gen_function_file)

adeg_function_file <- 

### IMPORT DATA

#vcf_in_0 <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00'
#vcf_0 <- read.table(vcf_in_0, header = F, stringsAsFactors = F, sep = '\t')

#data_dir <- ags[1]
data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/'

data_dir <- add_slash(data_dir)

#vcf_search_string <- args[2]
vcf_search_string <- 'Chr01K.polyploid.CDS.geosamps.vcf_*'

vcf_search_command <- paste('ls ', data_dir, vcf_search_string, sep = '')

vcf_files <- system(vcf_search_command, inter = T)

#vcf_files <- system('ls /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_*', inter = T)

vcf_1 <- read.table(vcf_files[1], header = F, stringsAsFactors = F, 
    sep = '\t')

#vcf_header_file <- args[3]
vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt'

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf_1) <- vcf_header

tet_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs <- intersect(tet_libs_0, vcf_header)

oct_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs <- intersect(oct_libs_0, vcf_header)

# maf_cut <- args[4]
maf_cut <- 0.002

#######
test_po_0 <- gen_gl_preobj(vcf = vcf_0, oct_libs = oct_libs, 
  tet_libs = tet_libs, maf_cut = maf_cut)

test_po_1 <- gen_gl_preobj(vcf = vcf_1, oct_libs = oct_libs,
  tet_libs = tet_libs, maf_cut = maf_cut)

data_gl_0 <- new('genlight', gen = t(as.matrix(test_po_0[['geno_mat']])),
  ploidy = test_po_0[['ploidy']],
  loc.names = test_po_0[['loc.names']],
  chromosome = test_po_0[['chromosome']],
  position = test_po_0[['position']],
  loc.all = test_po_0[['loc.all']]
)

data_gl_1 <- new('genlight', gen = t(as.matrix(test_po_1[['geno_mat']])),
  ploidy = test_po_1[['ploidy']],
  loc.names = test_po_1[['loc.names']],
  chromosome = test_po_1[['chromosome']],
  position = test_po_1[['position']],
  loc.all = test_po_1[['loc.all']]
)

data_gl_tot <- cbind(data_gl_0, data_gl_1)

# on HA pants
vcf_test <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/test_vcfs/Chr01K.polyploid.CDS.geosamps.vcf_00'
vcf_1 <- read.table(vcf_test, header = F, stringsAsFactors = F,
    sep = '\t')

vcf_header_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/test_vcfs/CDS.geosamps.vcf.header.txt'

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

tet_lib_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/test_vcfs/tetraploid_lib_names_May2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs <- intersect(tet_libs_0, vcf_header)

oct_lib_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/test_vcfs/octoploid_lib_names_May2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs <- intersect(oct_libs_0, vcf_header)

maf_cut <- 0.002


###

colnames(vcf_1) <- vcf_header

preobj_1 <- gen_gl_preobj(vcf = vcf_1, oct_libs = oct_libs, 
  tet_libs = tet_libs, maf_cut = maf_cut)

gl_tot <- new('genlight', gen = t(as.matrix(preobj_1[['geno_mat']])),
  ploidy = preobj_1[['ploidy']],
  loc.names = preobj_1[['loc.names']],
  chromosome = preobj_1[['chromosome']],
  position = preobj_1[['position']],
  loc.all = preobj_1[['loc.all']]
)

#test_out <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.genlight.rds'

#saveRDS(gl_tot, test_out)

if(length(vcf_files) > 1){
  for(i in c(2:length(vcf_files))){
    vcf_tmp <- read.table(vcf_files[i], header = T, stringsAsFactors = F,
      sep = '\t')
    colnames(vcf_tmp) <- vcf_header
    preobj_tmp <- gen_gl_preobj(vcf = vcf_tmp, oct_libs = oct_libs, 
      tet_libs = tet_libs, maf_cut = maf_cut)
    gl_tmp <- new('genlight', gen = t(as.matrix(preobj_tmp[['geno_mat']])),
      ploidy = preobj_tmp[['ploidy']],
      loc.names = preobj_tmp[['loc.names']],
      chromosome = preobj_tmp[['chromosome']],
      position = preobj_tmp[['position']],
      loc.all = preobj_tmp[['loc.all']]
    )
    gl_tot <- cbind(gl_tot, gl_tmp)
    print(i)
  }
}

temp_sub <- gl_tot[c(1:500) , 
  sort(sample(seq(nLoc(gl_tot)), size = 1000))]

#temp_sub <- gl_tot[c(1:250), ]

pca_test <- glPca(temp_sub, nf = 30, loadings = F, alleleAsUnit = F,
  useC = F)




test_vec <- paste(preobj_1[['geno_mat']][ ,204], sep = '', collapse = '')

result_vec <- c()
for(i in seq(203)){
  tmp_vec <- paste(preobj_1[['geno_mat']][ ,i], sep = '', collapse = '')
  tmp_res <- test_vec == tmp_vec
  result_vec <- c(result_vec, tmp_res)
}



