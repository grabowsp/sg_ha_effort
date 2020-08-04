#module load python/3.7-anaconda-2019.07
#source activate R_analysis

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

### LOAD PACKAGES ###
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.basename <- dirname(script.name)

#polyploid_function_file <- file.path(script.basename, 'polyploid_functions.r')
#polyploid_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/polyploid_functions.r'
#source(polyploid_function_file)

gen_function_file <- file.path(script.basename, 'general_functions.r')
#gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
source(gen_function_file)

### LOAD INPUTS ###
vcf_in <- args[1]
#vcf_in <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_00'
vcf_1 <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

### SET OUTPUT ###
out_dir <- args[2]
#out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs'
out_dir <- add_slash(out_dir)

out_file_short <- gsub('vcf_', 'vcf_dipcode_', basename(vcf_in))
out_file <- paste(out_dir, out_file_short, sep = '')

### SET VARIABLES ###
# tetrasomic genotypes that need to be converted; note that 2/2 stays the same
tet_genos <- c('4/0', '3/1', '1/3', '0/4')
dip_genos <- c('2/0', '2/1', '1/2', '0/2')

##################

for(i in seq(length(tet_genos))){
  vcf_1[vcf_1 == tet_genos[i]] <- dip_genos[i]
}

write.table(vcf_1, file = out_file, quote = F, sep = '\t', row.names = F,
  col.names = F)

quit(save = 'no')

