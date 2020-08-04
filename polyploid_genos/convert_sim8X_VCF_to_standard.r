# Script to convert allele-count simulated 8X VCFs to the same format as
#  standard VCFs so VCFs can be merged

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/R_analysis

data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/sim8X_vcfs/'

in_file_short <- 'Chr01K.polyploid.CDS.geosamps.geo_samp_sim8X_AltDosage.vcf_00'

in_file <- paste(data_dir, in_file_short, sep = '')


out_file <- gsub('AltDosage', 'standard', in_file)


vcf_a <- read.table(in_file, header = T, sep = '\t', stringsAsFactors = F)

vcf_b <- vcf_a

vcf_b[is.na(vcf_a)] <- './.'

vcf_b[vcf_a == 0] <- '4/0'

vcf_b[vcf_a == 1] <- '3/1'

vcf_b[vcf_a == 2] <- '2/2'

vcf_b[vcf_a == 3] <- '1/3'

vcf_b[vcf_a == 4] <- '0/4'

write.table(vcf_b, file = out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)


