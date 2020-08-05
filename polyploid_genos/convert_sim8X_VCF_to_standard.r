# Script to convert allele-count simulated 8X VCFs to the same format as
#  standard VCFs so VCFs can be merged

#module load python/3.7-anaconda-2019.07
#source activate /global/homes/g/grabowsp/.conda/envs/R_analysis

args = commandArgs(trailingOnly = TRUE)


### INPUTS ###
vcf_in <- args[1]

vcf_a <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

### SET OUTPUTS ###
out_file <- gsub('AltDosage', 'standard', vcf_in)

#####

vcf_b <- vcf_a

vcf_b[is.na(vcf_a)] <- './.'

vcf_b[vcf_a == 0] <- '4/0'

vcf_b[vcf_a == 1] <- '3/1'

vcf_b[vcf_a == 2] <- '2/2'

vcf_b[vcf_a == 3] <- '1/3'

vcf_b[vcf_a == 4] <- '0/4'

write.table(vcf_b, file = out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)


