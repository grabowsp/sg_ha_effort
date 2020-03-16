# Generate allele ratio files from VCF

# LOAD PACKAGES #

# INPUTS #
args = commandArgs(trailingOnly=TRUE)

vcf_in <- as.character(args[1])
vcf <- read.table(vcf_in, header = F, stringsAsFactors = F, sep = '\t')

samp_in <- as.character(args[2])
samps <- unlist(read.table(samp_in, header = F, sep = '\t', 
  stringsAsFactors = F))
# OUTPUTS #
out_file <- paste(vcf_in, '_alleleratios.rds', sep = '')

# SET VARIABLES #
snp_info_names <- c('CHROM', 'POS', 'ID', 'REF', 'ALT')

################
snp_info <- vcf[, c(1:5)]

counts <- apply(vcf[, c(10:ncol(vcf))], 2, function(x)
  unlist(lapply(strsplit(x, split = ':'), function(y) y[[2]]))
)

ratios <- apply(counts, 2, function(x)
  unlist(lapply(strsplit(x, split = ','), function(y)
  as.numeric(y[1])/sum(as.numeric(y))))
)

ratio_df <- data.frame(snp_info, ratios)

colnames(ratio_df) <- c(snp_info_names, samps)

saveRDS(ratio_df, file = out_file)

quit(save = 'no')


