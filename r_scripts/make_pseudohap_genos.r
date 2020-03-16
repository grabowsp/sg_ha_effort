# Generate allele ratio files from VCF

# LOAD PACKAGES #

# INPUTS #
args = commandArgs(trailingOnly=TRUE)

ratios_in <- as.character(args[1])
ratios <- readRDS(ratios_in)
 
# OUTPUTS #
out_file <- paste(gsub('_alleleratios.rds', '', ratios_in), 
  '_pseudohapgenos_v01.rds', sep = '')

# SET VARIABLES #

################
snp_info <- ratios[, c(1:5)]

pseudohaps <- apply(ratios[, c(6:ncol(ratios))], 2, function(x)
  sapply(x, function(y) rbinom(n = 1, size = 1, prob = y))
)

pseudohap_df <- data.frame(snp_info, pseudohaps, stringAsFactors = F)

saveRDS(pseudohap_df, file = out_file)

quit(save = 'no')


