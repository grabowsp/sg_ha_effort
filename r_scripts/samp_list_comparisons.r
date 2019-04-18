# Comparison of Reseq and Exome-capture sample lists

reseq_samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/sample_list_1005_to_Paul.txt'
reseq_samp_meta <- read.table(reseq_samp_meta_file, header = T, 
  sep = '\t', stringsAsFactors = F)

reseq_samp_meta <- read.table(reseq_samp_meta_file, header = T, sep = '\t', 
  stringsAsFactors = F, comment.char = '@')

quit(save = 'no')
