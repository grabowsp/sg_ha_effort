# Generate Files that Contain Genic and exon Positions

## Directory with genic position files
* `/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions`
### Example file
* Genic positions:
  * `Chr01K_genic_positions.rds`

## Generate Genic Position Files
* in R
```
sg_gff_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ref_stuff/annotation/Pvirgatum_516_v5.1/Pvirgatum_516_v5.1.gene.gff3'

sg_gff <- read.table(sg_gff_file, header = F, stringsAsFactors = F, sep = '\t')

sg_genes <- sg_gff[which(sg_gff[,3] == 'gene'), ]

tmp_chr_names <- unique(sg_genes[,1])
scaf_inds <- grep('scaffold', tmp_chr_names)
chr_names <- tmp_chr_names[-scaf_inds]

out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/'

for(sg_chr in chr_names){
  tmp_chr_inds <- grep(sg_chr, sg_genes[,1])
  tmp_test_df <- sg_genes[tmp_chr_inds, ]
  #
  tmp_genic_inds_list <- apply(tmp_test_df, 1, function(x) c(x[4]:x[5]))
  tmp_genic_inds <- sort(unique(unlist(tmp_genic_inds_list)))
  tmp_out_file <- paste(out_dir, sg_chr, '_genic_positions.rds', sep = '')
  saveRDS(tmp_genic_inds, tmp_out_file)
  print(sg_chr)
}
```

## Generate Exonic Position Files
* in R
```
sg_gff_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ref_stuff/annotation/Pvirgatum_516_v5.1/Pvirgatum_516_v5.1.gene.gff3'

sg_gff <- read.table(sg_gff_file, header = F, stringsAsFactors = F, sep = '\t')

sg_CDS <- sg_gff[which(sg_gff[,3] == 'CDS'), ]

tmp_chr_names <- unique(sg_CDS[,1])
scaf_inds <- grep('scaffold', tmp_chr_names)
chr_names <- tmp_chr_names[-scaf_inds]

out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/'

for(sg_chr in chr_names){
  tmp_chr_inds <- grep(sg_chr, sg_CDS[,1])
  tmp_test_df <- sg_CDS[tmp_chr_inds, ]
  #
  tmp_CDS_inds_list <- apply(tmp_test_df, 1, function(x) c(x[4]:x[5]))
  tmp_CDS_inds <- sort(unique(unlist(tmp_CDS_inds_list)))
  tmp_out_file <- paste(out_dir, sg_chr, '_CDS_positions.rds', sep = '')
  saveRDS(tmp_CDS_inds, tmp_out_file)
  print(sg_chr)
}
```


