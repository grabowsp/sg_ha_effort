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

sg_gene_bed_df <- sg_genes[, c(1,4,5,9)]

gene_name_short <- gsub('Name=', '', unlist(lapply(
  strsplit(sg_gene_bed_df[,4], split = ';'), function(x) x[2])))

sg_gene_bed_df[,4] <- gene_name_short

out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/'

gene_bed_out <- paste(out_dir, 'sg_v5_genes.bed', sep = '')

write.table(sg_gene_bed_df, file = gene_bed_out, quote = F, sep = '\t',
  row.names = F, col.names = F)

```

## Generate Exonic Position Files
* in R
```
sg_gff_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ref_stuff/annotation/Pvirgatum_516_v5.1/Pvirgatum_516_v5.1.gene.gff3'

sg_gff <- read.table(sg_gff_file, header = F, stringsAsFactors = F, sep = '\t')

sg_CDS <- sg_gff[which(sg_gff[,3] == 'CDS'), ]

sg_CDS_bed_df <- sg_CDS[, c(1,4,5,9)]

CDS_name_short <- gsub('ID=', '', unlist(lapply(
  strsplit(sg_CDS_bed_df[,4], split = ';'), function(x) x[1])))

sg_CDS_bed_df[,4] <- CDS_name_short

out_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/genic_positions/'

CDS_bed_out <- paste(out_dir, 'sg_v5_CDS.bed', sep = '')

write.table(sg_CDS_bed_df, file = CDS_bed_out, quote = F, sep = '\t',
  row.names = F, col.names = F)



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


