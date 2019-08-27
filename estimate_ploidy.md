# Steps for estimating ploidy in switchgrass data

## Location of Data
### Multi-allelic SNP files
* VCF of genotypes
  * `/home/f1p1/tmp/sujan/Pvirgatum_V5/0.RAW_CALLS/MNP/Pvirgatum_1035g_MNP.bz2`
* Counts of each type of genotype
  * `/home/f1p1/tmp/sujan/Pvirgatum_V5/0.RAW_CALLS/MNP/Pvirgatum_1035g_MNP.counts.bz2`

## Idea for processing data:
1. Unzip vcf (bzcat or bunzip)
2. recursively, use 'cut' to extract info for a sample
3. use awk to choose lines that have a genotype
4. Use sed to remove genotype info from those lines
5. Tally number of SNPs with min count for all 3 allels of 2, 3, 4, 5, etc


## Test with top 500 lines
### Directory
* `/home/grabowsky_scratch/sg_ploidy/test`
### Make test file
```
bzcat \
/home/f1p1/tmp/sujan/Pvirgatum_V5/0.RAW_CALLS/MNP/Pvirgatum_1035g_MNP.bz2 | \
head -500 > /home/grabowsky_scratch/sg_ploidy/test/mnp_sg_snps_test500.vcf
```

```
awk '{if ($1 > 1) print $11}' mnp_sg_snps_test500.vcf | cut -d ':' -f 2 \
> test_samp1.txt
```
awk '{if ($11 > 1) print $11}' mnp_sg_snps_test500.vcf | cut -d ':' -f 2 | \
awk '{FS=","; OFS="\t"}{if ($1 > 0) print $1,$2,$3}'

awk 'NR>5 {print $11}' mnp_sg_snps_test500.vcf | cut -d ':' -f 2 | \
awk 'BEGIN{FS=","; OFS="\t"}{if (($1 > 0) && ($2 > 0) && ($3 > 0)) \
print $1,$2,$3}'
# This isn't working yet, but is the right track


awk 'NR>5 {print $11}' mnp_sg_snps_test500.vcf | head -50


* in R
```
test_file <- 'test_samp1.txt'
test <- read.table(test_file, sep = '\t', header = F, stringsAsFactors = F)
```


