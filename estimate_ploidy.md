# Steps for estimating ploidy in switchgrass data

## Location of Data
### Multi-allelic SNP files
* VCF of genotypes
  * `/home/f1p1/tmp/sujan/Pvirgatum_V5/0.RAW_CALLS/MNP/Pvirgatum_1035g_MNP.vcf.bz2`
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
/home/f1p1/tmp/sujan/Pvirgatum_V5/0.RAW_CALLS/MNP/Pvirgatum_1035g_MNP.vcf.bz2 \
| head -500 > /home/grabowsky_scratch/sg_ploidy/test/mnp_sg_snps_test500.vcf
```
*
```
bzcat \
/home/f1p1/tmp/sujan/Pvirgatum_V5/0.RAW_CALLS/MNP/Pvirgatum_1035g_MNP.vcf.bz2 \
| head -100000 > ./mnp_snps_100k.vcf

bzcat \
/home/f1p1/tmp/sujan/Pvirgatum_V5/0.RAW_CALLS/MNP/Pvirgatum_1035g_MNP.vcf.bz2 \
| head -1000000 > ./mnp_snps_1M.vcf


for j in {900..903};
  do
  LIB_NAME=`cut -f $j samp_name_line.txt`
  echo $LIB_NAME
  cut -f $j mnp_snps_1M.vcf | cut -d ':' -f 2 \
    > tmp.counts
  for i in {0..9}; 
    do
    echo $i
    awk 'BEGIN{FS=","; count = 0} NR<=5 { next} \
    (($1 > '"$i"') && ($2 > '"$i"') && ($3 > '"$i"')) \
    {count++} END {print count}' tmp.counts >> $LIB_NAME.counts
  done
done
```
*
```
for j in {10..1044};
  do
  LIB_NAME=`cut -f $j samp_name_line.txt`
  echo $LIB_NAME
  cut -f $j mnp_snps_100k.vcf | cut -d ':' -f 2 |
  awk 'BEGIN{FS=","; count = 0} NR<=5 { next} \
    (($1 > 0) && ($2 > 0) && ($3 > 0)) \
    {count++} END {print count}' tmp.counts
done
```





* in R
```
test_file <- 'test_samp1.txt'
test <- read.table(test_file, sep = '\t', header = F, stringsAsFactors = F)
```


