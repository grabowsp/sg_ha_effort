# sort the CDS VCFs so that SNPs are in correct order

```

module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs

FORMAT_HEAD=/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/vcf_format_header.txt

CHRNAME=Chr01N

gunzip -c $CHRNAME'.polyploid.CDS.vcf.gz' > tmp_$CHRNAME'.vcf'

tail -n +5 tmp_$CHRNAME'.vcf' > tmp_$CHRNAME'.vcf_1'
cat $FORMAT_HEAD tmp_$CHRNAME'.vcf_1' > tmp_$CHRNAME'.vcf_2'

bcftools sort tmp_$CHRNAME'.vcf_2' -Oz -o $CHRNAME'.polyploid.CDS.sorted.vcf.gz'

rm tmp_$CHRNAME*

for CHR_SUB in 02K 02N 03K 03N 04K 04N;
do
CHRNAME=Chr$CHR_SUB
gunzip -c $CHRNAME'.polyploid.CDS.vcf.gz' > tmp_$CHRNAME'.vcf';
tail -n +5 tmp_$CHRNAME'.vcf' > tmp_$CHRNAME'.vcf_1';
cat $FORMAT_HEAD tmp_$CHRNAME'.vcf_1' > tmp_$CHRNAME'.vcf_2';
bcftools sort tmp_$CHRNAME'.vcf_2' -Oz -o \
$CHRNAME'.polyploid.CDS.sorted.vcf.gz';
rm tmp_$CHRNAME*;
done

```
* try this script for remaining Chromosomes
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs
sbatch sort_remain_Chrs.sh
```


