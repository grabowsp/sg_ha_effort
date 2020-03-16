# Use nQuire to determine ploidy of samples

## Fetch archived BAMs with jamo
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy

module load jamo
module unload python/2.7-anaconda-2019.07
module load python/3.7-anaconda-2019.07

cat Pvirgatum_V5_bams_arcived.txt | \
awk '{print "jamo fetch all filename "$1}' | sh
```

## Generate soft link to the switchgrass bams
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy

for LIB in `cat /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/Pvirgatum_V5_libs.txt`; 
  do
  ln -s /global/dna/dm_archive/hudson_alpha/SNPS/Pvirgatum/$LIB'.Pvirgatum.V5.merged.gatk.bam' \
  /global/cscratch1/sd/grabowsp/sg_ploidy/sg_bam_links/$LIB'.bam';
  done
```

## Generate indices for each bam
### Make sublists of samples for submitting
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire
cp Pvirgatum_V5_libs.txt ./subsamp_lists/Pvirgatum_V5_libs_full.txt

cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/subsamp_lists
split -l 100 -d Pvirgatum_V5_libs_full.txt Pvirgatum_v5_libs_sub100_

split -l 6 -d -a 3 Pvirgatum_V5_libs_full.txt Pvirgatum_v5_libs_sub6_

```
### Run index job(s)
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_bam_links
sbatch make_bam_index_00.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_bam_links
for SUBSET in {01..10};
  do
  sed 's/sub100_00/'sub100_"$SUBSET"'/g' make_bam_index_00.sh > \
make_bam_index_$SUBSET'.sh';
  done

cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_bam_links
for SUBSET in {01..10};
  do
  sbatch make_bam_index_$SUBSET'.sh';
  done

```

## Run for all libraries
### Test with one subset of libraries
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire

cp /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test/test_IENS_workflow.sh \
CDS_nquire_000.sh

cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire
sbatch CDS_nquire_000.sh

for SUBSET in {001..172};
  do
  sed 's/sub6_000/'sub6_"$SUBSET"'/g' CDS_nquire_000.sh > \
CDS_nquire_$SUBSET'.sh';
  done


cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire
for SUBSET in {001..172};
  do
  sbatch CDS_nquire_$SUBSET'.sh';
  done
```

## Re-run failed jobs
### List of failed jobs
* list of subset_6 id's that had node failures on cori
  * /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire/node_fail_list.txt
### Remove results that will be re-made
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire

for SUBSET in `cat node_fail_list.txt`;
#for SUBSET in `head -3 node_fail_list.txt`;
  do
  for TEST_LIB in `cat /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/subsamp_lists/Pvirgatum_v5_libs_sub6_$SUBSET`;
    do
#    echo $TEST_LIB;
    rm $TEST_LIB*;
    done;
  done
```
### Re-run jobs
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire

for SUBSET in `cat node_fail_list.txt`;
  do
  sbatch CDS_nquire_$SUBSET'.sh';
  done
```

## Combine results
* Combined result file on NERSC
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire/sg_reseq_nQuire_results_summary_total.txt`
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire

cat `head -1 IICY.nQuire_res_summary.txt` 

head -1 IICY.nQuire_res_summary.txt | cat - *nQuire_res_summary.no_head.txt \
> sg_reseq_nQuire_results_summary_total.txt
```

## Extract results for c30 and c40
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/sg_nquire/sg_CDS_nquire

for RES in `ls *lrd_results_processed.txt`;
do
head -3 $RES | tail -n 1 >> CDS_nquire_c30_results_total.tmp;
done

head -1 XXBC.lrd_results_processed.txt | cat - \
CDS_nquire_c30_results_total.tmp > CDS_nquire_c30_results_total.txt

for RES in `ls *lrd_results_processed.txt`;
do
head -4 $RES | tail -n 1 >> CDS_nquire_c40_results_total.tmp;
done

head -1 XXBC.lrd_results_processed.txt | cat - \
CDS_nquire_c40_results_total.tmp > CDS_nquire_c40_results_total.txt

```

## Analyze nQuire results
### Notes
* Want to see
  * c20_min vs c50_min
  * c20_min vs coverage
  * c50_min vs coverage
  * c20_min vs nSNPs
  * c50_min vs nSNPs
  * c20_min vs dip_slope
  * c50_min vs dip_slope
  * c20_min vs tet_slope
  * c50_min vs tet_slope
  * c20_min vs Sujan
  * c50_min vs Sujan
  * dip_slope vs Sujan
  * tet_slope vs Sujan
* 4X-vs-8X calls are the exact same for c20 proportions and slope_vs_c50 proportions
* Will include the nquire calls (which include 6X calls) and cluster-calls
* It looks like nQuire works best with a lot of SNPs
  * the slopes show a relationship to c50 nSNPs but not to c20 nSNPs
  * there are WAY more c20 SNPs
  * Might be worth looking at results for c30/40
    * Numbers of SNPs
    * How those results compare to c20 results
###
* scripts for exploring results
  * `~/sg_ha_effort/r_scripts/nquire_exploration.r`
  * `~/sg_ha_effort/r_scripts/nquire_exploration_2.r`
* script for calling ploidy
  * `~/sg_ha_effort/r_scripts/nquire_ploidy_calls.r`



## Test workflow
### Run nQuire and get lrd results
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test
sbatch test_ABHM_workflow.sh

sbatch test_IENS_workflow.sh


```
### Test R script for processing results
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/nquire_test

Rscript ~/tools/sg_ha_effort/r_scripts/process_nQuire_results.r \
ABHM_lrd_results_tot.txt ABHM_nSNPS.txt


```
