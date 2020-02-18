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
