# Use nQuire to determine ploidy of samples

## Fetch archived BAMs with jamo
```
/global/cscratch1/sd/grabowsp/sg_ploidy

module load jamo
module unload python/2.7-anaconda-2019.07
module load python/3.7-anaconda-2019.07

cat Pvirgatum_V5_bams_arcived.txt | \
awk '{print "jamo fetch all filename "$1}' | sh
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
