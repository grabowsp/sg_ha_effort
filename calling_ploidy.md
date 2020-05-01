# Process Used for Calling Ploidy in Switchgrass Resequencing Data

## Overview
* Use two independent methods to infer ploidy from Illumina Data
  * nQuire: ratio of reads for each allele compared to theoretical expectations
    * ratios should differ for disomic, trisomic, and tetrasomic individuals
  * Count of multi-allelic SNPs (MNP)
    * 8X samples should have positions with 3 or 4 alleles; 4X samples should not

## Output
* Ploidy call and other info
  * on NERSC
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v2.0.txt`
  * on HA
    * `/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/sg_ploidy_results_v2.0.txt`
