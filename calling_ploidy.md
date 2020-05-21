# Process Used for Calling Ploidy in Switchgrass Resequencing Data

## Overview
* Use two independent methods to infer ploidy from Illumina Data
  * nQuire: ratio of reads for each allele compared to theoretical expectations
    * ratios should differ for disomic, trisomic, and tetrasomic individuals
  * Count of multi-allelic SNPs (MNP)
    * 8X samples should have positions with 3 or 4 alleles; 4X samples should not
  * v3 has "total_ploidy_2" that has some adjustements to "total_ploidy"
    * Correct some nQuire-based calls that seem wrong
    * Adjust some samples to 6X
## Output
* Ploidy call and other info
  * on NERSC
    * `/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v3.0.txt`
  * on HA
    * `/home/t4c1/WORK/grabowsk/data/switchgrass/pseudohap/sg_ploidy_results_v3.0.txt`
