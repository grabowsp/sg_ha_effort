# Pop structure analysis using geographic samples for subset of data

## Generate `geo` Chr01K CDS Distance Matrix
* Only includes "Natural Collection" 8X samples
* Using polyploid genotype
  * 4X get disomic genotypes, 8X get tetrasomic genotypes
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/geo826_dists

sbatch gen_geo_poly_dists_test.sh
sbatch gen_geo_poly_dists_Chr01K.sh
```

## Generate `expand_geo` Chr01K CDS Distance Matrix
* Includes "Natural Collection" and "Cultivar" 8X samples
* Using polyploid genotypes
```
cd /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_genos_popstructure/expand_geo839_dists

sbatch gen_expandgeo_poly_dists_test.sh
sbatch gen_expandgeo_poly_dists_Chr01K.sh



```

