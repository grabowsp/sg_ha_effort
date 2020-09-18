# Steps for calculating delta-K statistices for `upexpand` structure results

## `upexpand` `all_tet`
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_tet

bash

RES_PRE=upexpand_alltet_
MAX_K=10
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_deltaK_steps.sh $RES_PRE $MAX_K $N_R

```

## `upexpand` `all_dip`
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/all_dip

bash

RES_PRE=upexpand_alldip_
MAX_K=10
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_deltaK_steps.sh $RES_PRE $MAX_K $N_R
```

## `upexpand` `partNA`
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/part_NA

bash

RES_PRE=upexpand_partNA_
MAX_K=10
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_deltaK_steps.sh $RES_PRE $MAX_K $N_R

```

## `upexpand` `pseudohap`
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/upexpand/pseudohap

bash

RES_PRE=upexpand_pseudohap_
MAX_K=10
N_R=3

/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/run_deltaK_steps.sh $RES_PRE $MAX_K $N_R

```

# Conclusions:
# k=2 is good for all but all_tet, and it's fine for all_tet
# k=5 is good for all 4 genotype styles
# k=3 is best for all_tet, but bad for the other 3 - weird



