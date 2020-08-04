# Steps and looking at deltaK and other metrics used to evaluate K as 
#    demonstrated in the Evanno paper

## Test robustness of results by  using different types of genotypes
* all samples with 4 lines
* all samples with 2 lines
* all samples with 4 lines but 4X (disomic) SG have NA (-9) in lines 3 and 4
* all samples with pseudhaploid genotypes

## `expand_geo` `all_tet` Results 
### Get Estimated Ln Prob of Data
#### Extract values from STRUCTURE output
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet

bash

OUT_PRE=expandgeo_alltet

for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' $OUT_PRE'_'$KT'.'$KR'_f' >> \
      $OUT_PRE'_LnProbData.txt';
    done;
  done;
```
#### Calculate delta-K and other Evanno et al metrics
```
bash
source activate R_analysis

# ADJUST THESE VARIABLES TO CORRECT NAMES AND VALUES
PARENT_DIR=/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/

FILE_SHORT=expandgeo_alltet_LnProbData.txt

NKS=10
NREPS=3

SAMP_SET_NAME=expandgeo_alltet

####
LN_DATA_FILE=$PARENT_DIR$FILE_SHORT

Rscript /home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/gen_deltaK_plot.r $LN_DATA_FILE $NKS $NREPS $SAMP_SET_NAME

```
#### Output files
* Interpretation:
  * K=2 is best; K=3 is okay, other K's don't show much support
* Figure
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_LnProbData_deltaK.pdf`
* Table
  * `/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_LnProbData_deltaK.txt`









## Run CLUMPP on results
### K=3
#### Concatenate results from the 3 replicates
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2

cat geo_v2_3.*_q > geo_v2_k3_combo_q_tmp
```
#### Format input file for CLUMPP
* in R
```
bash 
source activate R_analysis

in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo_q_tmp'

in_res <- read.table(in_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(in_res)/3

samp_int_vec <- rep(seq(nsamps), times = 3)

third_vec <- rep('(x)', times = nrow(in_res))
#third_vec <- rep(paste('(', in_res[,1], ')', sep = ''), times = 3)

pop_vec <- rep(NA, times = nrow(in_res))
fifth_vec <- rep(':', times = nrow(in_res))

out_res_tmp <- data.frame( C1=samp_int_vec, C2 = samp_int_vec, 
  C3 = third_vec, C4 = samp_int_vec, C5 = fifth_vec,
  stringsAsFactors = F)

out_res <- cbind(out_res_tmp, in_res[, c(2:ncol(in_res))])

out_file <- gsub('_q_tmp', '.clumpp.indfile', in_file)

write.table(out_res, out_file, quote=F, sep = '\t', row.names = F, 
  col.names = F)
```
#### Run CLUMPP
```
cd /home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2
qsub geo_v2_k3_clumpp.sh

```
#### Process output from CLUMPP
```
bash 
source activate R_analysis

pre_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo_q_tmp'
pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo.clumpp.out'
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

nreps <- 3
nsamps <- nrow(pre_info)/3

out_info <- data.frame(LIB = pre_info[seq(nsamps),1], 
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

out_file <- gsub('clumpp.out', 'clumpp.processed', res_file)
write.table(out_info, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
#### Make bargraph of CLUMPP results
```
bash
source activate R_analysis

library(ggplot2)

struc_func_file <- '/home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/struc_functions.r'
source(struc_func_file)

c_res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/geo_v2/geo_v2_k3_combo.clumpp.processed'

c_res <- read.table(c_res_file, header = F, sep = '\t', stringsAsFactors = F)

ploidy_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/sg_ploidy_results_v3.0.txt'
ploidy <- read.table(ploidy_file, header = T, stringsAsFactors = F, sep = '\t')

meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'
meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = '\t',
  quote = "", comment.char = '$')

### evaluate what the groups are

pop_info_list <- list()

test_order_mat <- apply(c_res[, -1], 1, order, decreasing = T)

for(i in seq(ncol(c_res)-1)){
  test_pop <- i
  pop_name <- paste('Population', i)
  test_libs <- c_res[which(apply(test_order_mat, 2, 
    function(x) x[1] == test_pop)), 1]
  ploidy_inds <- which(ploidy$lib %in% test_libs)
  meta_inds <- which(meta$LIBRARY %in% test_libs)
  tmp_list <- list()
  tmp_list[['Cloroplast Ecotype']] <- table(ploidy$ECOTYPE_SNP_CHLR
    [ploidy_inds])
  tmp_list[['Old Subpop']] <- table(ploidy$SUBPOP_SNP[ploidy_inds])
  tmp_list[['Ploidy']] <- table(ploidy$total_ploidy_2[ploidy_inds])
  tmp_list[['State']] <- table(meta$STATE[meta_inds])
  pop_info_list[[pop_name]] <- tmp_list
}

# pop 1 = EC lowland
# pop 2 = upland
# pop 3 = TX lowland

pop_order <- c(2,1,3)

######

test_samp_orders <- order_struc_samps(clumpp_result_df = c_res, 
  pop_order = pop_order, zero_cut = 0.1)

#####

res_ord <- c_res[test_samp_orders, ]

samp_vec <- rep(res_ord[,1], times = ncol(c_res)-1)
group_vec <- rep(seq(ncol(c_res)-1), each = nrow(res_ord))

res_vec <- c()
for(j in c(2:ncol(res_ord))){
  res_vec <- c(res_vec, res_ord[,j])
}

plot_df <- data.frame(LIB = samp_vec, GROUP = group_vec, MEMB = res_vec,
  stringsAsFactors = F)

plot_df$LIB <- factor(plot_df$LIB, levels = res_ord[,1])

plot_df$GROUP <- factor(plot_df$GROUP, levels = pop_order)

group_col_vec <- c()
group_col_vec['1'] <- 'yellow3'
group_col_vec['2'] <- 'blue2'
group_col_vec['3'] <- 'red2'

group_lab_vec <- c()
group_lab_vec['1'] <- 'EC'
group_lab_vec['2'] <- 'MW'
group_lab_vec['3'] <- 'TX'
group_lab_vec <- group_lab_vec[pop_order]

#group_palette <- scale_colour_manual(name = 'K_Group', values = group_col_vec)

gg_bar <- ggplot(plot_df, aes(y = MEMB, x = LIB, fill = GROUP)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual('Group', values = group_col_vec, 
    labels = group_lab_vec) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = 'Group Membership')

out_file <- gsub('clumpp.processed', 'STRUCTURE.memb.pdf', c_res_file)

pdf(out_file, width = 12*6, height = 4)
gg_bar
dev.off()


```
