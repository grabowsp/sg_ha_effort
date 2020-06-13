# Notes on selecting targets for additional sequencing

## Files for analysis
* Metadata file
  * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/Reseq_Metadata_Sept_2019_Edited_for_R.tsv`
* Ploidy file
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v3.0.txt`
* Accessions with multiple ploidy levels
  * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/accessions_with_multiple_ploidy.txt`
* Samples in accessions with multiple ploidy levels
  * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/samps_in_multiploidy_accessions.txt`
* Phenotype info
  * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/pheno_status_8X_lines_fv.csv`
* Samples in Genome Manuscript
  * `/global/cscratch1/sd/grabowsp/sg_ploidy/samp901_lib_names.txt`
* Reseq priority sample list from Joe
  * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/8X_Priority_List_FV.csv`

## Examine samples
```
module load python3/3.7-anaconda-2019.10
source activate r_adegenet_env

library(adegenet)

### LOAD DATA ###

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/Reseq_Metadata_Sept_2019_Edited_for_R.tsv'

meta <- read.table(meta_file, sep = '\t', header = T, stringsAsFactors = F, 
  quote = "", comment.char = '$')

ploidy_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/ploidy_calling/sg_ploidy_results_v3.0.txt'

ploidy <- read.table(ploidy_file, sep = '\t', header = T, stringsAsFactors = F)

pheno_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/pheno_status_8X_lines_fv.csv'

phenos <- read.table(pheno_file, sep = ',', header = T, stringsAsFactors = F)

multiploid_accession_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/accessions_with_multiple_ploidy.txt'

multi_acc <- read.table(multiploid_accession_file, sep = '\t', header = F, 
  stringsAsFactors = F)

multiploid_sample_file <- '/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/samps_in_multiploidy_accessions.txt'

multi_samp <- read.table(multiploid_sample_file, sep = '\t', header = T, 
  stringsAsFactors = F)

genome_lib_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/samp901_lib_names.txt'
genome_libs <- read.table(genome_lib_file, sep = '\t', header = F, 
  stringsAsFactors = F)

all_pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Combo.595K.polyploid.CDS.geosamps.genlight.PCAresults.rds'

all_pca_full <- readRDS(all_pca_res_file)
all_pca <- all_pca_full$scores

up_pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/up_geo_samps/Combo.sub.polyploid.CDS.upgeosamps.genlight.PCAresults.rds'
up_pca_full <- readRDS(up_pca_res_file)
up_pca <- up_pca_full$scores

low_pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/low_geo_samps/Combo.sub.polyploid.CDS.lowgeosamps.genlight.PCAresults.rds'
low_pca_full <- readRDS(low_pca_res_file)
low_pca <- low_pca_full$scores
####

# Look at Accessions without representation

well_rep_seq_inds <- intersect(
  which(phenos$Well_Replicated_Phenotypic_Data == 'Yes'), 
  which(is.na(phenos$LIBRARY)==F))

well_rep_seq_df <- phenos[well_rep_seq_inds, ]

well_rep_noseq_inds <- intersect(
  which(phenos$Well_Replicated_Phenotypic_Data == 'Yes'), 
  which(is.na(phenos$LIBRARY)))

well_rep_noseq_df <- phenos[well_rep_noseq_inds, ]

#

na_acc_inds <- which(is.na(meta$ACC))

length(intersect(meta$LIBRARY[na_acc_inds], genome_libs[,1]))
# 44 - so these ARE included in some of the analysis

table(meta$COLLECTION_TYPE[na_acc_inds])

na_acc_nat_inds <- intersect( na_acc_inds, 
  which(meta$COLLECTION_TYPE == 'Natural Collection'))
# CWAB = DAC6 which is a cultivar, J037
# ABHM = WBC3a - this is a weird sample from TX
# PZXS = PI 414065, BN-14668-65; this is also J008
# PYSX = PI 315723, BN-8358-62; this is also J003
# PZXO = same as PYSX, but this is 8X; J003
# PYSP = PI 315725; BN-14669-92; this is also J023
# OBWZ = same as PYSX and PZXO; J003
# OBXC = natural collection from TX; not much other info about the accession
# WZCY = Dacotah_W01; J037
# WZCW = plant specifically chosen from PI422001; this accession is also J020
# WZGC = Mexico

meta$ACC[which(meta$LIBRARY == 'CWAB')] <- 'J037'
meta$ACC[which(meta$LIBRARY == 'PZXS')] <- 'J008'
meta$ACC[which(meta$LIBRARY == 'PYSX')] <- 'J003'
meta$ACC[which(meta$LIBRARY == 'PZXO')] <- 'J003'
meta$ACC[which(meta$LIBRARY == 'PYSP')] <- 'J023'
meta$ACC[which(meta$LIBRARY == 'OBWZ')] <- 'J003'
meta$ACC[which(meta$LIBRARY == 'WZCY')] <- 'J037'
meta$ACC[which(meta$LIBRARY == 'WZCW')] <- 'J020'

meta$ACC[which(meta$LIBRARY == 'ABHM')] <- 'WBC3a'
meta$ACC[which(meta$LIBRARY == 'WZGC')] <- 'Mexico'

acc_vec <- setdiff(unique(meta$ACC), NA)

acc_n_nas <- sapply(acc_vec, function(x) 
  sum(is.na(meta$LIBRARY[which(meta$ACC == x)])))

acc_n_samps <- sapply(acc_vec, function(x)
  sum(meta$ACC == x, na.rm = T))

no_lib_acc_inds <- which(acc_n_nas == acc_n_samps)

no_lib_accs <- acc_vec[no_lib_acc_inds]

# 201 accessions have no libraries

no_lib_samp_inds <- c()
for(i in no_lib_accs){
  tmp_inds <- which(meta$ACC == i)
  no_lib_samp_inds <- c(no_lib_samp_inds, tmp_inds)
}

no_lib_acc_df <- meta[no_lib_samp_inds, ]

length(intersect(well_rep_noseq_df$PLANT_ID, no_lib_acc_df$PLANT_ID))
# 67

no_lib_acc_df$well_rep <- F

no_lib_acc_df$well_rep[
  (no_lib_acc_df$PLANT_ID %in% well_rep_noseq_df$PLANT_ID)] <- T

length(unique(no_lib_acc_df$ACC[no_lib_acc_df$well_rep == T]))
# 54

length(setdiff(unique(no_lib_acc_df$ACC), 
  unique(no_lib_acc_df$ACC[no_lib_acc_df$well_rep == T])))
# 147 Accessions without representation when include well repped samples

length(intersect(no_lib_acc_df$PLANT_ID, phenos$PLANT_ID[phenos$Any_Phenotypic_Data == 'Yes']))

no_lib_acc_df$any_pheno <- F

no_lib_acc_df$any_pheno[
  no_lib_acc_df$PLANT_ID %in% 
  phenos$PLANT_ID[phenos$Any_Phenotypic_Data == 'Yes']] <- T

length(unique(no_lib_acc_df$ACC[no_lib_acc_df$any_pheno == T]))
# 123; 69 that are not well replicated

length(setdiff(unique(no_lib_acc_df$ACC), 
  unique(no_lib_acc_df$ACC[no_lib_acc_df$any_pheno == T])))
# 78 accessions without representative that is phenotyped

not_well_inds <- setdiff(which(no_lib_acc_df$any_pheno == T),
  which(no_lib_acc_df$well_rep == T))

length(intersect(unique(no_lib_acc_df$ACC[not_well_inds]), 
  unique(no_lib_acc_df$ACC[no_lib_acc_df$well_rep == T])))
# 24 accessions have both well-replicated and poorly-replicated phenotyped 
#   samples

no_pheno_cult_accs <- unique(
  no_lib_acc_df$ACC[intersect(
  which(no_lib_acc_df$COLLECTION_TYPE == 'Cultivar'),
  which(no_lib_acc_df$any_pheno == F))])

length(setdiff(no_pheno_cult_accs, 
  unique(no_lib_acc_df$ACC[which(no_lib_acc_df$any_pheno == T)])))
# 7 non-phenotyped accessions are cultivars

no_pheno_breed_accs <- unique(
  no_lib_acc_df$ACC[intersect(
  which(no_lib_acc_df$COLLECTION_TYPE == 'Breeding Selection'),
  which(no_lib_acc_df$any_pheno == F))])

length(setdiff(no_pheno_breed_accs,
  unique(no_lib_acc_df$ACC[which(no_lib_acc_df$any_pheno == T)])))
# 35 non-phenotyped accessions are breeding selections

no_pheno_nat_accs <- unique(
  no_lib_acc_df$ACC[intersect(
  which(no_lib_acc_df$COLLECTION_TYPE == 'Natural Collection'),
  which(no_lib_acc_df$any_pheno == F))])

length(setdiff(no_pheno_nat_accs,
  unique(no_lib_acc_df$ACC[which(no_lib_acc_df$any_pheno == T)])))
# 36 non-phenotyped accessions are natural collections

no_pheno_nat_accs_2 <- setdiff(no_pheno_nat_accs,
  unique(no_lib_acc_df$ACC[which(no_lib_acc_df$any_pheno == T)]))

no_pheno_nat_samp_inds <- which(no_lib_acc_df$ACC %in% no_pheno_nat_accs_2)

no_lib_acc_df[no_pheno_nat_samp_inds, 
  c('ACC', 'PLANT_ID', 'PLOIDY', 'COLLECTION_TYPE', 
  'SOURCE_ID_1', 'SOURCE_ID_2', 'COUNTRY', 'STATE')]

# Look at mixed ploidy accessions
multi_samp[is.na(multi_samp$LIBRARY),]

multi_samp$LIBRARY[multi_samp$ACC == 'J179']

all_pca[which(rownames(all_pca) == 'IEXZ'), c(1:4)]
# PC1 = Lowland, PC2 = Eastcoast, PC3 = ambiguous (EC), PC4 = ambiguous (EC) 
all_pca[which(rownames(all_pca) == 'IEYA'), c(1:4)]
# PC1 = Lowland, PC2 = L_admix, GC?, PC3 = ?, PC4 = Not EC
all_pca[which(rownames(all_pca) == 'IEYB'), c(1:4)]

multi_samp$LIBRARY[multi_samp$ACC == 'J249']

all_pca[which(rownames(all_pca) == 'IELY'), c(1:4)]
# PC1 = Lowland, PC2 = TX, PC3 = ?, PC4 = ?
all_pca[which(rownames(all_pca) == 'IEYZ'), c(1:4)]
# PC1 = Lowland, PC2 = TX, PC3 = ?, PC4 = ?

multi_samp$LIBRARY[multi_samp$ACC == 'J288']

all_pca[which(rownames(all_pca) == 'IJAT'), c(1:4)]

multi_samp$LIBRARY[which(multi_samp$ACC == 'J447')]
multi_samp[which(multi_samp$ACC == 'J447'), ]
all_pca[which(rownames(all_pca) == 'IIKQ'), c(1:4)]
# not in PCA

multi_samp$LIBRARY[which(multi_samp$ACC == 'J464')]
multi_samp[which(multi_samp$ACC == 'J464'), ]
all_pca[which(rownames(all_pca) == 'ITAU'), c(1:4)]
# PC1 = lowland but closer to upland, PC2 = GC
all_pca[which(rownames(all_pca) == 'IILE'), c(1:4)]
# PC1 = lowland, PC2 = GC


# Check if the nQuire results identify additional mixed-ploidy samples

change_ploidy_inds <- which(ploidy$PLOIDY != ploidy$total_ploidy_2)
# IEYD - FL, shows mixed ploidy in population; all 4 samples are sequenced

for(i in seq(length(change_ploidy_inds))){
  print(i)
  tmp_lib <- ploidy$lib[change_ploidy_inds[i]]
  print(ploidy[change_ploidy_inds[i], c('lib', 'PLOIDY', 'total_ploidy_2')])
  tmp_acc <- meta$ACC[which(meta$LIBRARY == tmp_lib)]
  print(tmp_acc)
  print(meta[which(meta$ACC == tmp_acc), c('ACC', 'PLANT_ID', 'LIBRARY', 
    'SUBPOP_SNP', 'PLOIDY', 'COLLECTION_TYPE', 'STATE')])
}

all_pca[which(rownames(all_pca) == 'INGI'), c(1:4)]
# PC1 = upland, PC4 = interface of 4X and 8X

# HOW MANY SAMPLES PER ACCESSION?
table(table(meta$ACC[which(is.na(meta$LIBRARY) == F)]))
  1   2   3   4   5   6   7   8 
216  98 127  11   1   3   1   1

table(table(
  meta$ACC[intersect(which(is.na(meta$LIBRARY) == F), 
    which(meta$PLOIDY == '4X'))]))
  1   2   3   4   5   6   7   8 
218  80  99  11   2   2   1   1

table(table(
  meta$ACC[intersect(which(is.na(meta$LIBRARY) == F), 
    which(meta$PLOIDY == '8X'))]))
 1  2  3 
34 20 15

# most accessions, including most 4X accessions, represented by just
#  a single library

# SOUTHERN 8X LIBRARIES
phenos$STATE <- NA
for(i in seq(nrow(phenos))){
  tmp_m_ind <- which(meta$PLANT_ID == phenos$PLANT_ID[i])
  phenos$STATE[i] <- meta$STATE[tmp_m_ind]
}

table(phenos$STATE)
# FL = 10, GA = 9, LA = 6, MS = 14, NC = 12, TX = 21

table(phenos$STATE[which(phenos$Well_Replicated_Phenotypic_Data == 'No')])
# FL = 10, GA = 7, LA = 4, MS = 8, NC = 7, TX = 20


meta[intersect(which(is.na(meta$LIBRARY)), 
  which(meta$STATE == 'New Hampshire')), 
  c('ACC', 'PLANT_ID', 'LIBRARY', 'PLOIDY', 'COLLECTION_TYPE', 'LATITUDE', 
  'LONGITUDE')]

unique(meta$ACC[intersect(which(is.na(meta$LIBRARY)),
  which(meta$STATE == 'Texas'))])

meta[which(meta$STATE == 'Massachusetts'), 
  c('ACC', 'PLANT_ID', 'LIBRARY', 'PLOIDY', 'COLLECTION_TYPE', 'LATITUDE', 
  'LONGITUDE')]

meta[intersect(which(is.na(meta$LIBRARY) == F),
  intersect(which(meta$STATE == 'Texas'), which(meta$PLOIDY == '8X'))), 
  c('ACC', 'PLANT_ID', 'LIBRARY', 'PLOIDY', 'COLLECTION_TYPE', 'LATITUDE', 
  'LONGITUDE')]

unique(meta$ACC[intersect(which(is.na(meta$LIBRARY) == F),
  intersect(which(meta$STATE == 'Texas'), which(meta$PLOIDY == '4X')))])

all_pca[c('IHYU', 'IIBT', 'IELY', 'ASHOA', 'IRWY', 'IRWZ'), 
  c(1:4)]
up_pca[c('IHWQ', 'IRYM'), c(1:4)]
low_pca['IEXZ', c(1:4)]

### Generate List for Selecting Samples

phenos$ACC <- NA
for(i in seq(nrow(phenos))){
  tmp_ind <- which(meta$PLANT_ID == phenos$PLANT_ID[i])
  tmp_acc <- meta$ACC[tmp_ind]
  phenos$ACC[i] <- tmp_acc
}

for(naind in which(is.na(phenos$ACC))){
  phenos$ACC[naind] <- phenos$PLANT_ID[naind]
}


phenos$acc_pheno <- 'No'
for(sga in unique(phenos$ACC)){
  acc_inds <- which(phenos$ACC == sga)
  if('Yes' %in% phenos$Any_Phenotypic_Data[acc_inds]){
    phenos$acc_pheno[acc_inds] <- 'Yes'
  }
}

length(intersect(which(phenos$acc_pheno == 'Yes'), 
  intersect(which(phenos$Any_Phenotypic_Data == 'No'),
  which(is.na(phenos$LIBRARY)))))
# 63 unsequenced samples have no phenotypic data, but another sample from thier
#   accession has phenotypic data

phenos$mix_ploid_acc <- 'No'
for(mpx in seq(nrow(multi_samp))){
  tmp_ind <- which(phenos$PLANT_ID == multi_samp$PLANT_ID[mpx])
  phenos$mix_ploid_acc[tmp_ind] <- 'Yes'
}

j636A_df <- data.frame(PLANT_ID = 'J636.A', LIBRARY = NA, PLOIDY = '4X', 
  Any_Phenotypic_Data = 'No', Well_Replicated_Phenotypic_Data = 'No',
  ACC = 'J636', acc_pheno = 'No', mix_ploid_acc = 'Yes', stringsAsFactors = F)

phenos <- rbind(phenos, j636A_df)

phenos$STATE <- NA
phenos$COLLECTION_TYPE <- NA
phenos$LATITUDE <- NA
phenos$LONGITUDE <- NA
for(i in seq(nrow(phenos))){
  tmp_meta_ind <- which(meta$PLANT_ID == phenos$PLANT_ID[i])
  tmp_state <- meta$STATE[tmp_meta_ind]
  tmp_type <- meta$COLLECTION_TYPE[tmp_meta_ind]
  tmp_lat <- meta$LATITUDE[tmp_meta_ind]
  tmp_long <- meta$LONGITUDE[tmp_meta_ind]
  phenos$STATE[i] <- tmp_state
  phenos$COLLECTION_TYPE[i] <- tmp_type
  phenos$LATITUDE[i] <- tmp_lat
  phenos$LONGITUDE[i] <- tmp_long
}

phenos$south_8X <- 'No'
phenos$south_8X[which(phenos$STATE == 'Florida')] <- 'Yes'
phenos$south_8X[which(phenos$STATE == 'Georgia')] <- 'Yes'
phenos$south_8X[which(phenos$STATE == 'North Carolina')] <- 'Yes'
phenos$south_8X[which(phenos$STATE == 'Louisiana')] <- 'Yes'
phenos$south_8X[which(phenos$STATE == 'Mississippi')] <- 'Yes'
phenos$south_8X[which(phenos$STATE == 'Texas')] <- 'Yes'


ga_unseq_4X <- intersect(which(meta$PLOIDY == '4X'), 
  intersect(which(is.na(meta$LIBRARY)), which(meta$STATE == 'Georgia')))

ga_4X_df <- data.frame(
  PLANT_ID = meta$PLANT_ID[ga_unseq_4X],
  LIBRARY = rep(NA, times = 3),
  PLOIDY = rep('4X', times = 3),
  Any_Phenotypic_Data = rep('No', times = 3),
  Well_Replicated_Phenotypic_Data = rep('No', times = 3),
  ACC = meta$ACC[ga_unseq_4X],
  acc_pheno = rep('No', times = 3),
  mix_ploid_acc = rep('No', times = 3),
  STATE = meta$STATE[ga_unseq_4X],
  COLLECTION_TYPE = meta$COLLECTION_TYPE[ga_unseq_4X],
  LATITUDE = meta$LATITUDE[ga_unseq_4X],
  LONGITUDE = meta$LONGITUDE[ga_unseq_4X],
  south_8X = rep('Yes', times = 3),
  stringsAsFactors = F)

phenos <- rbind(phenos, ga_4X_df)

phenos$EC_8X <- 'No'
phenos$EC_8X[which(phenos$STATE == 'Virginia')] <- 'Yes'
phenos$EC_8X[which(phenos$STATE == 'Maryland')] <- 'Yes'
phenos$EC_8X[which(phenos$STATE == 'Massachusetts')] <- 'Yes'
phenos$EC_8X[which(phenos$STATE == 'New Hampshire')] <- 'Yes'
phenos$EC_8X[which(phenos$STATE == 'Maine')] <- 'Yes'
phenos$EC_8X[intersect(which(phenos$STATE == 'Maine'), 
  which(phenos$LATITUDE < 41.53))] <- 'Yes'

low_cov_libs <- ploidy$lib[which(ploidy$seq_cov < 1e10)]

phenos$low_cov <- 'No'
phenos$low_cov[which(phenos$LIBRARY %in% low_cov_libs)] <- 'Yes'

hi_rep_inds <- intersect(which(is.na(phenos$LIBRARY)), 
  which(phenos$Well_Replicated_Phenotypic_Data == 'Yes'))
# 76

south_8X_tmp_inds <- intersect(which(is.na(phenos$LIBRARY)),
  which(phenos$south_8X == 'Yes'))
south_8X_inds <- setdiff(south_8X_tmp_inds, hi_rep_inds)
length(south_8X_inds)
#[1] 45
length(unique(phenos$ACC[south_8X_inds]))
# 37 accessions

EC_8X_tmp_inds <- intersect(which(is.na(phenos$LIBRARY)),
  which(phenos$EC_8X == 'Yes'))
EC_8X_inds <- setdiff(EC_8X_tmp_inds, hi_rep_inds)
length(EC_8X_inds)
# 29
length(unique(phenos$ACC[EC_8X_inds]))
# 13




```



