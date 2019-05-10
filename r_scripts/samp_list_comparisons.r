# Comparison of Reseq and Exome-capture sample lists

reseq_samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/Pvirgatum_921samples_ploidy_FINAL.txt'

reseq_samp_meta <- read.table(reseq_samp_meta_file, header = T, sep = '\t', 
  stringsAsFactors = F, comment.char = '@')

j_res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/PVDIV_Master_Metadata_File_5-1-2019.txt'

j_res_meta <- read.table(j_res_file, header = T, sep = '\t', 
  stringsAsFactors = F, quote = '', comment.char = '$')

# I know that this sample name is incorrect
reseq_samp_meta$'Plant.ID'[
  which(reseq_samp_meta$'Plant.ID' == 'J664.T')] <- 'J664.E'

reseq_prob_inds <- which(
  (reseq_samp_meta$'Plant.ID' %in% j_res_meta$PLANT_ID) == F)

prob_name_list <- strsplit(reseq_samp_meta$'Plant.ID'[reseq_prob_inds], 
  split = '(', fixed = T)
prob_name_fixed_vec <- unlist(lapply(prob_name_list, 
  function(x) gsub(')', '', x[2], fixed = T)))

reseq_samp_meta$'Plant.ID'[reseq_prob_inds] <- prob_name_fixed_vec

sum((reseq_samp_meta$'Plant.ID' %in% j_res_meta$PLANT_ID) == F)
# 0

sum(duplicated(j_res_meta$PLANT_ID))
#0

reseq_samp_meta$j_source_1 <- NA
reseq_samp_meta$j_source_2 <- NA
reseq_samp_meta$j_lat <- NA
reseq_samp_meta$j_long <- NA

for(ri in seq(nrow(reseq_samp_meta))){
  j_ind <- which(j_res_meta$PLANT_ID == reseq_samp_meta$'Plant.ID'[ri])
  reseq_samp_meta[ri, c('j_source_1', 'j_source_2', 'j_lat', 'j_long')] <- (
    j_res_meta[j_ind, c('SOURCE_ID_1', 'SOURCE_ID_2', 'LATITUDE', 
    'LONGITUDE')])
}


#reseq_dup_samps <- which(duplicated(reseq_samp_meta_0[,5]))
#reseq_samp_meta <- reseq_samp_meta_0[-reseq_dup_samps, ]

reseq_samp_meta$LATITUDE[which(reseq_samp_meta$LATITUDE == '???')] <- NA
reseq_samp_meta$LONGITUDE[which(reseq_samp_meta$LONGITUDE == '???')] <- NA

reseq_samp_meta$LATITUDE <- as.numeric(reseq_samp_meta$LATITUDE)
reseq_samp_meta$LONGITUDE <- as.numeric(reseq_samp_meta$LONGITUDE)

reseq_samp_meta$exome_name_pop <- NA

exome_pop_meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/exome/metadata/combo_pop_metadata_v1.4.txt'
exome_pop_meta_0 <- read.table(exome_pop_meta_file, header = T, sep = '\t', 
  stringsAsFactors = F)

exome_pop_dup <- which(duplicated(exome_pop_meta_0[,1]))

exome_pop_meta <- exome_pop_meta_0[-exome_pop_dup, ]

exome_samp_meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/exome/metadata/combo_samp_metadata_v1.9.txt'
exome_samp_meta_0 <- read.table(exome_samp_meta_file, header = T, sep = '\t',
  stringsAsFactors = F)

exome_samp_meta <- exome_samp_meta_0[
  -c(which(exome_samp_meta_0$sample_status == 'bad'), which(exome_samp_meta_0$sample_status == 'clone')), ]

exome_samp_meta$v2_pop_name[
  which(exome_samp_meta$v2_pop_name == 'Summer')] <- 'Summer.1'
exome_samp_meta$v2_pop_name[
  which(exome_samp_meta$v2_pop_name == 'Ellsworth')] <- 'Ellsworth.1'


# Match exome population info to reseq sample names
match_list <- list()

for(i in seq(nrow(exome_pop_meta))){
  alias_vec <- gsub(' ', '',
    unlist(strsplit(exome_pop_meta$alias[i], split = ';')))
  tmp_match_vec <- c()
  for(j in seq(length(alias_vec))){
    tmp_matches <- c(
      grep(alias_vec[j], gsub(' ', '', reseq_samp_meta$'SampleName..onPlate.'),
        fixed = T), 
      grep(alias_vec[j], gsub(' ', '', reseq_samp_meta$'Plant.ID'), fixed = T),
      grep(alias_vec[j], gsub(' ', '', reseq_samp_meta$'j_source_1'), 
        fixed = T),
      grep(alias_vec[j], gsub(' ', '', reseq_samp_meta$'j_source_2'), 
        fixed = T))
    tmp_match_vec <- c(tmp_match_vec, tmp_matches)
  }
  match_list[[i]] <- unique(tmp_match_vec)
}

for(i in seq(length(match_list))){
  if(length(match_list[[i]]) > 0){
    for(j in seq(length(match_list[[i]]))){
      reseq_ind <- match_list[[i]][j]
      reseq_samp_meta$exome_name_pop[reseq_ind] <- (
        exome_pop_meta$info_pop_name[i])
    }
  }
}

# 170 samps
#

# I know that there is a Carthage sample...
reseq_samp_meta$exome_name_pop[grep('carthage', 
  reseq_samp_meta$'SampleName..onPlate.')] <- 'Carthage' 
# match reseq samples based on their sample names
reseq_samp_meta$samp_based_pop <- NA

j_inds <- grep('^J', reseq_samp_meta$'Plant.ID')

reseq_samp_meta$samp_based_pop[j_inds] <- substr(
  reseq_samp_meta$'Plant.ID'[j_inds], start = 1, stop = 4)

for(i in seq(nrow(reseq_samp_meta))){
  if(is.na(reseq_samp_meta$exome_name_pop[i]) != T){
    tmp_jpop <- reseq_samp_meta$samp_based_pop[i]
    tmp_matches <- which(reseq_samp_meta$samp_based_pop == tmp_jpop)
    for(j in tmp_matches){
      if(is.na(reseq_samp_meta$exome_name_pop[j])){
        reseq_samp_meta$exome_name_pop[j] <- reseq_samp_meta$exome_name_pop[i]
      }
    }
  }
}
# add pop names to an additional 1 samps
#########

# use geographic info to find matches
#  note: rounding to 1 decimal place seems to work the best
# plan:
# search for match
# if find match, then add as geo pop
# if find multiple matches, then paste them togehter

exome_geostring_1 <- paste(round(exome_pop_meta$lat, digits = 1), 
  round(exome_pop_meta$long, digits = 1), sep = '_')

reseq_geostring_1 <- paste(round(reseq_samp_meta$LATITUDE, digits = 1), 
  round(reseq_samp_meta$LONGITUDE, digits = 1), sep = '_')

reseq_geostring_j <- paste(round(reseq_samp_meta$j_lat, digits = 1),
  round(reseq_samp_meta$j_long, digits = 1), sep = '_')

reseq_samp_meta$exome_geo_pop <- NA
reseq_samp_meta$j_exome_geo_pop <- NA

for(i in seq(nrow(reseq_samp_meta))){
  tmp_geo_match <- which(exome_geostring_1 == reseq_geostring_1[i])
  if(reseq_geostring_1[i] == 'NA_NA'){tmp_geo_match <- c()}
  if(length(tmp_geo_match) > 0){
    reseq_samp_meta$exome_geo_pop[i] <- paste(
      exome_pop_meta$info_pop_name[tmp_geo_match], collapse = '_')
  }
}

for(i in seq(nrow(reseq_samp_meta))){
  tmp_geo_match <- which(exome_geostring_1 == reseq_geostring_j[i])
  if(reseq_geostring_j[i] == 'NA_NA'){tmp_geo_match <- c()}
  if(length(tmp_geo_match) > 0){
    reseq_samp_meta$j_exome_geo_pop[i] <- paste(
      exome_pop_meta$info_pop_name[tmp_geo_match], collapse = '_')
  }
}


sum(unique(exome_geostring_1) %in% reseq_geostring_1)
# 50 of 127 unique geographic locations

sum(unique(exome_geostring_1) %in% reseq_geostring_j)
# 50

sum(reseq_geostring_1 %in% exome_geostring_1)
# 229 of 921 total samples
# 50 of 306 unique geographic locations

sum(reseq_geostring_j %in% exome_geostring_1)
# 218

# CONTINUE FROM HERE INCORPORATING j_meta GEO info


reseq_samp_meta$exome_pop <- reseq_samp_meta$exome_name_pop
add_geo_pop_inds <- intersect(which(is.na(reseq_samp_meta$exome_pop)), 
  which(is.na(reseq_samp_meta$exome_geo_pop) == F))

reseq_samp_meta$exome_pop[add_geo_pop_inds] <- reseq_samp_meta$exome_geo_pop[
  add_geo_pop_inds]

nrow(reseq_samp_meta) - sum(is.na(reseq_samp_meta$exome_pop))
# 222 resequencing samples have Name or Geographic connection to exome data

add_j_geo_inds <- intersect(which(is.na(reseq_samp_meta$exome_pop)),
  which(is.na(reseq_samp_meta$j_exome_geo_pop) == F))

reseq_samp_meta$exome_pop[add_j_geo_inds] <- (
  reseq_samp_meta$j_exome_geo_pop[add_j_geo_inds])

nrow(reseq_samp_meta) - sum(is.na(reseq_samp_meta$exome_pop))
# 226 resequencing samples have Name or Geographic connection to exome data

#####
exome_pop_meta_0$reseq_name_match <- F

for(i in seq(nrow(exome_pop_meta_0))){
  tmp_matches <- grep(exome_pop_meta_0$info_pop_name[i], 
    reseq_samp_meta$exome_name_pop)
  if(length(tmp_matches) > 0){exome_pop_meta_0$reseq_name_match[i] <- T}
}

exome_pop_meta_0$reseq_geo_match <- F

for(i in seq(nrow(exome_pop_meta_0))){
  tmp_matches <- grep(exome_pop_meta_0$info_pop_name[i],
    reseq_samp_meta$exome_geo_pop)
  if(length(tmp_matches) > 0){exome_pop_meta_0$reseq_geo_match[i] <- T}
}

exome_pop_meta_0$reseq_any_match <- F

for(i in seq(nrow(exome_pop_meta_0))){
  tmp_matches <- grep(exome_pop_meta_0$info_pop_name[i],
    reseq_samp_meta$exome_pop)
  if(length(tmp_matches) > 0){exome_pop_meta_0$reseq_any_match[i] <- T}
}

######


exome_samp_meta$reseq_name_match <- F
exome_samp_meta$reseq_geo_match <- F
exome_samp_meta$reseq_any_match <- F

for(i in seq(nrow(exome_samp_meta))){
  tmp_pop_ind <- which(exome_pop_meta_0$v2_pop_name == 
    exome_samp_meta$v2_pop_name[i])
  exome_samp_meta[i, 
    c('reseq_name_match', 'reseq_geo_match', 'reseq_any_match')] <- (
    exome_pop_meta_0[tmp_pop_ind, 
    c('reseq_name_match', 'reseq_geo_match', 'reseq_any_match')])
}

#######

exome_pop_meta$reseq_name_match <- F
exome_pop_meta$reseq_geo_match <- F
exome_pop_meta$reseq_any_match <- F

exome_pop_meta$reseq_name_match[
  exome_pop_meta$info_pop_name %in% unique(
  exome_samp_meta$info_pop[exome_samp_meta$reseq_name_match])] <- T

exome_pop_meta$reseq_geo_match[
  exome_pop_meta$info_pop_name %in% unique(
  exome_samp_meta$info_pop[exome_samp_meta$reseq_geo_match])] <- T

exome_pop_meta$reseq_any_match[
  exome_pop_meta$info_pop_name %in% unique(
  exome_samp_meta$info_pop[exome_samp_meta$reseq_any_match])] <- T

######
# Generate figure showing exome populations

library(mapproj)

# exome populations
exome_pop_map_fig <- paste(
  '/home/t4c1/WORK/grabowsk/data/switchgrass/samp_comp_figs/',
  'exome_pop_comparison.pdf', sep = '')

pdf(exome_pop_map_fig)
map('usa', project = 'albers', par = c(39,45), col = 'gray85', fill = T,
  orientation = c(90, 0, 265), border = 'black')
exome_coordss <- mapproject(exome_pop_meta$long, exome_pop_meta$lat, 
  par = c(39,45), orientation = c(90,0,265))

exome_col_vec <- rep('red2', times = nrow(exome_pop_meta))
exome_col_vec[exome_pop_meta$reseq_any_match] <- 'gray15'
points(exome_coordss, pch = 19, cex = 1, col = exome_col_vec )

dev.off()

# reseq populations
reseq_pop_map_fig <- paste(
  '/home/t4c1/WORK/grabowsk/data/switchgrass/samp_comp_figs/',
  'reseq_pop_comparison.pdf', sep = '')

pdf(reseq_pop_map_fig)
map('usa', project = 'albers', par = c(39,45), col = 'gray85', fill = T,
  orientation = c(90, 0, 265), border = 'black')
reseq_coordss <- mapproject(reseq_samp_meta$LONGITUDE, 
  reseq_samp_meta$LATITUDE, par = c(39,45), orientation = c(90,0,265))

reseq_col_vec <- rep('gray15', times = nrow(reseq_samp_meta))
reseq_col_vec[is.na(reseq_samp_meta$exome_pop)] <- 'blue2'
points(reseq_coordss, pch = 19, cex = 1, col = reseq_col_vec )

dev.off()

# all populations (with GPS coordinates) from both sets
bothsets_pop_map_fig <- paste(
  '/home/t4c1/WORK/grabowsk/data/switchgrass/samp_comp_figs/',
  'bothsets_pop_comparison.pdf', sep = '')

pdf(bothsets_pop_map_fig)
map('usa', project = 'albers', par = c(39,45), col = 'gray85', fill = T,
  orientation = c(90, 0, 265), border = 'black')

reseq_only_coordss <- mapproject(
  reseq_samp_meta$LONGITUDE[is.na(reseq_samp_meta$exome_pop)], 
  reseq_samp_meta$LATITUDE[is.na(reseq_samp_meta$exome_pop)], 
  par = c(39,45), orientation = c(90,0,265))

points(reseq_only_coordss, pch = 19, cex = 1, col = 'blue2')
points(exome_coordss, pch = 19, cex = 1, col = exome_col_vec )

dev.off()


# Next: 

# include ploidy
# filled square (15) = shared tetraploid
# empty square (0) = not-shared tetraploid
# filled circle (19) = shared octoploid
# empty circle (1) = not-shared octoploid

reseq_tet_inds <- which(reseq_samp_meta$Consensus_ploidy == '4x')
reseq_oct_inds <- which(reseq_samp_meta$Consensus_ploidy == '8x')
reseq_na_ploid <- which(is.na(reseq_samp_meta$Consensus_ploidy))
reseq_tmp_tet <- intersect(reseq_na_ploid, 
  which(reseq_samp_meta$'Cyto.ploidy..TOM_SHEET.' == '4x'))
reseq_tet_tot <- c(reseq_tet_inds, reseq_tmp_tet)

reseq_notshared_inds <- which(is.na(reseq_samp_meta$exome_pop))
reseq_shared_inds <- setdiff(seq(nrow(reseq_samp_meta)), reseq_notshared_inds)

reseq_pch <- rep(NA, times = nrow(reseq_samp_meta))
reseq_pch[intersect(reseq_tet_tot, reseq_shared_inds)] <- 15
reseq_pch[intersect(reseq_tet_tot, reseq_notshared_inds)] <- 0
reseq_pch[intersect(reseq_oct_inds, reseq_shared_inds)] <- 19
reseq_pch[intersect(reseq_oct_inds, reseq_notshared_inds)] <- 1

exome_tet_inds <- which(exome_samp_meta$MASNP_3read_ploidy == 4)
exome_oct_inds <- which(exome_samp_meta$MASNP_3read_ploidy == 8)
exome_shared_inds <- which(exome_samp_meta$reseq_any_match)
exome_notshared_inds <- which(exome_samp_meta$reseq_any_match == F)

exome_pch <- rep(NA, times = nrow(exome_samp_meta))
exome_pch[intersect(exome_tet_inds, exome_shared_inds)] <- 15
exome_pch[intersect(exome_tet_inds, exome_notshared_inds)] <- 0
exome_pch[intersect(exome_oct_inds, exome_shared_inds)] <- 19
exome_pch[intersect(exome_oct_inds, exome_notshared_inds)] <- 1

bothsets_samp_map_fig <- paste(
  '/home/t4c1/WORK/grabowsk/data/switchgrass/samp_comp_figs/',
  'bothsets_samp_comparison_map_with_ploidy.pdf', sep = '')

pdf(bothsets_samp_map_fig)
map('usa', project = 'albers', par = c(39,45), col = 'gray85', fill = T,
  orientation = c(90, 0, 265), border = 'black')
exome_samp_coordss <- mapproject(exome_samp_meta$long, exome_samp_meta$lat,
  par = c(39,45), orientation = c(90,0,265))

reseq_coordss <- mapproject(reseq_samp_meta$LONGITUDE,
  reseq_samp_meta$LATITUDE, par = c(39,45), orientation = c(90,0,265))

points(reseq_coordss, pch = reseq_pch, cex = 1.25, col = 'blue2')
points(exome_samp_coordss, pch = exome_pch, cex = 1, col = 'red2')

dev.off()


# compare ploidy
# go by 'info_pop_name' at top of loop
# look for any reseq samples that match
# check reseq consensus ploidy
# get v2_pop_names in exome data
# for each v2_pop_name, get samp inds and thier inferred ploidy

for(ep in unique(exome_pop_meta_0$info_pop_name)){
  tmp_r_inds <- grep(ep, reseq_samp_meta$exome_pop)
  if(length(tmp_r_inds) > 0){
    print(ep)
    tmp_r_ploidy <- unique(reseq_samp_meta$Consensus_ploidy[tmp_r_inds])
    print(tmp_r_ploidy)
    tmp_e_ploidy <- c()
    tmp_e_pops <- exome_pop_meta_0$v2_pop_name[
      exome_pop_meta_0$info_pop_name == ep]
    for(mep in tmp_e_pops){
      tmp_e_inds <- which(exome_samp_meta$v2_pop_name == mep)
     tmp_e_ploidy <- c(tmp_e_ploidy, 
       exome_samp_meta$MASNP_3read_ploidy[tmp_e_inds])
    }
    tmp_e_ploidy <- unique(tmp_e_ploidy)
    print(tmp_e_ploidy)
  }  

}

# For the most part, the ploidy info matches between the two sets. Just
#   using the eye-ball test, I think some of the reseq ploidy determinations
#   are incorrect (for example, the Argentina and some Summer calls)

# Next:
### 1) add info to exome_pop_meta_0: name match, geo match, any match,
# 2) figure out how many exome samples overlap or don't overlap
# ...add name, geo, any-match columns to exome sample metadata
# 3) Make list of overlapping populations
# 4) Make list including geographic info
# 5) Try to match ploidy info between matches
# 6) Make maps showing overlap and non-overlap between sample sets

# Caluclate statistics:
length(unique(exome_pop_meta_0$info_pop_name))
# 140
length(unique(
  exome_pop_meta_0$info_pop_name[exome_pop_meta_0$reseq_name_match]))
# 66
length(unique(
  exome_pop_meta_0$info_pop_name[exome_pop_meta_0$reseq_geo_match]))
# 56
length(unique(
  exome_pop_meta_0$info_pop_name[exome_pop_meta_0$reseq_any_match]))
# 75

length(unique(reseq_geostring_1))
# 305 (actually 306, but one is NA_NA)
length(unique(
  reseq_geostring_1[which(is.na(reseq_samp_meta$exome_name_pop) == F)]))
# 59 (actually 60, but one is NA_NA) 
length(unique(
  reseq_geostring_1[which(is.na(reseq_samp_meta$exome_geo_pop) == F)]))
# 49
length(unique(
  reseq_geostring_1[which(is.na(reseq_samp_meta$exome_pop) == F)]))
# 68

nrow(exome_samp_meta)
# 1094
sum(exome_samp_meta$reseq_name_match)
# 670
sum(exome_samp_meta$reseq_geo_match)
# 507
sum(exome_samp_meta$reseq_any_match)
# 726

nrow(reseq_samp_meta)
# 921
length(which(is.na(reseq_samp_meta$exome_name_pop) == F))
# 84
length(which(is.na(reseq_samp_meta$exome_geo_pop) == F))
# 171
length(which(is.na(reseq_samp_meta$exome_pop) == F))
# 199

# Generate lists of populations:
exome_pop_meta$info_pop_name %in% unique(
  exome_samp_meta$info_pop[exome_samp_meta$reseq_name_match])

exome_v1name_share_pops <- unique(
  exome_samp_meta$info_pop[exome_samp_meta$reseq_name_match])
exome_v2name_share_pops <- unique(
  exome_samp_meta$v2_pop_name[exome_samp_meta$reseq_name_match])

# add exome ecotype and genepool info to reseq samples:
# plan:
# for pop match(es): 
#  1) find v2 pop names in exome_pop_meta_0
#  2) tally ecotype designations (upland, lowland, admixed) from 
#       exome_samp_meta that have the v2 pop name(s)
#  3) tally genepool designations from exome_samp_meta that have the v2 pop
#       names

reseq_samp_meta$exome_ecotype_info <- NA
reseq_samp_meta$exome_ecotype <- NA

reseq_samp_meta$exome_gp_info <- NA
reseq_samp_meta$exome_gp <- NA

reseq_samp_meta$exome_ploidy_info <- NA
reseq_samp_meta$exome_ploidy <- NA

rem_inds <- which(is.na(reseq_samp_meta$exome_pop) == F)

good_geo_exome_pops <- grep('core|permissive', exome_pop_meta_0$geo_use)

for(rem in rem_inds){
  info_m_name <- reseq_samp_meta$exome_pop[rem]
  ep_inds <- intersect(grep(info_m_name, exome_pop_meta_0$info_pop_name), 
    good_geo_exome_pops)
  es_inds <- c()
  for(ep in ep_inds){
    tmp_es_match <- which(exome_samp_meta$v2_pop_name == 
      exome_pop_meta_0$v2_pop_name[ep])
    es_inds <- c(es_inds, tmp_es_match)
  }
  # ecotype
  tmp_eco_tab <- table(exome_samp_meta$geo_ecotype[es_inds])
  if(length(tmp_eco_tab) == 1){
    reseq_samp_meta$exome_ecotype[rem] <- names(tmp_eco_tab)
  }
  if(length(tmp_eco_tab) > 1){
    reseq_samp_meta$exome_ecotype[rem] <- 'MIXED'
  }
  tmp_u <- tmp_eco_tab[names(tmp_eco_tab) == 'upland']
  if(length(tmp_u) == 0){tmp_u <- 0}
  tmp_l <- tmp_eco_tab[names(tmp_eco_tab) == 'lowland']
  if(length(tmp_l) == 0){tmp_l <- 0}
  tmp_a <- tmp_eco_tab[names(tmp_eco_tab) == 'ecotype_admixed']
  if(length(tmp_a) == 0){tmp_a <- 0}
  tmp_eco_info <- paste('L', tmp_l, '__U', tmp_u, '__A', tmp_a, sep = '.')
  reseq_samp_meta$exome_ecotype_info[rem] <- tmp_eco_info
  # gene pool
  tmp_gp_tab <- table(exome_samp_meta$geo_genepool[es_inds])
  if(length(tmp_gp_tab) == 1){
    reseq_samp_meta$exome_gp[rem] <- names(tmp_gp_tab)
  }
  if(length(tmp_gp_tab) > 1){
    reseq_samp_meta$exome_gp[rem] <- 'MIXED'
  }
  tmp_gp_info <- paste(names(tmp_gp_tab), tmp_gp_tab, sep = '.', 
    collapse = '__')
  reseq_samp_meta$exome_gp_info[rem] <- tmp_gp_info
  # ploidy
  tmp_ploidy_tab <- table(exome_samp_meta$MASNP_3read_ploidy[es_inds])
  if(length(tmp_ploidy_tab) == 1){
    reseq_samp_meta$exome_ploidy[rem] <- names(tmp_ploidy_tab)
  }
  if(length(tmp_ploidy_tab) > 1){
    reseq_samp_meta$exome_ploidy[rem] <- 'MIXED'
  }
  tmp_ploidy_info <- paste(names(tmp_ploidy_tab), tmp_ploidy_tab, sep = '.', 
    collapse = '__')
  reseq_samp_meta$exome_ploidy_info[rem] <- tmp_ploidy_info
}

reseq_with_exome_data_meta_out <- '/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/sg_reseq_metadata_with_exome_info.txt'

write.table(reseq_samp_meta, file = reseq_with_exome_data_meta_out, quote = F,
  sep = '\t', row.names = F, col.names = T)

quit(save = 'no')
