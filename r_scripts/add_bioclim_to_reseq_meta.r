# Script for adding bioclimate info to the resequencing data for Sujan

# LOAD PACKAGES
library(sp)
library(raster)

# LOAD DATA
reseq_meta_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/reseq/metadata/sg_reseq_metadata_with_exome_info.txt'
reseq_samp_meta <- read.table(reseq_meta_file, header = T, sep = '\t', 
  stringsAsFactors = F, comment.char = '@')

# For v1.4 bioclim data

bioclim_var <- c(
#  alt = 'Altitude',
  bio_1 = 'Annual Mean Temperature',
  bio_2 = 'Mean Diurnal Range',
  bio_3 = 'Isothermality',
  bio_4 = 'Temperature Seasonality',
  bio_5 = 'Max Temperature of Warmest Month',
  bio_6 = 'Min Temperature of Coldest Month',
  bio_7 = 'Temperature Annual Range',
  bio_8 = 'Mean Temperature of Wettest Quarter',
  bio_9 = 'Mean Temperature of Driest Quarter',
  bio_10 = 'Mean Temperature of Warmest Quarter',
  bio_11 = 'Mean Temperature of Coldest Quarter',
  bio_12 = 'Annual Precipitation',
  bio_13 = 'Precipitation of Wettest Month',
  bio_14 = 'Precipitation of Driest Month',
  bio_15 = 'Precipitation Seasonality',
  bio_16 = 'Precipitation of Wettest Quarter',
  bio_17 = 'Precipitation of Driest Quarter',
  bio_18 = 'Precipitation of Warmest Quarter',
  bio_19 = 'Precipitation of Coldest Quarter'
)

bioclim_var_full <- bioclim_var
bioclim_var <- gsub(' of', '', bioclim_var)
bioclim_var <- gsub(' ', '.', bioclim_var)
bioclim_var <- gsub('Month', 'M', bioclim_var)
bioclim_var <- gsub('Quarter', 'Q', bioclim_var)
bioclim_var <- gsub('Annual', 'A', bioclim_var)
bioclim_var <- gsub('Temperature', 'Temp', bioclim_var)
bioclim_var <- gsub('Precipitation', 'Prec', bioclim_var)
bioclim_var <- gsub('Wettest', 'Wet', bioclim_var)
bioclim_var <- gsub('Driest', 'Dry', bioclim_var)
bioclim_var <- gsub('Coldest', 'Cold', bioclim_var)
bioclim_var <- gsub('Warmest', 'Warm', bioclim_var)
bioclim_var <- gsub('Seasonality', 'Seas', bioclim_var)

data_dir <- '/home/t4c1/WORK/grabowsk/data/bioclim/bioclim_v1.4'

bioclim_data <- list()
for(bcv in names(bioclim_var) ){
  bioclim_data[[ bcv ]] <- raster( sprintf('%s/%s.bil', data_dir, bcv) )
}
names(bioclim_data) <- bioclim_var

for(bcd in names(bioclim_data)){
    dataType(bioclim_data[[bcd]]) <- 'INT2S'
}


# get bioclim data
for( bv in bioclim_var ){
  reseq_samp_meta[[ bv ]] <- extract( bioclim_data[[bv]],  
    reseq_samp_meta[, c('LONGITUDE','LATITUDE')] )
  ## Divide temperature by 10
  if(grepl('Temp', bv )) {
     reseq_samp_meta[[ bv ]] <-  reseq_samp_meta[[ bv ]] / 10
  }
}

reseq_samp_meta$Mean.Diurnal.Range <- reseq_samp_meta$Mean.Diurnal.Range / 10

# There are several spots where the lat/long is a spot without climate data, so need to adjust those measurements

ntbm_inds <- which(reseq_samp_meta$Prec.Wet.Q == -9999)

unique(paste(reseq_samp_meta$LATITUDE[ntbm_inds], 
  reseq_samp_meta$LONGITUDE[ntbm_inds], sep = '_')

# next steps: adjust lat/long so get true climate measures

quit(save = 'no')
