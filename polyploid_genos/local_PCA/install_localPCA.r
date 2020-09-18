# INSTALLING local PCA; library = lostruct

# cd /global/homes/g/grabowsp/tools
# module load python/3.7-anaconda-2019.07
# source activate local_PCA

#install.packages("data.table")
library(data.table)
Sys.setenv(TAR= '/bin/tar')
devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)

#############3

# Load local PCA

# cd /global/homes/g/grabowsp/tools
# module load python/3.7-anaconda-2019.07
# source activate local_PCA

library(data.table)
library(lostruct)



