# Generate some maps with pie charts to visualize structure results

# bash
# source activate R_maps

### LOAD PACKAGES ###
library(ggplot2)
library(mapdata)
library(scatterpie)

### LOAD INPUTS ###
res_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/up_geo/upgeo_structure_combined_results.txt'
res <- read.table(res_file, header = T, sep = '\t', stringsAsFactors = F)

### SET OUTPUTS ###
k_val = 3
out_string <- paste('k', k_val, '_pie_map.pdf', sep = '')
fig_out <- gsub('combined_results.txt', out_string, res_file)

### SET VARIABLES ###
## plot title variables
samp_set_name <- 'Upland Samples'
fig_title <- paste(samp_set_name, ', K=', k_val, sep = '')

## result variables
# result columns to include
column_vec <- c('K3_1', 'K3_2', 'K3_3')

# colors of groups
k_col_vec <- c()
k_col_vec['K3_1'] <- 'blue2'
k_col_vec['K3_2'] <- 'steelblue2'
k_col_vec['K3_3'] <- 'darkorchid2'

# labels of groups
k_lab_vec <- c()
k_lab_vec['K3_1'] <- '4X Dakota'
k_lab_vec['K3_2'] <- '4X Great Lakes'
k_lab_vec['K3_3'] <- 'Mainly 8X'

# line color for ploidy
p_col_vec <- c()
p_col_vec['4X'] <- 'black'
p_col_vec['8X'] <- 'gray80'

p_lab_vec <- c()
p_lab_vec['4X'] <- '4X'
p_lab_vec['8X'] <- '8X'

## map size variables
map_x_lims <- c(-105, -67)
map_y_lims <- c(25,50)

## plotting variables
fig_base_size = 8
fig_width_adj <- 1

x_dist <- map_x_lims[2]-map_x_lims[1]
y_dist <- map_y_lims[2]-map_y_lims[1]

tot_fig_width <- fig_base_size * (x_dist/y_dist) * fig_width_adj

#######################

# remove duplicate coordinates so all samples plot separately
res$coord <- paste(res$LATITUDE, res$LONGITUDE, sep = '')
dup_coord <- which(duplicated(res$coord))
res$LATITUDE[dup_coord] <- jitter(res$LATITUDE[dup_coord], 2)
res$LONGITUDE[dup_coord] <- jitter(res$LONGITUDE[dup_coord], 2)

#res_4X <- res[which(res$PLOIDY_v1 == '4X'), ]
res_8X <- res[which(res$PLOIDY_v1 == '8X'), ]

state_map <- map_data('state')
world_map <- map_data('worldHires')
mexico_map <- subset(world_map, region == 'Mexico')

gg_map <- ggplot(data = state_map) + 
  geom_polygon(data = state_map, aes(x = long, y = lat, group = group), 
    color = 'white', fill = 'gray50') +
  geom_polygon(data = mexico_map, aes(x = long, y = lat, group = group),
    color = 'white', fill = 'gray50') +
  geom_scatterpie(data = res, aes(x = LONGITUDE, y = LATITUDE, group = LIB, 
    color = PLOIDY_v1), cols = column_vec) +
  geom_scatterpie(data = res_8X, aes(x = LONGITUDE, y = LATITUDE, group = LIB),     cols = column_vec, color = 'gray80') +
  scale_fill_manual(name = 'Group', values = k_col_vec, labels = k_lab_vec) + 
  scale_colour_manual(name = 'Ploidy', values = p_col_vec, labels = p_lab_vec) +
  coord_fixed(xlim = map_x_lims, ylim = map_y_lims, ratio = 1.1) +
  ggtitle(fig_title)


pdf(fig_out, width = tot_fig_width, height = fig_base_size)
gg_map
dev.off()

quit(save = 'no')

