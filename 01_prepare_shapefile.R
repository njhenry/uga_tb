## #######################################################################################
##
## PREPARE SHAPEFILE FOR MERGE
##
## AUTHOR: Nat Henry
## CREATED: 20 September 2021
## PURPOSE: Prepare Uganda district-level shapefile for dissolve operation
##
## #######################################################################################

library(data.table); library(sf)

work_dir <- '{REDACTED}'
data_dir <- file.path(work_dir, 'raw_data')

shp <- sf::st_read(file.path(data_dir, 'UGA_consistent_2019_2020/DISTRICTS_2018_UTM_36N.shp'))
colnames(shp) <- tolower(colnames(shp))
colnames(shp)[colnames(shp)=='dnama2017'] <- 'dname2017'

# Add unique identifier column to create merged shapefile
uid_id_dt <- data.table(dname2017 = shp$dname2017, dname2018 = shp$dname2018)
collapse_fun <- function(vec) paste(sort(unique(vec)), collapse=', ')
uid_id_dt <- uid_id_dt[, .(merged2018 = collapse_fun(dname2018)), by=dname2017]
uid_id_dt[, uid := .I ]
shp_with_merge <- merge(x=shp, y=uid_id_dt, by='dname2017')

# Save for dissolve in QGIS
sf::st_write(
  shp_with_merge, dsn=file.path(data_dir, 'UGA_shp_consistent_uids.shp'), delete_layer=T
)
