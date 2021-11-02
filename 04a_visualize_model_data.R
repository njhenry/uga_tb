## #######################################################################################
##
## VISUALIZE MODEL INPUTS
##
## AUTHOR: Nat Henry
## CREATED: 1 September 2021
## PURPOSE: Visualize input regions, populations, and case notifications
##
## #######################################################################################

load_libs <- c('data.table', 'sf', 'ggplot2', 'glue')
invisible(lapply(load_libs, library, character.only=TRUE))

# Set filepaths
work_dir <- '{REDACTED}'
in_dir <- glue("{work_dir}/raw_data")
out_dir <- glue("{work_dir}/viz")
dir.create(out_dir, showWarnings = FALSE)

# Load datasets
notifs <- fread(file.path(in_dir, 'district_notifs_2019_2020_combined.csv'))
rhosp_matching <- fread(file.path(in_dir, 'referral_hospital_district_matching.csv'))
rhosp_gps <- fread(file.path(in_dir, 'referral_hospital_geolocation.csv'))
tb_prev_surv_data <- fread(
  file.path(in_dir, 'prev_survey/UGA_TB_prev_survey_appendices_all.csv')
)
# Load district shapefile
ad2_sf <- sf::st_read(file.path(in_dir, 'UGA_consistent_2019_2020/DISTRICTS_2018_UTM_36N.shp'))
ad2_sf$DName2017 <- ad2_sf$DNama2017
ad2_sf$DNama2017 <- NULL

# TEST SHAPEFILE MERGE
ad2_sf_data <- as.data.table(ad2_sf)[, .(F15Regions, DName2016, DName2017, DName2018)]
notifs[, DName2018 := toupper(district)]
test_merge <- merge(x=ad2_sf_data, y=notifs, by='DName2018', all=TRUE)

## MATCH DATA WITH SHAPEFILE
sf_distmatch <- x
sf_rrh <- merge(x=ad2_sf, y=rhosp_matching, by.x="DName2018", by.y='dname_upper', all=TRUE)
