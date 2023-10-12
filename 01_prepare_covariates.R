## #######################################################################################
##
## TB SPATIAL DATA PREP
##
## AUTHOR: Nat Henry
## CREATED: 2 September 2021
## PURPOSE: Load and prepare spatial data objects for TB modeling in Uganda
##
## #######################################################################################

# Prepared data version
versions <- list(prepped_data = '20231010')

# Load packages
load_pkgs <- c('sf', 'data.table', 'glue')
invisible(lapply(load_pkgs, library, character.only=TRUE))

# Load custom packages
repos_dir <- '~/repos/'
for(custom_pkg in c('versioning', 'pixel2poly', 'mbg')){
  devtools::load_all(file.path(repos_dir, custom_pkg))
}

config_fp <- file.path(repos_dir, 'uga_tb/config.yaml')
config <- versioning::Config$new(config_list = config_fp, versions = versions)

config$get_dir_path('prepped_data') |> dir.create()


## 01. Prepare district shapefile with covariates --------------------------------------->

model_shp <- config$read('shps', 'adm2')

# Create ID raster and aggregation table
id_raster <- pixel2poly::build_id_raster(polygons = terra::vect(model_shp))
aggregation_table <- pixel2poly::build_aggregation_table(
  polygons = terra::vect(model_shp),
  id_raster = id_raster,
  polygon_id_field = 'uid'
)

config$write(model_shp, 'prepped_data', 'model_shp')
config$write(id_raster, 'prepped_data', 'id_raster')
config$write(aggregation_table, 'prepped_data', 'link_table')

agg_cols <- config$get("shapefile_settings", "ids", "adm2")
# Load population and aggregate
population_raster <- lapply(2016:2022, function(year){
  mbg::load_covariates(
    directory = config$get_dir_path('covariates'),
    settings = config$get('pop_covariate'),
    id_raster = id_raster,
    year = year,
    file_format = config$get('covariate_settings', 'file_format'),
    add_intercept = FALSE
  )[[1]]
}) |> terra::rast()
admin_pop <- pixel2poly::aggregate_raster_to_polygons(
  data_raster = population_raster,
  aggregation_table = aggregation_table,
  aggregation_cols = agg_cols,
  method = 'sum',
  z_dimension_name = 'year',
  z_dimension = 2016:2022,
  aggregated_field = 'population'
)

# Load raster covariates and aggregate
covariates_list <- mbg::load_covariates(
  directory = config$get_dir_path('covariates'),
  settings = config$get('covariates'),
  id_raster = id_raster,
  year = config$get('year'),
  file_format = config$get('covariate_settings', 'file_format'),
  add_intercept = config$get('covariate_settings', 'add_intercept')
)
raster_cov_names <- names(config$get('covariates'))
admin_covs_list <- vector('list', length = length(covariates_list))
names(admin_covs_list) <- raster_cov_names
for(raster_cov_name in raster_cov_names){
  admin_covs_list[[raster_cov_name]] <- pixel2poly::aggregate_raster_to_polygons(
    data_raster = covariates_list[[raster_cov_name]],
    aggregation_table = aggregation_table,
    aggregation_cols = agg_cols,
    method = 'weighted.mean',
    weighting_raster = population_raster[[(2016:2022) == config$get('year')]],
    aggregated_field = raster_cov_name
  )
}
admin_covariates <- Reduce(
  f = function(x, y) merge(x=x, y=y, by = agg_cols),
  x = admin_covs_list
)
# Add intercept
admin_covariates[, intercept := 1 ]
# Log-transform nighttime lights
admin_covariates[, log_ntl := log1p(ntl) ]

# Merge on refugee data
refugee_dt <- config$read('raw_data', 'refugees')
admin_covariates[refugee_dt, refugees_pct := i.refugees_pct, on = 'uid']
admin_covariates[is.na(refugees_pct), refugees_pct := 0 ]

# Merge on HIV, household crowding, and cattle from an existing shapefile
old_shp <- config$read('old_data', 'shp')
old_shp$old_uid <- old_shp$uid
old_shp <- old_shp[, c('old_uid')]
old_covs <- config$read('old_data', 'covs')
setnames(old_covs, 'uid', 'old_uid')

centroids <- sf::st_centroid(model_shp[, c('uid')], of_largest_polygon = TRUE)
uid_match_tbl <- as.data.table(sf::st_join(x = centroids, y = old_shp))[, .(uid, old_uid)]
(admin_covariates
  [uid_match_tbl, old_uid := i.old_uid, on = 'uid']
  [
    old_covs,
    `:=` (cattle_pc = i.cattle_pc, hiv = i.hiv_test, hh_crowding = i.household_crowding),
    on = 'old_uid'
  ]
  [, old_uid := NULL ]
)

# Normalize all covariates
for(cov_name in c('tt_hcf', 'log_ntl', 'refugees_pct', 'cattle_pc', 'hiv', 'hh_crowding')){
  norm_cov <- paste0(cov_name, '_norm')
  admin_covariates[, (norm_cov) := (get(cov_name) - mean(get(cov_name))) / sd(get(cov_name))]
}

# Save admin covariates
admin_covariates <- admin_covariates[order(uid)]
config$write(admin_covariates, 'prepped_data', 'covariates')


## 02. Prepare prevalence data ---------------------------------------------------------->

# Merge on prevalence survey data by spatial location
prev_survey_raw <- config$read('raw_data', 'prev_survey')
# Convert to spatial object, then overlay with shapefile
prev_with_uids <- sf::st_as_sf(
  x = prev_survey_raw,
  coords = c('longitude', 'latitude'),
  crs = sf::st_crs(4326)
) |> 
  sf::st_join(y = model_shp) |>
  as.data.table()
# Aggregate by district
prev_agg <- (prev_with_uids
  [
  , .(sampsize = sum(num_partic), ptb_bc = sum(num_prevalent_cases)),
  by = agg_cols
  ]
  [, prev_per_100k_obs := ptb_bc / sampsize * 1e5 ]
  [admin_covariates, on = agg_cols, nomatch = NULL]
  [order(uid)]
)

config$write(prev_agg, 'prepped_data', 'prev_data')


## 03. Prepare case notifications ------------------------------------------------------->

# Load all notifications
notifs_list <- list(
  p1 = config$read('raw_data', 'notifs_p1'),
  p2 = config$read('raw_data', 'notifs_p2'),
  ped_p1 = config$read('raw_data', 'ped_notifs_p1'),
  ped_p2 = config$read('raw_data', 'ped_notifs_p2')
)
# Prepare old notifications
notifs_list <- lapply(notifs_list, melt, id.vars = 'organisationunitname')

first4 <- function(char_vec) substr(char_vec, 1, 4)
last4 <- function(char_vec) substr(char_vec, nchar(char_vec) - 3, nchar(char_vec))

# Format separately
for(grp in c('p1', 'p2')){
  dt <- copy(notifs_list[[grp]])
  dt$year <- dt$variable |> as.character() |> last4() |> as.integer()
  dt_agg <- (dt
    [, type := 'PCD' ]
    [grepl('P-BC', variable), type := 'PBC' ]
    [grepl('EPTB', variable), type := 'EPTB' ]
    [, value := nafill(value, fill = 0) ]
    [, .(value = sum(value)), by = .(year, type, organisationunitname)]
  )
  dt_wide <- dcast(dt_agg, organisationunitname + year ~ type)
  dt_wide[, value := PBC + PCD ]
  dt_wide[, pct_in_case_def := value / (PBC + PCD + EPTB) ]
  notifs_list[[grp]] <- copy(dt_wide)
}

notifs_list$ped_p1$year <- (notifs_list$ped_p1$variable |>
  as.character() |>
  first4() |>
  as.integer()
)
notifs_list$ped_p1[, value := nafill(value, fill = 0)]
notifs_list$ped_p1 <- notifs_list$ped_p1[
  , .(ped_value = sum(value)), by = .(organisationunitname, year)
]

notifs_list$ped_p2 <- (notifs_list$ped_p2
  [, year := as.integer(gsub('^y', '', variable)) ]
  [notifs_list$p2, pct_in_case_def := i.pct_in_case_def, on = c('year', 'organisationunitname')]
  [, ped_value := value * pct_in_case_def ]
  [, c('variable', 'pct_in_case_def', 'value') := NULL ]
)

notifs_combined <- merge(
  x = notifs_list[c('p1','p2')] |> rbindlist(use.names = T),
  y = notifs_list[c('ped_p1','ped_p2')] |> rbindlist(use.names = T),
  by = c('year', 'organisationunitname'),
  all = T
)[, notif_count := value - ped_value ][, c('value', 'ped_value') := NULL ]

# Standardize district names
notifs_combined$ADM2_EN <- (notifs_combined$organisationunitname |>
  gsub(pattern = ' District$', replacement = '') |>
  gsub(pattern = ' City$', replacement = '')
)
(notifs_combined
  [, organisationunitname := NULL ]
  [ADM2_EN == "Fort Portal", ADM2_EN := "Kabarole" ]
  [ADM2_EN == "Madi-Okollo", ADM2_EN := 'Madi Okollo' ]
  [ADM2_EN == "Sembabule", ADM2_EN := 'Ssembabule' ]
)
notifs_agg <- (notifs_combined
  [
    , lapply(.SD, sum, na.rm=T),
    .SDcols = c('notif_count', 'PBC', 'PCD', 'EPTB'),
    by = .(ADM2_EN, year)]
  [order(ADM2_EN, year)]
)
notifs_final <- merge(
  x = admin_covariates[, .(uid, ADM1_EN, ADM2_EN, tt_hcf, tt_hcf_norm)],
  y = notifs_agg,
  by = 'ADM2_EN',
  all = TRUE
)

# Merge on over-15 population as the denominator
over_15s <- config$read('raw_data', 'prop_over_15')
admin_pop[over_15s, prop_over_15 := i.prop_over_15, on = 'year']
admin_pop[, pop_over_15 := population * prop_over_15 ][, prop_over_15 := NULL ]
notifs_final[admin_pop, pop_over_15 := i.pop_over_15, on = c('uid', 'year')]

# Save notifications
config$write(notifs_final, 'prepped_data', 'notif_data')
