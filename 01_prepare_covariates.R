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
population_raster <- mbg::load_covariates(
  directory = config$get_dir_path('covariates'),
  settings = config$get('pop_covariate'),
  id_raster = id_raster,
  year = config$get('year'),
  file_format = config$get('covariate_settings', 'file_format'),
  add_intercept = FALSE
)[[1]]
admin_pop <- pixel2poly::aggregate_raster_to_polygons(
  data_raster = population_raster,
  aggregation_table = aggregation_table,
  aggregation_cols = agg_cols,
  method = 'sum',
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
    weighting_raster = population_raster,
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

# TODO