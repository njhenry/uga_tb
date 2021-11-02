## #######################################################################################
##
## TB SPATIAL DATA PREP
##
## AUTHOR: Nat Henry
## CREATED: 2 September 2021
## PURPOSE: Load and prepare spatial data objects for TB modeling in Uganda
##
## #######################################################################################

# Load libraries
load_libs <- c(
  'sf','sp','rgdal','raster','data.table','glue','lbd.mbg','ggplot2','grid','gridExtra'
)
invisible(lapply(load_libs, library, character.only=TRUE))

# Settings
data_version <- '20210930'

# Set input and output filepaths
work_dir <- '{REDACTED}'
data_dir <- file.path(work_dir, 'raw_data')
prepped_dir <- file.path(work_dir, 'prepped_data', data_version)
viz_dir <- file.path(work_dir, 'viz', data_version)
dir.create(prepped_dir, showWarnings=FALSE)
dir.create(viz_dir, showWarnings=FALSE)


## PREP SHAPEFILE AND RASTER DATA ------------------------------------------------------->

## Load and prep shapefile
shp <- sf::st_read(file.path(data_dir, 'UGA_shp_consistent_uids.shp'))
colnames(shp)[colnames(shp)=='f15regions'] <- 'subregion'
shp <- shp[, c('uid', 'subregion', paste0('dname',2016:2018), 'geometry')]
# Load merged shapefile
shp_merged <- sf::st_read(file.path(data_dir, 'UGA_shp_stable.shp'))

# Merge on un-capitalized district names from notifications
n_old <- fread(file.path(data_dir, 'district_notifs_2017.csv'))
n_old[, dname2017 := toupper(district) ]
shp <- merge(x=shp, y=n_old[, .(district, dname2017)])
shp_merged <- merge(x=shp_merged, y=n_old[, .(district, dname2017)])
# Convert to lat-long
shp_merged <- sf::st_transform(shp_merged, crs = sf::st_crs(4326))
# Save prepped shapefile
shp_out_fp <- glue('{prepped_dir}/model_shp.shp')
sf::st_write(shp_merged, dsn = shp_out_fp, delete_dsn = TRUE)

## Set up fractional link table
fp_list <- list(covariate_root = '/home/j/WORK/11_geospatial/01_covariates')
link_list <- suppressMessages(suppressWarnings(lbd.mbg::build_link_table(
  shapefile_version = NULL, cores = 1, region = NULL,
  custom_shapefile_path = shp_out_fp, custom_shapefile_field = 'uid'
)))
id_raster <- link_list$id_raster
link_table <- link_list$link_table
# Save to file
fwrite(link_table, file = file.path(prepped_dir, 'link_table.csv'))
writeRaster(id_raster, filename = file.path(prepped_dir, 'id_raster.tif'))
sf::st_write(link_list$id_poly, dsn = file.path(prepped_dir, 'id_poly.shp'), delete_dsn=T)

## Extract covariates - population, population density, GDP PC, proximity to HFs
# Population data by district-year
pop_measures <- c('total','a0004t', 'a0514t')
plist <- vector('list', length=length(pop_measures))
names(plist) <- pop_measures
for(pm in pop_measures){
  plist[[pm]] <- frac_agg_covs(
    cov_config = data.table(covariate='worldpop', measure=pm, agg_method = 'sum'),
    years = c(2016:2019),
    shapefile_path = shp_out_fp, shapefile_field = 'uid',
    core_repo = '/share/code/geospatial/lbd_core/',
    link_table = link_table,
    id_raster = id_raster
  )
}
pop_dt <- Reduce(f = merge, x = plist)
setnames(pop_dt, pop_measures, paste0('pop_',pop_measures))
pop_dt[, pop_under15 := pop_a0004t + pop_a0514t ]
pop_dt[, pop_over15 := pop_total - pop_under15 ]
pop_dt[, pct_over15 := pop_over15 / pop_total ]
pop_dt[, c('pop_a0004t', 'pop_a0514t') := NULL ]
# Covariate data by district-year
cov_config <- data.table(
  covariate = c('tt_hcf', 'gdp', 'ntl_harm', 'ghslurbanicity', 'cattle', 'hdi'),
  measure = c('motorized_tt', 'percapita_log10', 'mean', 'mean', 'da_mean', 'mean'),
  agg_method = 'pop_weight'
)
covs_dt <- frac_agg_covs(
  cov_config = cov_config, years = 2016:2019, shapefile_path = shp_out_fp,
  shapefile_field = 'uid',
  core_repo = '/share/code/geospatial/lbd_core/',
  worldpop_age_sex_release = "2020_03_20",
  link_table = link_table,
  id_raster = id_raster
)
# Normalize covariates
for(cov in c(cov_config$covariate, 'year')){
  norm_cov <- paste0(cov,'_norm')
  covs_dt[, (norm_cov) := (get(cov) - mean(get(cov))) / sd(get(cov))]
}
# Combine district-level tabular information
dist_idx <- as.data.table(shp_merged)[, geometry := NULL]
dist_dt <- (CJ(uid = shp_merged$uid, year = 2016:2019)
  [dist_idx[, .(uid, district)], on='uid']
  [pop_dt, on = c('uid','year')]
  [covs_dt, on=c('uid','year')]
)
dist_dt <- dist_dt[order(year, uid)]
fwrite(dist_dt, file = file.path(prepped_dir, 'dist_dt.csv'))

## Create spatial adjacency matrix
adjmat <- create_adj_matrix(shp_out_fp)
saveRDS(adjmat, file = file.path(prepped_dir, 'adjmat.RDS'))


## PREP PREVALENCE SURVEY --------------------------------------------------------------->

prev_surv <- fread(file.path(data_dir, 'prev_survey_summary.csv'))
setnames(prev_surv, c('num_partic','num_prevalent_cases'), c('n','tb_cases'))
prev_pts <- sf::st_as_sf(prev_surv, coords = c('longitude','latitude'), crs = st_crs(4326))
# Aggregate by district
prev_agg <- as.data.table(sf::st_join(
  x = prev_pts, y=shp_merged[, c('uid','district')]
))[, geometry := NULL ][, lapply(.SD, sum), .SDcols=c('tb_cases','n'), by=.(uid, district)]
# Merge on covariate data, then save
prev_prepped <- merge(x=prev_agg, y=covs_dt[year==2016,], by='uid')
fwrite(prev_prepped, file = file.path(prepped_dir, 'prev_data.csv'))


## PREP NOTIFICATIONS ------------------------------------------------------------------->

dname_merge_dt <- as.data.table(shp)[, geometry := NULL ]

notifs_2017 <- fread(file.path(data_dir, 'district_notifs_2017.csv'))
notifs_2017[, year := 2017][, dname2017 := toupper(district)]
notifs_2017 <- merge(notifs_2017, unique(dname_merge_dt[, .(uid, dname2017)]), by='dname2017')

notifs_2019 <- fread(file.path(data_dir, 'district_notifs_2019.csv'))
notifs_2019[, year := 2019][, dname2018 := toupper(district)]
notifs_2019 <- merge(notifs_2019, unique(dname_merge_dt[, .(uid, dname2018)]), by='dname2018')

notifs_all <- rbindlist(list(notifs_2017, notifs_2019), use.names=T, fill=T)
notifs_agg <- notifs_all[, .(notif_all = sum(notif_all)), by=.(year, uid)]

# Merge onto population and covariate data
notifs_merged <- merge(x=notifs_agg, y=dist_dt, by = c('year','uid'))

# Adjust for extrapulmonary and under-15s
age_adj <- 0.885
ep_adj <- 0.927
notifs_merged[, notif_adj := notif_all * age_adj * ep_adj]
notifs_merged[, pop_adj := pop_over15 ]

# Save to file
fwrite(notifs_merged, file = file.path(prepped_dir, 'notif_data.csv'))


## PREP HIV PREVALENCE BY DISTRICT ------------------------------------------------------>

covs_dir <- '/home/j/WORK/11_geospatial/01_covariates/00_MBG_STANDARD'
pop_rast <- raster::raster(
  file.path(covs_dir, 'worldpop/total/2020_03_20/1y/worldpop_total_1y_2017_00_00.tif')
)
hiv_rast <- raster::raster(
  file.path(covs_dir, 'hiv_test/mean/2019_06_10/1y/hiv_test_mean_1y_2017_00_00.tif')
)
pop_sub <- raster::crop(x=pop_rast, y=id_raster)
hiv_sub <- raster::crop(x=hiv_rast, y=id_raster)
hiv_for_agg <- data.table(
  hiv = as.vector(hiv_sub),
  pop = as.vector(pop_sub),
  pixel_id = as.vector(id_raster)
)
hiv_for_agg <- merge(
  x = hiv_for_agg,
  y = link_table[, .(pixel_id, uid, area_fraction)],
  by = 'pixel_id'
)
hiv_for_agg[, pop_adj := pop * area_fraction ]

# Fractional aggregation
hiv_agg <- hiv_for_agg[, .(hiv = weighted.mean(hiv, w=pop_adj, na.rm=T)), by=uid][order(uid)]

# Save to file
fwrite(hiv_agg, file = file.path(prepped_dir, 'hiv_prev.csv'))
