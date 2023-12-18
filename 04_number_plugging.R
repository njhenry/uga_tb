## #######################################################################################
##
## NUMBER PLUGGING SCRIPT FOR UGANDA TB MANUSCRIPT
##
## Author: Nat Henry
## Created: 3 September 2021
## Purpose: Format outputs in a manner that is convenient for number-plugging
##  the paper
##
## #######################################################################################

# Settings
MODEL_VERSION <- '20231218_full'
PREV_VERSION <- '20231218_prev'
SUMMARY_YEAR <- 2019

load_libs <- c('data.table','glue','matrixStats')
invisible(lapply(load_libs, library, character.only=TRUE))

# Load custom packages
repos_dir <- '~/repos/'
devtools::load_all(file.path(repos_dir, 'versioning'))

# Load configs
config <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/model_results', MODEL_VERSION, 'config.yaml')
)
config_p <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/model_results', PREV_VERSION, 'config.yaml')
)

# Load input data
prev_data <- config$read('prepped_data', 'prev_data')
notif_data <- config$read('prepped_data', 'notif_data')
dist_dt <- config$read('prepped_data', 'covariates')
prev_summ <- config$read('model_results', 'prevalence_summary')
comp_summ <- config$read('model_results', 'completeness_observations')
prev_summ_no_notifs <- config_p$read('model_results', 'prevalence_summary')

model_preds <- config$read('model_results', 'model_preds')
prev_draws <- model_preds$pred_draws_prevalence
comp_draws <- model_preds$pred_draws_completeness

# Merge population onto the district table
dist_dt[notif_data, pop_over_15 := i.pop_over_15, on = c('uid', 'year')]

## Incidence and completeness summary --------------------------------------------------->

ival <- function(inc) round(inc * 1e5)
comp_reporting <- (copy(comp_summ)
  [ year == SUMMARY_YEAR, ]
  [, inc_ciwidth := ival(inc_upper - inc_lower)]
  [, `:=` (inc_mean = ival(inc_mean), inc_upper = ival(inc_upper), inc_lower = ival(inc_lower))]
  [(order(-inc_mean))]
)

# Completeness summary at the district level
start_year <- comp_summ[, min(year)]

c_to_d <- function(comp_val) (
  config$get('duration', 'intercept') - config$get('duration', 'slope') * comp_val
)
d_summ <- (
  copy(comp_summ)
  [year %in% c(start_year, SUMMARY_YEAR), ]
  [, `:=` (d_mean = c_to_d(mean), d_lower = c_to_d(upper), d_upper = c_to_d(lower))]
)


## NUMBER PLUGGING ---------------------------------------------------------------------->

## Count change in case notification rate from 2016 to 2019

## SECTION: Spatial variation in tuberculosis incidence across Uganda
maincols <- c('ADM1_EN','ADM2_EN','inc_mean','inc_lower','inc_upper')
head(comp_reporting[, ..maincols])
tail(comp_reporting[, ..maincols])

# Change in district completeness over time
comp_change <- merge(
  x = comp_summ[year == start_year, .(ADM1_EN, ADM2_EN, uid, mean)],
  y = comp_summ[year == SUMMARY_YEAR, .(ADM1_EN, ADM2_EN, uid, mean)],
  by = c("ADM1_EN", "ADM2_EN", "uid"),
  suffixes = c('_start', '_end')
)
comp_change[, .(mean(mean_end > mean_start))]
comp_change[, .(sum(mean_end > mean_start))]


# Comparison of completeness in start year versus comparison year
comp_summ[(year == start_year) , .(sum(mean > .7) / .N)]
comp_summ[(year == start_year), .(sum(mean < .5) / .N)]
comp_summ[(year == SUMMARY_YEAR), .(sum(mean > .7) / .N)]
comp_summ[(year == SUMMARY_YEAR), .(sum(mean < .5) / .N)]

# Comparison of average durations at the start versus in the summary year
d_summ[, .(d_mean = weighted.mean(d_mean, w = pop_over_15)), by = year]
(d_summ
  [year == SUMMARY_YEAR, ]
  [order(-d_mean)]
  [c(1, .N), .(ADM1_EN, ADM2_EN, comp = mean, d_mean)]
)

## SECTION: Completeness of case notifications increases over time
# Case notifs, 2017
comp_cols <- c('ADM1_EN', 'ADM2_EN', 'mean', 'lower', 'upper')
head(comp_summ[year==start_year, ..comp_cols][order(-mean)])
tail(comp_summ[year==SUMMARY_YEAR, ..comp_cols][order(-mean)])

## SECTION: Effect of including notifications on estimates of TB burden

# UI widths of both results sets
prev_summ[, range(ival(upper-lower))]
prev_summ[, .(ival(mean(upper-lower)))]

prev_summ_no_notifs[, tt := ival(upper-lower)][order(tt)][c(1, 2, .N-3, .N-2, .N-1, .N), tt]
prev_summ_no_notifs[, .(ival(mean(upper-lower)))]

# How many districts fall outside of 'low' or 'high' prevalence bounds?
LOW_PREV <- 300/1e5
HIGH_PREV <- 600/1e5
calc_bounds <- function(mean_vec, lower_vec, upper_vec){
  output <- rep('Neither', length(mean_vec))
  output[mean_vec > HIGH_PREV] <- 'High burden, low confidence'
  output[lower_vec > HIGH_PREV] <- 'High burden, high confidence'
  output[mean_vec < LOW_PREV] <- 'Low burden, low confidence'
  output[upper_vec < LOW_PREV] <- 'Low burden, high confidence'
  return(output)
}

prev_summ$bounds_lab <- calc_bounds(prev_summ$mean, prev_summ$lower, prev_summ$upper)
prev_summ[, .N, by=bounds_lab][order(bounds_lab)]

prev_summ_no_notifs$bounds_lab <- calc_bounds(
  prev_summ_no_notifs$mean, prev_summ_no_notifs$lower, prev_summ_no_notifs$upper
)
prev_summ_no_notifs[, .N, by=bounds_lab][order(bounds_lab)]
