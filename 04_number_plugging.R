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

load_libs <- c('data.table','glue','matrixStats')
invisible(lapply(load_libs, library, character.only=TRUE))

# Settings
data_version <- '20211003'
prev_data_version <- '20210920_prev_only'
duration_intercept <- 1.62
duration_slope <- .56

# Input and output filepaths
work_dir <- 'REDACTED'
prepped_dir <- file.path(work_dir, 'prepped_data', data_version)
model_dir <- file.path(work_dir, 'model_results', data_version)

# Load input data
prev_data <- fread(file.path(prepped_dir, 'prev_data.csv'))
notif_data <- fread(file.path(prepped_dir, 'notif_data.csv'))
dist_dt <- fread(file.path(prepped_dir, 'dist_dt.csv'))
prev_summ <- fread(file.path(model_dir, 'prev_summ.csv'))
comp_summ <- fread(file.path(model_dir, 'comp_summ.csv'))
prev_summ_no_notifs <- fread(file.path(
  work_dir, 'model_results', prev_data_version, 'prev_summ.csv'
))

model_preds <- readRDS(file.path(model_dir, 'model_preds.RDS'))
prev_draws <- model_preds$pred_draws_prevalence
comp_draws <- model_preds$pred_draws_completeness

## Prepare prevalence and completeness outputs ------------------------------------------>

# Prevalence summary at the district level
pval <- function(prev) round(prev*1E5, 1)
prev_summ[, ciwidth := pval(upper - lower)]
prev_summ[, `:=` (mean = pval(mean), lower=pval(lower), upper=pval(upper))]
prev_summ[dist_dt[year==2016], pop := i.pop_over15, on='uid']
prev_summ <- prev_summ[order(-mean)]
# Same for no-notifications results
prev_summ_no_notifs[, ciwidth := pval(upper - lower)]
prev_summ_no_notifs[, `:=` (mean = pval(mean), lower=pval(lower), upper=pval(upper))]
prev_summ_no_notifs <- prev_summ_no_notifs[order(-mean)]

# Completeness summary at the district level
c_to_d <- function(comp_val) duration_intercept - duration_slope * comp_val
comp_summ[, `:=` (d_mean = c_to_d(mean), d_lower = c_to_d(upper), d_upper = c_to_d(lower))]

agg_draws <- function(draws_mat, weights){
  draws <- sapply(1:1000, function(dd) weighted.mean(draws_mat[, dd], w=weights))
  summ <- c(mean(draws), quantile(draws, probs=c(0.025, 0.975)))
  return(summ)
}
natl_prev_summ <- agg_draws(prev_draws, w=dist_dt[year==2016, pop_over15]) * 1E5
natl_comp_summ <- list(
  'yr2017' = agg_draws(comp_draws[1:122,], w=dist_dt[year==2017, pop_over15]),
  'yr2019' = agg_draws(comp_draws[123:244,], w=dist_dt[year==2019, pop_over15])
)

# Get change at the district level, 2017 to 2019
dist_comp_change <- comp_draws[123:244, ] / comp_draws[1:122, ]
comp_change_summ <- dist_dt[year==2016, .(uid, district)]
comp_change_summ$change_mean <- rowMeans(dist_comp_change)
comp_change_summ$change_lower <- rowQuantiles(dist_comp_change, probs=0.025)
comp_change_summ$change_upper <- rowQuantiles(dist_comp_change, probs=0.975)
comp_change_summ <- comp_change_summ[order(-change_mean)]

## NUMBER PLUGGING ---------------------------------------------------------------------->

## SECTION: Spatial variation in tuberculosis prevalence across Uganda
# Prevalence
maincols <- c('district','mean','lower','upper')
natl_prev_summ
head(prev_summ[, ..maincols])
tail(prev_summ[, ..maincols])


## SECTION: Completeness of case notifications increases over time
# Case notifs, 2017
comp_summ <- comp_summ[order(year, -mean)]
head(comp_summ[year==2017, ..maincols])
tail(comp_summ[year==2017, ..maincols])

# Case notifs, 2019
head(comp_summ[year==2019, ..maincols])
tail(comp_summ[year==2019, ..maincols])

# Change in completeness over time
comp_change_summ[, sum(change_mean > 1) ]
comp_change_summ[, sum(change_mean > 1) / .N ]
comp_change_summ[, mean(change_mean)]

# Average duration, 2017
c_to_d(natl_comp_summ$yr2017)
dcols <- c('district','d_mean','d_lower','d_upper')
d_2017 <- comp_summ[year==2017, ..dcols][order(-d_mean)]
head(d_2017)
tail(d_2017)

# Average duration, 2019
c_to_d(natl_comp_summ$yr2019)
d_2019 <- comp_summ[year==2019, ..dcols][order(-d_mean)]
head(d_2019)
tail(d_2019)


## SECTION: Effect of including notifications on estimates of TB burden

# UI widths of both results sets
mean(prev_summ$ciwidth)
quantile(prev_summ$ciwidth)

mean(prev_summ_no_notifs$ciwidth)
quantile(prev_summ_no_notifs$ciwidth)

# How many districts fall outside of the prevalence bounds?
prev_low <- 200
prev_high <- 759
calc_bounds <- function(mean_vec, lower_vec, upper_vec){
  output <- rep('Neither', length(mean_vec))
  output[mean_vec > 759] <- 'High burden, low confidence'
  output[lower_vec > 759] <- 'High burden, high confidence'
  output[mean_vec < 253] <- 'Low burden, low confidence'
  output[upper_vec < 253] <- 'Low burden, high confidence'
  return(output)
}

prev_summ$bounds_lab <- calc_bounds(prev_summ$mean, prev_summ$lower, prev_summ$upper)
prev_summ[, .N, by=bounds_lab][order(bounds_lab)]

prev_summ_no_notifs$bounds_lab <- calc_bounds(
  prev_summ_no_notifs$mean, prev_summ_no_notifs$lower, prev_summ_no_notifs$upper
)
prev_summ_no_notifs[, .N, by=bounds_lab][order(bounds_lab)]
