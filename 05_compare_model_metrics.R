## #######################################################################################
##
## VISUALIZE OUT OF SAMPLE
##
## Author: Nat Henry
## Created: 25 October 2023
## Purpose: Results from prevalence survey OOS runs
##
## #######################################################################################

model_versions <- list(
  is = '20231218_full',
  oos = '20231218_full_oos',
  prev_is = '20231218_prev',
  prev_oos = '20231218_prev_oos'
)

load_libs <- c('data.table','sf','ggplot2','grid','gridExtra','glue','terra','versioning')
invisible(lapply(load_libs, library, character.only=TRUE))

configs <- lapply(model_versions, function(vv) versioning::Config$new(
  file.path('~/temp_data/uga/model_results', vv, 'config.yaml')
))
names(configs) <- names(model_versions)

viz_dir <- file.path(configs$is$get_dir_path('model_results'), 'viz_comparison')
dir.create(viz_dir, showWarnings = FALSE)

rmse <- function(x, y) sqrt(mean((y - x)**2))

## Get predictive metrics for the "national average" model ------------------------------>

# Load prevalence survey data
prev_data <- configs$is$read('prepped_data', 'prev_data')
notif_data <- configs$is$read('prepped_data', 'notif_data')

ps_is_est <- prev_data[, sum(ptb_bc) / sum(sampsize)]
ps_obs <- prev_data[, ptb_bc / sampsize ]
ps_oos_est <- sapply(1:nrow(prev_data), function(row_i){
  prev_data[-row_i, sum(ptb_bc) / sum(sampsize)]
})

ps_is_rmse <- rmse(ps_is_est, ps_obs) * 1e5
ps_oos_rmse <- rmse(ps_oos_est, ps_obs) * 1e5
ps_oos_cor <- cor(ps_oos_est, ps_obs)


## Get predictive metrics for the in-sample models
prev_summaries <- lapply(
  configs[c('is', 'prev_is')],
  function(cc) cc$read('model_results', 'prevalence_summary')
)
fm_is_est <- prev_data[prev_summaries$is, fm_mean := i.mean, on = 'uid'][, fm_mean ]
fm_is_rmse <- rmse(fm_is_est, ps_obs) * 1e5
fm_is_cor <- cor(fm_is_est, ps_obs)

pm_is_est <- prev_data[prev_summaries$prev_is, pm_mean := i.mean, on = 'uid'][, pm_mean]
pm_is_rmse <- rmse(pm_is_est, ps_obs) * 1e5
pm_is_cor <- cor(pm_is_est, ps_obs)

## Get predictive metrics for the out-of-sample models
prev_summaries_oos <- lapply(
  configs[c('oos', 'prev_oos')],
  function(cc) cc$read('model_results', 'oos_summary')
)
fm_oos_est <- prev_data[prev_summaries_oos$oos, fm_oos_est := i.mean, on = 'uid'][, fm_oos_est]
fm_oos_rmse <- rmse(fm_oos_est, ps_obs) * 1e5
fm_oos_cor <- cor(fm_oos_est, ps_obs)

pm_oos_est <- prev_data[prev_summaries_oos$prev_oos, pm_oos_est := i.mean, on = 'uid'][, pm_oos_est ]
pm_oos_rmse <- rmse(pm_oos_est, ps_obs) * 1e5
pm_oos_cor <- cor(pm_oos_est, ps_obs)


## Get validation metrics against the adjusted validation dataset ----------------------->

# Get the "prevalence floor": minimum prevalence consistent with case notifications
# Minimum Inc:Prev ratio = adjustment rate when completeness is 1
adjust_rate <- configs$is$get('duration','intercept') - configs$is$get('duration','slope')
notifs_agg <- notif_data[
  year %in% 2016:2019,
  .(prev_floor = mean(notif_count / pop_over_15) * adjust_rate),
  by = uid
]

adj_data <- (
  copy(prev_data)
  [notifs_agg, prev_floor := i.prev_floor, on = 'uid']
  [, prev_floor_per_100k := prev_floor * 1e5 ]
  [, adj_rate := pmax(ptb_bc / sampsize, prev_floor) ]
  [, adj_per_100k := adj_rate * 1e5 ]
)

# Run validation metrics
adj_obs <- adj_data$adj_rate

# National average in-sample
ps_is_rmse_adj <- rmse(ps_is_est, adj_obs) * 1e5
# National average out-of-sample
ps_oos_rmse_adj <- rmse(ps_oos_est, adj_obs) * 1e5
ps_oos_cor_adj <- cor(ps_oos_est, adj_obs)

# Prevalence only in-sample
pm_is_rmse_adj <- rmse(pm_is_est, adj_obs) * 1e5
pm_is_cor_adj <- cor(pm_is_est, adj_obs)
# Prevalence only out-of-sample
pm_oos_rmse_adj <- rmse(pm_oos_est, adj_obs) * 1e5
pm_oos_cor_adj <- cor(pm_oos_est, adj_obs)

# Full model in-sample
fm_is_rmse_adj <- rmse(fm_is_est, adj_obs) * 1e5
fm_is_cor_adj <- cor(fm_is_est, adj_obs)
# Full model out-of-sample
fm_oos_rmse_adj <- rmse(fm_oos_est, adj_obs) * 1e5
fm_oos_cor_adj <- cor(fm_oos_est, adj_obs)


## Plot original data points and validation dataset ------------------------------------->

validation_fig <- ggplot() + 
  geom_abline(intercept=0, slope=1, linetype=2, lwd=.5, color='#888888') +
  geom_segment(
    data = adj_data[prev_floor_per_100k == adj_per_100k, ],
    aes(
      x = prev_floor_per_100k, xend = prev_floor_per_100k,
      y = prev_per_100k_obs, yend = adj_per_100k - 17
    ),
    arrow = arrow(length = unit(0.1, 'cm'))
  ) +
  geom_point(
    data = adj_data,
    aes(x = prev_floor_per_100k, y = prev_per_100k_obs),
    size = 1.5, pch = 16, color = 'dodgerblue1'
  ) + 
  geom_point(
    data = adj_data,
    aes(x = prev_floor_per_100k, y = adj_per_100k),
    size = 1.5, pch = 1, fill = NA, color = 'black'
  ) + 
  annotate('text', x = 1000 - 68, y = 1000 + 30, label = "Y = X line", color = '#444444') +
  labs(
    title = 'Adjustment of prevalence survey data',
    x = "Minimum 'prevalence floor' consistent with case notifications",
    y = 'Prevalence survey data by district\n(blue: original, black: adjusted)'
  ) + 
  scale_x_continuous(labels = scales::comma, limits = c(0, 1000)) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 2500)) +
  theme_bw()

png(file.path(viz_dir, 'validation_data.png'), height=7, width=7, units='in', res=300)
print(validation_fig)
dev.off()
