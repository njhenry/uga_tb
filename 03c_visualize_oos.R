## #######################################################################################
##
## VISUALIZE OUT OF SAMPLE
##
## Author: Nat Henry
## Created: 25 October 2023
## Purpose: Results from prevalence survey OOS runs
##
## #######################################################################################

MODEL_VERSION_OOS <- '20231218_full_oos'
PREV_VERSION_OOS <- '20231218_prev_oos'

load_libs <- c('data.table','sf','ggplot2','grid','gridExtra','glue','terra')
invisible(lapply(load_libs, library, character.only=TRUE))

# Load custom packages
repos_dir <- '~/repos/'
devtools::load_all(file.path(repos_dir, 'versioning'))

config <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/model_results', MODEL_VERSION_OOS, 'config.yaml')
)
config_p <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/model_results', PREV_VERSION_OOS, 'config.yaml')
)

viz_dir <- file.path(config$get_dir_path('model_results'), 'viz')
dir.create(viz_dir, showWarnings = FALSE)


## Data loading and prep ---------------------------------------------------------------->

# Load input data
oos_full <- config$read('model_results', 'oos_summary')
oos_full[, plot_label := "Joint model"]
oos_prev_only <- config_p$read('model_results', 'oos_summary')
oos_prev_only[, plot_label := "Prevalence-only model"]
measure_cols <- c('mean', 'median', 'lower', 'upper')
for(measure_col in measure_cols){
  oos_full[[measure_col]] <- oos_full[[measure_col]] * 1e5
  oos_prev_only[[measure_col]] <- oos_prev_only[[measure_col]] * 1e5
}

oos_avg <- (
  copy(oos_full)[, .(uid, sampsize, ptb_bc, prev_per_100k_obs)]
  [, c('lower', 'upper', 'median') := NA_real_ ]
  [, plot_label := "National average model"]
)
oos_avg$mean <- sapply(oos_avg[, .I], function(row_i) oos_avg[-row_i, sum(ptb_bc)/sum(sampsize)*1e5])

# Estimate duration when completeness is 100%
min_duration <- config$get('duration', 'intercept') - config$get('duration', 'slope')
inc_data <- config$read("prepped_data", "notif_data")[year %in% config$get('model_years'), ]
inc_floor <- (inc_data
  [ year %in% config$get('model_years'), ]
  [, inc := notif_count / pop_over_15 * 1e5 ]
  [, .(inc = mean(inc, na.rm = T)), by = .(ADM1_EN, ADM2_EN, uid)]
  [, prev_floor := inc * min_duration ]
)

## Prepare data for side-by-side scatters
scatter_keep_cols <- c('uid', 'sampsize', 'prev_per_100k_obs', 'plot_label', measure_cols)
scatter_dt <- rbindlist(list(
  oos_full[, ..scatter_keep_cols], oos_prev_only[, ..scatter_keep_cols],
  oos_avg[, ..scatter_keep_cols]
))

# Merge on "prevalence floor" data
scatter_dt[inc_floor, prev_floor := i.prev_floor, on = 'uid']

rmse <- function(est, obs) sqrt(mean((est - obs)**2))
scatter_dt[, .(rmse(mean, prev_per_100k_obs)), by = plot_label]
scatter_dt[, .(cor(mean, prev_per_100k_obs)), by = plot_label ]
scatter_dt[, .(mean(upper - lower)), by = plot_label]

## Side-by-side scatter plots of prevalence survey estimates vs. out-of-sample preds ---->

scatter_fig <- ggplot(
  data = scatter_dt,
  aes(x = prev_per_100k_obs, y = mean, ymin = lower, ymax = upper)
) +
  facet_wrap('plot_label', ncol = 2) + 
  geom_abline(intercept=0, slope=1, linetype=2, lwd=.5, color='#888888') +
  geom_crossbar(alpha=.5) +
  geom_point(size=2, aes(fill = sampsize), shape=21) +
  viridis::scale_fill_viridis() +
  scale_x_continuous(oob=scales::squish, limits=c(0, 2500), labels=scales::comma) +
  scale_y_continuous(oob=scales::squish, limits=c(0, 2500), labels=scales::comma) +
  labs(
    title='Out-of-sample comparison to TB Prevalence Survey',
    fill='Prevalence survey\nsample size',
    x = 'Unadjusted survey prevalence', y = 'Out-of-sample modeled TB prevalence'
  ) +
  theme_bw() +
  theme(legend.position = c(.75, .25))
png(file.path(viz_dir, 'oos_data_scatter.png'), height=7.5, width=7.5, units='in', res=300)
print(scatter_fig)
dev.off()


## Plot results against minimum possible prevalence given case notifications ------------>

prev_floor_labels <- c(
  '95% UI bounds above prevalence floor',
  'Lower bound of 95% UI below prevalence floor',
  'Mean estimate below prevalence floor'
)
prev_floor_colors <- c(
  rgb(0, 0.6, 0.2, alpha = 0.5),
  rgb(1, 0.6, 0, alpha = 0.5),
  rgb(0.8, 0, 0, alpha = 0.5)
)
names(prev_floor_colors) <- prev_floor_labels

(scatter_dt
  [, pf_label := prev_floor_labels[1] ]
  [ lower < prev_floor, pf_label := prev_floor_labels[2] ]
  [ mean < prev_floor, pf_label := prev_floor_labels[3] ]
)

scatter_dt[, .N, by = .(plot_label, pf_label)]


scatter_by_prev_floor <- ggplot(
  data = scatter_dt,
  aes(x = prev_floor, y = mean, ymin = lower, ymax = upper, color = pf_label)
) +
  facet_wrap('plot_label', ncol = 2) + 
  geom_abline(intercept=0, slope=1, linetype=2, lwd=.5, color='#888888') +
  geom_crossbar(show.legend = FALSE) +
  geom_point(size = 2) + 
  scale_color_manual(values = prev_floor_colors) +
  scale_x_continuous(limits = c(0, 1500), oob = scales::squish, labels = scales::comma) + 
  scale_y_continuous(limits = c(0, 2000), oob = scales::squish, labels = scales::comma) +
  labs(
    title = "Out-of-sample comparison to 'prevalence floor' from case notifications",
    color = "Estimate compared to prevalence floor",
    x = "Minimum reasonable prevalence given case notifications",
    y = 'Out-of-sample modeled TB prevalence'
  ) + 
  theme_bw() + 
  theme(legend.position = c(.75, .25))
png(file.path(viz_dir, 'prev_floor_scatter.png'), height=7.5, width=7.5, units='in', res=300)
print(scatter_by_prev_floor)
dev.off()
