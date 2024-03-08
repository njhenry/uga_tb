## #######################################################################################
##
## VISUALIZE MODEL RESULTS
##
## Author: Nat Henry
## Created: 3 September 2021
## Purpose: Show Uganda TB model results in relation to the underlying data
##
## #######################################################################################

MODEL_VERSION <- '20231218_full'
PREV_VERSION <- '20231218_prev'

# Load packages
load_libs <- c(
  'data.table','sf','ggplot2','ggspatial','grid','gridExtra','glue','terra','versioning'
)
invisible(lapply(load_libs, library, character.only=TRUE))

# Load configuration files for the default and prevalence-only models
config <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/model_results', MODEL_VERSION, 'config.yaml')
)
config_p <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/model_results', PREV_VERSION, 'config.yaml')
)

viz_dir <- file.path(config$get_dir_path('model_results'), 'viz')
dir.create(viz_dir, showWarnings = FALSE)

# Load spatial data
uga <- config$read('shps', 'viz_adm2')
outline_sf <- config$read('shps', 'viz_adm0')

# Load input data
prev_data <- config$read("prepped_data", "prev_data")
notif_data <- config$read("prepped_data", "notif_data")
notif_data[, notifs_per_100k := notif_count / pop_over_15 * 1e5 ]
dist_dt <- config$read("prepped_data", "covariates")
prev_summ <- config$read("model_results", "prevalence_summary")
comp_summ <- config$read("model_results", "completeness_summary")
prev_summ_no_notifs <- config_p$read("model_results", "prevalence_summary")
inc_summ <- config$read('model_results', 'incidence_summary')
comp_obs <- config$read('model_results', 'completeness_observations')

model_years <- config$get('model_years')
notif_data <- copy(notif_data[year %in% model_years, ])
dist_dt <- copy(dist_dt[year %in% model_years, ])


## DATA PREP ---------------------------------------------------------------------------->

# Merge prevalence data, estimates, and 2019 completeness onto spatial data
prev_sf <- merge(
  x=uga[, c('uid')],
  y=prev_data[, .(uid, ADM2_EN, ptb_bc, sampsize)],
  all.x=T
)
measures <- c('mean','lower','upper')
setnames(prev_summ, measures, paste0('prev_',measures))
prev_sf <- merge(x=prev_sf, y=prev_summ[, .(uid,prev_mean,prev_lower,prev_upper)], all.x=T)
setnames(comp_obs, measures, paste0('comp_',measures))
prev_sf <- merge(
  x=prev_sf, y=comp_obs[ year==2019, .(uid,comp_mean,comp_lower,comp_upper)], all.x=TRUE
)
prev_sf$raw_prev <- prev_sf$ptb_bc / prev_sf$sampsize
setnames(inc_summ, measures, paste0('inc_', measures))
prev_sf <- merge(
  x = prev_sf, y = inc_summ[year == 2019, .(uid, inc_mean, inc_lower, inc_upper)],
  all.x = TRUE
)
# Add on no-notification prevalence data if it exists
if(exists('prev_summ_no_notifs')){
  psnn_cols <- paste0('nn_prev_',measures)
  setnames(prev_summ_no_notifs, measures, psnn_cols)
  psnn_keep_cols <- c('uid',psnn_cols)
  prev_sf <- merge(x=prev_sf, y=prev_summ_no_notifs[, ..psnn_keep_cols], by='uid', all.x=T)
}
# Set prevalence to rates per 100,000
for(prev_col in grep('prev|inc', colnames(prev_sf), value=T)){
  prev_sf[[prev_col]] <- prev_sf[[prev_col]] * 1E5
}


## MAKE PLOTS --------------------------------------------------------------------------->

## Plot: Raw prevalence versus modeled prevalence, and scatter
prev_cols <- RColorBrewer::brewer.pal(n=9, name='YlOrBr')
cn_cols <- RColorBrewer::brewer.pal(n=9, name='Blues')
prev_breaks <- seq(0, 1000, by=250)
prev_lims <- range(prev_breaks)
country_outline <- geom_sf(data=outline_sf, fill=NA, lwd=.25, linetype=2, color='#222222')


## FIG: INPUT DATA PLOT - TB PREVALENCE

fig_ps_data <- ggplot() +
  geom_sf(data=prev_sf, aes(fill=raw_prev), lwd=0.25, color='#222222') +
  country_outline +
  ggspatial::annotation_scale(location='br', width_hint = 0.5) +
  scale_fill_gradientn(
    colors = prev_cols, breaks=prev_breaks, limits=prev_lims, oob = scales::squish,
    na.value = '#AAAAAA', labels=scales::comma
  ) +
  labs(title=NULL, fill='Observed\nTB prevalence\nper 100,000', x=NULL, y=NULL) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent')
  )

png(file.path(viz_dir, 'data_prev_survey.png'), height=3.5, width=4.5, units='in', res=300)
plot(fig_ps_data)
dev.off()
pdf(file.path(viz_dir, 'Figure_1-top.pdf'), height=3.5, width=4.5)
plot(fig_ps_data)
dev.off()

## FIG: INPUT DATA PLOT - CASE NOTIFICATIONS IN ALL MODEL YEARS

uid_year_square <- data.table::CJ(uid = unique(uga$uid), year = unique(notif_data$year))
cn_data_sf <- merge(x = uga[, c('uid')], y = uid_year_square, by = 'uid')
cn_data_sf <- merge(
  x = cn_data_sf,
  y = notif_data[, .(uid, ADM2_EN, year, notifs_per_100k)],
  by = c('uid', 'year'),
  all.x = TRUE
)

fig_cn_data <- ggplot() +
  facet_wrap('year', nrow=ceiling(sqrt(length(model_years)))) +
  geom_sf(data=cn_data_sf, aes(fill=notifs_per_100k), lwd=0.25, color='#222222') +
  country_outline +
  scale_fill_gradientn(
    colors = cn_cols, breaks=prev_breaks, limits=prev_lims, oob = scales::squish,
    na.value = '#AAAAAA', labels=scales::comma
  ) +
  labs(title=NULL, fill='Case notifications\nper 100,000\n(ages 15+)', x=NULL, y=NULL) +
  theme_bw() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent')
  )
png(file.path(viz_dir, 'data_case_notifs.png'), height=7.5, width=7.5, units='in', res=300)
plot(fig_cn_data)
dev.off()
pdf(file.path(viz_dir, 'Figure_1-bottom.pdf'), height=7.5, width=7.5)
plot(fig_cn_data)
dev.off()


## FIG: MAIN MODELED RESULTS - TB INCIDENCE PER 100k

main_inc <- ggplot() + 
  geom_sf(data = prev_sf, aes(fill = inc_mean), lwd = 0.25, color = '#222222') +
  country_outline +
  ggspatial::annotation_scale(location='br', width_hint = 0.5) +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(n = 9, name = 'Greens'), breaks = prev_breaks,
    limits = prev_lims, oob = scales::squish, labels = scales::comma
  ) + 
  labs(
    title = 'Modeled TB Incidence (Mean Estimate, 2019)',
    fill = "TB incidence\nper 100,000",
    x = NULL, y = NULL
  ) + 
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )
png(file.path(viz_dir, 'incidence_map.png'), height=6, width=6, units='in', res=300)
print(main_inc)
dev.off()
pdf(file.path(viz_dir, 'Figure_3.pdf'), height=6, width=6)
print(main_inc)
dev.off()

## FIG: MAIN MODELED RESULTS - TB PREVALENCE PER 100k

main_prev_a <- ggplot() +
  geom_sf(data=prev_sf, aes(fill=prev_mean), lwd=0.25, color='#222222') +
  country_outline +
  scale_fill_gradientn(
    colors = prev_cols, breaks=prev_breaks, limits=prev_lims,
    oob = scales::squish, labels=scales::comma
  ) +
  labs(
    title='Modeled TB Prevalence (Mean Estimate)', fill='TB prevalence\nper 100,000',
    x=NULL, y=NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )
png(file.path(viz_dir, 'prev_map.png'), height=6, width=6, units='in', res=300)
print(main_prev_a)
dev.off()


scatter_dt <- as.data.table(prev_sf)[, geometry := NULL ][!is.na(raw_prev)]
col_breaks <- seq(0, 10, by=2)
col_labs <- c(as.character(col_breaks[1:5]), '10+')
main_prev_b <- ggplot(
    data=na.omit(scatter_dt),
    aes(x=raw_prev, y=prev_mean, ymin=prev_lower, ymax=prev_upper)
  ) +
  geom_abline(intercept=0, slope=1, linetype=2, lwd=.5, color='#888888') +
  geom_crossbar(alpha=.5) +
  geom_point(size=2, aes(fill=ptb_bc), shape=21) +
  viridis::scale_fill_viridis(
    limits = range(col_breaks), breaks=col_breaks, labels=col_labs, oob = scales::squish
  ) +
  scale_x_continuous(oob=scales::squish, limits=c(0, 2500), labels=scales::comma) +
  scale_y_continuous(oob=scales::squish, limits=c(0, 1500), labels=scales::comma) +
  labs(
    title='Comparison to TB Prevalence Survey',
    fill='Number of\nobserved\nTB cases',
    x='Unadjusted survey prevalence', y='Modeled TB prevalence'
  ) +
  theme_bw()
png(file.path(viz_dir, 'prev_data_scatter.png'), height=6, width=6, units='in', res=300)
print(main_prev_b)
dev.off()


notif_scatter_dt <- merge(
  x = inc_summ[, .(uid, year, inc_mean, inc_lower, inc_upper)],
  y = notif_data[, .(uid, year, notif_count, pop_over_15)],
  by = c('uid', 'year')
)
notif_scatter_dt[, cnr_adj := notif_count / pop_over_15 * 1E5 ]
notif_scatter_dt[, `:=` (
  inc_mean=inc_mean*1E5, inc_lower=inc_lower*1E5, inc_upper=inc_upper*1E5
)]
notif_scatter_dt[inc_upper >= 1200, inc_upper := 1200 - 1E-5]

n_breaks <- seq(0, 1000, by=200)
n_labs <- c(as.character(n_breaks[1:5]), '1,000+')
inc_scatter <- ggplot(
  data=na.omit(notif_scatter_dt),
  aes(x=cnr_adj, y=inc_mean, ymin=inc_lower, ymax=inc_upper)
) +
  facet_wrap('year', ncol = ceiling(sqrt(length(model_years)))) +
  geom_abline(intercept=0, slope=1, linetype=2, lwd=.5, color='#888888') +
  geom_crossbar(lwd=.25, color='#AAAAAA') +
  geom_point(size=1.5, aes(fill=notif_count), shape=21, alpha=.5) +
  viridis::scale_fill_viridis(
    limits = range(n_breaks), breaks=n_breaks, labels=n_labs, oob = scales::squish
  ) +
  scale_x_continuous(oob=scales::squish, limits=c(0, 1200), labels=scales::comma) +
  scale_y_continuous(oob=scales::squish, limits=c(0, 1200), labels=scales::comma) +
  labs(
    title='Comparison to TB Case Notification Rate',
    fill='Annual\ncase\nnotifications',
    x='Average CNR for pulmonary TB, ages 15+', y='Modeled TB incidence'
  ) +
  theme_bw()
png(file.path(viz_dir, 'notifs_data_scatter.png'), height=6, width=6, units='in', res=300)
print(inc_scatter)
dev.off()

png(file.path(viz_dir, 'prev_combined_plot.png'), height=9, width=6, units='in', res=300)
grid.arrange(
  ggplotGrob(main_prev_a), ggplotGrob(main_prev_b), ggplotGrob(inc_scatter),
  layout_matrix=matrix(1:3, ncol=1)
)
dev.off()


## Plot: Completeness in first and last model years
comp_sf <- merge(
  x = uga[, c('uid')],
  y = comp_obs[, .(uid, year, comp_mean, comp_lower, comp_upper)],
  by = 'uid'
)
comp_sf <- comp_sf[comp_sf$year %in% range(comp_sf$year), ]

comp_cols <- RColorBrewer::brewer.pal(name='BuPu', n=9)
comp_breaks <- seq(0.2, 0.8, by=.1)
comp_lims <- range(comp_breaks)

comp_fig <- ggplot(data=comp_sf) +
  facet_wrap('year', nrow = 1) +
  country_outline +
  geom_sf(data=comp_sf, aes(fill=comp_mean), lwd=.15, color='#222222') +
  scale_fill_gradientn(
    colors = comp_cols, limits = comp_lims, breaks = comp_breaks, oob = scales::squish,
    labels = function(x) scales::percent(x, accuracy=1)
  ) +
  theme_bw() +
  labs(fill='TB notification\ncompleteness', title=NULL, x=NULL, y=NULL) +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )

png(file.path(viz_dir, 'notif_completeness.png'), height=4, width=7.5, units='in', res=300)
plot(comp_fig)
dev.off()
pdf(file.path(viz_dir, 'Figure_4.pdf'), height = 4, width = 7.5)
plot(comp_fig)
dev.off()

## Plot: Average duration in all model years

comp_sf$duration <- (
  config$get('duration', 'intercept') - config$get('duration', 'slope') * comp_sf$comp_mean
)
dur_cols <- viridis::magma(n=10)
dur_breaks <- seq(1.3, 1.7, by=.1)
dur_lims <- range(dur_breaks)

duration_fig <- ggplot(data=comp_sf) +
  facet_wrap('year', ncol = ceiling(sqrt(length(model_years)))) +
  country_outline +
  geom_sf(data=comp_sf, aes(fill=duration), lwd=.15, color='#222222') +
  scale_fill_gradientn(
    colors = dur_cols, limits = dur_lims, breaks = dur_breaks, oob = scales::squish
  ) +
  theme_bw() +
  labs(fill='Average\nDuration\n(years)', title=NULL, x=NULL, y=NULL) +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )

png(file.path(viz_dir, 'average_duration.png'), height=3.5, width=7.5, units='in', res=300)
plot(duration_fig)
dev.off()



## Plot: Prevalence-incidence ratio compared to completeness
ip_ratio_dt <- merge(
  x = prev_data[, .(uid, ptb_bc, sampsize)],
  y = notif_data[, .(uid, year, notif_count, pop_over_15)]
)
ip_ratio_dt <- merge(
  x = ip_ratio_dt,
  y = comp_summ[, .(uid, year, comp_mean = mean, comp_lower = lower, comp_upper = upper)],
  by = c('uid','year')
)
ip_ratio_dt[, inc_prev_ratio := (notif_count/pop_over_15)/(ptb_bc/sampsize) * .82 ]
ip_ratio_dt <- ip_ratio_dt[(inc_prev_ratio < Inf) & !is.na(inc_prev_ratio),  ]
ip_ratio_dt[, yr_label := as.character(year) ]

ip_ratio_fig <- ggplot(
    data=ip_ratio_dt,
    aes(x=inc_prev_ratio, y=comp_mean, ymin=comp_lower, ymax=comp_upper)
  ) +
  facet_wrap('year', nrow=1) +
  geom_crossbar(alpha=.3) +
  geom_point(aes(fill=yr_label), shape=21, size=1.5) +
  geom_smooth(method = 'lm', linetype = 3) +
  labs(
    x='Notification-prevalence ratio (observed)', y='Estimated completeness', fill='Year') +
  scale_x_continuous(limits=c(0, 1.5)) +
  theme_bw()

png(file.path(viz_dir, 'notif_prev_ratio.png'), height=4, width=8, units='in', res=300)
plot(ip_ratio_fig)
dev.off()


## Plot: Differences with prevalence-only results
prev_sf$overlap <- 'Not significant'
prev_sf[prev_sf$prev_lower > prev_sf$nn_prev_upper, 'overlap'] <- 'Increased prevalence'
prev_sf[prev_sf$prev_upper < prev_sf$nn_prev_lower, 'overlap'] <- 'Decreased prevalence'
prev_sf$nn_diff <- prev_sf$prev_mean - prev_sf$nn_prev_mean

nn_fig_a <- ggplot() +
  country_outline +
  geom_sf(data=prev_sf, aes(fill=nn_prev_mean), lwd=0.15, color='#222222') +
  scale_fill_gradientn(
    colors = prev_cols, breaks=prev_breaks, limits=prev_lims,
    oob = scales::squish, labels = scales::comma
  ) +
  labs(
    title='A. Modeled TB Prevalence, No Notifications',
    fill='Mean estimated\nTB prevalence\nper 100,000',
    x=NULL, y=NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )

diff_cols <- rev(RColorBrewer::brewer.pal(name='RdYlBu', n=9))
diff_breaks <- seq(-300, 300, by=150)
diff_labs <- c('-300', '-150', 'No change', '+150', '+300')
diff_limits <- range(diff_breaks)

nn_fig_b <- ggplot() +
  country_outline +
  geom_sf(data=prev_sf, aes(fill=nn_diff), lwd=0.15, color='#222222') +
  scale_fill_gradientn(
    colors = diff_cols, breaks=diff_breaks, limits=diff_limits, labels = diff_labs,
    oob = scales::squish
  ) +
  labs(
    title='B. Effect of Notifications on Mean Estimate',
    fill='Difference in mean\nprevalence estimate', x=NULL, y=NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )

diff_ui_cols <- c(
  'red', '#EEEEEE', 'blue'
)
names(diff_ui_cols) <- c('Increased prevalence', 'Not significant', 'Decreased prevalence')

nn_fig_c <- ggplot() +
  country_outline +
  geom_sf(data=prev_sf, aes(fill=overlap), lwd=.15, color='#222222') +
  scale_fill_manual(values = diff_ui_cols) +
  labs(
    title='C. Effect of Notifications on Uncertainty Intervals',
    fill = 'Effect of adding\nnotifications', x=NULL, y=NULL) +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )

png(file.path(viz_dir, 'prev_notifs_diff.png'), height=8, width=5, units='in', res=300)
grid.arrange(
  ggplotGrob(nn_fig_a), ggplotGrob(nn_fig_b), ggplotGrob(nn_fig_c),
  layout_matrix=matrix(1:3, ncol=1)
)
dev.off()


## FIG: Comparing ability to distinguish between low and high-burden districts

thresh_low <- 300
thresh_high <- 600
tlabs <- c(
  'High burden (high confidence)',
  'High burden (low confidence)',
  'Low burden (high confidence)',
  'Low burden (low confidence)',
  'Neither'
)
tcols <- c('#800000', '#f9cac9', '#000080', '#cddff6', '#EEEEEE')
names(tcols) <- tlabs
set_thresh <- function(mean_vec, low_vec, high_vec){
  tvec <- rep(tlabs[5], length(mean_vec))
  tvec[mean_vec > thresh_high] <- tlabs[2]
  tvec[low_vec > thresh_high] <- tlabs[1]
  tvec[mean_vec < thresh_low] <- tlabs[4]
  tvec[high_vec < thresh_low] <- tlabs[3]
  return(tvec)
}

prev_sf$thresh_full_model <- set_thresh(
  prev_sf$prev_mean, prev_sf$prev_lower, prev_sf$prev_upper
)
prev_sf_t1 <- prev_sf[, c('thresh_full_model')]
prev_sf_t1$model_type <- "Full model"
prev_sf$thresh_no_notifs <- set_thresh(
  prev_sf$nn_prev_mean, prev_sf$nn_prev_lower, prev_sf$nn_prev_upper
)
prev_sf_t2 <- prev_sf[, c('thresh_no_notifs')]
prev_sf_t2$model_type <- "Survey-only model"
colnames(prev_sf_t1) <- colnames(prev_sf_t2) <- c('thresh', 'geometry', 'model_type')
thresh_plot_sf <- rbind(prev_sf_t1, prev_sf_t2)

thresh_fig <- ggplot() +
  facet_wrap('model_type', nrow=1) +
  country_outline +
  geom_sf(data=thresh_plot_sf, aes(fill=thresh), lwd=.15, color='#222222') +
  scale_fill_manual(values = tcols) +
  labs(title=NULL, fill = 'Predicted burden category', x=NULL, y=NULL) +
  theme_bw() +
  theme(
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.text.y=element_blank(), axis.ticks.y=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )

png(file.path(viz_dir, 'thresh_comparison.png'), height=3.5, width=7.5, units='in', res=300)
print(thresh_fig)
dev.off()
pdf(file.path(viz_dir, 'Figure_5.pdf'), height = 3.5, width = 7.5)
print(thresh_fig)
dev.off()
