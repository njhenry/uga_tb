## #######################################################################################
##
## 03d: Compare incidence and completeness results across two model runs
##
## AUTHOR: Nat Henry, nat.henry7@gmail.com
## CREATED: 1 November 2023
## PURPOSE: Test the inclusion of 2020-2022 data on model estimates
## 
## #######################################################################################

versions <- list(A = '20231024_v2', B = '20231101')
version_labels <- list(A = 'Notifs 2016 to 2019', B = 'Notifs 2016 to 2022')
comparison_year <- 2019

load_libs <- c('data.table','sf','ggplot2','grid','gridExtra','glue','terra')
invisible(lapply(load_libs, library, character.only=TRUE))

# Load custom packages
repos_dir <- '~/repos/'
devtools::load_all(file.path(repos_dir, 'versioning'))

configs <- lapply(versions, function(vv) versioning::Config$new(
  config_list = file.path("~/temp_data/uga/model_results", vv, "config.yaml")
))
names(configs) <- names(versions)

district_sf <- configs$A$read('shps', 'viz_adm2')
outline_sf <- configs$A$read('shps', 'viz_adm0')

comparison_dir <- file.path(
  configs$A$config_list$directories$model_results$path,
  'comparisons',
  paste0(versions$B, '_vs_', versions$A)
)
dir.create(comparison_dir, recursive = T, showWarnings = F)


## Load observed completeness/incidence data table -------------------------------------->

notifs_a <- (
  configs$A$read('model_results', 'completeness_observations')
  [year == comparison_year, ]
  [, model := paste("Model:", version_labels$A) ]
)
notifs_b <- (
  configs$B$read('model_results', 'completeness_observations')
  [year == comparison_year, ]
  [, model := paste("Model:", version_labels$B) ]
)

notifs_long <- rbindlist(list(notifs_a, notifs_b), use.names = T, fill = T)
notifs_long_sf <- merge(
  x = district_sf, y = notifs_long, by = c('ADM1_EN', 'ADM2_EN', 'uid')
)

notifs_wide <- merge(
  x = notifs_a, y = notifs_b,
  by = c('ADM1_EN','ADM2_EN','uid','year'),
  suffixes = c('_a','_b')
)
(notifs_wide
  [, inc_mean_diff := inc_mean_b - inc_mean_a ]
  [, comp_mean_diff := mean_b - mean_a ]
  [, inc_ui_diff := 'Not significant']
  [ inc_lower_b > inc_upper_a, inc_ui_diff := 'Increased incidence']
  [ inc_upper_b < inc_lower_a, inc_ui_diff := 'Decreased incidence']
  [, comp_ui_diff := 'Not significant']
  [ lower_b > upper_a, comp_ui_diff := 'Increased completeness']
  [ upper_b < lower_a, comp_ui_diff := 'Decreased completeness']
)
notifs_wide_sf <- merge(
  x = district_sf, y = notifs_wide, by = c('ADM1_EN', 'ADM2_EN', 'uid')
)

## Map relative incidence and compare --------------------------------------------------->

incidence_breaks <- seq(0, 1000, by = 250)
incidence_limits <- range(incidence_breaks)
incidence_colors <- RColorBrewer::brewer.pal(name = 'Greens', n = 9)

country_outline <- geom_sf(data=outline_sf, fill=NA, lwd=.25, linetype=2, color='#222222')

incidence_duo <- ggplot() + 
  facet_wrap('model', ncol = 2) +
  geom_sf(data = notifs_long_sf, aes(fill = inc_mean * 1e5), lwd = 0.25, color = '#222222') +
  country_outline +
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(n = 9, name = 'Greens'), breaks = incidence_breaks,
    limits = incidence_limits, oob = scales::squish, labels = scales::comma
  ) + 
  labs(
    title = glue::glue('Modeled TB Incidence (Mean Estimate, {comparison_year})'),
    fill = "TB incidence\nper 100,000",
    x = NULL, y = NULL
  ) + 
  theme_minimal() +
  theme(
    axis.text=element_blank(), axis.ticks=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )
png(file.path(comparison_dir, 'incidence_maps.png'), height=6, width=10, units='in', res=300)
print(incidence_duo)
dev.off()


mean_diff_breaks <- seq(-500, 500, by = 250)
mean_diff_limits <- range(mean_diff_breaks)
mean_diff_labels <- sapply(mean_diff_breaks, function(x) if(x>0) paste0('+',x) else x)

incidence_comparison_mean <- ggplot() + 
  geom_sf(data = notifs_wide_sf, aes(fill = inc_mean_diff * 1e5), lwd = 0.25, color = '#222222') +
  country_outline + 
  scale_fill_gradientn(
    colors = RColorBrewer::brewer.pal(name = 'RdBu', n = 9) |> rev(),
    breaks = mean_diff_breaks, limits = mean_diff_limits, labels = mean_diff_labels,
    oob = scales::squish
  ) + 
  labs(
    title = NULL, fill = 'Difference in mean\nincidence estimate\n(per 100k)',
    x = NULL, y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text=element_blank(), axis.ticks=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )

ui_breaks <- c('blue', 'white', 'red')
names(ui_breaks) <- c('Decreased incidence', 'Not significant', 'Increased incidence')

incidence_comparison_uis <- ggplot() + 
  geom_sf(data = notifs_wide_sf, aes(fill = inc_ui_diff), lwd = 0.25, color = '#222222') +
  country_outline + 
  scale_fill_manual(values = ui_breaks) +
  labs(title = NULL, fill = 'Significance\nof differences', x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text=element_blank(), axis.ticks=element_blank(),
    panel.grid.major = element_line(colour = 'transparent'),
  )  

# Combine plots and save
png(
  file.path(comparison_dir, 'incidence_comparisons.png'),
  height=9, width=6, units='in', res=300
)
overall_title <- glue(
  'Difference in estimated incidence between models in {comparison_year}',
  ':\n({version_labels$B}) - ({version_labels$A})'
)
grid.arrange(
  ggplotGrob(incidence_comparison_mean), ggplotGrob(incidence_comparison_uis),
  layout_matrix=matrix(1:2, ncol=1), top = overall_title
)
dev.off()
