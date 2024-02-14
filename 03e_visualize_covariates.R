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

load_libs <- c('data.table','sf','ggplot2','grid','gridExtra','glue','terra')
invisible(lapply(load_libs, library, character.only=TRUE))

# Load custom packages
repos_dir <- '~/repos/'
devtools::load_all(file.path(repos_dir, 'versioning'))

config <- versioning::Config$new(
  config_list = file.path('~/temp_data/uga/model_results', MODEL_VERSION, 'config.yaml')
)

viz_dir <- file.path(config$get_dir_path('model_results'), 'viz')
dir.create(viz_dir, showWarnings = FALSE)

# Load input data
model_years <- config$get('model_years')
final_year <- max(model_years)

covs_dt <- config$read("prepped_data", "covariates")
inc_summ <- config$read('model_results', 'incidence_summary')[year == final_year, ]
comp_summ <- config$read("model_results", "completeness_summary")[year == final_year, ]

model_predictions <- config$read('model_results', 'model_preds')


## PREPARE FIXED EFFECTS TABLE ---------------------------------------------------------->

fe_names <- c(config$get('covs', 'prev'), config$get('covs','comp'), 'comp_time')
fe_matrix <- model_predictions$param_draws[
  which(model_predictions$param_names %in% c('beta_prev', 'beta_comp', 'beta_comp_t')),  
] |> as('matrix')
fe_table <- data.table::data.table(
  fe_name = fe_names,
  mean = rowMeans(fe_matrix),
  lower = matrixStats::rowQuantiles(fe_matrix, probs = 0.025),
  upper = matrixStats::rowQuantiles(fe_matrix, probs = 0.975)
)
data.table::fwrite(fe_table, file = file.path(viz_dir, 'fes_summary_table.csv'))


## Prepare incidence covariate scatter plot --------------------------------------------->

inc_cov_labels <- c(
  hh_crowding = 'Household crowding',
  log_ntl = 'Nighttime light brightness,\nlog-transformed',
  hiv = 'HIV prevalence',
  refugees_pct = 'Refugees per capita',
  cattle_pc = 'Cattle per capita'
)
inc_covs <- names(inc_cov_labels)

inc_covs_table <- copy(covs_dt)[year == final_year, c('uid', inc_covs), with = F] |>
  melt(id.var = 'uid') |>
  merge(y = inc_summ[, .(uid, mean, lower, upper)], by = 'uid')

inc_covs_table$var_label <- inc_cov_labels[inc_covs_table$variable]

inc_cov_scatter <- ggplot(
  data = inc_covs_table,
  aes(x = value, y = mean*1e5, ymin = lower*1e5, ymax = upper*1e5)
) + 
  facet_wrap('var_label', ncol = 3, scales = 'free_x') + 
  geom_crossbar(color = rgb(.2, .2, .2, 0.4)) + 
  geom_point(aes(color = var_label), alpha = .5) + 
  labs(
    x = 'Covariate value',
    y = glue::glue('Estimated TB incidence per 100,000 ({final_year})')
  ) +
  theme_bw() + 
  theme(legend.position = 'none')
png(file.path(viz_dir, 'inc_covs_scatter.png'), height = 6, width = 8, units = 'in', res = 300)
print(inc_cov_scatter)
dev.off()
