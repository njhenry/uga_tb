## #######################################################################################
##
## RUN TMB MODEL FOR UGANDA TB JOINT ESTIMATION
##
## AUTHOR: Nat Henry
## CREATED: 2 September 2021
## PURPOSE: Execute TMB model and fit outputs
##
## #######################################################################################

versions <- list(prepped_data = '20231024', model_results = '20231024_v2_oos')

# Load packages
load_libs <- c('data.table','Matrix','TMB','glue','matrixStats','tictoc','optimx')
invisible(lapply(load_libs, library, character.only=TRUE))

# Load custom packages
repos_dir <- '~/repos/'
devtools::load_all(file.path(repos_dir, 'versioning'))
source(file.path(repos_dir, 'uga_tb/model_functions.R'))

# Load config file and save
config_fp <- file.path(repos_dir, 'uga_tb/config.yaml')
config <- versioning::Config$new(config_list = config_fp, versions = versions)

config$get_dir_path('model_results') |> dir.create()
config$write_self('model_results')


## 01. Load model inputs ---------------------------------------------------------------->

dist_dt <- config$read('prepped_data', 'covariates')
prev_data <- config$read('prepped_data', 'prev_data')
notif_data <- config$read('prepped_data', 'notif_data')
adjmat <- config$read('prepped_data', 'adjmat')
model_years <- config$get('model_years')
prev_cov_names <- config$get('covs', 'prev')
comp_cov_names <- config$get('covs', 'comp')
# Different execution flows for in-sample vs. out-of-sample and full vs. prevalence only
prev_data_only <- config$get('prev_data_only')
out_of_sample <- config$get('out_of_sample')

# Order the districts DT correctly
dist_dt <- copy(dist_dt[year %in% model_years, ][order(year, uid)])
dist_dt$intercept <- 1
# Fixed effects - TB prevalence
covs_prev <- as.matrix(dist_dt[year == 2016, ..prev_cov_names])
# Fixed effects - notifications completeness
covs_comp <- as.matrix(dist_dt[year %in% model_years, ..comp_cov_names])
notif_data <- notif_data[year %in% model_years, ][order(year, uid)]

# OPTION - TRIM COMPLETENESS COVARIATES
if(config$get('covariate_settings', 'trim_comp_covariates')){
  n_yr <- length(unique(notif_data$year))
  n_d <- length(unique(notif_data$uid))
  ps_overlap_rows <- rep(prev_data$uid, n_yr) + rep(0:(n_yr-1), each=nrow(prev_data))*n_d
  # Subset to only values observed overlapping with a prevalence survey
  for(i_col in 1:ncol(covs_comp)){
    col_name <- comp_cov_names[i_col]
    cov_vals <- covs_comp[, i_col]
    trim_min <- min(cov_vals[ps_overlap_rows])
    trim_max <- max(cov_vals[ps_overlap_rows])
    covs_comp[cov_vals < trim_min ,i_col] <- trim_min
    covs_comp[cov_vals > trim_max ,i_col] <- trim_max
    dist_dt[[col_name]][ dist_dt[[col_name]] < trim_min] <- trim_min
    dist_dt[[col_name]][ dist_dt[[col_name]] > trim_max] <- trim_max
  }
}
dist_dt <- copy(dist_dt)

# Zero-center years
z_c_year_index <- data.table::data.table(
  year = model_years,
  zero_centered_year = seq(-0.5, 0.5, length.out = length(model_years))
)
dist_dt[z_c_year_index, zero_centered_year := i.zero_centered_year, on = 'year']

# DATA STACK FOR TMB
tmb_data_stack <- list(
  holdout = 0,
  # Prevalence data
  prev_y = prev_data$ptb_bc,
  prev_n = prev_data$sampsize,
  idx_loc_prev = prev_data$uid - 1,
  idx_holdout_prev = seq_len(nrow(prev_data)),
  # Notifications data
  notif_y = notif_data$notif_count,
  notif_n = notif_data$pop_over_15,
  idx_loc_notif = notif_data$uid - 1,
  idx_year_notif = notif_data$year - min(notif_data$year),
  zero_centered_years = z_c_year_index$zero_centered_year,
  # Function for TB duration, determining the incidence-to-prevalence ratio
  # Calculated as: Duration = <intercept> - <slope> * completeness
  duration_intercept = config$get("duration", "intercept"),
  duration_slope = config$get("duration", "slope"),
  # TB prevalence covariates
  X_prev = covs_prev,
  # Notification completeness covariates
  X_comp = covs_comp,
  # Adjacency matrix
  adjmat = adjmat
)
# Add priors from config
priors_list <- config$get('priors')
tmb_data_stack <- c(tmb_data_stack, priors_list)

# PARAMETERS FOR TMB
num_locs <- nrow(covs_prev)
num_years <- length(model_years)
params_list <- list(
  # Defining prevalence surface
  beta_prev = c(-10, rep(0, ncol(covs_prev)-1)),
  beta_comp = rep(0, ncol(covs_comp)),
  beta_comp_t = 1,
  # Hyperparameters
  logit_rho_Z_prev = 0,
  logit_rho_Z1_comp = 0,
  logit_rho_Z2_comp = 0,
  log_tau_Z_prev = 0,
  log_tau_Z1_comp = 0,
  log_tau_Z2_comp = 0,
  log_tau_e_prev = 0,
  log_tau_e_comp = 0,
  # Random effects
  Z_prev = rep(0, num_locs),
  Z1_comp = rep(0, num_locs),
  Z2_comp = rep(0, num_locs),
  e_prev = rep(0, num_locs),
  e_comp = rep(0, num_locs * num_years)
)

# DROP NOTIFICATIONS DATA IF SPECIFIED
map_list <- list()
tmb_random <- c('Z_prev', 'Z1_comp', 'Z2_comp', 'e_prev', 'e_comp')
if(prev_data_only){
  # Drop input data
  null_fields <- c('notif_y','notif_n','idx_loc_notif','idx_year_notif')
  for(null_field in null_fields) tmb_data_stack[[null_field]] <- numeric(0)
  # Fix parameters related to completeness
  fix_params <- c(
    'beta_comp', 'beta_comp_t', 'logit_rho_Z1_comp', 'logit_rho_Z2_comp', 
    'log_tau_Z1_comp', 'log_tau_Z2_comp', 'log_tau_e_comp', 'Z1_comp', 'Z2_comp',
    'e_comp'
  )
  for(fix_param in fix_params){
    map_list[[fix_param]] <- as.factor(rep(NA, times = length(params_list[[fix_param]])))
  }
  tmb_random <- setdiff(tmb_random, fix_params)
}


## FIT TMB MODEL ------------------------------------------------------------------------>

# Prediction templates
prev_template <- copy(dist_dt[year == 2016, ])
comp_template <- copy(dist_dt)

if(out_of_sample){
  # OUT OF SAMPLE EXECUTION
  holdout_ids <- tmb_data_stack$idx_holdout_prev
  prev_data[, holdout_idx := .I ]
  oos_comparison_list <- vector('list', length = length(holdout_ids))

  # Run holding out every prevalence data point, then compare data to OOS prediction
  for(holdout_id in holdout_ids) try({
    tmb_data_stack$holdout <- holdout_id
    model_fit <- setup_run_tmb(
      tmb_data_stack = tmb_data_stack, params_list = params_list,
      tmb_random = tmb_random, tmb_map = map_list,
      template_fp = file.path(repos_dir, 'uga_tb/joint_completeness_tmb_model.cpp'),
      tmb_outer_maxsteps = 1.5E5, tmb_inner_maxsteps = 1E5, parallel_model = FALSE,
      optimization_method = 'L-BFGS-B',
      model_name="UGA joint model", verbose=TRUE, inner_verbose=FALSE
    )
    sdrep <- sdreport(model_fit$obj, bias.correct = TRUE, getJointPrecision = TRUE)
    model_preds <- generate_draws(
      tmb_sdreport = sdrep, prev_template = prev_template, comp_template = comp_template,
      prev_cov_names = prev_cov_names, comp_cov_names = comp_cov_names,
      prev_data_only = prev_data_only, num_draws = 250
    )
    # Compare to held-out prevalence survey data
    merge_cols <- c(config$get('shapefile_settings','ids','adm2'), 'year')
    prev_summ <- cbind(
      prev_template[, ..merge_cols],
      summarize_draws(model_preds$pred_draws_prevalence)
    )
    oos_keep_cols <- c(merge_cols, 'mean', 'median', 'lower', 'upper')
    oos_comparison_list[[holdout_id]] <- merge(
      x = prev_data[holdout_idx == holdout_id,],
      y = prev_summ[, ..oos_keep_cols],
      by = merge_cols
    )
  })
  # Combine list and save
  oos_comparison_dt <- rbindlist(oos_comparison_list)
  config$write(oos_comparison_dt, 'model_results', 'oos_summary')

} else {

  # IN-SAMPLE MODEL EXECUTION
  model_fit <- setup_run_tmb(
    tmb_data_stack = tmb_data_stack, params_list = params_list,
    tmb_random = tmb_random, tmb_map = map_list,
    template_fp = file.path(repos_dir, 'uga_tb/joint_completeness_tmb_model.cpp'),
    tmb_outer_maxsteps = 1.5E5, tmb_inner_maxsteps = 1E5, parallel_model = FALSE,
    optimization_method = 'L-BFGS-B', model_name="UGA joint model",
    verbose = FALSE, inner_verbose = FALSE
  )
  sdrep <- sdreport(model_fit$obj, bias.correct = TRUE, getJointPrecision = TRUE)
  model_preds <- generate_draws(
    tmb_sdreport = sdrep, prev_template = prev_template, comp_template = comp_template,
    prev_cov_names = prev_cov_names, comp_cov_names = comp_cov_names,
    prev_data_only = prev_data_only, num_draws = 250
  )
  # Summarize draws
  keep_cols_prev <- c(config$get('shapefile_settings','ids','adm2'), 'year', prev_cov_names)
  prev_summ <- cbind(
    prev_template[, ..keep_cols_prev],
    summarize_draws(model_preds$pred_draws_prevalence)
  )
  if(prev_data_only == FALSE){
    keep_cols_comp <- c(config$get('shapefile_settings','ids','adm2'), 'year', comp_cov_names)
    comp_summ <- cbind(
      comp_template[, ..keep_cols_comp],
      summarize_draws(model_preds$pred_draws_completeness)
    )
  }
  # Save all to file
  config$write(model_fit, 'model_results', 'model_fit')
  config$write(model_preds, 'model_results', 'model_preds')
  config$write(prev_summ, 'model_results', 'prevalence_summary')
  if(prev_data_only == FALSE) config$write(comp_summ, 'model_results', 'completeness_summary')
}
