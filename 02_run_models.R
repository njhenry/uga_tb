## #######################################################################################
##
## RUN TMB MODEL FOR UGANDA TB JOINT ESTIMATION
##
## AUTHOR: Nat Henry
## CREATED: 2 September 2021
## PURPOSE: Execute TMB model and fit outputs
##
## #######################################################################################

load_libs <- c('data.table','Matrix','TMB','glue','INLA','matrixStats','tictoc','optimx')
invisible(lapply(load_libs, library, character.only=TRUE))

# Settings
data_version <- '20211003'
prev_cov_names <- c('intercept','hdi_norm','cattle_norm')
comp_cov_names <- c('intercept','tt_hcf_norm','ntl_harm_norm','ghslurbanicity_norm')
trim_comp_covariates <- TRUE
prev_data_only <- FALSE
duration_intercept <- 1.88
duration_slope <- 0.67

## Set input and output paths
work_dir <- '{REDACTED}'
prepped_dir <- file.path(work_dir, 'prepped_data', data_version)
model_dir <- file.path(work_dir, 'model_results', data_version)
dir.create(model_dir, showWarnings=FALSE)
repo_fp <- '{REDACTED}'
source("model_functions.R")


## Load prepped data; create TMB inputs ------------------------------------------------->

prev_data <- fread(file.path(prepped_dir, 'prev_data.csv'))
notif_data <- fread(file.path(prepped_dir, 'notif_data.csv'))
adjmat <- readRDS(file.path(prepped_dir, 'adjmat.RDS'))
dist_dt <- fread(file.path(prepped_dir, 'dist_dt.csv'))

# Create ICAR precision matrix from adjacency matrix
q_icar <- icar_precision_from_adjacency(W=adjmat)

# Order the districts DT correctly
dist_dt <- dist_dt[order(year, uid)]
dist_dt[, intercept := 1 ]
# Fixed effects - TB prevalence
covs_prev <- as.matrix(dist_dt[year == 2016, ..prev_cov_names])
# Fixed effects - notifications completeness
covs_comp <- as.matrix(dist_dt[year %in% c(2017, 2019), ..comp_cov_names])

# OPTION - TRIM COMPLETENESS COVARIATES
if(trim_comp_covariates){
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

# DATA STACK FOR TMB
tmb_data_stack <- list(
  holdout = 0,
  # Prevalence data
  prev_y = prev_data$tb_cases,
  prev_n = prev_data$n,
  idx_loc_prev = prev_data$uid - 1,
  idx_holdout_prev = rep(1, nrow(prev_data)),
  # Notifications data
  notif_y = notif_data$notif_adj,
  notif_n = notif_data$pop_adj,
  idx_loc_notif = notif_data$uid - 1,
  idx_year_notif = (notif_data$year - 2017)/2,
  idx_holdout_notif = rep(1, nrow(notif_data)),
  zero_centered_years = c(-0.5, 0.5),
  # Function for TB duration, determining the incidence-to-prevalence ratio
  # Calculated as: Duration = <intercept> - <slope> * completeness
  duration_intercept = duration_intercept,
  duration_slope = duration_slope,
  # TB prevalence covariates
  X_prev = covs_prev,
  # Notification completeness covariates
  X_comp = covs_comp,
  # ICAR precision matrix for Ugandan districts
  Q_icar = q_icar
)

# PARAMETERS FOR TMB
num_locs <- nrow(covs_prev)
params_list <- list(
  # Defining prevalence surface
  beta_prev = c(-10, rep(0, ncol(covs_prev)-1)),
  log_tau_loc_prev = 0,
  logit_phi_loc_prev = 1,
  # Defining notification completeness surface
  beta_comp = rep(0, ncol(covs_comp)),
  log_tau_loc_comp = 0,
  logit_phi_loc_comp = 1,
  beta_comp_t = 1,
  log_tau_nu_comp = 0,
  # Random effects
  res_prev = rep(0, num_locs),
  res_comp = rep(0, num_locs),
  nu_comp = rep(0, num_locs)
)

# DROP NOTIFICATIONS DATA IF SPECIFIED
map_list <- list()
tmb_random <- c('res_prev')
if(prev_data_only == TRUE){
  # Drop input data
  null_fields <- c('notif_y','notif_n','idx_loc_notif','idx_year_notif','idx_holdout_notif')
  for(null_field in null_fields) tmb_data_stack[[null_field]] <- numeric(0)
  # Fix parameters related to completeness
  fix_params <- c(
    'beta_comp', 'log_tau_loc_comp', 'logit_phi_loc_comp', 'beta_comp_t',
    'log_tau_nu_comp', 'res_comp', 'nu_comp'
  )
  for(fix_param in fix_params){
    map_list[[fix_param]] <- as.factor(rep(NA, times = length(params_list[[fix_param]])))
  }
  tmb_random <- setdiff(tmb_random, fix_params)
}


## FIT TMB MODEL ------------------------------------------------------------------------>

model_fit <- setup_run_tmb(
  tmb_data_stack = tmb_data_stack, params_list = params_list,
  tmb_random = tmb_random, tmb_map = map_list,
  template_fp = 'joint_completeness_tmb_model.cpp',
  tmb_outer_maxsteps = 1.5E5, tmb_inner_maxsteps = 1E5, parallel_model = TRUE,
  optimization_method = 'L-BFGS-B',
  model_name="UGA joint model", verbose=TRUE, inner_verbose=FALSE
)
sdrep <- sdreport(model_fit$obj, bias.correct = TRUE, getJointPrecision = TRUE)

dist_dt[, zero_centered_year := (year - 2018) / 2]
prev_template <- dist_dt[year == 2016, ]
comp_template <- dist_dt[year %in% c(2017, 2019), ]
model_preds <- generate_draws(
  tmb_sdreport = sdrep, prev_template = prev_template, comp_template = comp_template,
  prev_cov_names = prev_cov_names, comp_cov_names = comp_cov_names,
  prev_data_only = prev_data_only, num_draws = 1000
)

# Summarize draws
keep_cols_prev <- c('uid','year','district',prev_cov_names)
prev_summ <- cbind(
  prev_template[, ..keep_cols_prev],
  summarize_draws(model_preds$pred_draws_prevalence)
)
if(prev_data_only == FALSE){
  keep_cols_comp <- c('uid','year','district',comp_cov_names)
  comp_summ <- cbind(
    comp_template[, ..keep_cols_comp],
    summarize_draws(model_preds$pred_draws_completeness)
  )
}

# Save all to file
saveRDS(model_preds, file=file.path(model_dir, 'model_preds.RDS'))
saveRDS(model_fit, file=file.path(model_dir, 'model_fit.RDS'))
saveRDS(sdrep, file=file.path(model_dir, 'sdrep.RDS'))
fwrite(prev_summ, file=file.path(model_dir, 'prev_summ.csv'))
if(prev_data_only == FALSE) fwrite(comp_summ, file=file.path(model_dir, 'comp_summ.csv'))
