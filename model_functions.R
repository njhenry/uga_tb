## MODEL RUNNING AND POSTESTIMATION FUNCTIONS ------------------------------------------->


#' Set up and run TMB
#'
#' @description Generic TMB model run handler. Sets up the ADFun object, applies
#'   model speedups and fixes as specified, and optimizes using `nlminb`. This
#'   is meant to be a helper function run by more specific model functions.
#'
#' @param tmb_data_stack List containing all data inputs expected in the TMB
#'   CPP file
#' @param params_list List containing all parameters expected in the TMB CPP file
#' @param tmb_random Character vector containing all random effects that will be
#'   optimized in the inner optimizer
#' @param tmb_map Named list containing parameters that will be treated in a
#'   particular way by the optimizer
#' @param template_fp Filepath containing the TMB objective function
#' @param tmb_outer_maxsteps Max number of steps taken by the outer optimizer
#' @param tmb_inner_maxsteps Max number of steps taken by the inner optimizer
#'   in a single outer optimizer step
#' @param parallel_model [bool, default FALSE] Is the model implemented in parallel? If
#'   TRUE, opens multiple OMP threads before fitting
#' @param optimization_method [char, default 'nlminb'] Outer optimization method to use
#'   for fitting, implemented in the optimx library. Recommended options include 'nlminb'
#'   and 'L-BFGS-B'
#' @param model_name [char, default "model"] name of the model
#' @param verbose [boolean, default FALSE] Should this function return logging
#'   information about the stage of model fitting, including the outer optizimer
#'   sampling? This will be somewhat long (>100 lines)
#' @param inner_verbose [boolean, default FALSE] Should this function return
#'   logging information about inner optimizer sampling? This can be useful for
#'   debugging, but will show very verbose (>10k lines) output.
#'
#' @return list of two objects: obj (ADFunction object), and opt (optimized
#'   nlminb object)
#'
#' @import TMB glue tictoc optimx
#' @export
setup_run_tmb <- function(
  tmb_data_stack, params_list, tmb_random, tmb_map, template_fp, tmb_outer_maxsteps = 1E3,
  tmb_inner_maxsteps = 1E3, parallel_model = FALSE, optimization_method = 'nlminb',
  model_name="model", verbose=FALSE, inner_verbose=FALSE
){
  # Helper function to send a message only if verbose
  vbmsg <- function(x) if(verbose) message(x)
  # Setup
  vbmsg(paste0(c("\n",rep("*",nchar(model_name)+14)),collapse=''))
  vbmsg(glue::glue("***  {model_name} RUN  ***"))
  vbmsg(paste0(c(rep("*",nchar(model_name)+14),"\n"),collapse=''))

  # Compile TMB C++ template
  TMB::compile(template_fp)
  current_dir <- getwd()
  setwd(dirname(template_fp))
  compiled_path <- tools::file_path_sans_ext(basename(template_fp))
  dyn.load(TMB::dynlib(compiled_path))

  if(parallel_model){
    # Set up openmp threads
    threads <- system('echo $OMP_NUM_THREADS', intern = TRUE)
    if(threads != '') {
      vbmsg(sprintf('Detected %s threads in OMP environmental variable.',threads))
      TMB::openmp(as.numeric(threads))
    } else {
      vbmsg("Did not detect environmental OMP variable, defaulting to 4 cores. \n
             You can set this using OMP_NUM_THREADS.")
      TMB::openmp(4)
    }
  }

  # Make Autodiff function
  vbmsg("Constructing ADFunction...")
  tictoc::tic("  Making Model ADFun")
  obj <- TMB::MakeADFun(
    data = tmb_data_stack, parameters = params_list, random = tmb_random,
    map = tmb_map, DLL = compiled_path, silent = !inner_verbose,
    inner.control = list(trace=inner_verbose) # tol=1E-11
  )
  obj$env$tracemgc <- as.integer(verbose)
  tictoc::toc()

  # Optimize using the specified outer optimizer, implemented in optimx
  tictoc::tic("  Optimization")
  vbmsg(glue("\n** OPTIMIZING USING METHOD {optimization_method} **"))
  opt <- optimx::optimx(
    par = obj$par, fn = function(x) as.numeric(obj$fn(x)), gr = obj$gr,
    method = optimization_method,
    itnmax = tmb_outer_maxsteps,
    hessian = FALSE,
    control = list(
      trace = as.integer(verbose), follow.on = TRUE,
      dowarn = as.integer(verbose), maxit = tmb_inner_maxsteps,
      factr = 1E-10
    )
  )
  conv_code <- opt$convcode
  vbmsg(glue::glue(
    "{model_name} optimization finished with convergence code {conv_code}.\n"
  ))
  # Clean up
  tictoc::toc()
  setwd(current_dir)
  vbmsg(glue::glue("*** {model_name} RUN COMPLETE **************\n\n"))
  return(list(obj=obj, opt=opt))
}


#' Take multivariate normal draws given a mean vector and precision matrix
#'
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#'
#' @return length(mu) by n.sims matrix of parameter draws
#'
#' @import Matrix
#' @export
rmvnorm_prec <- function(mu, prec, n.sims) {
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv = Matrix::Cholesky(prec, super=TRUE)
  return(mu + solve(as(L_inv, 'pMatrix'), solve(t(as.matrix(as(L_inv, 'Matrix'))), z)))
}


#' Helper function: rbind a matrix N times
replicate_mat_n_times <- function(m, n_times = 1){
  matrix( rep(t(m), n_times), ncol = ncol(m), byrow = TRUE )
}


#' Create post-estimation predictive draws
#'
#' @description Given the outputs from a fitted TMB model object, create
#'   an object with posterior predictive draws for all groupings specified by a
#'   template data.table
#'
#' @param tmb_sdreport output of `TMB::sdreport()` on the fitted model object.
#'   Should include a joint precision matrix (by specifying
#'   `getJointPrecision = TRUE` in the call to `sdreport()`). This object will be
#'   parsed to check for fixed effects, random effects, and the Fourier time
#'   series terms.
#' @param data_template Prepped data with covariates, in random effects order
#' @param num_draws [int] How many posterior predictive samples to take?
#' @param covariate_names [char] All covariate field names, including 'intercept'
#'
#' @return A named list with three items:
#'    - 'param_names': Vector of parameter names in the order they have been
#'         extracted
#'    - 'param_draws': Matrix of parameter draws
#'    - 'pred_draws': Matrix of mortality predictive draws, taken at the
#'         observation points specified in the `template_dt`
#'
#' @import data.table
#' @export
generate_draws <- function(
  tmb_sdreport, prev_template, comp_template, prev_cov_names, comp_cov_names,
  prev_data_only = FALSE, num_draws = 250
){
  # Copy input data
  p_templ <- data.table::copy(prev_template)
  c_templ <- data.table::copy(comp_template)
  num_locs <- nrow(p_templ)
  num_years <- nrow(c_templ) / nrow(p_templ)

  # Get parameter names
  mu <- c(tmb_sdreport$par.fixed, tmb_sdreport$par.random)
  parnames <- names(mu)

  ## Input data checks
  # Check that joint precision matrix was created
  if(!"jointPrecision" %in% names(tmb_sdreport)) stop("Missing joint precision matrix")
  # Check that all the covariate names are there
  missing_covs_prev <- setdiff(prev_cov_names, colnames(p_templ))
  if(length(missing_covs_prev) > 0){
    stop("Missing prevalence covariates: ", paste0(missing_covs_prev, collapse=', '))
  }
  if(length(prev_cov_names) != sum(parnames == 'beta_prev')){
    stop("Wrong number of prevalence covariates in model fit")
  }
  if(prev_data_only == FALSE){
    missing_covs_comp <- setdiff(comp_cov_names, colnames(c_templ))
    if(length(missing_covs_comp) > 0){
      stop("Missing completeness covariates: ", paste0(missing_covs_comp, collapse=', '))
    }
    if(length(comp_cov_names) != sum(parnames == 'beta_comp')){
      stop("Wrong number of completeness covariates in model fit")
    }
  }

  ## Get parameter draws
  message(sprintf(" - Generating %i parameter draws...", num_draws))
  prec_mat <- tmb_sdreport$jointPrecision
  if(any(colnames(prec_mat) != parnames )) stop("Issue with parameter ordering")
  # If the matrix is not initially positive definite, try dropping hyperparameters
  tryCatch({
    keep_pars <<- rep(TRUE, length(mu))
    param_draws <<- rmvnorm_prec(
      mu = mu[keep_pars],
      prec = prec_mat[keep_pars, keep_pars],
      n.sims = num_draws
    )
  }, error = function(e){
    keep_pars <<- grepl('^beta|^Z|^e', parnames)
    param_draws <<- rmvnorm_prec(
      mu = mu[keep_pars],
      prec = prec_mat[keep_pars, keep_pars],
      n.sims = num_draws
    )
  })
  parnames <- parnames[keep_pars]
  rownames(param_draws) <- parnames

  ## Generate predictive draws from parameter draws
  ## 1) Prevalence
  p_fes <- as.matrix(p_templ[, ..prev_cov_names]) %*% param_draws[parnames=='beta_prev',]
  p_spatial_res <- param_draws[parnames == 'Z_prev', ]
  p_iid_res <- param_draws[parnames == 'e_prev', ]
  p_preds <- (p_fes + p_spatial_res + p_iid_res) |> as.matrix() |> exp()
  ## 2) Completeness
  if(prev_data_only == TRUE){
    c_preds <- NULL
  } else {
    c_fes <- as.matrix(c_templ[, ..comp_cov_names]) %*% param_draws[parnames=='beta_comp',]
    c_spatial_res <- replicate_mat_n_times(param_draws[parnames == 'Z1_comp', ], num_years)
    c_iid_res <- param_draws[parnames == 'e_comp', ]
    c_r_slopes <- replicate_mat_n_times(
      replicate_mat_n_times(t(param_draws[parnames == 'beta_comp_t',]), num_locs) +
        param_draws[parnames == 'Z2_comp' ],
      n_times = num_years
    )
    c_preds <- (c_fes + c_spatial_res + c_iid_res + c_templ$zero_centered_year * c_r_slopes) |>
      as.matrix() |>
      plogis()
  }

  # Return parameter draws and predictive draws
  return(list(
    param_names = parnames,
    param_draws = param_draws,
    pred_draws_prevalence = p_preds,
    pred_draws_completeness = c_preds
  ))
}


#' Summarize predictive draws
#'
#' @description Summarize the mean and select quantiles of a matrix of posterior
#'   draws, where draws are stored in columns
#'
#' @param draws [matrix] matrix of dimensions (num obs) by (num draws)
#'
#' @return data.table with columns 'mean','median','upper','lower' (of 95% UI)
#'
#' @import data.table matrixStats
#' @export
summarize_draws <- function(draws){
  if('data.table' %in% class(draws)) draws <- as.matrix(draws)
  summs <- cbind(
    rowMeans(draws), matrixStats::rowQuantiles(draws, probs=c(0.5, 0.025, 0.975))
  )
  colnames(summs) <- c('mean','median','lower','upper')
  return(as.data.table(summs))
}

