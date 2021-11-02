// ///////////////////////////////////////////////////////////////////////////////////////
//
// Author: Nat Henry
// Created: 2 September 2021
// Purpose: TMB objective function for Uganda TB notifications completeness model
//
// ///////////////////////////////////////////////////////////////////////////////////////


#include <TMB.hpp>
using namespace density;
using Eigen::SparseMatrix;


// HELPER FUNCTIONS --------------------------------------------------------------------->

// Function for rescaling a precision matrix to have standard deviation sigma
//
// Parameter Q: Unscaled precision matrix
// Parameter sigma: Standard deviation to scale to
//
template<class Type>
SparseMatrix<Type> scale_precision(SparseMatrix<Type> Q, Type sigma){
  SparseMatrix<Type> Q_scaled = Q / (sigma * sigma);
  return Q_scaled;
}

// Function to create an IID precision matrix (AKA a scaled identity matrix)
//
// Parameter dim: Number of rows (and columns) for the precision matrix
// Parameter sigma: Standard deviation of the iid process
//
template<class Type>
SparseMatrix<Type> iid_precision(int dim, Type sigma = 1.0){
  SparseMatrix<Type> I(dim, dim);
  for(int ii=0; ii < dim; ii++){
    I.insert(ii, ii) = 1.0;
  }
  SparseMatrix<Type> I_scaled = scale_precision(I, sigma);
  return I_scaled;
}

// Function for preparing a precision matrix corresponding to a BYM2 spatial model,
//   based on a scaled ICAR precision matrix. For more details, see:
//   Riebler et al. (2016). An intuitive Bayesian sptial model for disease mapping that
//     accounts for scaling. Statistical methods in medical research, 25(4):1145-65.
//
// Parameter Q_icar: A precision matrix corresponding to an intrinsic correlated
//   autoregressive (ICAR) model in space, scaled to have generalized variance 1
// Parameter phi: A mixing parameter indicating the relative contribution of spatial and
//   IID variation, strictly between 0 and 1.
// Parameter sigma: Standard deviation of the LCAR process
//
template<class Type>
SparseMatrix<Type> bym2_precision(SparseMatrix<Type> Q_icar, Type phi, Type sigma = 1.0){
  SparseMatrix<Type> I = iid_precision(Q_icar.rows(), Type(1.0));
  SparseMatrix<Type> Q_bym2 = phi * Q_icar + (1 - phi) * I;
  SparseMatrix<Type> Q_bym2_scaled = scale_precision(Q_bym2, sigma);
  return Q_bym2_scaled;
}

// Function for evaluating a penalized complexity (PC) prior, defined by parameters mu and
//   alpha. The interpretation of a PC prior for some term tau is P(tau > mu) = alpha,
//   mu > 0, 0 < alpha < 1. Typically used to interpret precision parameters:
//   tau = 1/sigma**2. For more documentation, see the INLA manual:
//   https://inla.r-inla-download.org/r-inla.org/doc/prior/pc.prec.pdf
//
// Parameter tau: Value to evaluate against a penalized complexity prior
// Parameter mu: Defines the PC prior - evaluated compared to tau. Tau and mu both > 0
// Parameter alpha: Defines the PC prior - probability P(tau > mu). 0 < alpha < 1 - often
//   0.05.
//
template<class Type>
Type pc_prior_log_density(Type tau, Type mu, Type alpha){
  Type lambda = -log(alpha) / mu;
  // Evaluate prior
  // Density = lambda / 2 * tau**(-3/2) * exp(-lambda / sqrt(tau))
  Type log_dens = log(lambda/2.) - 3./2. * log(tau) - lambda / sqrt(tau);
  return log_dens;
}

// Robust inverse logit that sets min and max values to avoid numerical instability
template<class Type>
Type invlogit_robust(Type x){
  if (x < -20.723){
    x = -20.723; // corresponds to p=1e-9
  } else if ( x > 20.723 ){
    x = 20.723;  // cooresponds to p=1-1e-9
  }
  return 1 / (1 + exp( -1.0 * x ));
}


// OBJECTIVE FUNCTION ------------------------------------------------------------------->

template<class Type>
Type objective_function<Type>::operator() () {

  // INPUT DATA ------------------------------------------------------------------------->

  // OPTION: Holdout number
  // Any observation where `idx_holdout_vr` or `idx_holdout_bh` is equal to `holdout`
  //   will have the corresponding data type excluded from this model fit
  // Holdouts are 1-indexed, so if holdout==0 all data will be used
  DATA_INTEGER(holdout);

  // Input prevalence data
  DATA_VECTOR(prev_y);
  DATA_VECTOR(prev_n);
  // Index: district associated with each observation
  DATA_IVECTOR(idx_loc_prev);
  // Index: holdout associates with each observation
  DATA_IVECTOR(idx_holdout_prev);

  // Input notifications data
  DATA_VECTOR(notif_y);
  DATA_VECTOR(notif_n);
  // Index: district associated with each notification
  DATA_IVECTOR(idx_loc_notif);
  // Index: zero-indexed year associated with each notification
  DATA_IVECTOR(idx_year_notif);
  // Index: holdout associated with each notification
  DATA_IVECTOR(idx_holdout_notif);

  // Zero-centered vector of years, normalized to have range 1
  // Used for estimating random slopes
  // Length: number of years
  DATA_VECTOR(zero_centered_years);

  // TB duration as a function of notification completeness
  // Formula is:
  //   Duration = <Intercept> - <slope> * completeness
  //   0 < slope < intercept, 0 < completeness < 1
  DATA_SCALAR(duration_intercept);
  DATA_SCALAR(duration_slope);
  // TB prevalence fixed effects: dim (num locs) by (num covs)
  DATA_MATRIX(X_prev);
  // Motification completeness fixed effects: dim (num locs x num years) by (num covs)
  DATA_MATRIX(X_comp);

  // Precision matrix for an ICAR spatial model - used for both spatial surfaces
  DATA_SPARSE_MATRIX(Q_icar);


  // INPUT PARAMETERS ------------------------------------------------------------------->

  // PARAMETERS DETERMINING PREVALENCE SURFACE
  // Fixed effects
  PARAMETER_VECTOR(beta_prev);
  // Log precision of the spatial random effects
  PARAMETER(log_tau_loc_prev);
  // Mixing parameter controlling spatial vs. nonspatial correlation
  PARAMETER(logit_phi_loc_prev);

  // PARAMETERS DETERMINING NOTIFICATION COMPLETENESS SURFACE
  // Fixed effects for notification completeness, NOT including year
  PARAMETER_VECTOR(beta_comp);
  // Log precision of the spatial random effects
  PARAMETER(log_tau_loc_comp);
  // Mixing parameter controlling spatial vs. nonspatial correlation
  PARAMETER(logit_phi_loc_comp);
  // Overall slope on completeness change over time, beta_year
  PARAMETER(beta_comp_t);
  // Hyperparameter: log precision of the completeness random slopes
  PARAMETER(log_tau_nu_comp);

  // REs across prevalence surface - 1D, length (num locs)
  PARAMETER_VECTOR(res_prev);
  // REs across completeness surface - 1D, length (num locs)
  PARAMETER_VECTOR(res_comp);
  // Random slopes by district on completeness change over time, nu
  PARAMETER_VECTOR(nu_comp);

  // Useful indices
  int num_locs = res_prev.size();
  int num_years = zero_centered_years.size();
  int num_obs_prev = prev_y.size();
  int num_notifs = notif_y.size();
  // Standard deviation for completeness random slopes
  Type sigma_loc_prev = exp(log_tau_loc_prev * -0.5);
  Type sigma_loc_comp = exp(log_tau_loc_comp * -0.5);
  Type sigma_nu_comp = exp(log_tau_nu_comp * -0.5);


  // Instantiate joint negative log-likelihood (JNLL) ----------------------------------->

  parallel_accumulator<Type> jnll(this);


  // JNLL CONTRIBUTION FROM PRIORS ------------------------------------------------------>

  // PC hyperpriors for standard deviations of spatial latent surfaces and random slope
  jnll -= pc_prior_log_density(sigma_loc_prev, Type(2.0), Type(0.05));
  jnll -= pc_prior_log_density(sigma_loc_comp, Type(2.0), Type(0.05));
  jnll -= pc_prior_log_density(sigma_nu_comp, Type(0.5), Type(0.05));

  // Spatial effect = BYM2 (scaled CAR) model using municipal neighborhood structure
  // Prevalence
  SparseMatrix<Type> Q_loc_prev = bym2_precision(
    Q_icar, invlogit(logit_phi_loc_prev), sigma_loc_prev
  );
  jnll += GMRF(Q_loc_prev)(res_prev);
  // Notification completeness random effects = BYM2 in space
  SparseMatrix<Type> Q_loc_comp = bym2_precision(
    Q_icar, invlogit(logit_phi_loc_comp), sigma_loc_comp
  );
  jnll += GMRF(Q_loc_comp)(res_comp);
  // Notification completeness random slopes: IID by district
  jnll -= dnorm(nu_comp, Type(0.0), sigma_nu_comp, true).sum();

  // N(mean=0, sd=3) prior for fixed effects
  // Skip the intercept (index 0)
  for(int cov_j_prev = 1; cov_j_prev < beta_prev.size(); cov_j_prev++){
    jnll -= dnorm(beta_prev(cov_j_prev), Type(0.0), Type(3.0), true);
  }
  for(int cov_j_comp = 1; cov_j_comp < beta_comp.size(); cov_j_comp++){
    jnll -= dnorm(beta_comp(cov_j_comp), Type(0.0), Type(3.0), true);
  }
  jnll -= dnorm(beta_comp_t, Type(0.0), Type(3.0), true);


  // JNLL CONTRIBUTION FROM DATA -------------------------------------------------------->

  // Log(True prevalence) = fixed effects + random effects
  // Length: num locs
  vector<Type> fes_prev = X_prev * beta_prev.matrix();
  vector<Type> prev_est(num_locs);
  for(int loc_i = 0; loc_i < num_locs; loc_i++){
    prev_est[loc_i] = exp(fes_prev(loc_i) + res_prev(loc_i));
  }

  // True completeness = inverse.logit(fixed effects + random effects)
  // Length: (num locs) by (num years)
  // Index: all ordered locs from year 1, then year 2, etc
  vector<Type> fes_comp = X_comp * beta_comp.matrix();
  vector<Type> comp_est(num_locs * num_years);
  // Duration: estimated as a function of completeness
  // Same indexing as comp_est
  vector<Type> duration(num_locs * num_years);
  // Populate estimates for completeness and duration
  for(int year_i = 0; year_i < num_years; year_i++){
    for(int loc_i = 0; loc_i < num_locs; loc_i++){
      comp_est[year_i * num_locs + loc_i] = invlogit_robust(
        fes_comp[year_i * num_locs + loc_i] + res_comp[loc_i] +
          (beta_comp_t + nu_comp[loc_i]) * zero_centered_years[year_i]
      );
      duration[year_i * num_locs + loc_i] = (
        duration_intercept - duration_slope * comp_est[year_i * num_locs + loc_i]
      );
    }
  }

  // Incorporate prevalence data
  // Probability density function:
  //   Y^{prev}_i = Poisson(N^{prev}_i * p_{s(i)})
  for(int p_obs_i = 0; p_obs_i < num_obs_prev; p_obs_i++){
    if((idx_holdout_prev(p_obs_i) != holdout) & (prev_n(p_obs_i) > 0.)){
      jnll -= dpois(
        prev_y(p_obs_i),
        prev_est(idx_loc_prev(p_obs_i)) * prev_n(p_obs_i),
        true
      );
    }
  }
  // Incorporate incidence data
  // Probability density function:
  //   Y^{notif}_i = Poisson(p_{s(i)} / D(\pi_{s(i),t(i)}) * N^{notif}_i * \pi_{s(i),t(i)})
  for(int nn_i = 0; nn_i < num_notifs; nn_i++){
    if((idx_holdout_notif(nn_i) != holdout) & (notif_n(nn_i) > 0.)){
      jnll -= dpois(
        notif_y(nn_i),
        prev_est(idx_loc_notif(nn_i)) /
          duration(idx_year_notif(nn_i)*num_locs+idx_loc_notif(nn_i)) *
          notif_n(nn_i) *
          comp_est(idx_year_notif(nn_i)*num_locs+idx_loc_notif(nn_i)),
        true
      );
    }
  }


  // RETURN JNLL ------------------------------------------------------------------------>

  return jnll;

} // END objective function
