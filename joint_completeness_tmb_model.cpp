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
#include "matrix_helpers.hpp"

// PRIOR EVALUATION --------------------------------------------------------------------->

// Read in an R list specifying a prior
// List should contain three objects: "name", "par1", "par2"
template<class Type>
struct prior_type {
  std::string name;
  Type par1;
  Type par2;

  prior_type(SEXP x){
    name = CHAR(STRING_ELT(getListElement(x,"name"), 0));
    par1 = asVector<float>(getListElement(x,"par1"))[0];
    par2 = asVector<float>(getListElement(x,"par2"))[0];
  }
};

// evaluate a prior for an object
template<class Type>
Type evaluate_prior_density(prior_type<Type> prior, Type param, bool log_density = true){
  Type density;

  if(prior.name == "gaussian"){
    density = dnorm(param, prior.par1, prior.par2, log_density);
  } else if(prior.name == "gamma"){
    density = dgamma(param, prior.par1, prior.par2, log_density);
  } else if(prior.name == "beta"){
    density = dbeta(param, prior.par1, prior.par2, log_density);
  } else {
    // Prior name must be one of "gaussian", "gamma", or "beta"
    exit(1);
  }
  return density;
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

  // Adjacency matrix for Ugandan districts
  DATA_SPARSE_MATRIX(adjmat);

  // PRIORS
  DATA_STRUCT(prior_beta_covs, prior_type);
  DATA_STRUCT(prior_beta_comp_t, prior_type);
  DATA_STRUCT(prior_rho_Z_prev, prior_type);
  DATA_STRUCT(prior_rho_Z1_comp, prior_type);
  DATA_STRUCT(prior_rho_Z2_comp, prior_type);

  DATA_STRUCT(prior_sigma_Z_prev, prior_type);
  DATA_STRUCT(prior_sigma_Z1_comp, prior_type);
  DATA_STRUCT(prior_sigma_Z2_comp, prior_type);
  DATA_STRUCT(prior_sigma_e_prev, prior_type);
  DATA_STRUCT(prior_sigma_e_comp, prior_type);


  // INPUT PARAMETERS ------------------------------------------------------------------->

  // PARAMETERS DETERMINING PREVALENCE SURFACE
  // Fixed effects
  PARAMETER_VECTOR(beta_prev);
  PARAMETER_VECTOR(beta_comp);
  PARAMETER(beta_comp_t);

  // Spatial autocorrelation hyperparameters
  PARAMETER(logit_rho_Z_prev);
  PARAMETER(logit_rho_Z1_comp);
  PARAMETER(logit_rho_Z2_comp);

  // Variance hyperparameters
  PARAMETER(log_tau_Z_prev);
  PARAMETER(log_tau_Z1_comp);
  PARAMETER(log_tau_Z2_comp);
  PARAMETER(log_tau_e_prev);
  PARAMETER(log_tau_e_comp);

  // Spatial random intercepts across prevalence surface
  PARAMETER_VECTOR(Z_prev);
  // Spatial random intercepts across completeness surface
  PARAMETER_VECTOR(Z1_comp);
  // Spatial random time slopes across completeness surface
  PARAMETER_VECTOR(Z2_comp);

  // Non-spatial IID random effect on prevalence surface
  // Length: number of locations
  PARAMETER_VECTOR(e_prev);
  // Non-spatial IID random effect on completeness surface
  // Length: number of locations * number of years
  PARAMETER_VECTOR(e_comp);


  // Useful indices --------------------------------------------------------------------->

  int num_locs = Z_prev.size();
  int num_years = zero_centered_years.size();
  int num_obs_prev = prev_y.size();
  int num_notifs = notif_y.size();
  int num_prev_covs = beta_prev.size();
  int num_comp_covs = beta_comp.size();


  // Transform some terms --------------------------------------------------------------->

  Type rho_Z_prev = invlogit(logit_rho_Z_prev);
  Type rho_Z1_comp = invlogit(logit_rho_Z1_comp);
  Type rho_Z2_comp = invlogit(logit_rho_Z2_comp);

  Type sigma_Z_prev = exp(log_tau_Z_prev * -0.5);
  Type sigma_Z1_comp = exp(log_tau_Z1_comp * -0.5);
  Type sigma_Z2_comp = exp(log_tau_Z2_comp * -0.5);
  Type sigma_e_prev = exp(log_tau_e_prev * -0.5);
  Type sigma_e_comp = exp(log_tau_e_comp * -0.5);


  // Instantiate joint negative log-likelihood (JNLL) ----------------------------------->

  parallel_accumulator<Type> jnll(this);


  // JNLL CONTRIBUTION FROM PRIORS ------------------------------------------------------>

  // Priors on fixed effects, skipping intercepts
  for(int prev_cov_i = 1; prev_cov_i < num_prev_covs; prev_cov_i++){
    jnll -= evaluate_prior_density(prior_beta_covs, beta_prev(prev_cov_i));
  }
  for(int comp_cov_i = 1; comp_cov_i < num_comp_covs; comp_cov_i++){
    jnll -= evaluate_prior_density(prior_beta_covs, beta_comp(comp_cov_i));
  }
  jnll -= evaluate_prior_density(prior_beta_comp_t, beta_comp_t);

  // Priors on hyperparameters
  jnll -= evaluate_prior_density(prior_rho_Z_prev, rho_Z_prev);
  jnll -= evaluate_prior_density(prior_rho_Z1_comp, rho_Z1_comp);
  jnll -= evaluate_prior_density(prior_rho_Z2_comp, rho_Z2_comp);

  jnll -= evaluate_prior_density(prior_sigma_Z_prev, sigma_Z_prev);
  jnll -= evaluate_prior_density(prior_sigma_Z1_comp, sigma_Z1_comp);
  jnll -= evaluate_prior_density(prior_sigma_Z2_comp, sigma_Z2_comp);
  jnll -= evaluate_prior_density(prior_sigma_e_prev, sigma_e_prev);
  jnll -= evaluate_prior_density(prior_sigma_e_comp, sigma_e_comp);

  // Spatial random effects
  SparseMatrix<Type> Q_prev = car_precision(adjmat, rho_Z_prev);
  SparseMatrix<Type> Q1_comp = car_precision(adjmat, rho_Z1_comp);
  SparseMatrix<Type> Q2_comp = car_precision(adjmat, rho_Z2_comp);

  jnll += SCALE(GMRF(Q_prev), sigma_Z_prev)(Z_prev);
  jnll += SCALE(GMRF(Q1_comp), sigma_Z1_comp)(Z1_comp);
  jnll += SCALE(GMRF(Q2_comp), sigma_Z2_comp)(Z2_comp);

  // Non-spatial random effects
  jnll -= dnorm(e_prev, Type(0.0), sigma_e_prev, true).sum();
  jnll -= dnorm(e_comp, Type(0.0), sigma_e_comp, true).sum();


  // ESTIMATE PREVALENCE, COMPLETENESS, DURATION, AND INCIDENCE ------------------------->

  // True prevalence = Exp(covariate effects + spatial REs + non-spatial REs)
  // Length: num locs
  vector<Type> fes_prev = X_prev * beta_prev.matrix();
  vector<Type> prev_est(num_locs);
  for(int loc_i = 0; loc_i < num_locs; loc_i++){
    prev_est(loc_i) = exp(fes_prev(loc_i) + Z_prev(loc_i) + e_prev(loc_i));
  } 

  // True completeness
  vector<Type> fes_comp = X_comp * beta_comp.matrix();
  vector<Type> comp_est(num_locs * num_years);

  // Duration: estimated as a function of completeness
  vector<Type> duration(num_locs * num_years);

  // Incidence = prevalence / duration
  vector<Type> incidence_est(num_locs * num_years);

  // Populate estimates for completeness, duration, and incidence
  for(int year_i = 0; year_i < num_years; year_i++){
    for(int loc_i = 0; loc_i < num_locs; loc_i++){
      // Completeness
      comp_est(year_i * num_locs + loc_i) = invlogit(
        fes_comp(year_i * num_locs + loc_i) + 
        Z1_comp(loc_i) + 
        (beta_comp_t + Z2_comp(loc_i)) * zero_centered_years(year_i) +
        e_comp(year_i * num_locs + loc_i) 
      );
      // Duration
      duration(year_i * num_locs + loc_i) = (duration_intercept - duration_slope * comp_est(year_i * num_locs + loc_i));
      // Incidence
      incidence_est(year_i * num_locs + loc_i) = prev_est(loc_i) / duration(year_i * num_locs + loc_i);
    }
  }


  // EVALUATE PROBABILITY OF DATA GIVEN ESTIMATES --------------------------------------->

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
  //   Y^{notif}_i = Poisson(Incidence_{s(i), t(i)} * N^{notif}_i * Completeness_{s(i),t(i)})
  vector<int> notifs_long_index = idx_year_notif * num_locs + idx_loc_notif;
  for(int nn_i = 0; nn_i < num_notifs; nn_i++){
    if(notif_n(nn_i) > 0.){
      jnll -= dpois(
        notif_y(nn_i),
        incidence_est(notifs_long_index(nn_i)) * 
          notif_n(nn_i) *
          comp_est(notifs_long_index(nn_i)),
        true
      );
    }
  }


  // RETURN JNLL ------------------------------------------------------------------------>

  return jnll;

} // END objective function
