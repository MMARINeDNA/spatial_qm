//Quant MB model without qPCR

data {
  // ===========================================
  // DIMENSIONS
  // ===========================================
  int<lower=1> N_species;
  int<lower=1> N_locations;
  int<lower=0> N_edna_obs;
  // int<lower=0> N_qpcr_reps;

  // ===========================================
  // REFERENCE SPECIES
  // ===========================================
  int<lower=1, upper=N_species> reference_species;

  // ===========================================
  // eDNA METABARCODING DATA
  // ===========================================
  array[N_edna_obs, N_species] int<lower=0> edna_obs_data;
  array[N_edna_obs] int<lower=1, upper=N_locations> edna_obs_location;

  // ===========================================
  // qPCR REPLICATE-LEVEL DATA
  // ===========================================
  // array[N_qpcr_reps] int<lower=1, upper=N_locations> qpcr_rep_location;
  // array[N_qpcr_reps] int<lower=0, upper=1>            qpcr_rep_detected;
  // vector[N_qpcr_reps]                                  qpcr_rep_log_conc;  // sentinel 0.0 when detected==0

  // ===========================================
  // KNOWN AMPLIFICATION EFFICIENCIES
  // ===========================================
  vector[N_species] alpha_known;

  // ===========================================
  // PCR PARAMETERS
  // ===========================================
  real<lower=1> N_PCR;

  // ===========================================
  // LOCATION-LEVEL AVAILABILITY
  // ===========================================
  array[N_locations] int<lower=0, upper=1> has_edna;
  // array[N_locations] int<lower=0, upper=1> has_qpcr;

  // ===========================================
  // PRIORS — concentration
  // ===========================================
  real log_conc_prior_mean;
  real<lower=0> log_conc_prior_sd;

  // ===========================================
  // PRIORS — qPCR hurdle
  // ===========================================
  // real<lower=0> sigma_qpcr_prior_sd;
  // real<lower=0> qpcr_detection_intercept_prior_sd;
  // real          qpcr_detection_slope_prior_mean;
  // real<lower=0> qpcr_detection_slope_prior_sd;

  // ===========================================
  // PRIORS — eDNA dispersion (Beta-Binomial)
  // ===========================================
  real<lower=0> beta0_prior_sd;
  real<lower=0> gamma0_prior_mean;
  real<lower=0> gamma0_prior_sd;
  real          log_gamma1_prior_mean;
  real<lower=0> log_gamma1_prior_sd;
}

transformed data {
  array[N_edna_obs] int total_reads_obs;
  array[N_locations] int n_edna_per_loc;
  array[N_locations] int n_qpcr_per_loc;
  array[N_locations] int n_detected_per_loc;

  // Precompute alpha × N_PCR and amplification deviation from mean
  vector[N_species] alpha_pcr = N_PCR * alpha_known;
  real mean_alpha = mean(alpha_known);
  vector[N_species] kappa;
  for (sp in 1:N_species)
    kappa[sp] = alpha_known[sp] - mean_alpha;

  for (loc in 1:N_locations) {
    n_edna_per_loc[loc]     = 0;
    n_qpcr_per_loc[loc]     = 0;
    n_detected_per_loc[loc] = 0;
  }
  for (obs in 1:N_edna_obs) {
    total_reads_obs[obs] = sum(edna_obs_data[obs,]);
    n_edna_per_loc[edna_obs_location[obs]] += 1;
  }
  // for (rep in 1:N_qpcr_reps) {
  //   n_qpcr_per_loc[qpcr_rep_location[rep]]     += 1;
  //   n_detected_per_loc[qpcr_rep_location[rep]] += qpcr_rep_detected[rep];
  // }
}

parameters {
  // ===========================================
  // LOG CONCENTRATION (directly estimated)
  // ===========================================
  matrix[N_locations, N_species] log_concentration;

  // ===========================================
  // qPCR OBSERVATION ERROR
  // ===========================================
  // real log_sigma_qpcr;

  // ===========================================
  // qPCR HURDLE PARAMETERS
  // ===========================================
  // real qpcr_detection_intercept;
  // real qpcr_detection_slope;

  // ===========================================
  // eDNA DISPERSION PARAMETERS (Beta-Binomial)
  // beta1<upper=0>: larger alpha → slightly less overdispersed
  // gamma0<lower=0>: extra overdispersion at low concentration
  // ===========================================
  real          beta0;
  real<upper=0> beta1;
  real<lower=0> gamma0;
  real          log_gamma1;
}

transformed parameters {
  // real<lower=0> sigma_qpcr = exp(log_sigma_qpcr);

  // lambda = exp(log_concentration); lambda_K = ZTP mean
  matrix[N_locations, N_species] lambda;
  matrix[N_locations, N_species] lambda_K;
  matrix[N_locations, N_species] log_lambda_K;
  vector[N_locations]            total_lambda;
  vector[N_locations]            log_total_lambda;
  matrix[N_locations, N_species] log_phi_mat;

  {
    real gamma1 = exp(log_gamma1);
    for (loc in 1:N_locations) {
      for (sp in 1:N_species) {
        lambda[loc, sp]   = exp(log_concentration[loc, sp]);
        log_lambda_K[loc, sp] = log_concentration[loc, sp] - log1m_exp(-log_concentration[loc, sp]);
      }
      log_total_lambda[loc] = log_sum_exp(log_concentration[loc, ]);
      for (sp in 1:N_species)
        log_phi_mat[loc, sp] = beta0 + log_total_lambda[loc]
                               + gamma0 * exp(-gamma1 * log_lambda_K[loc, sp])
                               + beta1 * kappa[sp];
    }
  }

  // Metabarcoding proportions per PCR replicate via ZTP correction
  matrix<lower=0, upper=1>[N_edna_obs, N_species] mu_metabar_obs;
  for (obs in 1:N_edna_obs) {
    int loc = edna_obs_location[obs];
    real log_lK_ref = log_lambda_K[loc, reference_species];
    vector[N_species] log_nu;
    for (sp in 1:N_species)
      log_nu[sp] = log_lambda_K[loc, sp] - log_lK_ref + alpha_pcr[sp];
    mu_metabar_obs[obs,] = to_row_vector(softmax(log_nu));
  }
}

model {
  // ============================================
  // LOG CONCENTRATION PRIOR (wide, independent)
  // ============================================
  to_vector(log_concentration) ~ normal(log_conc_prior_mean, log_conc_prior_sd);

  // ============================================
  // qPCR OBSERVATION ERROR PRIOR
  // Half-normal on sigma_qpcr; Jacobian for log parameterization
  // ============================================
  // sigma_qpcr ~ normal(0, sigma_qpcr_prior_sd);
  // target += log_sigma_qpcr;   // Jacobian: log|d(exp(x))/dx| = x

  // ============================================
  // qPCR HURDLE PRIORS
  // ============================================
  // qpcr_detection_intercept ~ normal(0, qpcr_detection_intercept_prior_sd);
  // qpcr_detection_slope     ~ normal(qpcr_detection_slope_prior_mean,
  //                                   qpcr_detection_slope_prior_sd);

  // ============================================
  // eDNA DISPERSION PRIORS
  // beta1 prior centered far below 0 (weakly informative of negative effect)
  // ============================================
  beta0      ~ normal(0, beta0_prior_sd);
  beta1      ~ normal(-10, 5);
  gamma0     ~ normal(gamma0_prior_mean, gamma0_prior_sd);
  log_gamma1 ~ normal(log_gamma1_prior_mean, log_gamma1_prior_sd);

  // ============================================
  // eDNA LIKELIHOOD — Poisson ZI + Beta-Binomial
  // Applied independently per species (factored approximation)
  // ============================================
  for (obs in 1:N_edna_obs) {
    if (total_reads_obs[obs] > 0) {
      int loc = edna_obs_location[obs];
      for (sp in 1:N_species) {
        real lam       = lambda[loc, sp];
        real log_pzero = -lam;
        real log_ppos  = log1m_exp(-lam);
        real phi_sp    = exp(log_phi_mat[loc, sp]);
        real bb_a      = mu_metabar_obs[obs, sp] * phi_sp;
        real bb_b      = (1.0 - mu_metabar_obs[obs, sp]) * phi_sp;

        if (edna_obs_data[obs, sp] > 0) {// If there are 1 or more reads observed
          target += log_ppos
                    + beta_binomial_lpmf(edna_obs_data[obs, sp]
                                         | total_reads_obs[obs], bb_a, bb_b);
        } else { // If there are exactly zero reads observed
          target += log_sum_exp(
            log_pzero,
            log_ppos + beta_binomial_lpmf(0 | total_reads_obs[obs], bb_a, bb_b)
          );
        }
      }
    }
  }

  // ============================================
  // qPCR HURDLE LIKELIHOOD
  // ============================================
  // for (rep in 1:N_qpcr_reps) {
  //   int loc = qpcr_rep_location[rep];
  //   real log_conc_ref = log_concentration[loc, reference_species];
  //   real logit_p = qpcr_detection_intercept + qpcr_detection_slope * log_conc_ref;
  //   target += bernoulli_logit_lpmf(qpcr_rep_detected[rep] | logit_p);
  //   if (qpcr_rep_detected[rep] == 1)
  //     target += normal_lpdf(qpcr_rep_log_conc[rep] | log_conc_ref, sigma_qpcr);
  // }
}

generated quantities {
  // ============================================
  // CONCENTRATION AND PROPORTIONS
  // ============================================
  matrix<lower=0>[N_locations, N_species] concentration;
  vector<lower=0>[N_locations] total_concentration;
  matrix<lower=0, upper=1>[N_locations, N_species] true_proportions;

  for (loc in 1:N_locations) {
    concentration[loc,] = exp(log_concentration[loc,]);
    total_concentration[loc] = sum(concentration[loc,]);
    for (sp in 1:N_species)
      true_proportions[loc, sp] = concentration[loc, sp] / total_concentration[loc];
  }

  // ============================================
  // LOCATION-LEVEL METABARCODING PROPORTIONS
  // ============================================
  matrix<lower=0, upper=1>[N_locations, N_species] mu_metabar_loc;
  for (loc in 1:N_locations)
    for (i in 1:N_species) mu_metabar_loc[loc, i] = 0;

  for (obs in 1:N_edna_obs) {
    int loc = edna_obs_location[obs];
    for (i in 1:N_species)
      mu_metabar_loc[loc, i] += mu_metabar_obs[obs, i];
  }
  for (loc in 1:N_locations)
    if (n_edna_per_loc[loc] > 0)
      for (i in 1:N_species)
        mu_metabar_loc[loc, i] /= n_edna_per_loc[loc];

  // ============================================
  // DERIVED: PSI AND SPECIES OFFSETS
  // ============================================
  real psi = mean(to_vector(log_concentration));

  vector[N_species] species_mean_log_conc;
  vector[N_species] species_offset;
  for (sp in 1:N_species) {
    species_mean_log_conc[sp] = mean(log_concentration[, sp]);
    species_offset[sp] = species_mean_log_conc[sp] - psi;
  }

  // ============================================
  // DISPERSION SUMMARY
  // ============================================
  real gamma1   = exp(log_gamma1);
  real mean_phi = mean(exp(log_phi_mat));

  // ============================================
  // LOCATION-LEVEL EXPECTED qPCR
  // ============================================
  vector[N_locations] expected_qpcr_loc;
  for (loc in 1:N_locations)
    expected_qpcr_loc[loc] = log_concentration[loc, reference_species];

  // ============================================
  // qPCR HURDLE: POSTERIOR DETECTION PROBABILITY
  // ============================================
  // vector[N_locations] p_detect_loc;
  // for (loc in 1:N_locations)
  //   p_detect_loc[loc] = inv_logit(
  //     qpcr_detection_intercept +
  //     qpcr_detection_slope * log_concentration[loc, reference_species]
  //   );

  // ============================================
  // SITE TOTALS
  // ============================================
  vector[N_locations] total_conc_loc;
  vector[N_locations] log_total;
  for (loc in 1:N_locations) {
    total_conc_loc[loc] = sum(concentration[loc,]);
    log_total[loc]      = log(total_concentration[loc]);
  }

  // ============================================
  // SPECIES-LEVEL SUMMARIES
  // ============================================
  vector[N_species] mean_concentration_species;
  vector[N_species] mean_true_proportion_species;
  vector[N_species] mean_metabar_proportion_species;

  int n_edna_locations = 0;
  int n_qpcr_locations = 0;
  for (loc in 1:N_locations) {
    if (has_edna[loc] == 1) n_edna_locations += 1;
    // if (has_qpcr[loc] == 1) n_qpcr_locations += 1;
  }

  for (i in 1:N_species) {
    mean_concentration_species[i]   = mean(concentration[, i]);
    mean_true_proportion_species[i] = mean(true_proportions[, i]);
    real sum_metabar = 0;
    for (loc in 1:N_locations)
      if (has_edna[loc] == 1) sum_metabar += mu_metabar_loc[loc, i];
    mean_metabar_proportion_species[i] =
      n_edna_locations > 0 ? sum_metabar / n_edna_locations : 0;
  }

  // ============================================
  // AMPLIFICATION BIAS SUMMARY
  // ============================================
  vector[N_species] amplification_bias_ratio;
  for (i in 1:N_species) {
    if (mean_true_proportion_species[i] > 1e-10)
      amplification_bias_ratio[i] =
        mean_metabar_proportion_species[i] / mean_true_proportion_species[i];
    else
      amplification_bias_ratio[i] = 1.0;
  }

  // ============================================
  // eDNA PRESENCE PROBABILITY
  // ============================================
  matrix<lower=0, upper=1>[N_edna_obs, N_species] p_present_pred;
  for (obs in 1:N_edna_obs)
    for (sp in 1:N_species)
      p_present_pred[obs, sp] = 1.0 - exp(-lambda[edna_obs_location[obs], sp]);

  // ============================================
  // POSTERIOR PREDICTIVE — eDNA (ZI-masked multinomial draw)
  // ============================================
  array[N_edna_obs, N_species] int edna_pred;
  for (obs in 1:N_edna_obs) {
    if (total_reads_obs[obs] > 0) {
      int loc = edna_obs_location[obs];
      vector[N_species] pi_masked;
      for (sp in 1:N_species) {
        int sp_present = bernoulli_rng(1.0 - exp(-lambda[loc, sp]));
        pi_masked[sp] = sp_present > 0 ? mu_metabar_obs[obs, sp] : 0.0;
      }
      real sum_pi = sum(pi_masked);
      if (sum_pi > 0) {
        pi_masked /= sum_pi;
        edna_pred[obs,] = multinomial_rng(pi_masked, total_reads_obs[obs]);
      } else {
        for (sp in 1:N_species) edna_pred[obs, sp] = 0;
      }
    } else {
      for (sp in 1:N_species) edna_pred[obs, sp] = 0;
    }
  }

  // ============================================
  // POSTERIOR PREDICTIVE — qPCR (hurdle)
  // ============================================
  // array[N_qpcr_reps] int qpcr_pred_detected;
  // vector[N_qpcr_reps] qpcr_pred_log_conc;
  // 
  // for (rep in 1:N_qpcr_reps) {
  //   int loc = qpcr_rep_location[rep];
  //   real log_conc_ref = log_concentration[loc, reference_species];
  //   real logit_p = qpcr_detection_intercept + qpcr_detection_slope * log_conc_ref;
  //   qpcr_pred_detected[rep] = bernoulli_logit_rng(logit_p);
  //   if (qpcr_pred_detected[rep] == 1)
  //     qpcr_pred_log_conc[rep] = normal_rng(log_conc_ref, sigma_qpcr);
  //   else
  //     qpcr_pred_log_conc[rep] = 0.0;
  // }

  // ============================================
  // LOCATION-LEVEL POSTERIOR PREDICTIVE — eDNA
  // ============================================
  array[N_locations, N_species] int edna_pred_loc;
  for (loc in 1:N_locations)
    for (i in 1:N_species) edna_pred_loc[loc, i] = 0;
  for (obs in 1:N_edna_obs) {
    int loc = edna_obs_location[obs];
    for (i in 1:N_species) edna_pred_loc[loc, i] += edna_pred[obs, i];
  }

  // ============================================
  // LOG-LIKELIHOOD — eDNA (ZI + BB)
  // ============================================
  vector[N_edna_obs] log_lik_edna;
  for (obs in 1:N_edna_obs) {
    if (total_reads_obs[obs] > 0) {
      int loc = edna_obs_location[obs];
      real ll = 0;
      for (sp in 1:N_species) {
        real lam       = lambda[loc, sp];
        real log_pzero = -lam;
        real log_ppos  = log1m_exp(-lam);
        real phi_sp    = exp(log_phi_mat[loc, sp]);
        real bb_a      = mu_metabar_obs[obs, sp] * phi_sp;
        real bb_b      = (1.0 - mu_metabar_obs[obs, sp]) * phi_sp;
        if (edna_obs_data[obs, sp] > 0) {
          ll += log_ppos
                + beta_binomial_lpmf(edna_obs_data[obs, sp]
                                      | total_reads_obs[obs], bb_a, bb_b);
        } else {
          ll += log_sum_exp(
            log_pzero,
            log_ppos + beta_binomial_lpmf(0 | total_reads_obs[obs], bb_a, bb_b)
          );
        }
      }
      log_lik_edna[obs] = ll;
    } else {
      log_lik_edna[obs] = 0;
    }
  }

  // ============================================
  // LOG-LIKELIHOOD — qPCR (hurdle, per replicate)
  // ============================================
  // vector[N_qpcr_reps] log_lik_qpcr;
  // for (rep in 1:N_qpcr_reps) {
  //   int loc = qpcr_rep_location[rep];
  //   real log_conc_ref = log_concentration[loc, reference_species];
  //   real logit_p = qpcr_detection_intercept + qpcr_detection_slope * log_conc_ref;
  //   real ll = bernoulli_logit_lpmf(qpcr_rep_detected[rep] | logit_p);
  //   if (qpcr_rep_detected[rep] == 1)
  //     ll += normal_lpdf(qpcr_rep_log_conc[rep] | log_conc_ref, sigma_qpcr);
  //   log_lik_qpcr[rep] = ll;
  // }

  // ============================================
  // LOCATION-LEVEL LOG-LIKELIHOOD
  // ============================================
  vector[N_locations] log_lik_edna_loc;
  vector[N_locations] log_lik_qpcr_loc;
  for (loc in 1:N_locations) {
    log_lik_edna_loc[loc] = 0;
    log_lik_qpcr_loc[loc] = 0;
  }
  for (obs in 1:N_edna_obs)
    log_lik_edna_loc[edna_obs_location[obs]] += log_lik_edna[obs];
  // for (rep in 1:N_qpcr_reps)
  //   log_lik_qpcr_loc[qpcr_rep_location[rep]] += log_lik_qpcr[rep];
}
