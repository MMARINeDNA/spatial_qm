/**
 * Joint Model for Metabarcoding and qPCR Data
 *
 * VERSION 22: ISOTROPIC HSGP (SINGLE rho PER SPECIES)
 *
 * Changes from v21:
 *
 * 1. Replace anisotropic length scales (rho_x[sp], rho_y[sp]) with a single
 *    isotropic length scale rho[sp].
 *    v21: log_rho_x[N_species] + log_rho_y[N_species]  (2 * N_species params)
 *    v22: rho[N_species] (truncated-normal prior, km scale) (N_species params)
 *
 *    Motivation: simulations with rho_x = rho_y (isotropic truth) gave the
 *    anisotropic model two equivalent solutions (large rho_x + small rho_y
 *    and small rho_x + large rho_y), creating a bimodal or diffuse joint
 *    posterior even when concentrations are estimated well. Removing this
 *    symmetry is the first diagnostic step before attributing mis-estimation
 *    to the rho–gp_alpha–nugget competition.
 *
 * 2. Replace hsgp_2d_aniso_spectral_weights with hsgp_2d_iso_spectral_weights.
 *    Isotropic 2D SE spectral density:
 *      S = alpha^2 * 2*pi * rho^2 * exp(-0.5 * rho^2 * (lambda_x[j1] + lambda_y[j2]))
 *    (equivalent to anisotropic formula with rho_x = rho_y = rho)
 *
 * Everything else is identical to v21:
 *   - Poisson ZI + Beta-Binomial eDNA model
 *   - Per-replicate hurdle qPCR model
 *   - Decoupled sigma_nugget (own hierarchy, independent of gp_alpha)
 *   - No nugget mean-centering
 *   - Independent (non-hierarchical) length scales per species
 *   - Dispersion model (beta0/beta1/gamma0/gamma1)
 */

functions {
  vector hsgp_eigenvalues(int m, real L) {
    vector[m] lambda;
    for (j in 1:m)
      lambda[j] = square(j * pi() / (2 * L));
    return lambda;
  }

  matrix hsgp_basis(vector x, int m, real L) {
    int N = rows(x);
    matrix[N, m] PHI;
    real sqrt_L_inv = 1.0 / sqrt(L);
    for (j in 1:m) {
      real freq = j * pi() / (2 * L);
      for (i in 1:N)
        PHI[i, j] = sqrt_L_inv * sin(freq * (x[i] + L));
    }
    return PHI;
  }

  /**
   * Isotropic 2D squared-exponential spectral density weights.
   *
   * S[idx] = alpha^2 * 2*pi * rho^2 *
   *          exp(-0.5 * rho^2 * (lambda_x[j1] + lambda_y[j2]))
   *
   * This is the special case of the anisotropic formula with rho_x = rho_y = rho.
   */
  vector hsgp_2d_iso_spectral_weights(int m_x, int m_y,
                                      real L_x, real L_y,
                                      real alpha, real rho) {
    int M = m_x * m_y;
    vector[M] weights;
    vector[m_x] lambda_x = hsgp_eigenvalues(m_x, L_x);
    vector[m_y] lambda_y = hsgp_eigenvalues(m_y, L_y);
    int idx = 1;
    for (j1 in 1:m_x) {
      for (j2 in 1:m_y) {
        weights[idx] = square(alpha) * 2 * pi() * square(rho) *
                       exp(-0.5 * square(rho) * (lambda_x[j1] + lambda_y[j2]));
        idx += 1;
      }
    }
    return sqrt(weights);
  }

  matrix hsgp_2d_basis(matrix PHI_x, matrix PHI_y) {
    int N = rows(PHI_x);
    int m_x = cols(PHI_x);
    int m_y = cols(PHI_y);
    int M = m_x * m_y;
    matrix[N, M] PHI;
    int idx = 1;
    for (j1 in 1:m_x) {
      for (j2 in 1:m_y) {
        PHI[, idx] = PHI_x[, j1] .* PHI_y[, j2];
        idx += 1;
      }
    }
    return PHI;
  }
}

data {
  int<lower=1> N_species;
  int<lower=1> N_locations;
  int<lower=0> N_edna_obs;
  int<lower=0> N_qpcr_reps;

  int<lower=1, upper=N_species> reference_species;

  int<lower=1> m_x;
  int<lower=1> m_y;
  real<lower=1> c_x;
  real<lower=1> c_y;

  vector[N_locations] x_utm;
  vector[N_locations] y_utm;

  array[N_edna_obs, N_species] int<lower=0> edna_obs_data;
  array[N_edna_obs] int<lower=1, upper=N_locations> edna_obs_location;

  array[N_qpcr_reps] int<lower=1, upper=N_locations> qpcr_rep_location;
  array[N_qpcr_reps] int<lower=0, upper=1>            qpcr_rep_detected;
  vector[N_qpcr_reps]                                  qpcr_rep_log_conc;

  vector[N_species] alpha_known;
  real<lower=1> N_PCR;

  array[N_locations] int<lower=0, upper=1> has_edna;
  array[N_locations] int<lower=0, upper=1> has_qpcr;

  // Concentration field priors
  real log_conc_prior_mean;
  real<lower=0> log_conc_prior_sd;

  // GP amplitude priors
  real gp_alpha_prior_mean;
  real<lower=0> gp_alpha_prior_sd;

  // Isotropic length scale prior (single rho per species; v22)
  real<lower=0> rho_prior_mean;
  real<lower=0> rho_prior_sd;

  // qPCR observation error prior
  real<lower=0> sigma_qpcr_prior_sd;

  // qPCR hurdle priors
  real<lower=0> qpcr_detection_intercept_prior_sd;
  real          qpcr_detection_slope_prior_mean;
  real<lower=0> qpcr_detection_slope_prior_sd;

  // Nugget SD priors (hierarchical across species, independent of gp_alpha)
  real          sigma_nugget_prior_mean;   // prior center for mu_log_sigma_nugget (log scale)
  real<lower=0> sigma_nugget_prior_sd;     // prior SD for mu_log_sigma_nugget
  real<lower=0> sigma_nugget_sigma_rate;   // exponential rate for sigma_sigma_nugget

  // eDNA dispersion priors (unchanged)
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

  int M = m_x * m_y;
  real L_x = c_x * max(fabs(x_utm)) + 1.0;
  real L_y = c_y * max(fabs(y_utm)) + 1.0;
  matrix[N_locations, m_x] PHI_x = hsgp_basis(x_utm, m_x, L_x);
  matrix[N_locations, m_y] PHI_y = hsgp_basis(y_utm, m_y, L_y);
  matrix[N_locations, M]   PHI   = hsgp_2d_basis(PHI_x, PHI_y);

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
  for (rep in 1:N_qpcr_reps) {
    n_qpcr_per_loc[qpcr_rep_location[rep]]     += 1;
    n_detected_per_loc[qpcr_rep_location[rep]] += qpcr_rep_detected[rep];
  }
}

parameters {
  vector[N_species] mu_species;

  real log_sigma_qpcr;

  real qpcr_detection_intercept;
  real qpcr_detection_slope;

  real              beta0;
  real<upper=0>     beta1;
  real<lower=0>     gamma0;
  real              log_gamma1;

  // Isotropic length scale per species (v22: single rho, replaces rho_x + rho_y)
  // Truncated-normal prior directly on rho (km); lower=0 enforced by Stan
  vector<lower=0>[N_species] rho;

  vector[N_species] log_gp_alpha;
  matrix[M, N_species] z_basis;

  // Nugget SD: hierarchical log-normal, INDEPENDENT of gp_alpha (from v21)
  real              mu_log_sigma_nugget;
  real<lower=0>     sigma_sigma_nugget;
  vector[N_species] log_sigma_nugget_raw;
  matrix[N_locations, N_species] nugget_raw;
}

transformed parameters {
  // rho: truncated-normal primary parameter; gp_alpha: lognormal via log_gp_alpha
  vector<lower=0>[N_species] gp_alpha = exp(log_gp_alpha);

  real<lower=0> sigma_qpcr = exp(log_sigma_qpcr);

  // Nugget SD: own hierarchy, decoupled from gp_alpha (from v21)
  vector[N_species] log_sigma_nugget =
      mu_log_sigma_nugget + sigma_sigma_nugget * log_sigma_nugget_raw;
  vector<lower=0>[N_species] sigma_nugget = exp(log_sigma_nugget);

  matrix[N_locations, N_species] spatial_effect;
  for (sp in 1:N_species) {
    vector[M] spd_weights_sp = hsgp_2d_iso_spectral_weights(
      m_x, m_y, L_x, L_y, gp_alpha[sp], rho[sp]);
    vector[N_locations] raw_effect = PHI * (spd_weights_sp .* z_basis[, sp]);
    spatial_effect[, sp] = raw_effect - mean(raw_effect);
  }

  // Nugget: iid N(0, sigma_nugget^2) — no mean-centering (from v21)
  matrix[N_locations, N_species] nugget;
  for (sp in 1:N_species)
    nugget[, sp] = sigma_nugget[sp] * nugget_raw[, sp];

  matrix[N_locations, N_species] log_concentration;
  for (loc in 1:N_locations)
    for (sp in 1:N_species)
      log_concentration[loc, sp] = mu_species[sp]
                                   + spatial_effect[loc, sp]
                                   + nugget[loc, sp];

  matrix[N_locations, N_species] lambda;
  matrix[N_locations, N_species] lambda_K;
  vector[N_locations]            total_lambda;
  matrix[N_locations, N_species] log_phi_mat;

  {
    real gamma1 = exp(log_gamma1);
    for (loc in 1:N_locations) {
      for (sp in 1:N_species) {
        lambda[loc, sp]   = exp(log_concentration[loc, sp]);
        lambda_K[loc, sp] = lambda[loc, sp] / (-expm1(-lambda[loc, sp]));
      }
      total_lambda[loc] = sum(lambda[loc,]);
      for (sp in 1:N_species)
        log_phi_mat[loc, sp] = beta0 + log(total_lambda[loc])
                               + gamma0 * exp(-gamma1 * log(lambda_K[loc, sp]))
                               + beta1 * kappa[sp];
    }
  }

  matrix<lower=0, upper=1>[N_edna_obs, N_species] mu_metabar_obs;
  for (obs in 1:N_edna_obs) {
    int loc = edna_obs_location[obs];
    real log_lK_ref = log(lambda_K[loc, reference_species]);
    vector[N_species] log_nu;
    for (sp in 1:N_species)
      log_nu[sp] = log(lambda_K[loc, sp]) - log_lK_ref + alpha_pcr[sp];
    mu_metabar_obs[obs,] = to_row_vector(softmax(log_nu));
  }
}

model {
  mu_species ~ normal(log_conc_prior_mean, log_conc_prior_sd);

  sigma_qpcr ~ normal(0, sigma_qpcr_prior_sd);
  target += log_sigma_qpcr;

  qpcr_detection_intercept ~ normal(0, qpcr_detection_intercept_prior_sd);
  qpcr_detection_slope     ~ normal(qpcr_detection_slope_prior_mean,
                                    qpcr_detection_slope_prior_sd);

  beta0      ~ normal(0, beta0_prior_sd);
  beta1      ~ normal(-10, 5);
  gamma0     ~ normal(gamma0_prior_mean, gamma0_prior_sd);
  log_gamma1 ~ normal(log_gamma1_prior_mean, log_gamma1_prior_sd);

  // Truncated-normal prior on rho (km): N+(rho_prior_mean, rho_prior_sd^2)
  // Stan handles truncation automatically given lower=0 constraint on rho
  rho ~ normal(rho_prior_mean, rho_prior_sd);

  log_gp_alpha     ~ normal(gp_alpha_prior_mean, gp_alpha_prior_sd);
  to_vector(z_basis) ~ std_normal();

  // Nugget SD hierarchy — decoupled from gp_alpha (from v21)
  mu_log_sigma_nugget  ~ normal(sigma_nugget_prior_mean, sigma_nugget_prior_sd);
  sigma_sigma_nugget   ~ exponential(sigma_nugget_sigma_rate);
  log_sigma_nugget_raw ~ std_normal();
  to_vector(nugget_raw) ~ std_normal();

  // eDNA likelihood — Poisson ZI + Beta-Binomial
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

        if (edna_obs_data[obs, sp] > 0) {
          target += log_ppos
                    + beta_binomial_lpmf(edna_obs_data[obs, sp]
                                         | total_reads_obs[obs], bb_a, bb_b);
        } else {
          target += log_sum_exp(
            log_pzero,
            log_ppos + beta_binomial_lpmf(0 | total_reads_obs[obs], bb_a, bb_b)
          );
        }
      }
    }
  }

  // qPCR hurdle likelihood
  for (rep in 1:N_qpcr_reps) {
    int loc = qpcr_rep_location[rep];
    real log_conc_ref = log_concentration[loc, reference_species];
    real logit_p = qpcr_detection_intercept + qpcr_detection_slope * log_conc_ref;
    target += bernoulli_logit_lpmf(qpcr_rep_detected[rep] | logit_p);
    if (qpcr_rep_detected[rep] == 1)
      target += normal_lpdf(qpcr_rep_log_conc[rep] | log_conc_ref, sigma_qpcr);
  }
}

generated quantities {
  matrix<lower=0>[N_locations, N_species] concentration;
  vector<lower=0>[N_locations] total_concentration;
  matrix<lower=0, upper=1>[N_locations, N_species] true_proportions;

  for (loc in 1:N_locations) {
    concentration[loc,] = exp(log_concentration[loc,]);
    total_concentration[loc] = sum(concentration[loc,]);
    for (sp in 1:N_species)
      true_proportions[loc, sp] = concentration[loc, sp] / total_concentration[loc];
  }

  matrix<lower=0, upper=1>[N_locations, N_species] mu_metabar_loc;
  for (loc in 1:N_locations)
    for (i in 1:N_species)
      mu_metabar_loc[loc, i] = 0;

  for (obs in 1:N_edna_obs) {
    int loc = edna_obs_location[obs];
    for (i in 1:N_species)
      mu_metabar_loc[loc, i] += mu_metabar_obs[obs, i];
  }
  for (loc in 1:N_locations)
    if (n_edna_per_loc[loc] > 0)
      for (i in 1:N_species)
        mu_metabar_loc[loc, i] /= n_edna_per_loc[loc];

  real psi = mean(mu_species);
  vector[N_species] species_offset = mu_species - psi;

  real mean_gp_alpha = mean(gp_alpha);

  // Cross-species mean of isotropic length scale (v22)
  real rho_mean = mean(rho);

  // Nugget SD population summaries (lognormal moments of hierarchical prior)
  real sigma_nugget_pop_mean   = exp(mu_log_sigma_nugget + 0.5 * square(sigma_sigma_nugget));
  real sigma_nugget_pop_median = exp(mu_log_sigma_nugget);
  real sigma_nugget_mean       = mean(sigma_nugget);

  vector[N_species] prop_var_gp;
  vector[N_species] prop_var_nugget;
  for (sp in 1:N_species) {
    real var_gp_sp     = variance(spatial_effect[, sp]);
    real var_nugget_sp = square(sigma_nugget[sp]);
    real total_var     = var_gp_sp + var_nugget_sp;
    if (total_var > 1e-10) {
      prop_var_gp[sp]     = var_gp_sp / total_var;
      prop_var_nugget[sp] = var_nugget_sp / total_var;
    } else {
      prop_var_gp[sp]     = 0.5;
      prop_var_nugget[sp] = 0.5;
    }
  }

  vector[N_locations] expected_qpcr_loc;
  for (loc in 1:N_locations)
    expected_qpcr_loc[loc] = log_concentration[loc, reference_species];

  vector[N_locations] total_conc_loc;
  vector[N_locations] log_total;
  for (loc in 1:N_locations) {
    total_conc_loc[loc] = sum(concentration[loc,]);
    log_total[loc]      = log(total_concentration[loc]);
  }

  vector[N_species] mean_concentration_species;
  vector[N_species] mean_true_proportion_species;
  vector[N_species] mean_metabar_proportion_species;

  int n_edna_locations = 0;
  int n_qpcr_locations = 0;
  for (loc in 1:N_locations) {
    if (has_edna[loc] == 1) n_edna_locations += 1;
    if (has_qpcr[loc] == 1) n_qpcr_locations += 1;
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

  vector[N_species] amplification_bias_ratio;
  for (i in 1:N_species) {
    if (mean_true_proportion_species[i] > 1e-10)
      amplification_bias_ratio[i] =
        mean_metabar_proportion_species[i] / mean_true_proportion_species[i];
    else
      amplification_bias_ratio[i] = 1.0;
  }

  vector[N_species] spatial_effect_range_species;
  vector[N_species] spatial_effect_sd_species;
  for (sp in 1:N_species) {
    spatial_effect_range_species[sp] = max(spatial_effect[, sp]) - min(spatial_effect[, sp]);
    spatial_effect_sd_species[sp]    = sd(spatial_effect[, sp]);
  }
  real mean_spatial_effect_sd = mean(spatial_effect_sd_species);

  vector[N_species] mean_spatial_effect_by_species;
  for (sp in 1:N_species)
    mean_spatial_effect_by_species[sp] = mean(spatial_effect[, sp]);

  vector[N_locations] p_detect_loc;
  for (loc in 1:N_locations)
    p_detect_loc[loc] = inv_logit(
      qpcr_detection_intercept +
      qpcr_detection_slope * log_concentration[loc, reference_species]);

  matrix<lower=0, upper=1>[N_edna_obs, N_species] p_present_pred;
  for (obs in 1:N_edna_obs)
    for (sp in 1:N_species)
      p_present_pred[obs, sp] = 1.0 - exp(-lambda[edna_obs_location[obs], sp]);

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

  array[N_qpcr_reps] int qpcr_pred_detected;
  vector[N_qpcr_reps] qpcr_pred_log_conc;
  for (rep in 1:N_qpcr_reps) {
    int loc = qpcr_rep_location[rep];
    real log_conc_ref = log_concentration[loc, reference_species];
    real logit_p = qpcr_detection_intercept + qpcr_detection_slope * log_conc_ref;
    qpcr_pred_detected[rep] = bernoulli_logit_rng(logit_p);
    if (qpcr_pred_detected[rep] == 1)
      qpcr_pred_log_conc[rep] = normal_rng(log_conc_ref, sigma_qpcr);
    else
      qpcr_pred_log_conc[rep] = 0.0;
  }

  array[N_locations, N_species] int edna_pred_loc;
  for (loc in 1:N_locations)
    for (i in 1:N_species) edna_pred_loc[loc, i] = 0;
  for (obs in 1:N_edna_obs) {
    int loc = edna_obs_location[obs];
    for (i in 1:N_species) edna_pred_loc[loc, i] += edna_pred[obs, i];
  }

  // Log-likelihood — eDNA
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

  // Log-likelihood — qPCR
  vector[N_qpcr_reps] log_lik_qpcr;
  for (rep in 1:N_qpcr_reps) {
    int loc = qpcr_rep_location[rep];
    real log_conc_ref = log_concentration[loc, reference_species];
    real logit_p = qpcr_detection_intercept + qpcr_detection_slope * log_conc_ref;
    real ll = bernoulli_logit_lpmf(qpcr_rep_detected[rep] | logit_p);
    if (qpcr_rep_detected[rep] == 1)
      ll += normal_lpdf(qpcr_rep_log_conc[rep] | log_conc_ref, sigma_qpcr);
    log_lik_qpcr[rep] = ll;
  }

  vector[N_locations] log_lik_edna_loc;
  vector[N_locations] log_lik_qpcr_loc;
  for (loc in 1:N_locations) {
    log_lik_edna_loc[loc] = 0;
    log_lik_qpcr_loc[loc] = 0;
  }
  for (obs in 1:N_edna_obs)
    log_lik_edna_loc[edna_obs_location[obs]] += log_lik_edna[obs];
  for (rep in 1:N_qpcr_reps)
    log_lik_qpcr_loc[qpcr_rep_location[rep]] += log_lik_qpcr[rep];

  real gamma1   = exp(log_gamma1);
  real mean_phi = mean(exp(log_phi_mat));
}
