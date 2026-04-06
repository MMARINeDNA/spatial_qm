#!/usr/bin/env Rscript
# sim_example_v27_mwe.R
#
# Minimal working example for testing the saveRDS pipeline.
# Small but viable problem: 30 locations, 4 species, K_f=2, M=12 basis functions.
# 2 chains x 400 iter (200 warmup) — should finish in ~20 minutes on HYAK.
# Sources run_qpcr_metabar_workflow_v27.R for all utilities.
#
# Usage (run from project root):
#   Rscript code/spatial_qPCR_mb/sim_example_v27_mwe.R [output_dir]

args       <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else
  "/mmfs1/gscratch/coenv/rpkelly/output/sim_example_v27_mwe"

cat("=== sim_example_v27_mwe.R ===\n")
cat(sprintf("output_dir: %s\n\n", output_dir))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

# --------------------------------------------------------------------------
# Source workflow
# --------------------------------------------------------------------------
workflow_file <- file.path("code", "spatial_qPCR_mb", "run_qpcr_metabar_workflow_v27.R")
if (!file.exists(workflow_file))
  stop("Cannot find workflow file: ", workflow_file,
       "\nRun this script from the project root directory.")

source(workflow_file)

# --------------------------------------------------------------------------
# Run — small but uses K_f=2 / 4 species to match verified working profile
# --------------------------------------------------------------------------
result <- run_workflow_v27(
  stan_file           = "code/spatial_qPCR_mb/quant_metabar_qpcr_v27.stan",
  n_locations         = 30,
  n_species           = 4,
  K_f                 = 2,
  reference_species   = 1,
  n_qpcr_per_location = 3,
  chains              = 2,
  iter                = 400,
  warmup              = 200,
  seed                = 42L,
  # HSGP basis dimensions (M = 3*4 = 12; M < 30/2 = 15)
  m_x                 = 3,
  m_y                 = 4,
  max_rho             = 300,
  rho_min             = 20,
  # True simulation parameters
  rho_factor_km       = 50,
  gp_alpha_factor_true = NULL,
  loading_true        = NULL,
  sigma_nugget_true   = 0.3,
  log_lambda_mean     = log(5),
  mu_species_true     = NULL,
  sigma_qpcr_true     = 0.3,
  qpcr_detection_intercept_true = -2.0,
  qpcr_detection_slope_true     = 1.5,
  total_reads_per_rep = 10000,
  beta0_true          = 5.0,
  beta1_true          = -2.0,
  gamma0_true         = 2.0,
  gamma1_true         = exp(1.0),
  # Priors
  rho_factor_prior_mean = 100.0,
  rho_factor_prior_sd   = 200.0,
  gp_alpha_factor       = 1.0,
  loading_prior_sd      = 1.0,
  sigma_nugget_prior_mean = log(0.1),
  sigma_nugget_prior_sd   = 1e-5,
  sigma_nugget_sigma_rate = 1.0,
  sigma_qpcr_prior_sd     = 1.0,
  qpcr_detection_intercept_prior_sd = 5.0,
  qpcr_detection_slope_prior_mean   = 2.0,
  qpcr_detection_slope_prior_sd     = 2.0,
  log_conc_prior_sd   = 5.0,
  beta0_prior_sd      = 1.0,
  gamma0_prior_mean   = 2.0,
  gamma0_prior_sd     = 1.0,
  log_gamma1_prior_mean = 1.0,
  log_gamma1_prior_sd   = 0.75,
  # Relaxed HMC tuning for speed
  adapt_delta         = 0.90,
  max_treedepth       = 10,
  # Output
  save_outputs        = TRUE,
  output_dir          = output_dir
)

# --------------------------------------------------------------------------
# Report + save
# --------------------------------------------------------------------------
cat("\n=== Done ===\n")
cat(sprintf("Fit successful: %s\n", result$fit_successful))

diag <- result$diagnostics
cat(sprintf("Divergences: %s | Max treedepth hits: %s\n",
            diag$n_divergent, diag$max_treedepth_hits))
cat(sprintf("Max Rhat: %.3f | # Rhat > 1.01: %d\n",
            diag$max_rhat, diag$n_high_rhat))

fit_obj <- result[["fit"]]
if (!is.null(fit_obj) && inherits(fit_obj, "stanfit")) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  rds_path <- file.path(output_dir, sprintf("fit_v27_mwe_%s.rds", ts))
  posterior_draws <- rstan::extract(fit_obj, permuted = FALSE)
  saveRDS(posterior_draws, rds_path)
  rds_size_mb <- file.info(rds_path)$size / 1e6
  cat(sprintf("Stan draws saved to: %s (%.1f MB)\n", rds_path, rds_size_mb))
} else {
  cat(sprintf("WARNING: result[[\"fit\"]] is class '%s' — RDS not saved\n",
              paste(class(fit_obj), collapse = "/")))
}

cat(sprintf("All outputs written to: %s\n", output_dir))
