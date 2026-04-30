#!/usr/bin/env Rscript
# sim_example_v22_mwe.R
#
# MINIMAL WORKING EXAMPLE of sim_example_v22.R.
# Reduced dimensions for fast testing / debugging:
#   - n_locations   : 80  -> 15
#   - n_species     : 5   -> 2
#   - chains        : 3   -> 1
#   - iter          : 2000 -> 500
#   - warmup        : 1000 -> 200
#   - m_x / m_y     : 4/8 -> 2/4  (HSGP basis)
#   - adapt_delta   : 0.99 -> 0.90
#   - max_treedepth : 14   -> 10
#
# Sources run_qpcr_metabar_workflow_v22.R for all utilities.
#
# Usage (run from project root):
#   Rscript code/spatial_qPCR_mb/sim_example_v22_mwe.R [options]
#
# Optional command-line arguments (positional, in order):
#   1  seed        (integer, default 42)
#   2  chains      (integer, default 1)
#   3  iter        (integer, default 500)
#   4  warmup      (integer, default 200)
#   5  output_dir  (string,  default "/mmfs1/gscratch/coenv/rpkelly/output/sim_example_v22_mwe")

args      <- commandArgs(trailingOnly = TRUE)
`%||%`    <- function(x, y) if (is.null(x)) y else x

seed       <- if (length(args) >= 1) as.integer(args[1]) else 42L
chains     <- if (length(args) >= 2) as.integer(args[2]) else 1L
iter       <- if (length(args) >= 3) as.integer(args[3]) else 500L
warmup     <- if (length(args) >= 4) as.integer(args[4]) else 200L
output_dir <- if (length(args) >= 5) args[5] else
  "/mmfs1/gscratch/coenv/rpkelly/output/sim_example_v22_mwe"

cat("=== sim_example_v22_mwe.R ===\n")
cat(sprintf("seed=%d | chains=%d | iter=%d | warmup=%d\n",
            seed, chains, iter, warmup))
cat(sprintf("output_dir: %s\n\n", output_dir))

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

# --------------------------------------------------------------------------
# Source workflow
# --------------------------------------------------------------------------
workflow_file <- file.path("code", "spatial_qPCR_mb", "run_qpcr_metabar_workflow_v22.R")
if (!file.exists(workflow_file))
  stop("Cannot find workflow file: ", workflow_file,
       "\nRun this script from the project root directory.")

source(workflow_file)

# --------------------------------------------------------------------------
# Run (minimal settings)
# --------------------------------------------------------------------------
result <- run_workflow_v22(
  stan_file           = "code/spatial_qPCR_mb/quant_metabar_qpcr_v22.stan",
  n_locations         = 15,
  n_species           = 2,
  reference_species   = 1,
  n_qpcr_per_location = 3,
  chains              = chains,
  iter                = iter,
  warmup              = warmup,
  seed                = seed,
  # HSGP basis dimensions (reduced)
  m_x                 = 2,
  m_y                 = 4,
  max_rho             = 300,
  # True simulation parameters (unchanged)
  rho_km              = 50,
  gp_alpha_true       = NULL,
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
  # Priors (unchanged)
  rho_prior_mean      = 100.0,
  rho_prior_sd        = 200.0,
  gp_alpha_prior_mean = 1.0,
  gp_alpha_prior_sd   = 0.3,
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
  # HMC tuning (relaxed for speed)
  adapt_delta         = 0.90,
  max_treedepth       = 10,
  # Output
  save_outputs        = TRUE,
  output_dir          = output_dir
)

cat("\n=== Done ===\n")
cat(sprintf("Fit successful: %s\n", result$fit_successful))
if (result$fit_successful) {
  diag <- result$diagnostics
  cat(sprintf("Divergences: %s | Max treedepth hits: %s\n",
              diag$n_divergent, diag$max_treedepth_hits))
  cat(sprintf("Max Rhat: %.3f | # Rhat > 1.01: %d\n",
              diag$max_rhat, diag$n_high_rhat))
  cat(sprintf("Outputs written to: %s\n", output_dir))
}
