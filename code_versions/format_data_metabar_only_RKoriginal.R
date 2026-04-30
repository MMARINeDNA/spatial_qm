# =============================================================================
# Data Formatter + Fitter: quant_metabar_only.stan
# Poisson ZI + Beta-Binomial eDNA-only model (no qPCR)
# =============================================================================
#
# Functions:
#   format_data_metabar_only()    — reshape raw eDNA data to indexed arrays
#   make_stan_data_metabar_only() — assemble the Stan data list
#   fit_metabar_only()            — compile, sample, print diagnostics
#   plot_bb_params_metabar()      — BB dispersion parameter posteriors
#   plot_edna_fit_metabar()       — observed vs predicted eDNA proportions
#   plot_species_offsets_metabar()— species offset posteriors
#   plot_concentration_metabar()  — true vs predicted concentration (optional)
#   fit_and_plot_metabar_only()   — convenience wrapper: fit + all plots
#
# Expected input (edna_data):
#   Long-format data frame with one row per (PCR replicate, species) combination.
#   Required columns (names configurable via arguments):
#     - pcr_id_col     : PCR replicate identifier (character/factor)
#     - location_id_col: Location identifier (character/factor)
#     - species_col    : Species name (character/factor)
#     - reads_col      : Integer read count
#
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

`%||%` <- function(x, y) if (is.null(x)) y else x


# =============================================================================
# SECTION 1: IMPORT RAW DATA
# =============================================================================
eDNA_data<-read.csv("Mod2_data/RL2501_3cetacean_3fish_combined_otu.csv")


# =============================================================================
# SECTION 1: FORMAT RAW DATA
# =============================================================================

#' Format eDNA metabarcoding data for quant_metabar_only.stan
#'
#' @param edna_data         Long-format data frame (one row per replicate × species)
#' @param alpha_known       Numeric vector of amplification efficiencies, one per species,
#'                          in the same order as species_names
#' @param reference_species Integer index (1-based) of the reference species
#' @param species_names     Character vector of species to include, in model order.
#'                          If NULL, uses all species found in edna_data (sorted).
#' @param pcr_id_col        Column name for PCR replicate IDs
#' @param location_id_col   Column name for location IDs
#' @param species_col       Column name for species
#' @param reads_col         Column name for read counts
#'
#' @return Named list consumed by make_stan_data_metabar_only()
#'
format_data_metabar_only <- function(
    edna_data,
    alpha_known,
    reference_species  = 1,
    species_names      = NULL,
    pcr_id_col         = "pcr_replicate_id",
    location_id_col    = "location_id",
    species_col        = "species",
    reads_col          = "Nreads"
) {
  # --- Step 1: Standardise column names ---
  edna_work <- edna_data %>%
    rename(
      pcr_replicate_id = !!sym(pcr_id_col),
      location_id      = !!sym(location_id_col),
      species          = !!sym(species_col),
      Nreads           = !!sym(reads_col)
    )

  # --- Step 2: Filter to target species ---
  if (!is.null(species_names)) {
    edna_work   <- edna_work %>% filter(species %in% species_names)
    all_species <- species_names
  } else {
    all_species <- sort(unique(edna_work$species))
    if (length(all_species) == 0)
      stop("No species found in edna_data. Supply species_names or check data.")
  }

  # --- Step 3: Collapse haplotypes — sum reads within (species, location, replicate) ---
  edna_work <- edna_work %>%
    group_by(species, location_id, pcr_replicate_id) %>%
    summarise(Nreads = sum(Nreads, na.rm = TRUE), .groups = "drop")

  # --- Step 4: Build location index ---
  location_list <- edna_work %>%
    distinct(location_id) %>%
    arrange(location_id) %>%
    mutate(location_idx = row_number())

  N_locations <- nrow(location_list)

  # --- Step 5: Species index ---
  N_species <- length(all_species)

  if (length(alpha_known) != N_species)
    stop(sprintf("alpha_known length (%d) must match N_species (%d)",
                 length(alpha_known), N_species))

  sp_list <- data.frame(
    species     = all_species,
    species_idx = seq_len(N_species),
    stringsAsFactors = FALSE
  )

  # --- Step 6: Build eDNA observation (PCR-replicate) index ---
  edna_work <- edna_work %>%
    left_join(location_list, by = "location_id")

  pcr_info <- edna_work %>%
    group_by(pcr_replicate_id) %>%
    summarise(location_idx = first(location_idx), .groups = "drop") %>%
    arrange(location_idx, pcr_replicate_id) %>%
    mutate(edna_obs_idx = row_number())

  N_edna_obs <- nrow(pcr_info)

  edna_work <- edna_work %>%
    left_join(sp_list, by = "species") %>%
    left_join(pcr_info %>% select(pcr_replicate_id, edna_obs_idx),
              by = "pcr_replicate_id")

  # --- Step 7: Fill edna_obs_data matrix [N_edna_obs, N_species] ---
  edna_obs_data <- matrix(0L, nrow = N_edna_obs, ncol = N_species)
  for (i in seq_len(nrow(edna_work))) {
    obs_idx <- edna_work$edna_obs_idx[i]
    sp_idx  <- edna_work$species_idx[i]
    reads   <- edna_work$Nreads[i]
    if (!is.na(obs_idx) && !is.na(sp_idx) && !is.na(reads))
      edna_obs_data[obs_idx, sp_idx] <- as.integer(reads)
  }

  edna_obs_location <- as.integer(pcr_info$location_idx)

  has_edna <- rep(0L, N_locations)
  has_edna[unique(pcr_info$location_idx)] <- 1L

  # --- Summary ---
  cat("\n=== Formatted Data Summary (metabar_only) ===\n")
  cat(sprintf("Locations  : %d\n", N_locations))
  cat(sprintf("Species    : %d (%s)\n", N_species, paste(all_species, collapse = ", ")))
  cat(sprintf("Ref species: %s (index %d)\n", all_species[reference_species], reference_species))
  cat(sprintf("eDNA obs   : %d PCR replicates\n", N_edna_obs))
  cat(sprintf("Total reads: %d\n", sum(edna_obs_data)))

  list(
    N_species         = N_species,
    N_locations       = N_locations,
    N_edna_obs        = N_edna_obs,
    reference_species = as.integer(reference_species),
    edna_obs_data     = edna_obs_data,
    edna_obs_location = edna_obs_location,
    alpha_known       = alpha_known,
    has_edna          = has_edna,
    sp_list           = sp_list,
    location_list     = location_list
  )
}


# =============================================================================
# SECTION 2: ASSEMBLE STAN DATA LIST
# =============================================================================

#' Build Stan data list for quant_metabar_only.stan
#'
#' @param formatted_data      Output from format_data_metabar_only()
#' @param N_pcr               Number of PCR cycles (N_PCR in the model)
#' @param log_conc_prior_mean Prior mean for log_concentration. If NULL, defaults to 2.0.
#' @param log_conc_prior_sd   Prior SD for log_concentration (wide, e.g. 5.0)
#' @param beta0_prior_sd      Normal SD for BB dispersion intercept
#' @param gamma0_prior_mean   Normal mean for gamma0 (low-concentration overdispersion)
#' @param gamma0_prior_sd     Normal SD for gamma0
#' @param log_gamma1_prior_mean  Normal mean for log_gamma1
#' @param log_gamma1_prior_sd    Normal SD for log_gamma1
#'
#' @return Named list ready to pass to Stan
#'
make_stan_data_metabar_only <- function(
    formatted_data,
    N_pcr                 = 35,
    log_conc_prior_mean   = NULL,
    log_conc_prior_sd     = 5.0,
    beta0_prior_sd        = 5.0,
    gamma0_prior_mean     = 2.0,
    gamma0_prior_sd       = 2.0,
    log_gamma1_prior_mean = 0.0,
    log_gamma1_prior_sd   = 1.0
) {
  log_conc_prior_mean <- log_conc_prior_mean %||% 2.0

  cat(sprintf("log_conc_prior_mean = %.2f (~%.0f copies/uL)\n",
              log_conc_prior_mean, exp(log_conc_prior_mean)))
  cat(sprintf("log_conc_prior_sd   = %.1f\n", log_conc_prior_sd))
  cat(sprintf("BB dispersion: beta0_prior_sd=%.1f, gamma0~N(%.1f,%.1f), log_gamma1~N(%.1f,%.1f)\n",
              beta0_prior_sd, gamma0_prior_mean, gamma0_prior_sd,
              log_gamma1_prior_mean, log_gamma1_prior_sd))

  list(
    # Dimensions
    N_species         = formatted_data$N_species,
    N_locations       = formatted_data$N_locations,
    N_edna_obs        = formatted_data$N_edna_obs,

    # Reference species
    reference_species = formatted_data$reference_species,

    # eDNA data
    edna_obs_data     = formatted_data$edna_obs_data,
    edna_obs_location = formatted_data$edna_obs_location,

    # Amplification efficiencies
    alpha_known       = formatted_data$alpha_known,

    # PCR cycles
    N_PCR             = as.numeric(N_pcr),

    # Availability flags
    has_edna          = formatted_data$has_edna,

    # Concentration priors
    log_conc_prior_mean = log_conc_prior_mean,
    log_conc_prior_sd   = log_conc_prior_sd,

    # BB dispersion priors
    beta0_prior_sd        = beta0_prior_sd,
    gamma0_prior_mean     = gamma0_prior_mean,
    gamma0_prior_sd       = gamma0_prior_sd,
    log_gamma1_prior_mean = log_gamma1_prior_mean,
    log_gamma1_prior_sd   = log_gamma1_prior_sd
  )
}


# =============================================================================
# SECTION 3: THEME
# =============================================================================

.theme_mb <- function() {
  theme_bw() +
    theme(
      panel.background  = element_rect(fill = "gray98", color = NA),
      plot.background   = element_rect(fill = "white",  color = NA),
      panel.grid.major  = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.minor  = element_line(color = "gray92", linewidth = 0.2),
      panel.border      = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(color = "gray40", size = 10),
      strip.background  = element_rect(fill = "gray90", color = "gray70"),
      strip.text        = element_text(face = "bold", size = 10),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA)
    )
}


# =============================================================================
# SECTION 4: FIT MODEL
# =============================================================================

#' Compile and sample quant_metabar_only.stan
#'
#' @param formatted_data  Output from format_data_metabar_only()
#' @param stan_file       Path to quant_metabar_only.stan
#' @param stan_data       Pre-built Stan data list. If NULL, built from formatted_data
#'                        using make_stan_data_metabar_only() defaults.
#' @param chains          Number of MCMC chains
#' @param iter            Total iterations per chain (including warmup)
#' @param warmup          Warmup iterations per chain
#' @param seed            Random seed
#' @param adapt_delta     HMC adapt_delta (increase toward 1.0 for difficult posteriors)
#' @param max_treedepth   HMC max_treedepth
#'
#' @return rstan stanfit object
#'
fit_metabar_only <- function(
    formatted_data,
    stan_file     = "quant_metabar_only.stan",
    stan_data     = NULL,
    chains        = 3,
    iter          = 2000,
    warmup        = 1000,
    seed          = 42,
    adapt_delta   = 0.95,
    max_treedepth = 12
) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("rstan is required. Install with: install.packages('rstan')")

  if (is.null(stan_data)) {
    cat("\n=== Building Stan data ===\n")
    stan_data <- make_stan_data_metabar_only(formatted_data)
  }

  cat("\n=== Compiling Stan model ===\n")
  stan_mod <- rstan::stan_model(stan_file)

  cat(sprintf("\n=== Sampling (chains=%d, iter=%d, warmup=%d) ===\n",
              chains, iter, warmup))
  cat(sprintf("adapt_delta = %.3f | max_treedepth = %d\n", adapt_delta, max_treedepth))

  fit <- rstan::sampling(
    stan_mod, data = stan_data,
    chains = chains, iter = iter, warmup = warmup, seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  )

  # --- Diagnostics ---
  cat("\n=== Sampling Diagnostics ===\n")
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_div  <- tryCatch(
    sum(sapply(sampler_params, function(x) if (is.matrix(x)) sum(x[, "divergent__"])           else 0)),
    error = function(e) NA)
  max_td <- tryCatch(
    sum(sapply(sampler_params, function(x) if (is.matrix(x)) sum(x[, "treedepth__"] >= max_treedepth) else 0)),
    error = function(e) NA)

  fit_summary <- summary(fit)$summary
  max_rhat    <- max(fit_summary[, "Rhat"], na.rm = TRUE)
  n_high_rhat <- sum(fit_summary[, "Rhat"] > 1.01, na.rm = TRUE)

  cat(sprintf("Divergences      : %s\n", ifelse(is.na(n_div), "NA", n_div)))
  cat(sprintf("Max treedepth hit: %s\n", ifelse(is.na(max_td), "NA", max_td)))
  cat(sprintf("Max Rhat         : %.3f\n", max_rhat))
  cat(sprintf("# Rhat > 1.01    : %d\n", n_high_rhat))

  # Top Rhat parameters (flagging worst convergence)
  top_rhat <- data.frame(
    parameter = rownames(fit_summary),
    Rhat  = fit_summary[, "Rhat"],
    n_eff = fit_summary[, "n_eff"],
    sd    = fit_summary[, "sd"],
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(Rhat), sd > 1e-10) %>%
    arrange(desc(Rhat)) %>%
    head(10)

  cat("\n=== Top 10 parameters by Rhat ===\n")
  print(top_rhat, row.names = FALSE)

  # Key scalar parameter summaries
  cat("\n=== Derived: psi, species_offset ===\n")
  print(fit, pars = c("psi", "species_offset"))

  cat("\n=== BB Dispersion Parameters ===\n")
  print(fit, pars = c("beta0", "beta1", "gamma0", "log_gamma1", "mean_phi"))

  attr(fit, "top_rhat")   <- top_rhat
  attr(fit, "stan_data")  <- stan_data
  fit
}


# =============================================================================
# SECTION 5: PLOT FUNCTIONS
# =============================================================================

#' Plot BB dispersion parameter posteriors
#'
#' @param fit            stanfit from fit_metabar_only()
#' @param true_values    Optional named list with beta0_true, beta1_true,
#'                       gamma0_true, log_gamma1_true for simulation checks
#'
plot_bb_params_metabar <- function(fit, true_values = NULL) {
  posterior <- rstan::extract(fit)

  param_data <- bind_rows(
    data.frame(parameter = "beta0",      value = posterior$beta0),
    data.frame(parameter = "beta1",      value = posterior$beta1),
    data.frame(parameter = "gamma0",     value = posterior$gamma0),
    data.frame(parameter = "log_gamma1", value = posterior$log_gamma1),
    data.frame(parameter = "mean_phi",   value = posterior$mean_phi)
  )

  true_df <- NULL
  if (!is.null(true_values)) {
    true_df <- bind_rows(
      data.frame(parameter = "beta0",      true_value = true_values$beta0_true      %||% NA),
      data.frame(parameter = "beta1",      true_value = true_values$beta1_true      %||% NA),
      data.frame(parameter = "gamma0",     true_value = true_values$gamma0_true     %||% NA),
      data.frame(parameter = "log_gamma1", true_value = true_values$log_gamma1_true %||% NA)
    ) %>% filter(!is.na(true_value))
    if (nrow(true_df) == 0) true_df <- NULL
  }

  p <- ggplot(param_data, aes(x = value)) +
    geom_density(fill = "darkorange3", alpha = 0.5, color = "darkorange4") +
    facet_wrap(~ parameter, scales = "free", ncol = 3) +
    labs(
      title    = "Beta-Binomial Dispersion Parameters",
      subtitle = if (!is.null(true_df)) "Red dashed = true values" else NULL,
      x = "Value", y = "Density"
    ) +
    .theme_mb()

  if (!is.null(true_df))
    p <- p + geom_vline(data = true_df, aes(xintercept = true_value),
                        color = "red", linetype = "dashed", linewidth = 1)
  p
}


#' Plot observed vs posterior eDNA proportions per species
#'
#' @param fit            stanfit from fit_metabar_only()
#' @param formatted_data Output from format_data_metabar_only()
#'
plot_edna_fit_metabar <- function(fit, formatted_data) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species
  N_species     <- length(species_names)

  plots <- list()
  for (sp_idx in seq_len(N_species)) {
    sp <- species_names[sp_idx]

    p <- tryCatch({
      obs_df <- data.frame(
        edna_obs_idx = seq_len(formatted_data$N_edna_obs),
        location_idx = formatted_data$edna_obs_location,
        obs_reads    = formatted_data$edna_obs_data[, sp_idx],
        total_reads  = rowSums(formatted_data$edna_obs_data)
      ) %>%
        mutate(obs_prop = obs_reads / total_reads) %>%
        group_by(location_idx) %>%
        summarise(obs_prop_mean = mean(obs_prop, na.rm = TRUE),
                  n_reps = n(), .groups = "drop")

      mu_samples <- posterior$mu_metabar_loc[, , sp_idx]
      pred_df <- data.frame(
        location_idx = seq_len(ncol(mu_samples)),
        pred_mean    = colMeans(mu_samples),
        pred_q025    = apply(mu_samples, 2, quantile, 0.025),
        pred_q975    = apply(mu_samples, 2, quantile, 0.975)
      )

      compare_df <- obs_df %>% left_join(pred_df, by = "location_idx")

      ggplot(compare_df, aes(x = obs_prop_mean, y = pred_mean)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
        geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975),
                      width = 0, alpha = 0.3, color = "steelblue") +
        geom_point(size = 2, color = "steelblue") +
        coord_equal() +
        labs(
          title    = paste("eDNA Proportions:", gsub("_", " ", sp)),
          subtitle = "Observed vs. Posterior mean ± 95% CI",
          x = "Observed Proportion", y = "Posterior Proportion"
        ) +
        .theme_mb()
    }, error = function(e) {
      warning(paste("Error in eDNA plot for", sp, ":", e$message)); NULL
    })

    if (!is.null(p)) plots[[sp]] <- p
  }
  plots
}


#' Plot species offset posteriors
#'
#' @param fit            stanfit from fit_metabar_only()
#' @param formatted_data Output from format_data_metabar_only()
#' @param true_values    Optional named list with species_offsets (numeric vector)
#'
plot_species_offsets_metabar <- function(fit, formatted_data, true_values = NULL) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species

  if (!"species_offset" %in% names(posterior)) {
    warning("species_offset not found in posterior."); return(NULL)
  }

  offset_df <- do.call(bind_rows, lapply(seq_along(species_names), function(i) {
    data.frame(species = species_names[i], value = posterior$species_offset[, i])
  }))

  true_df <- NULL
  if (!is.null(true_values) && !is.null(true_values$species_offsets)) {
    true_df <- data.frame(
      species    = species_names,
      true_value = true_values$species_offsets
    )
  }

  p <- ggplot(offset_df, aes(x = value, fill = species)) +
    geom_density(alpha = 0.6, color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray60") +
    facet_wrap(~ species, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = "Species Offset Posteriors",
      subtitle = if (!is.null(true_df)) "Red dashed = true values" else "Relative to grand mean",
      x = "Offset (log scale)", y = "Density"
    ) +
    .theme_mb() +
    theme(legend.position = "none")

  if (!is.null(true_df))
    p <- p + geom_vline(data = true_df, aes(xintercept = true_value),
                        color = "red", linetype = "dashed", linewidth = 1)
  p
}


#' Plot posterior concentration vs true concentration (simulation check)
#'
#' @param fit              stanfit from fit_metabar_only()
#' @param formatted_data   Output from format_data_metabar_only()
#' @param true_concentration Matrix [N_locations, N_species] of true concentrations
#'
plot_concentration_metabar <- function(fit, formatted_data, true_concentration) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species
  N_species     <- length(species_names)
  conc_samples  <- posterior$concentration

  plots <- list()
  for (sp_idx in seq_len(N_species)) {
    sp <- species_names[sp_idx]

    p <- tryCatch({
      sp_conc <- conc_samples[, , sp_idx]
      df <- data.frame(
        location_idx = seq_len(ncol(sp_conc)),
        true_conc    = true_concentration[, sp_idx],
        pred_mean    = colMeans(sp_conc),
        pred_q025    = apply(sp_conc, 2, quantile, 0.025),
        pred_q975    = apply(sp_conc, 2, quantile, 0.975)
      ) %>%
        mutate(covered = true_conc >= pred_q025 & true_conc <= pred_q975)

      coverage <- mean(df$covered) * 100
      cor_val  <- cor(log(df$true_conc), log(df$pred_mean))

      ggplot(df, aes(x = true_conc, y = pred_mean)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
        geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975),
                      width = 0, alpha = 0.4, color = "darkorange3") +
        geom_point(aes(color = covered), size = 2.5) +
        scale_color_manual(values = c("FALSE" = "red", "TRUE" = "darkorange"),
                           name = "In 95% CI") +
        scale_x_log10() + scale_y_log10() +
        labs(
          title    = paste("True vs Predicted:", gsub("_", " ", sp)),
          subtitle = sprintf("95%% CI coverage: %.0f%%  |  r(log) = %.2f", coverage, cor_val),
          x = "True Concentration (copies/uL)",
          y = "Posterior Concentration (copies/uL)"
        ) +
        .theme_mb() +
        theme(legend.position = "bottom")
    }, error = function(e) {
      warning(paste("Error in concentration plot for", sp, ":", e$message)); NULL
    })

    if (!is.null(p)) plots[[sp]] <- p
  }
  plots
}


# =============================================================================
# SECTION 6: FIT + PLOT WRAPPER
# =============================================================================

#' Fit quant_metabar_only.stan and produce all diagnostic plots
#'
#' @param formatted_data  Output from format_data_metabar_only()
#' @param stan_file       Path to quant_metabar_only.stan
#' @param stan_data       Pre-built Stan data list (NULL = build with defaults)
#' @param chains,iter,warmup,seed,adapt_delta,max_treedepth  Sampling controls
#' @param true_values     Optional list for simulation checks. Named fields:
#'                          beta0_true, beta1_true, gamma0_true, log_gamma1_true,
#'                          species_offsets, true_concentration [N_loc × N_sp matrix]
#' @param save_outputs    Write PNG files to output_dir
#' @param output_dir      Directory for saved plots
#'
#' @return List with fit, stan_data, all plot objects, and diagnostics
#'
fit_and_plot_metabar_only <- function(
    formatted_data,
    stan_file     = "code/qPCR_mb/quant_metabar_only.stan",
    stan_data     = NULL,
    chains        = 3,
    iter          = 2000,
    warmup        = 1000,
    seed          = 42,
    adapt_delta   = 0.95,
    max_treedepth = 12,
    true_values   = NULL,
    save_outputs  = FALSE,
    output_dir    = "code/qPCR_mb/outputs_metabar_only"
) {
  if (save_outputs && !dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

  # Fit
  fit <- fit_metabar_only(
    formatted_data = formatted_data,
    stan_file      = stan_file,
    stan_data      = stan_data,
    chains         = chains,
    iter           = iter,
    warmup         = warmup,
    seed           = seed,
    adapt_delta    = adapt_delta,
    max_treedepth  = max_treedepth
  )

  # Retrieve stan_data used (may have been built inside fit_metabar_only)
  used_stan_data <- attr(fit, "stan_data")

  # Plots
  cat("\n=== Generating plots ===\n")

  bb_plot      <- plot_bb_params_metabar(fit, true_values)
  edna_plots   <- plot_edna_fit_metabar(fit, formatted_data)
  offset_plot  <- plot_species_offsets_metabar(fit, formatted_data, true_values)

  # Combine eDNA per-species plots into one panel
  edna_plots_clean <- Filter(Negate(is.null), edna_plots)
  combined_edna <- if (length(edna_plots_clean) > 0) {
    n_cols <- ceiling(sqrt(length(edna_plots_clean)))
    wrap_plots(edna_plots_clean, ncol = n_cols) +
      plot_annotation(title = "eDNA Model Fit (Proportions)",
                      theme = theme(plot.title = element_text(size = 14, face = "bold")))
  } else NULL

  # Optional: true vs predicted concentration
  conc_plots <- NULL
  combined_conc <- NULL
  if (!is.null(true_values) && !is.null(true_values$true_concentration)) {
    conc_plots <- plot_concentration_metabar(fit, formatted_data,
                                             true_values$true_concentration)
    conc_plots_clean <- Filter(Negate(is.null), conc_plots)
    if (length(conc_plots_clean) > 0) {
      n_cols <- ceiling(sqrt(length(conc_plots_clean)))
      combined_conc <- wrap_plots(conc_plots_clean, ncol = n_cols) +
        plot_annotation(title = "True vs Predicted Concentration",
                        theme = theme(plot.title = element_text(size = 14, face = "bold")))
    }
  }

  # Print plots to active device
  if (!is.null(bb_plot))      print(bb_plot)
  if (!is.null(combined_edna)) print(combined_edna)
  if (!is.null(offset_plot))  print(offset_plot)
  if (!is.null(combined_conc)) print(combined_conc)

  # Save
  if (save_outputs) {
    if (!is.null(bb_plot))
      ggsave(file.path(output_dir, "bb_params.png"),     bb_plot,       width = 12, height = 6,  dpi = 150)
    if (!is.null(combined_edna))
      ggsave(file.path(output_dir, "edna_fit.png"),      combined_edna, width = 10, height = 8,  dpi = 150)
    if (!is.null(offset_plot))
      ggsave(file.path(output_dir, "species_offsets.png"), offset_plot, width = 10, height = 8,  dpi = 150)
    if (!is.null(combined_conc))
      ggsave(file.path(output_dir, "concentration.png"), combined_conc, width = 12, height = 10, dpi = 150)
    top_rhat <- attr(fit, "top_rhat")
    if (!is.null(top_rhat))
      write.csv(top_rhat, file.path(output_dir, "top_rhat.csv"), row.names = FALSE)
    cat("Saved outputs to:", output_dir, "\n")
  }

  list(
    fit            = fit,
    stan_data      = used_stan_data,
    formatted_data = formatted_data,
    bb_plot        = bb_plot,
    edna_plots     = edna_plots,
    combined_edna  = combined_edna,
    offset_plot    = offset_plot,
    conc_plots     = conc_plots,
    combined_conc  = combined_conc,
    diagnostics    = list(
      top_rhat = attr(fit, "top_rhat")
    )
  )
}

