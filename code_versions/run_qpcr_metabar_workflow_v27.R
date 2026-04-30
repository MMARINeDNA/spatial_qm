# run_qpcr_metabar_workflow_v27.R — SELF-CONTAINED
# Companion workflow for quant_metabar_qpcr_v27.stan
#
# v27 model: Factor-Analytic HSGP — K_f shared latent fields + lower-triangular loading.
#   Each species' spatial effect is a weighted sum of K_f latent fields:
#     spatial_effect[loc, sp] = sum_k loading[sp,k] * latent_field[loc,k]
#   Lower-triangular loading with positive diagonal resolves rotational and
#   sign-flip non-identifiabilities. gp_alpha_factor is fixed data (no prior),
#   removing the multiplicative ridge between alpha, rho, and loading_diag.
#   No depth strata (contrast with v24/v25, which were abandoned).
#   Everything else is identical to v22:
#     - Poisson ZI + Beta-Binomial eDNA observation model
#     - Per-replicate hurdle qPCR model
#     - Decoupled sigma_nugget hierarchy (independent of loading amplitude)
#     - No nugget mean-centering
#     - Dispersion model (beta0 / beta1 / gamma0 / gamma1)
#
# This file has NO source() dependencies on earlier workflow files.
# All required utilities, formatting, extraction, and plotting functions are
# included directly.
#
# Stan model: code/spatial_qPCR_mb/quant_metabar_qpcr_v27.stan

library(tidyverse)
library(rstan)
library(sf)
library(akima)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# =============================================================================
# SECTION 1: UTILITIES
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

theme_light_custom <- function() {
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

get_utm_zone <- function(lon) floor((lon + 180) / 6) + 1

convert_to_utm <- function(lon, lat, target_zone = NULL) {
  if (is.null(target_zone)) target_zone <- get_utm_zone(mean(lon))
  hemisphere <- ifelse(mean(lat) >= 0, "north", "south")
  epsg_code  <- if (hemisphere == "north") 32600 + target_zone else 32700 + target_zone
  pts        <- st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  pts_utm    <- st_transform(pts, crs = epsg_code)
  coords_utm <- st_coordinates(pts_utm)
  data.frame(x_utm = coords_utm[, 1] / 1000, y_utm = coords_utm[, 2] / 1000)
}

cluster_coordinates <- function(coords, tolerance = 0.001) {
  if (nrow(coords) == 0) return(coords %>% mutate(cluster_id = integer(0)))
  if (nrow(coords) == 1) return(coords %>% mutate(cluster_id = 1L))
  coords$cluster_id <- NA_integer_
  current_cluster   <- 0L
  for (i in seq_len(nrow(coords))) {
    if (is.na(coords$cluster_id[i])) {
      current_cluster           <- current_cluster + 1L
      coords$cluster_id[i]      <- current_cluster
      distances                 <- sqrt((coords$lon - coords$lon[i])^2 +
                                        (coords$lat - coords$lat[i])^2)
      nearby                    <- which(distances <= tolerance & is.na(coords$cluster_id))
      coords$cluster_id[nearby] <- current_cluster
    }
  }
  coords
}

# NA-aware separable Gaussian blur for interpolated grids.
# Operates on a matrix; NAs are excluded from the kernel sum (border/hole safe).
smooth_grid_gaussian <- function(mat, sigma = 2) {
  r  <- max(1L, round(3 * sigma))
  ks <- seq(-r, r)
  w  <- dnorm(ks, sd = sigma); w <- w / sum(w)

  blur_1d <- function(v) {
    n      <- length(v)
    padded <- c(rep(NA_real_, r), v, rep(NA_real_, r))
    num    <- rep(0, n); wt <- rep(0, n)
    for (k in seq_along(w)) {
      vals        <- padded[k:(k + n - 1)]
      ok          <- !is.na(vals)
      num[ok]     <- num[ok] + w[k] * vals[ok]
      wt[ok]      <- wt[ok]  + w[k]
    }
    ifelse(wt > 0, num / wt, NA_real_)
  }

  # Row-wise then column-wise (separable)
  row_sm <- t(apply(mat, 1, blur_1d))
  apply(row_sm, 2, blur_1d) |> matrix(nrow = nrow(mat))
}

interpolate_surface <- function(x, y, z, resolution = 100, log_transform = TRUE,
                                smooth_sigma = 2) {
  x_range  <- range(x, na.rm = TRUE)
  y_range  <- range(y, na.rm = TRUE)
  x_expand <- diff(x_range) * 0.05
  y_expand <- diff(y_range) * 0.05
  x_grid   <- seq(x_range[1] - x_expand, x_range[2] + x_expand, length.out = resolution)
  y_grid   <- seq(y_range[1] - y_expand, y_range[2] + y_expand, length.out = resolution)
  if (log_transform && all(z > 0, na.rm = TRUE)) {
    z_interp       <- log(z)
    back_transform <- TRUE
  } else {
    z_interp       <- z
    back_transform <- FALSE
  }
  interp_result <- akima::interp(x = x, y = y, z = z_interp,
                                  xo = x_grid, yo = y_grid,
                                  linear = TRUE, extrap = FALSE, remove = FALSE,
                                  duplicate = "mean")
  if (!is.null(smooth_sigma) && smooth_sigma > 0)
    interp_result$z <- smooth_grid_gaussian(interp_result$z, sigma = smooth_sigma)
  interp_df <- expand.grid(x = interp_result$x, y = interp_result$y) %>%
    mutate(value = as.vector(interp_result$z))
  if (back_transform) interp_df <- interp_df %>% mutate(value = exp(value))
  filter(interp_df, !is.na(value))
}

#' Extract posterior concentration samples into a tidy data frame
extract_concentration <- function(fit, formatted_data, parameter = "concentration") {
  posterior     <- rstan::extract(fit)
  locations     <- formatted_data$location_list
  species_names <- formatted_data$sp_list$species

  if (parameter == "concentration") {
    samples     <- posterior$concentration
    result_list <- list()
    for (sp_idx in seq_along(species_names)) {
      sp_samples <- samples[, , sp_idx]
      summary_df <- data.frame(
        location_idx = 1:ncol(sp_samples),
        species      = species_names[sp_idx],
        mean         = colMeans(sp_samples),
        median       = apply(sp_samples, 2, median),
        sd           = apply(sp_samples, 2, sd),
        q025         = apply(sp_samples, 2, quantile, 0.025),
        q975         = apply(sp_samples, 2, quantile, 0.975)
      )
      result_list[[sp_idx]] <- locations %>% left_join(summary_df, by = "location_idx")
    }
    result <- bind_rows(result_list)
  } else if (parameter == "total_concentration") {
    samples <- posterior$total_concentration
    result  <- data.frame(
      location_idx = 1:ncol(samples),
      mean         = colMeans(samples),
      median       = apply(samples, 2, median),
      sd           = apply(samples, 2, sd),
      q025         = apply(samples, 2, quantile, 0.025),
      q975         = apply(samples, 2, quantile, 0.975)
    ) %>% left_join(locations, by = "location_idx")
  }
  attr(result, "parameter") <- parameter
  result
}

#' Re-index true_concentration rows to match formatted_location_list order.
#' Handles cases where all-zero-qPCR locations were excluded during formatting.
match_true_to_formatted <- function(true_concentration, sim_locations, formatted_location_list) {
  n_fmt   <- nrow(formatted_location_list)
  n_sp    <- ncol(true_concentration)
  matched <- matrix(NA_real_, nrow = n_fmt, ncol = n_sp)
  for (i in seq_len(n_fmt)) {
    dists        <- sqrt((sim_locations$lon - formatted_location_list$lon[i])^2 +
                         (sim_locations$lat - formatted_location_list$lat[i])^2)
    matched[i, ] <- true_concentration[which.min(dists), ]
  }
  matched
}

#' Return the n parameters with the highest Rhat for convergence diagnosis
get_top_rhat <- function(fit, n = 10) {
  fit_summary <- summary(fit)$summary
  df <- data.frame(parameter = rownames(fit_summary),
                   Rhat  = fit_summary[, "Rhat"],
                   n_eff = fit_summary[, "n_eff"],
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$Rhat), , drop = FALSE]
  df[order(df$Rhat, decreasing = TRUE), , drop = FALSE][seq_len(min(n, nrow(df))), , drop = FALSE]
}

#' Arrange a list of ggplot objects into a single patchwork figure
combine_plots <- function(plot_list, title = "") {
  plot_list <- plot_list[!sapply(plot_list, is.null)]
  n_plots   <- length(plot_list)
  if (n_plots == 0)
    return(ggplot() + theme_void() + labs(title = paste(title, "(no plots available)")))
  if (n_plots == 1)
    return(plot_list[[1]] + plot_annotation(title = title))
  wrap_plots(plot_list, ncol = ceiling(sqrt(n_plots))) +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
}


# =============================================================================
# SECTION 2: DATA FORMATTING
# =============================================================================

#' Format raw eDNA + qPCR data into the list expected by make_stan_data_v27
#'
#' Locations with all-zero qPCR detection are excluded before fitting (the
#' hurdle qPCR likelihood requires at least one detected replicate per location).
#'
#' @param edna_data      Data frame: one row per (pcr_replicate, species).
#' @param qpcr_data      Data frame: one row per qPCR replicate.
#' @param alpha_known    Numeric vector length N_species. PCR amplification
#'   efficiencies on log scale relative to mean; reference species = 0.
#' @param reference_species Integer index (1-based) of the qPCR target species.
#' @param species_names  Optional character vector of species names in order.
#' @param pcr_id_col     Column name in edna_data carrying PCR replicate ID.
#' @param lon_col, lat_col Column names for longitude and latitude.
#' @param species_col    Column name for species identifier.
#' @param reads_col      Column name for read counts.
#' @param qpcr_detected_col Column name for binary detection (0/1).
#' @param qpcr_log_conc_col Column name for log-concentration when detected.
#' @param coord_tolerance Decimal-degree tolerance for clustering nearby coords.
#' @return Named list for use in make_stan_data_v27.
format_qpcr_metabar_data_v19 <- function(
    edna_data,
    qpcr_data,
    alpha_known,
    reference_species = 1,
    species_names     = NULL,
    pcr_id_col        = "pcr_replicate_id",
    location_id_col   = "location_id",
    lon_col           = "lon",
    lat_col           = "lat",
    species_col       = "species",
    reads_col         = "Nreads",
    qpcr_detected_col = "detected",
    qpcr_log_conc_col = "qpcr_log_conc"
) {
  edna_work <- edna_data %>%
    rename(pcr_replicate_id = !!sym(pcr_id_col),
           location_id      = !!sym(location_id_col),
           lon = !!sym(lon_col), lat = !!sym(lat_col),
           species = !!sym(species_col), Nreads = !!sym(reads_col))
  qpcr_work_raw <- qpcr_data %>%
    rename(location_id   = !!sym(location_id_col),
           lon = !!sym(lon_col), lat = !!sym(lat_col),
           detected      = !!sym(qpcr_detected_col),
           qpcr_log_conc = !!sym(qpcr_log_conc_col))

  # Exclude locations where all qPCR reps are non-detects
  zero_locs_qpcr <- qpcr_work_raw %>%
    group_by(location_id) %>%
    summarise(any_detected = any(detected == 1), .groups = "drop") %>%
    filter(!any_detected) %>%
    pull(location_id)

  if (length(zero_locs_qpcr) > 0) {
    cat(sprintf("Excluding %d location(s) with all-zero qPCR detection: %s\n",
                length(zero_locs_qpcr), paste(as.character(zero_locs_qpcr), collapse = ", ")))
    edna_work     <- edna_work     %>% filter(!location_id %in% zero_locs_qpcr)
    qpcr_work_raw <- qpcr_work_raw %>% filter(!location_id %in% zero_locs_qpcr)
  } else {
    cat("No all-zero qPCR locations to exclude.\n")
  }

  # Build location list from unique location_ids; mean lat/lon per location_id
  location_coords <- bind_rows(
    edna_work     %>% distinct(location_id, lon, lat),
    qpcr_work_raw %>% distinct(location_id, lon, lat)
  ) %>%
    group_by(location_id) %>%
    summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")

  location_list <- location_coords %>%
    arrange(location_id) %>%
    mutate(location_idx = row_number())
  N_locations <- nrow(location_list)

  utm_coords <- convert_to_utm(location_list$lon, location_list$lat)
  location_list$x_utm          <- utm_coords$x_utm
  location_list$y_utm          <- utm_coords$y_utm
  utm_center                   <- list(x_mean = mean(location_list$x_utm),
                                       y_mean = mean(location_list$y_utm))
  location_list$x_utm_centered <- location_list$x_utm - utm_center$x_mean
  location_list$y_utm_centered <- location_list$y_utm - utm_center$y_mean

  if (is.null(species_names)) {
    all_species <- unique(edna_work$species)
    if (length(all_species) == 0) stop("No species found. Provide species_names.")
  } else {
    all_species <- species_names
  }
  sp_list   <- data.frame(species = all_species, species_idx = seq_along(all_species),
                           stringsAsFactors = FALSE)
  N_species <- length(all_species)
  if (length(alpha_known) != N_species)
    stop(sprintf("alpha_known length (%d) must match N_species (%d)",
                 length(alpha_known), N_species))

  # Assign location_idx to PCR replicates via location_id join (no distance arithmetic)
  loc_idx_map <- location_list %>% dplyr::select(location_id, location_idx)

  edna_obs <- edna_work %>%
    group_by(pcr_replicate_id) %>%
    summarise(location_id = first(location_id), .groups = "drop") %>%
    left_join(loc_idx_map, by = "location_id") %>%
    mutate(edna_obs_idx = row_number())

  N_edna_obs <- nrow(edna_obs)

  edna_work <- edna_work %>%
    left_join(sp_list, by = "species") %>%
    left_join(edna_obs %>% dplyr::select(pcr_replicate_id, edna_obs_idx, location_idx),
              by = "pcr_replicate_id")

  edna_obs_data <- matrix(0L, nrow = N_edna_obs, ncol = N_species)
  for (i in seq_len(nrow(edna_work))) {
    obs_idx <- edna_work$edna_obs_idx[i]
    sp_idx  <- edna_work$species_idx[i]
    reads   <- edna_work$Nreads[i]
    if (!is.na(obs_idx) && !is.na(sp_idx) && !is.na(reads))
      edna_obs_data[obs_idx, sp_idx] <- as.integer(reads)
  }
  edna_obs_location               <- as.integer(edna_obs$location_idx)
  has_edna                        <- rep(0L, N_locations)
  has_edna[unique(edna_obs$location_idx)] <- 1L

  # Assign location_idx to qPCR replicates via location_id join
  qpcr_work <- qpcr_work_raw %>%
    left_join(loc_idx_map, by = "location_id")
  N_qpcr_reps       <- nrow(qpcr_work)
  qpcr_rep_location <- as.integer(qpcr_work$location_idx)
  qpcr_rep_detected <- as.integer(qpcr_work$detected)
  qpcr_rep_log_conc <- ifelse(qpcr_work$detected == 1, qpcr_work$qpcr_log_conc, 0.0)
  has_qpcr          <- rep(0L, N_locations)
  has_qpcr[unique(qpcr_work$location_idx)] <- 1L
  n_detected        <- sum(qpcr_rep_detected)

  pcr_per_location <- edna_obs %>%
    group_by(location_idx) %>%
    summarise(n_pcr_replicates = n(), .groups = "drop")

  cat(sprintf("Locations: %d | Species: %d | Ref: %s (idx %d)\n",
              N_locations, N_species, all_species[reference_species], reference_species))
  cat(sprintf("eDNA obs: %d | qPCR reps: %d total, %d detected (%.0f%%)\n",
              N_edna_obs, N_qpcr_reps, n_detected, 100 * n_detected / max(N_qpcr_reps, 1)))
  cat(sprintf("UTM extent: %.1f km (x) x %.1f km (y)\n",
              diff(range(location_list$x_utm)), diff(range(location_list$y_utm))))

  list(
    N_species          = N_species,
    N_locations        = N_locations,
    N_edna_obs         = N_edna_obs,
    N_qpcr_reps        = N_qpcr_reps,
    reference_species  = as.integer(reference_species),
    x_utm              = location_list$x_utm_centered,
    y_utm              = location_list$y_utm_centered,
    edna_obs_data      = edna_obs_data,
    edna_obs_location  = edna_obs_location,
    qpcr_rep_location  = qpcr_rep_location,
    qpcr_rep_detected  = qpcr_rep_detected,
    qpcr_rep_log_conc  = qpcr_rep_log_conc,
    alpha_known        = alpha_known,
    has_edna           = has_edna,
    has_qpcr           = has_qpcr,
    sp_list            = sp_list,
    location_list      = location_list,
    edna_obs           = edna_obs,
    pcr_per_location   = pcr_per_location,
    utm_center         = utm_center
  )
}

# v27 uses the same formatter as v19 (no structural changes needed)
format_qpcr_metabar_data_v27 <- format_qpcr_metabar_data_v19


# =============================================================================
# SECTION 3: SIMULATION
# =============================================================================

#' Simulate joint eDNA + qPCR data for v27 (factor-analytic GP)
#'
#' Generates spatially correlated log-concentrations via K_f shared isotropic SE
#' GP latent fields combined through a lower-triangular loading matrix.
#' Observation models:
#'   eDNA  — Poisson copies -> zero-truncated multinomial proportions -> Beta-Binomial reads
#'   qPCR  — Bernoulli detection + Normal log-conc when detected
#'
#' @param n_locations       Number of sampling locations.
#' @param n_species         Number of species.
#' @param K_f               Number of shared latent GP fields.
#' @param n_qpcr_per_location qPCR replicates per location.
#' @param reference_species Index (1-based) of the qPCR target species.
#' @param rho_factor_km     GP length scale (km); scalar or K_f-length vector.
#' @param gp_alpha_factor_true GP marginal SD per factor; NULL = rep(1, K_f).
#' @param loading_true      True loading matrix [n_species, K_f] lower-triangular;
#'   NULL = random lower-triangular with positive diagonal.
#' @param sigma_nugget_true Independent noise SD; scalar or n_species vector.
#' @param log_lambda_mean   Mean log-copies per aliquot; lower = more zeros.
#' @param mu_species_true   Per-species log-concentration mean; NULL = random.
#' @param sigma_qpcr_true   qPCR measurement SD (log scale).
#' @param qpcr_detection_intercept_true Logistic detection intercept.
#' @param qpcr_detection_slope_true     Logistic detection slope (log-conc).
#' @param N_PCR             Total PCR cycles.
#' @param total_reads_per_rep eDNA sequencing depth per PCR replicate.
#' @param beta0_true,beta1_true,gamma0_true,gamma1_true Dispersion parameters.
#' @param lat_range,lon_range Coordinate bounding box for simulated locations.
#' @param alpha_known PCR efficiency offsets (NULL = small random values).
#' @param seed Random seed.
#' @return List: edna_data, qpcr_data, locations, species_names, reference_species,
#'   alpha_known, true_values.
simulate_qpcr_metabar_data_v27 <- function(
    n_locations = 25,
    n_species = 4,
    K_f = 3,
    n_qpcr_per_location = 3,
    reference_species = 1,
    rho_factor_km = 50,
    gp_alpha_factor_true = NULL,
    loading_true = NULL,
    sigma_nugget_true = 0.3,
    log_lambda_mean = log(5),
    mu_species_true = NULL,
    sigma_qpcr_true = 0.3,
    qpcr_detection_intercept_true = -2.0,
    qpcr_detection_slope_true = 1.5,
    N_PCR = 35,
    total_reads_per_rep = 10000,
    beta0_true = 5.0,
    beta1_true = -2.0,
    gamma0_true = 2.0,
    gamma1_true = exp(1.0),
    lat_range = c(42, 48),
    lon_range = c(-126, -124),
    alpha_known = NULL,
    seed = 42
) {
  set.seed(seed)
  if (beta1_true > 0) stop("beta1_true must be <= 0")
  if (K_f > n_species) stop("K_f must be <= n_species")

  species_names   <- paste("species", 1:n_species, sep = "_")
  sigma_nugget_sp <- if (length(sigma_nugget_true) == 1) rep(sigma_nugget_true, n_species) else sigma_nugget_true
  if (length(sigma_nugget_sp) != n_species) stop("sigma_nugget_true must be scalar or length n_species")

  gp_alpha_factor_true <- if (is.null(gp_alpha_factor_true)) rep(1.0, K_f) else gp_alpha_factor_true
  rho_k <- if (length(rho_factor_km) == 1) rep(rho_factor_km, K_f) else rho_factor_km
  if (length(rho_k) != K_f) stop("rho_factor_km must be scalar or length K_f")
  if (length(gp_alpha_factor_true) != K_f) stop("gp_alpha_factor_true must be scalar or length K_f")

  if (is.null(alpha_known)) {
    alpha_known <- rnorm(n_species, 0, 0.02)
    alpha_known[reference_species] <- 0
  }

  lat <- runif(n_locations, lat_range[1], lat_range[2])
  lon <- runif(n_locations, lon_range[1], lon_range[2])
  locations <- data.frame(
    location_id = paste0("LOC_", sprintf("%02d", 1:n_locations)),
    lat = lat, lon = lon
  )
  utm_coords       <- convert_to_utm(lon, lat)
  locations$x_utm  <- utm_coords$x_utm
  locations$y_utm  <- utm_coords$y_utm
  x_km             <- locations$x_utm - mean(locations$x_utm)
  y_km             <- locations$y_utm - mean(locations$y_utm)

  if (!is.null(mu_species_true)) {
    if (length(mu_species_true) != n_species) stop("mu_species_true must be length n_species")
    species_offsets <- mu_species_true - mean(mu_species_true)
  } else {
    species_offsets <- numeric(n_species)
    for (sp in 1:n_species)
      if (sp != reference_species) species_offsets[sp] <- rnorm(1, 0, 1.0)
    mu_species_true <- log_lambda_mean + species_offsets
  }

  # Build true lower-triangular loading matrix [n_species, K_f]
  if (is.null(loading_true)) {
    loading_diag_true    <- abs(rnorm(K_f, 1.0, 0.3))
    n_offdiag            <- K_f * n_species - K_f * (K_f + 1) / 2
    loading_offdiag_true <- rnorm(n_offdiag, 0, 0.5)
    loading_true <- matrix(0, n_species, K_f)
    offdiag_idx  <- 1
    for (j in 1:K_f) {
      for (i in 1:n_species) {
        if (i < j)
          loading_true[i, j] <- 0
        else if (i == j)
          loading_true[i, j] <- loading_diag_true[j]
        else {
          loading_true[i, j] <- loading_offdiag_true[offdiag_idx]
          offdiag_idx        <- offdiag_idx + 1
        }
      }
    }
  }

  # Draw K_f latent fields using isotropic SE GP, mean-centered
  latent_fields <- matrix(0, n_locations, K_f)
  for (k in 1:K_f) {
    K_mat <- matrix(0, n_locations, n_locations)
    for (i in 1:n_locations)
      for (j in 1:n_locations) {
        d2 <- (x_km[i] - x_km[j])^2 + (y_km[i] - y_km[j])^2
        K_mat[i, j] <- gp_alpha_factor_true[k]^2 * exp(-0.5 * d2 / rho_k[k]^2)
      }
    diag(K_mat) <- diag(K_mat) + 1e-9
    L_k         <- t(chol(K_mat))
    field_k     <- as.vector(L_k %*% rnorm(n_locations))
    latent_fields[, k] <- field_k - mean(field_k)
  }

  # True log-concentration: mu + loading * latent_fields + nugget
  true_log_concentration <- matrix(0, n_locations, n_species)
  nugget_effects         <- matrix(0, n_locations, n_species)
  for (sp in 1:n_species) {
    nugget_effects[, sp]         <- rnorm(n_locations, 0, sigma_nugget_sp[sp])
    true_log_concentration[, sp] <- mu_species_true[sp] +
      as.vector(latent_fields %*% loading_true[sp, ]) +
      nugget_effects[, sp]
  }

  true_lambda <- exp(true_log_concentration)
  lambda_K    <- true_lambda / (1 - exp(-true_lambda))
  lambda_K[true_lambda < 1e-10] <- 1.0

  total_lambda <- rowSums(true_lambda)
  mean_alpha   <- mean(alpha_known)
  kappa        <- alpha_known - mean_alpha

  # ZTP proportions: pi[loc, sp] = softmax(log(lambda_K) + N_PCR * alpha)
  alpha_pcr <- N_PCR * alpha_known
  true_pi   <- matrix(0, n_locations, n_species)
  for (loc in 1:n_locations) {
    log_nu         <- log(lambda_K[loc, ]) + alpha_pcr
    log_nu         <- log_nu - max(log_nu)
    true_pi[loc, ] <- exp(log_nu) / sum(exp(log_nu))
  }

  # Beta-Binomial dispersion: log_phi = beta0 + log(Λ) + gamma0*exp(-gamma1*log(λ_K)) + beta1*kappa
  true_phi <- matrix(0, n_locations, n_species)
  for (loc in 1:n_locations) {
    log_total <- log(max(total_lambda[loc], 1e-10))
    for (sp in 1:n_species) {
      log_phi_val       <- beta0_true + log_total +
        gamma0_true * exp(-gamma1_true * log(max(lambda_K[loc, sp], 1.0))) +
        beta1_true * kappa[sp]
      true_phi[loc, sp] <- exp(log_phi_val)
    }
  }

  true_concentration  <- exp(true_log_concentration)
  total_concentration <- rowSums(true_concentration)
  true_proportions    <- true_concentration / total_concentration

  # eDNA: Poisson copies -> Beta-Binomial reads
  edna_data_list <- list()
  for (loc_i in 1:n_locations) {
    pcr_id <- paste0("LOC", sprintf("%02d", loc_i), "_B1_PCR1")
    copies  <- rpois(n_species, true_lambda[loc_i, ])
    present <- copies > 0
    if (any(present)) {
      pi_draw <- numeric(n_species)
      for (sp in 1:n_species) {
        if (present[sp]) {
          bb_a        <- max(true_pi[loc_i, sp] * true_phi[loc_i, sp], 1e-3)
          bb_b        <- max((1 - true_pi[loc_i, sp]) * true_phi[loc_i, sp], 1e-3)
          pi_draw[sp] <- rbeta(1, bb_a, bb_b)
        } else {
          pi_draw[sp] <- 0
        }
      }
      sum_pi  <- sum(pi_draw)
      pi_draw <- if (sum_pi > 0) pi_draw / sum_pi else {
        pi_draw[present] <- 1 / sum(present); pi_draw
      }
      reads <- as.vector(rmultinom(1, total_reads_per_rep, pi_draw))
    } else {
      reads <- rep(0L, n_species)
    }
    for (sp_i in 1:n_species) {
      edna_data_list[[length(edna_data_list) + 1]] <- data.frame(
        pcr_replicate_id = pcr_id,
        bottle_id        = paste0("LOC", sprintf("%02d", loc_i), "_B1"),
        location_id      = locations$location_id[loc_i],
        lat = locations$lat[loc_i], lon = locations$lon[loc_i],
        species = species_names[sp_i], Nreads = reads[sp_i],
        stringsAsFactors = FALSE
      )
    }
  }
  edna_data <- dplyr::bind_rows(edna_data_list)

  # qPCR: per-replicate hurdle model (Bernoulli detection + Normal log-conc)
  qpcr_data_list <- list()
  for (loc_i in 1:n_locations) {
    true_log_conc <- true_log_concentration[loc_i, reference_species]
    p_detect      <- plogis(qpcr_detection_intercept_true + qpcr_detection_slope_true * true_log_conc)
    for (qpcr_rep in 1:n_qpcr_per_location) {
      detected     <- rbinom(1, 1, p_detect)
      log_conc_obs <- if (detected == 1) rnorm(1, true_log_conc, sigma_qpcr_true) else NA_real_
      qpcr_data_list[[length(qpcr_data_list) + 1]] <- data.frame(
        qpcr_id       = paste0("qPCR_LOC", sprintf("%02d", loc_i), "_R", qpcr_rep),
        location_id   = locations$location_id[loc_i],
        lat = locations$lat[loc_i], lon = locations$lon[loc_i],
        species       = species_names[reference_species],
        detected      = as.integer(detected),
        qpcr_log_conc = log_conc_obs,
        stringsAsFactors = FALSE
      )
    }
  }
  qpcr_data <- dplyr::bind_rows(qpcr_data_list)

  total_obs    <- n_locations * n_species
  zero_obs     <- sum(edna_data$Nreads == 0)
  zero_loc_ids <- qpcr_data %>%
    dplyr::group_by(location_id) %>%
    dplyr::summarise(any_detected = any(detected == 1), .groups = "drop") %>%
    dplyr::filter(!any_detected) %>%
    dplyr::pull(location_id)

  cat("=== Simulated Data Summary (v27 factor-analytic) ===\n")
  cat(sprintf("Locations: %d | Species: %d | K_f: %d | Ref: %s\n",
              n_locations, n_species, K_f, species_names[reference_species]))
  cat(sprintf("log_lambda_mean: %.2f | mean lambda: %.3f | mean lambda_K: %.3f\n",
              log_lambda_mean, mean(true_lambda), mean(lambda_K)))
  cat(sprintf("Dispersion: beta0=%.1f, beta1=%.1f, gamma0=%.1f, gamma1=%.2f\n",
              beta0_true, beta1_true, gamma0_true, gamma1_true))
  cat(sprintf("eDNA zeros: %d / %d (%.1f%%) — Poisson sampling\n",
              zero_obs, total_obs, 100 * zero_obs / total_obs))
  cat(sprintf("qPCR detected: %d / %d reps (%.0f%%)\n",
              sum(qpcr_data$detected), nrow(qpcr_data), 100 * mean(qpcr_data$detected)))
  if (length(zero_loc_ids) > 0)
    cat(sprintf("WARNING: %d location(s) with all-zero qPCR: %s\n",
                length(zero_loc_ids), paste(zero_loc_ids, collapse = ", ")))

  list(
    edna_data         = edna_data,
    qpcr_data         = qpcr_data,
    locations         = locations,
    species_names     = species_names,
    reference_species = reference_species,
    alpha_known       = alpha_known,
    true_values = list(
      true_concentration            = true_concentration,
      true_proportions              = true_proportions,
      true_log_concentration        = true_log_concentration,
      total_concentration           = total_concentration,
      true_lambda                   = true_lambda,
      lambda_K                      = lambda_K,
      true_phi                      = true_phi,
      true_pi                       = true_pi,
      alpha_known                   = alpha_known,
      psi_mean                      = mean(mu_species_true),
      mu_species                    = mu_species_true,
      species_offsets               = species_offsets,
      sigma_qpcr_true               = sigma_qpcr_true,
      qpcr_detection_intercept_true = qpcr_detection_intercept_true,
      qpcr_detection_slope_true     = qpcr_detection_slope_true,
      latent_fields                 = latent_fields,
      loading_true                  = loading_true,
      rho_factor_km                 = rho_k,
      gp_alpha_factor_true          = gp_alpha_factor_true,
      sigma_nugget                  = sigma_nugget_sp,
      nugget_effects                = nugget_effects,
      beta0_true                    = beta0_true,
      beta1_true                    = beta1_true,
      gamma0_true                   = gamma0_true,
      gamma1_true                   = gamma1_true
    )
  )
}


# =============================================================================
# SECTION 4: STAN DATA
# =============================================================================

#' Build the Stan data list for quant_metabar_qpcr_v27.stan
#'
#' Converts formatted data (from format_qpcr_metabar_data_v27) into the named
#' list required by the Stan model.  Key v27 parameters:
#'   K_f              — number of shared latent GP fields
#'   gp_alpha_factor  — fixed GP amplitude per factor (scalar, converted to vector)
#'   rho_min          — hard lower bound on rho_factor (km)
#'   rho_factor_prior_mean/sd — truncated-normal prior on rho_factor in km
#'   loading_prior_sd — regularizing prior SD on loading elements
#'
#' @param formatted_data Output of format_qpcr_metabar_data_v27.
#' @param N_pcr          Total PCR cycles.
#' @param m_x,m_y        HSGP basis counts per spatial axis.
#' @param max_rho        Maximum rho (km); used to set HSGP domain half-width.
#' @param K_f            Number of shared latent fields.
#' @param gp_alpha_factor Scalar GP amplitude per factor (default 1.0).
#' @param rho_min        Hard lower bound on rho_factor (km).
#' @param rho_factor_prior_mean,rho_factor_prior_sd Truncated-normal prior on rho_factor (km).
#' @param loading_prior_sd Regularizing prior SD on loading elements.
#' @param log_conc_prior_mean Mean of log-concentration prior; auto-computed if NULL.
#' @param log_conc_prior_sd  SD of log-concentration prior.
#' @param sigma_qpcr_prior_sd  Half-normal prior SD for sigma_qpcr.
#' @param qpcr_detection_intercept_prior_sd  Normal prior SD for detection intercept.
#' @param qpcr_detection_slope_prior_mean,sd  Normal prior for detection slope.
#' @param sigma_nugget_prior_mean  Log-scale center of hierarchical nugget prior.
#' @param sigma_nugget_prior_sd    SD of hierarchical nugget prior; near-zero pins nugget.
#' @param sigma_nugget_sigma_rate  Exponential rate for nugget hyper-SD.
#' @param beta0_prior_sd,gamma0_prior_mean,gamma0_prior_sd,
#'   log_gamma1_prior_mean,log_gamma1_prior_sd Dispersion priors.
#' @return Named list for rstan::sampling().
make_stan_data_v27 <- function(
    formatted_data,
    N_pcr = 35,
    m_x = 10, m_y = 8,
    max_rho               = 300,
    K_f                   = 3,
    gp_alpha_factor       = 1.0,
    rho_min               = 20.0,
    rho_factor_prior_mean = 100.0,
    rho_factor_prior_sd   = 200.0,
    loading_prior_sd      = 1.0,
    log_conc_prior_mean   = NULL,
    log_conc_prior_sd     = 5.0,
    sigma_qpcr_prior_sd   = 1.0,
    qpcr_detection_intercept_prior_sd = 5.0,
    qpcr_detection_slope_prior_mean   = 2.0,
    qpcr_detection_slope_prior_sd     = 2.0,
    sigma_nugget_prior_mean  = log(0.1),
    sigma_nugget_prior_sd    = 1e-5,
    sigma_nugget_sigma_rate  = 1.0,
    beta0_prior_sd          = 1.0,
    gamma0_prior_mean       = 2.0,
    gamma0_prior_sd         = 1.0,
    log_gamma1_prior_mean   = 1.0,
    log_gamma1_prior_sd     = 0.75
) {
  if (K_f > formatted_data$N_species)
    stop(sprintf("K_f (%d) must be <= N_species (%d)", K_f, formatted_data$N_species))

  M <- m_x * m_y
  if (M >= formatted_data$N_locations)
    stop(sprintf("m_x*m_y (%d) must be < N_locations (%d). Reduce m_x or m_y.",
                 M, formatted_data$N_locations))
  if (M > formatted_data$N_locations / 2)
    warning(sprintf("m_x*m_y (%d) > N_locations/2 (%d). Consider reducing basis size.",
                    M, as.integer(formatted_data$N_locations / 2)))

  if (is.null(log_conc_prior_mean)) {
    detected_concs      <- formatted_data$qpcr_rep_log_conc[formatted_data$qpcr_rep_detected == 1]
    log_conc_prior_mean <- if (length(detected_concs) > 0) mean(detected_concs) else 2.0
  }

  max_abs_x  <- max(abs(formatted_data$x_utm))
  max_abs_y  <- max(abs(formatted_data$y_utm))
  L_target   <- max_rho * pi / 4
  c_x        <- max(1.0, (L_target - 1) / max_abs_x)
  c_y        <- max(1.0, (L_target - 1) / max_abs_y)

  if ((L_target - 1) / max_abs_x < 1.0)
    warning(sprintf("max_rho=%.0f km smaller than x extent; c_x clamped to 1.0.", max_rho))
  if ((L_target - 1) / max_abs_y < 1.0)
    warning(sprintf("max_rho=%.0f km smaller than y extent; c_y clamped to 1.0.", max_rho))

  n_zero_cells  <- sum(formatted_data$edna_obs_data == 0)
  n_total_cells <- formatted_data$N_edna_obs * formatted_data$N_species
  cat(sprintf("log_conc_prior_mean=%.2f | HSGP: c_x=%.2f, c_y=%.2f | K_f=%d\n",
              log_conc_prior_mean, c_x, c_y, K_f))
  cat(sprintf("eDNA zero cells: %d / %d (%.1f%%)\n",
              n_zero_cells, n_total_cells, 100 * n_zero_cells / n_total_cells))

  list(
    N_species          = formatted_data$N_species,
    N_locations        = formatted_data$N_locations,
    N_edna_obs         = formatted_data$N_edna_obs,
    N_qpcr_reps        = formatted_data$N_qpcr_reps,
    reference_species  = formatted_data$reference_species,
    m_x = as.integer(m_x), m_y = as.integer(m_y),
    c_x = c_x, c_y = c_y,
    x_utm              = formatted_data$x_utm,
    y_utm              = formatted_data$y_utm,
    edna_obs_data      = formatted_data$edna_obs_data,
    edna_obs_location  = formatted_data$edna_obs_location,
    qpcr_rep_location  = formatted_data$qpcr_rep_location,
    qpcr_rep_detected  = formatted_data$qpcr_rep_detected,
    qpcr_rep_log_conc  = formatted_data$qpcr_rep_log_conc,
    alpha_known        = formatted_data$alpha_known,
    N_PCR              = N_pcr,
    has_edna           = formatted_data$has_edna,
    has_qpcr           = formatted_data$has_qpcr,
    log_conc_prior_mean               = log_conc_prior_mean,
    log_conc_prior_sd                 = log_conc_prior_sd,
    K_f                               = as.integer(K_f),
    gp_alpha_factor                   = rep(gp_alpha_factor, K_f),
    rho_min                           = rho_min,
    rho_factor_prior_mean             = rho_factor_prior_mean,
    rho_factor_prior_sd               = rho_factor_prior_sd,
    loading_prior_sd                  = loading_prior_sd,
    sigma_qpcr_prior_sd               = sigma_qpcr_prior_sd,
    qpcr_detection_intercept_prior_sd = qpcr_detection_intercept_prior_sd,
    qpcr_detection_slope_prior_mean   = qpcr_detection_slope_prior_mean,
    qpcr_detection_slope_prior_sd     = qpcr_detection_slope_prior_sd,
    sigma_nugget_prior_mean           = sigma_nugget_prior_mean,
    sigma_nugget_prior_sd             = sigma_nugget_prior_sd,
    sigma_nugget_sigma_rate           = sigma_nugget_sigma_rate,
    beta0_prior_sd                    = beta0_prior_sd,
    gamma0_prior_mean                 = gamma0_prior_mean,
    gamma0_prior_sd                   = gamma0_prior_sd,
    log_gamma1_prior_mean             = log_gamma1_prior_mean,
    log_gamma1_prior_sd               = log_gamma1_prior_sd
  )
}


# =============================================================================
# SECTION 5: EXTRACTION
# =============================================================================

#' Extract qPCR observation-model parameters from a fitted Stan model
extract_qpcr_params_v15 <- function(fit, formatted_data = NULL) {
  posterior          <- rstan::extract(fit)
  psi_samples        <- posterior$psi
  sigma_qpcr_samples <- posterior$sigma_qpcr
  intercept_samples  <- posterior$qpcr_detection_intercept
  slope_samples      <- posterior$qpcr_detection_slope
  p_detect_samples   <- posterior$p_detect_loc

  p_detect_summary <- data.frame(
    location_idx  = 1:ncol(p_detect_samples),
    p_detect_mean = colMeans(p_detect_samples),
    p_detect_q025 = apply(p_detect_samples, 2, quantile, 0.025),
    p_detect_q975 = apply(p_detect_samples, 2, quantile, 0.975)
  )
  if (!is.null(formatted_data))
    p_detect_summary <- p_detect_summary %>%
      left_join(formatted_data$location_list, by = "location_idx")

  result <- list(
    psi = list(mean = mean(psi_samples), median = median(psi_samples),
               sd = sd(psi_samples), q025 = quantile(psi_samples, 0.025),
               q975 = quantile(psi_samples, 0.975)),
    sigma_qpcr = list(mean = mean(sigma_qpcr_samples), median = median(sigma_qpcr_samples),
                      sd = sd(sigma_qpcr_samples), q025 = quantile(sigma_qpcr_samples, 0.025),
                      q975 = quantile(sigma_qpcr_samples, 0.975)),
    qpcr_detection_intercept = list(mean = mean(intercept_samples), sd = sd(intercept_samples),
                                    q025 = quantile(intercept_samples, 0.025),
                                    q975 = quantile(intercept_samples, 0.975)),
    qpcr_detection_slope     = list(mean = mean(slope_samples), sd = sd(slope_samples),
                                    q025 = quantile(slope_samples, 0.025),
                                    q975 = quantile(slope_samples, 0.975)),
    p_detect_loc = p_detect_summary
  )
  if ("mu_species" %in% names(posterior)) {
    mu_species_samples <- posterior$mu_species
    sp_names           <- if (!is.null(formatted_data)) formatted_data$sp_list$species else
      paste0("sp_", 1:ncol(mu_species_samples))
    result$mu_species <- purrr::map_dfr(seq_len(ncol(mu_species_samples)), function(i) {
      vals <- mu_species_samples[, i]
      data.frame(species = sp_names[i], mean = mean(vals), median = median(vals),
                 sd = sd(vals), q025 = quantile(vals, 0.025), q975 = quantile(vals, 0.975))
    })
  }
  result
}

#' Extract nugget variance parameters (v21/v22 parameterization)
extract_nugget_params_v21 <- function(fit, formatted_data) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species
  N_species     <- length(species_names)

  species_summary <- purrr::map_dfr(seq_len(N_species), function(i) {
    vals <- posterior$sigma_nugget[, i]
    data.frame(
      species = species_names[i],
      mean    = mean(vals),
      median  = median(vals),
      sd      = sd(vals),
      q025    = quantile(vals, 0.025),
      q975    = quantile(vals, 0.975)
    )
  })

  var_decomp <- purrr::map_dfr(seq_len(N_species), function(i) {
    data.frame(
      species         = species_names[i],
      prop_var_gp     = mean(posterior$prop_var_gp[, i]),
      prop_var_nugget = mean(posterior$prop_var_nugget[, i])
    )
  })

  list(
    species_summary     = species_summary,
    var_decomp          = var_decomp,
    mu_log_sigma_nugget = list(
      mean   = mean(posterior$mu_log_sigma_nugget),
      median = median(posterior$mu_log_sigma_nugget),
      q025   = quantile(posterior$mu_log_sigma_nugget, 0.025),
      q975   = quantile(posterior$mu_log_sigma_nugget, 0.975)
    ),
    sigma_sigma_nugget  = list(
      mean   = mean(posterior$sigma_sigma_nugget),
      median = median(posterior$sigma_sigma_nugget),
      q025   = quantile(posterior$sigma_sigma_nugget, 0.025),
      q975   = quantile(posterior$sigma_sigma_nugget, 0.975)
    ),
    sigma_nugget_pop_mean   = list(
      mean = mean(posterior$sigma_nugget_pop_mean),
      q025 = quantile(posterior$sigma_nugget_pop_mean, 0.025),
      q975 = quantile(posterior$sigma_nugget_pop_mean, 0.975)
    ),
    sigma_nugget_pop_median = list(
      mean = mean(posterior$sigma_nugget_pop_median)
    ),
    sigma_nugget_mean       = list(
      mean = mean(posterior$sigma_nugget_mean),
      q025 = quantile(posterior$sigma_nugget_mean, 0.025),
      q975 = quantile(posterior$sigma_nugget_mean, 0.975)
    )
  )
}

#' Extract log-concentration posterior summaries per species
extract_concentration_summaries <- function(fit, formatted_data) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species
  N_species     <- length(species_names)
  log_conc      <- posterior$log_concentration
  species_summary <- purrr::map_dfr(seq_len(N_species), function(i) {
    sp_samples   <- log_conc[, , i]
    empirical_sd <- apply(sp_samples, 1, sd)
    sp_mean      <- apply(sp_samples, 1, mean)
    data.frame(
      species             = species_names[i],
      mean_log_conc_mean  = mean(sp_mean),
      mean_log_conc_q025  = quantile(sp_mean, 0.025),
      mean_log_conc_q975  = quantile(sp_mean, 0.975),
      empirical_sd_mean   = mean(empirical_sd),
      empirical_sd_median = median(empirical_sd),
      empirical_sd_q025   = quantile(empirical_sd, 0.025),
      empirical_sd_q975   = quantile(empirical_sd, 0.975)
    )
  })
  list(species_summary = species_summary)
}

#' Extract factor-analytic GP parameters (v27)
#'
#' Returns rho_factor, loading_diag, and full loading matrix summaries.
extract_factor_params_v27 <- function(fit, formatted_data) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species
  N_species     <- length(species_names)
  K_f           <- ncol(posterior$rho_factor)  # dim: [n_samp, K_f]

  rho_factor_summary <- purrr::map_dfr(seq_len(K_f), function(k) {
    vals <- posterior$rho_factor[, k]
    data.frame(factor = k, mean = mean(vals), median = median(vals), sd = sd(vals),
               q025 = quantile(vals, 0.025), q975 = quantile(vals, 0.975),
               eff_range = 2.15 * mean(vals))
  })

  loading_diag_summary <- purrr::map_dfr(seq_len(K_f), function(k) {
    vals <- posterior$loading_diag[, k]
    data.frame(factor = k, mean = mean(vals), median = median(vals), sd = sd(vals),
               q025 = quantile(vals, 0.025), q975 = quantile(vals, 0.975))
  })

  # loading is [n_samples, N_species, K_f]
  loading_post_mean <- apply(posterior$loading, c(2, 3), mean)
  loading_post_sd   <- apply(posterior$loading, c(2, 3), sd)
  loading_df <- expand.grid(species = species_names, factor = seq_len(K_f)) %>%
    mutate(post_mean = as.vector(loading_post_mean),
           post_sd   = as.vector(loading_post_sd))

  list(
    rho_factor_summary   = rho_factor_summary,
    loading_diag_summary = loading_diag_summary,
    loading_post_mean    = loading_post_mean,
    loading_df           = loading_df,
    rho_mean = list(mean = mean(posterior$rho_mean),
                    q025 = quantile(posterior$rho_mean, 0.025),
                    q975 = quantile(posterior$rho_mean, 0.975))
  )
}


# =============================================================================
# SECTION 6: VISUALIZATION
# =============================================================================

#' Plot observed vs posterior predicted qPCR log-concentrations
plot_qpcr_fit_v15 <- function(fit, formatted_data) {
  posterior        <- rstan::extract(fit)
  qpcr_params      <- extract_qpcr_params_v15(fit, formatted_data)
  expected_samples <- posterior$expected_qpcr_loc
  pred_df <- data.frame(
    location_idx = 1:ncol(expected_samples),
    pred_mean    = colMeans(expected_samples),
    pred_q025    = apply(expected_samples, 2, quantile, 0.025),
    pred_q975    = apply(expected_samples, 2, quantile, 0.975)
  )
  rep_df <- data.frame(
    rep_idx      = seq_len(formatted_data$N_qpcr_reps),
    location_idx = formatted_data$qpcr_rep_location,
    detected     = formatted_data$qpcr_rep_detected,
    obs_log_conc = formatted_data$qpcr_rep_log_conc
  ) %>% left_join(pred_df, by = "location_idx")

  rep_detected     <- rep_df %>% filter(detected == 1)
  rep_not_detected <- rep_df %>% filter(detected == 0)
  subtitle_text    <- sprintf(
    "sigma_qpcr = %.2f (%.2f-%.2f) | intercept = %.2f | slope = %.2f",
    qpcr_params$sigma_qpcr$mean, qpcr_params$sigma_qpcr$q025, qpcr_params$sigma_qpcr$q975,
    qpcr_params$qpcr_detection_intercept$mean, qpcr_params$qpcr_detection_slope$mean)

  p <- ggplot(rep_detected, aes(x = obs_log_conc, y = pred_mean)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975), width = 0, alpha = 0.3, color = "darkred") +
    geom_point(size = 2, color = "darkred", alpha = 0.8) +
    labs(title    = "qPCR: Observed vs. Predicted (detected replicates)",
         subtitle = subtitle_text,
         x = "Observed log(copies/uL)", y = "Posterior expected log(copies/uL)") +
    theme_light_custom()
  if (nrow(rep_not_detected) > 0)
    p <- p +
      geom_point(data = rep_not_detected, aes(x = pred_mean, y = pred_mean),
                 shape = 4, size = 2, color = "gray50", alpha = 0.5) +
      labs(caption = sprintf("X marks: %d non-detected reps (plotted at pred_mean diagonal)",
                             nrow(rep_not_detected)))
  p
}

#' Plot qPCR detection probability curve vs posterior log-concentration
plot_detection_curve_v15 <- function(fit, formatted_data) {
  posterior            <- rstan::extract(fit)
  ref_sp               <- formatted_data$reference_species
  log_conc_ref_mean    <- colMeans(posterior$log_concentration[, , ref_sp])
  p_detect_samples     <- posterior$p_detect_loc
  p_detect_mean        <- colMeans(p_detect_samples)
  p_detect_q025        <- apply(p_detect_samples, 2, quantile, 0.025)
  p_detect_q975        <- apply(p_detect_samples, 2, quantile, 0.975)

  obs_detect_df <- data.frame(
    location_idx = formatted_data$qpcr_rep_location,
    detected     = formatted_data$qpcr_rep_detected
  ) %>%
    group_by(location_idx) %>%
    summarise(n_reps = n(), n_detected = sum(detected),
              obs_prop = n_detected / n_reps, .groups = "drop")

  detect_curve_df <- data.frame(
    location_idx      = 1:formatted_data$N_locations,
    log_conc_ref_mean = log_conc_ref_mean,
    p_detect_mean     = p_detect_mean,
    p_detect_q025     = p_detect_q025,
    p_detect_q975     = p_detect_q975
  ) %>% left_join(obs_detect_df, by = "location_idx")

  qpcr_params   <- extract_qpcr_params_v15(fit, formatted_data)
  subtitle_text <- sprintf(
    "intercept = %.2f (%.2f-%.2f) | slope = %.2f (%.2f-%.2f)",
    qpcr_params$qpcr_detection_intercept$mean,
    qpcr_params$qpcr_detection_intercept$q025, qpcr_params$qpcr_detection_intercept$q975,
    qpcr_params$qpcr_detection_slope$mean,
    qpcr_params$qpcr_detection_slope$q025, qpcr_params$qpcr_detection_slope$q975)

  ggplot(detect_curve_df, aes(x = log_conc_ref_mean)) +
    geom_ribbon(aes(ymin = p_detect_q025, ymax = p_detect_q975), alpha = 0.2, fill = "steelblue") +
    geom_line(aes(y = p_detect_mean), color = "steelblue4", linewidth = 1) +
    geom_point(aes(y = obs_prop, size = n_reps), color = "darkred", alpha = 0.7,
               position = position_jitter(height = 0.02, width = 0)) +
    scale_size_continuous(name = "# reps", range = c(1.5, 5)) +
    ylim(0, 1) +
    labs(title    = "qPCR Detection Probability vs Log Concentration",
         subtitle = subtitle_text,
         x = "Posterior Mean log(Concentration) — Reference Species",
         y = "Detection Probability") +
    theme_light_custom()
}

#' Plot observed vs posterior eDNA proportions per species
plot_edna_fit <- function(fit, formatted_data) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species
  N_species     <- length(species_names)
  plots         <- list()
  for (sp_idx in seq_len(N_species)) {
    sp <- species_names[sp_idx]
    p  <- tryCatch({
      obs_edna <- data.frame(
        edna_obs_idx = 1:formatted_data$N_edna_obs,
        location_idx = formatted_data$edna_obs_location,
        obs_reads    = formatted_data$edna_obs_data[, sp_idx],
        total_reads  = rowSums(formatted_data$edna_obs_data)
      ) %>%
        mutate(obs_prop = obs_reads / total_reads) %>%
        group_by(location_idx) %>%
        summarise(obs_prop_mean = mean(obs_prop), obs_prop_sd = sd(obs_prop),
                  n_reps = n(), .groups = "drop")

      mu_samples <- posterior$mu_metabar_loc[, , sp_idx]
      pred_edna  <- data.frame(
        location_idx = 1:ncol(mu_samples),
        pred_mean    = colMeans(mu_samples),
        pred_q025    = apply(mu_samples, 2, quantile, 0.025),
        pred_q975    = apply(mu_samples, 2, quantile, 0.975)
      )
      edna_compare <- obs_edna %>% left_join(pred_edna, by = "location_idx")

      ggplot(edna_compare, aes(x = obs_prop_mean, y = pred_mean)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
        geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975),
                      width = 0, alpha = 0.3, color = "steelblue") +
        geom_point(size = 2, color = "steelblue") +
        labs(title    = paste("eDNA Proportions:", gsub("_", " ", sp)),
             subtitle = "Observed vs. Posterior (mean +/- 95% CI)",
             x = "Observed Proportion", y = "Posterior Proportion") +
        coord_equal() + theme_minimal() + theme(plot.title = element_text(face = "bold"))
    }, error = function(e) { warning(paste("eDNA plot failed for", sp, ":", e$message)); NULL })
    if (!is.null(p)) plots[[paste0(sp, "_edna")]] <- p
  }
  plots
}

#' Plot variance decomposition: proportion of log-conc variance from GP vs nugget (v22)
plot_variance_decomposition_v22 <- function(fit, formatted_data) {
  nugget_params <- extract_nugget_params_v21(fit, formatted_data)
  var_decomp    <- nugget_params$var_decomp

  vd_long <- dplyr::bind_rows(
    data.frame(species = var_decomp$species, component = "GP",     mean = var_decomp$prop_var_gp),
    data.frame(species = var_decomp$species, component = "Nugget", mean = var_decomp$prop_var_nugget)
  )
  vd_long$component <- factor(vd_long$component, levels = c("Nugget", "GP"))

  ggplot(vd_long, aes(x = species, y = mean, fill = component)) +
    geom_col(position = "stack", alpha = 0.85) +
    scale_fill_manual(values = c("GP" = "steelblue3", "Nugget" = "coral2"),
                      name = "Variance\nComponent") +
    labs(title    = "Variance Decomposition by Species",
         subtitle = "Proportion of log-concentration variance: GP + Nugget",
         x = "Species", y = "Proportion of Variance") +
    coord_flip() + theme_light_custom()
}

#' Plot true vs posterior-predicted concentrations (simulation validation)
plot_predicted_vs_true_concentration <- function(fit, formatted_data, true_concentration) {
  posterior     <- rstan::extract(fit)
  species_names <- formatted_data$sp_list$species
  N_species     <- length(species_names)
  conc_samples  <- posterior$concentration
  plots         <- list()
  for (sp_idx in seq_len(N_species)) {
    sp <- species_names[sp_idx]
    p  <- tryCatch({
      sp_conc    <- conc_samples[, , sp_idx]
      compare_df <- data.frame(
        location_idx = 1:ncol(sp_conc),
        true_conc    = true_concentration[, sp_idx],
        pred_mean    = colMeans(sp_conc),
        pred_q025    = apply(sp_conc, 2, quantile, 0.025),
        pred_q975    = apply(sp_conc, 2, quantile, 0.975)
      )
      compare_df$covered <- compare_df$true_conc >= compare_df$pred_q025 &
                            compare_df$true_conc <= compare_df$pred_q975
      coverage <- mean(compare_df$covered) * 100
      if (nrow(compare_df) != formatted_data$N_locations)
        warning(sprintf("true_concentration rows (%d) != N_locations (%d); use match_true_to_formatted()",
                        nrow(compare_df), formatted_data$N_locations))
      cor_val  <- tryCatch(
        cor(log(pmax(compare_df$true_conc, 1e-10)), log(pmax(compare_df$pred_mean, 1e-10))),
        error = function(e) NA_real_)
      ggplot(compare_df, aes(x = true_conc, y = pred_mean)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
        geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975),
                      width = 0, alpha = 0.4, color = "darkorange3") +
        geom_point(aes(color = covered), size = 2.5, shape = 16) +
        scale_color_manual(values = c("FALSE" = "red", "TRUE" = "darkorange"), name = "In 95% CI") +
        scale_x_log10() + scale_y_log10() +
        labs(title    = paste("True vs Predicted:", gsub("_", " ", sp)),
             subtitle = sprintf("Coverage: %.0f%% | Correlation (log): r = %.2f", coverage, cor_val),
             x = "True Concentration (copies/uL)", y = "Posterior Concentration (copies/uL)") +
        theme_light_custom() + theme(legend.position = "bottom")
    }, error = function(e) { warning(paste("True vs predicted failed for", sp, ":", e$message)); NULL })
    if (!is.null(p)) plots[[paste0(sp, "_true_vs_pred")]] <- p
  }
  plots
}

#' Plot K_f latent fields as spatial heatmaps (v27)
#'
#' One panel per factor showing posterior mean latent field values.
#' If true_latent_fields provided, overlays them as points.
plot_latent_fields_v27 <- function(fit, formatted_data, true_latent_fields = NULL) {
  posterior  <- rstan::extract(fit)
  K_f        <- dim(posterior$latent_field)[3]
  locations  <- formatted_data$location_list
  plot_list  <- list()

  for (k in seq_len(K_f)) {
    p <- tryCatch({
      field_samples <- posterior$latent_field[, , k]  # [n_samp, N_locations]
      field_mean    <- colMeans(field_samples)
      loc_df        <- data.frame(
        x_utm     = locations$x_utm_centered,
        y_utm     = locations$y_utm_centered,
        field_val = field_mean
      )
      interp_df <- interpolate_surface(x = loc_df$x_utm, y = loc_df$y_utm,
                                       z = loc_df$field_val, resolution = 80,
                                       log_transform = FALSE)
      p_base <- ggplot2::ggplot() +
        ggplot2::geom_raster(data = interp_df, ggplot2::aes(x = x, y = y, fill = value)) +
        ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "coral2",
                                      midpoint = 0, name = "Latent\nField") +
        ggplot2::labs(title = sprintf("Latent Field %d", k),
                      x = "UTM Easting (km)", y = "UTM Northing (km)") +
        ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_minimal()
      if (!is.null(true_latent_fields)) {
        true_df <- data.frame(x_utm = loc_df$x_utm, y_utm = loc_df$y_utm,
                              true_val = true_latent_fields[, k])
        p_base <- p_base +
          ggplot2::geom_point(data = true_df,
                              ggplot2::aes(x = x_utm, y = y_utm, fill = true_val),
                              shape = 21, size = 2.5, color = "white", stroke = 0.5)
      } else {
        p_base <- p_base +
          ggplot2::geom_point(data = loc_df, ggplot2::aes(x = x_utm, y = y_utm),
                              shape = 21, fill = "gray50", color = "white", size = 1.5)
      }
      p_base
    }, error = function(e) { warning(paste("Latent field plot failed for factor", k, ":", e$message)); NULL })
    if (!is.null(p)) plot_list[[paste0("factor_", k)]] <- p
  }
  combine_plots(plot_list, "Posterior Mean Latent Fields (v27)")
}

#' Plot the loading matrix as a species x factor heatmap (v27)
#'
#' Shows posterior mean loading values. Blue = negative, white = zero, red = positive.
#' If true_loading provided, overlays true values as text labels.
plot_loading_matrix_v27 <- function(fit, formatted_data, true_loading = NULL) {
  factor_params <- extract_factor_params_v27(fit, formatted_data)
  loading_df    <- factor_params$loading_df
  loading_df$factor_label <- paste0("Factor ", loading_df$factor)

  p <- ggplot2::ggplot(loading_df, ggplot2::aes(x = factor_label, y = species, fill = post_mean)) +
    ggplot2::geom_tile(color = "gray70") +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "coral2",
                                  midpoint = 0, name = "Loading\n(post. mean)") +
    ggplot2::labs(title    = "Loading Matrix (v27)",
                  subtitle = "Posterior mean; lower-triangular; positive diagonal",
                  x = "Factor", y = "Species") +
    theme_light_custom()

  if (!is.null(true_loading)) {
    species_names <- formatted_data$sp_list$species
    K_f           <- ncol(true_loading)
    true_df <- expand.grid(species = species_names, factor = seq_len(K_f)) %>%
      mutate(factor_label = paste0("Factor ", factor),
             true_val     = as.vector(true_loading),
             label        = sprintf("%.2f", true_val))
    p <- p + ggplot2::geom_text(data = true_df,
                                ggplot2::aes(label = label), size = 3, color = "black",
                                inherit.aes = FALSE)
  }
  p
}

#' Plot posterior densities of rho_factor per factor (v27)
plot_rho_factor_posteriors_v27 <- function(fit, formatted_data, true_rho_factor = NULL) {
  posterior <- rstan::extract(fit)
  K_f       <- ncol(posterior$rho_factor)

  plot_list <- lapply(seq_len(K_f), function(k) {
    df <- data.frame(rho = posterior$rho_factor[, k])
    p  <- ggplot2::ggplot(df, ggplot2::aes(x = rho)) +
      ggplot2::geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
      ggplot2::labs(title = paste0("Factor ", k), x = "rho_factor (km)", y = "Count") +
      ggplot2::theme_minimal()
    if (!is.null(true_rho_factor)) {
      true_v <- if (length(true_rho_factor) == 1) rep(true_rho_factor, K_f) else true_rho_factor
      p <- p + ggplot2::geom_vline(xintercept = true_v[k], color = "red", linewidth = 1)
    }
    p
  })
  patchwork::wrap_plots(plot_list, nrow = 1) +
    patchwork::plot_annotation(title = "Posterior: rho_factor per shared latent field (v27)")
}


# =============================================================================
# SECTION 7: FIT AND PLOT
# =============================================================================

#' Compile and sample the v27 Stan model, then extract results and generate plots
#'
#' @param formatted_data  Output of format_qpcr_metabar_data_v27.
#' @param stan_file       Path to quant_metabar_qpcr_v27.stan.
#' @param stan_data       Pre-built Stan data list; built from formatted_data if NULL.
#' @param chains,iter,warmup,seed  Sampling settings.
#' @param m_x,m_y        HSGP basis function counts per axis.
#' @param max_rho        Maximum rho (km) for HSGP domain.
#' @param K_f            Number of shared latent GP fields.
#' @param gp_alpha_factor Fixed GP amplitude per factor.
#' @param rho_min        Hard lower bound on rho_factor (km).
#' @param rho_factor_prior_mean,rho_factor_prior_sd  Truncated-normal prior on rho_factor.
#' @param loading_prior_sd Regularizing prior on loading elements.
#' @param sigma_qpcr_prior_sd,qpcr_detection_intercept_prior_sd,
#'   qpcr_detection_slope_prior_mean,qpcr_detection_slope_prior_sd  qPCR priors.
#' @param sigma_nugget_prior_mean,sigma_nugget_prior_sd,sigma_nugget_sigma_rate  Nugget priors.
#' @param log_conc_prior_sd  SD of log-concentration prior.
#' @param beta0_prior_sd,gamma0_prior_mean,gamma0_prior_sd,
#'   log_gamma1_prior_mean,log_gamma1_prior_sd  Dispersion priors.
#' @param adapt_delta,max_treedepth  HMC tuning parameters.
#' @param true_values    List of true simulation parameters for validation; NULL for real data.
#' @param save_outputs   Write plots and diagnostics to output_dir.
#' @param output_dir     Directory for saved outputs.
#' @return Named list containing fit, stan_data, formatted_data, parameter summaries,
#'   diagnostic statistics, and ggplot objects.
fit_and_plot_v27 <- function(
    formatted_data,
    stan_file  = "code/spatial_qPCR_mb/quant_metabar_qpcr_v27.stan",
    stan_data  = NULL,
    chains     = 3,
    iter       = 2000,
    warmup     = 1000,
    seed       = 42,
    m_x = 10, m_y = 8,
    max_rho               = 300,
    K_f                   = 3,
    gp_alpha_factor       = 1.0,
    rho_min               = 20.0,
    rho_factor_prior_mean = 100.0,
    rho_factor_prior_sd   = 200.0,
    loading_prior_sd      = 1.0,
    sigma_qpcr_prior_sd   = 1.0,
    qpcr_detection_intercept_prior_sd = 5.0,
    qpcr_detection_slope_prior_mean   = 2.0,
    qpcr_detection_slope_prior_sd     = 2.0,
    sigma_nugget_prior_mean  = log(0.1),
    sigma_nugget_prior_sd    = 1e-5,
    sigma_nugget_sigma_rate  = 2.0,
    log_conc_prior_sd        = 5.0,
    beta0_prior_sd           = 1.0,
    gamma0_prior_mean        = 2.0,
    gamma0_prior_sd          = 1.0,
    log_gamma1_prior_mean    = 1.0,
    log_gamma1_prior_sd      = 0.75,
    adapt_delta              = 0.99,
    max_treedepth            = 14,
    true_values              = NULL,
    save_outputs             = TRUE,
    output_dir               = "code/spatial_qPCR_mb/outputs_v27"
) {
  if (save_outputs && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (is.null(stan_data)) {
    cat("\n=== Preparing Stan data (v27) ===\n")
    stan_data <- make_stan_data_v27(
      formatted_data,
      m_x = m_x, m_y = m_y,
      max_rho = max_rho,
      K_f = K_f,
      gp_alpha_factor = gp_alpha_factor,
      rho_min = rho_min,
      rho_factor_prior_mean = rho_factor_prior_mean,
      rho_factor_prior_sd   = rho_factor_prior_sd,
      loading_prior_sd      = loading_prior_sd,
      log_conc_prior_sd     = log_conc_prior_sd,
      sigma_qpcr_prior_sd   = sigma_qpcr_prior_sd,
      qpcr_detection_intercept_prior_sd = qpcr_detection_intercept_prior_sd,
      qpcr_detection_slope_prior_mean   = qpcr_detection_slope_prior_mean,
      qpcr_detection_slope_prior_sd     = qpcr_detection_slope_prior_sd,
      sigma_nugget_prior_mean = sigma_nugget_prior_mean,
      sigma_nugget_prior_sd   = sigma_nugget_prior_sd,
      sigma_nugget_sigma_rate = sigma_nugget_sigma_rate,
      beta0_prior_sd          = beta0_prior_sd,
      gamma0_prior_mean       = gamma0_prior_mean,
      gamma0_prior_sd         = gamma0_prior_sd,
      log_gamma1_prior_mean   = log_gamma1_prior_mean,
      log_gamma1_prior_sd     = log_gamma1_prior_sd
    )
  }

  cat("\n=== Compiling Stan model (v27) ===\n")
  stan_mod <- rstan::stan_model(stan_file)

  cat("\n=== Sampling (v27) ===\n")
  cat(sprintf("adapt_delta=%.3f | max_treedepth=%d\n", adapt_delta, max_treedepth))
  fit <- tryCatch(
    rstan::sampling(
      stan_mod, data = stan_data,
      chains = chains, iter = iter, warmup = warmup, seed = seed,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
    ),
    error = function(e) {
      cat("\nSampling ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(fit)) stop("rstan::sampling() failed — see error above.")

  n_samples <- tryCatch(sum(sapply(fit@sim$samples, length)), error = function(e) 0L)
  if (n_samples == 0L)
    stop("Stan sampler returned 0 samples. Check stderr (SLURM error log) for chain-level errors.")

  cat("\n=== Diagnostics ===\n")
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  n_div  <- tryCatch(sum(sapply(sampler_params, function(x) if (is.matrix(x)) sum(x[, "divergent__"])  else 0)), error = function(e) NA)
  max_td <- tryCatch(sum(sapply(sampler_params, function(x) if (is.matrix(x)) sum(x[, "treedepth__"] >= max_treedepth) else 0)), error = function(e) NA)
  fit_summary <- summary(fit)$summary
  max_rhat    <- if (!is.null(fit_summary) && "Rhat" %in% colnames(fit_summary))
                   max(fit_summary[, "Rhat"], na.rm = TRUE) else NA_real_
  n_high_rhat <- if (!is.na(max_rhat))
                   sum(fit_summary[, "Rhat"] > 1.01, na.rm = TRUE) else 0L
  cat(sprintf("Divergences: %s | Max treedepth hits: %s\n",
              ifelse(is.na(n_div), "NA", n_div), ifelse(is.na(max_td), "NA", max_td)))
  cat(sprintf("Max Rhat: %.3f | # Rhat > 1.01: %d\n", max_rhat, n_high_rhat))
  top_rhat <- get_top_rhat(fit, n = 10)

  cat("\n=== Parameter Summaries ===\n")
  print(fit, pars = c("mu_species"))
  print(fit, pars = c("psi", "species_offset", "sigma_qpcr"))
  print(fit, pars = c("qpcr_detection_intercept", "qpcr_detection_slope"))
  print(fit, pars = c("beta0", "beta1", "gamma0", "log_gamma1", "gamma1", "mean_phi"))
  print(fit, pars = c("rho_factor", "rho_mean", "loading_diag"))
  print(fit, pars = c("mu_log_sigma_nugget", "sigma_sigma_nugget",
                      "sigma_nugget", "sigma_nugget_pop_mean", "sigma_nugget_pop_median",
                      "sigma_nugget_mean"))
  print(fit, pars = c("prop_var_gp", "prop_var_nugget"))

  species_names  <- formatted_data$sp_list$species
  qpcr_params    <- extract_qpcr_params_v15(fit, formatted_data)
  factor_params  <- extract_factor_params_v27(fit, formatted_data)
  nugget_summary <- extract_nugget_params_v21(fit, formatted_data)
  conc_summary   <- extract_concentration_summaries(fit, formatted_data)
  true_rho_f     <- if (!is.null(true_values)) true_values$rho_factor_km else NULL
  true_loading   <- if (!is.null(true_values)) true_values$loading_true  else NULL

  if (!is.null(true_values)) {
    cat("\n=== True vs Estimated ===\n")
    if (!is.null(true_values$mu_species)) {
      for (i in seq_along(species_names)) {
        mu_sp <- qpcr_params$mu_species[i, ]
        cat(sprintf("  %s: true=%.2f | est=%.2f (95%% CI: %.2f-%.2f)\n",
                    mu_sp$species, true_values$mu_species[i], mu_sp$mean, mu_sp$q025, mu_sp$q975))
      }
    }
    if (!is.null(true_values$qpcr_detection_intercept_true))
      cat(sprintf("qPCR detection_intercept: true=%.2f | est=%.2f (%.2f-%.2f)\n",
                  true_values$qpcr_detection_intercept_true,
                  qpcr_params$qpcr_detection_intercept$mean,
                  qpcr_params$qpcr_detection_intercept$q025,
                  qpcr_params$qpcr_detection_intercept$q975))
    if (!is.null(true_values$sigma_qpcr_true))
      cat(sprintf("sigma_qpcr: true=%.2f | est=%.2f (%.2f-%.2f)\n",
                  true_values$sigma_qpcr_true, qpcr_params$sigma_qpcr$mean,
                  qpcr_params$sigma_qpcr$q025, qpcr_params$sigma_qpcr$q975))
    posterior <- rstan::extract(fit)
    for (par_nm in c("beta0", "beta1", "gamma0", "log_gamma1")) {
      tv_nm <- paste0(par_nm, "_true")
      if (!is.null(true_values[[tv_nm]]) && !is.null(posterior[[par_nm]])) {
        samp <- posterior[[par_nm]]
        cat(sprintf("%-12s true=%.2f | est=%.2f (%.2f-%.2f)\n",
                    par_nm, true_values[[tv_nm]], mean(samp),
                    quantile(samp, 0.025), quantile(samp, 0.975)))
      }
    }
    if (!is.null(true_values$gamma1_true) && !is.null(posterior$gamma1)) {
      samp <- posterior$gamma1
      cat(sprintf("%-12s true=%.2f | est=%.2f (%.2f-%.2f)\n",
                  "gamma1", true_values$gamma1_true, mean(samp),
                  quantile(samp, 0.025), quantile(samp, 0.975)))
    }
    if (!is.null(true_values$rho_factor_km)) {
      true_rho_vec <- rep(true_values$rho_factor_km, length.out = nrow(factor_params$rho_factor_summary))
      for (i in seq_len(nrow(factor_params$rho_factor_summary))) {
        row <- factor_params$rho_factor_summary[i, ]
        cat(sprintf("  rho_factor[%d]: true=%.1f | est=%.1f (%.1f-%.1f)\n",
                    row$factor, true_rho_vec[i], row$mean, row$q025, row$q975))
      }
    }
    if (!is.null(true_values$loading_true)) {
      cat("  Loading matrix (true vs posterior mean):\n")
      for (sp_i in seq_along(species_names)) {
        K_f_val <- ncol(true_values$loading_true)
        for (k in seq_len(K_f_val)) {
          true_l <- true_values$loading_true[sp_i, k]
          est_l  <- factor_params$loading_post_mean[sp_i, k]
          cat(sprintf("    loading[%s, %d]: true=%.2f | est=%.2f\n",
                      species_names[sp_i], k, true_l, est_l))
        }
      }
    }
    if (!is.null(true_values$sigma_nugget)) {
      for (i in seq_len(nrow(nugget_summary$species_summary))) {
        row <- nugget_summary$species_summary[i, ]
        cat(sprintf("  sigma_nugget %s: true=%.2f | est=%.2f (%.2f-%.2f)\n",
                    row$species, true_values$sigma_nugget[i], row$mean, row$q025, row$q975))
      }
    }
  }

  if (save_outputs) {
    # Save immediately after sampling — before any plotting that could fail.
    # Note: the @stanmodel@dso slot (compiled C++ model pointer) is not
    # preserved across sessions, but all sample data in @sim is fully intact,
    # so rstan::extract(), summary(), print() etc. all work on the loaded object.
    saveRDS(fit, file.path(output_dir, "fit_v27.rds"))
    cat(sprintf("Stan fit saved to: %s\n", file.path(output_dir, "fit_v27.rds")))
  }

  cat("\n=== Generating Plots ===\n")
  true_conc_for_plot <- if (!is.null(true_values) && !is.null(true_values$true_concentration))
    true_values$true_concentration else NULL

  # Concentration heatmaps
  heatmap_plots <- tryCatch({
    conc_df  <- extract_concentration(fit, formatted_data, "concentration")
    sp_names <- unique(conc_df$species)
    plots    <- list()
    for (sp in sp_names) {
      p <- tryCatch({
        sp_data <- dplyr::filter(conc_df, species == sp)
        sp_idx  <- which(formatted_data$sp_list$species == sp)
        if (!is.null(true_conc_for_plot) && length(sp_idx) == 1)
          sp_data$true_conc <- true_conc_for_plot[sp_data$location_idx, sp_idx]
        interp_df <- interpolate_surface(x = sp_data$x_utm, y = sp_data$y_utm,
                                         z = sp_data$mean, resolution = 80, log_transform = TRUE)
        p_base <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = interp_df, ggplot2::aes(x = x, y = y, fill = value)) +
          ggplot2::labs(title = paste("DNA Concentration:", gsub("_", " ", sp)),
                        subtitle = "Posterior mean (v27 factor-analytic HSGP)",
                        x = "UTM Easting (km)", y = "UTM Northing (km)") +
          ggplot2::coord_fixed(ratio = 1) + ggplot2::theme_minimal()
        if (!is.null(true_conc_for_plot) && "true_conc" %in% names(sp_data)) {
          p_base + ggplot2::geom_point(data = sp_data,
                                       ggplot2::aes(x = x_utm, y = y_utm, fill = true_conc),
                                       color = "white", size = 2.5, shape = 21, stroke = 0.5) +
            ggplot2::scale_fill_viridis_c(option = "plasma", name = "Conc.\nDots=true", trans = "log10")
        } else {
          p_base +
            ggplot2::scale_fill_viridis_c(option = "plasma", name = "Conc.", trans = "log10") +
            ggplot2::geom_point(data = sp_data,
                                ggplot2::aes(x = x_utm, y = y_utm),
                                color = "white", fill = "gray30", size = 2.0,
                                shape = 21, stroke = 0.5)
        }
      }, error = function(e) { warning(paste("Heatmap failed for", sp, ":", e$message)); NULL })
      if (!is.null(p)) plots[[sp]] <- p
    }
    plots
  }, error = function(e) { warning(paste("Heatmap loop failed:", e$message)); list() })

  combined_heatmaps <- combine_plots(heatmap_plots, "Species Concentration Heatmaps (copies/uL)")
  qpcr_plot         <- plot_qpcr_fit_v15(fit, formatted_data)
  detection_plot    <- plot_detection_curve_v15(fit, formatted_data)
  edna_plots        <- plot_edna_fit(fit, formatted_data)
  combined_edna     <- combine_plots(edna_plots, "eDNA Model Fit (Proportions)")
  var_decomp_plot   <- plot_variance_decomposition_v22(fit, formatted_data)

  true_vs_pred_plots    <- NULL
  combined_true_vs_pred <- NULL
  if (!is.null(true_conc_for_plot)) {
    true_vs_pred_plots    <- plot_predicted_vs_true_concentration(fit, formatted_data, true_conc_for_plot)
    combined_true_vs_pred <- combine_plots(true_vs_pred_plots, "True vs Predicted Concentration")
  }

  true_latent_fields_for_plot <- if (!is.null(true_values)) true_values$latent_fields else NULL
  latent_field_plot <- tryCatch(
    plot_latent_fields_v27(fit, formatted_data, true_latent_fields_for_plot),
    error = function(e) { warning(paste("Latent field plot failed:", e$message)); NULL })

  loading_plot <- tryCatch(
    plot_loading_matrix_v27(fit, formatted_data, true_loading),
    error = function(e) { warning(paste("Loading matrix plot failed:", e$message)); NULL })

  rho_factor_plot <- tryCatch(
    plot_rho_factor_posteriors_v27(fit, formatted_data, true_rho_f),
    error = function(e) { warning(paste("rho_factor plot failed:", e$message)); NULL })

  if (save_outputs) {
    safe_ggsave <- function(path, plot, ...) {
      tryCatch(ggplot2::ggsave(path, plot, ...),
               error = function(e) warning(sprintf("ggsave failed for %s: %s", basename(path), e$message)))
    }
    safe_ggsave(file.path(output_dir, "concentration_heatmaps_v27.png"), combined_heatmaps,
                width = 12, height = 10, dpi = 150)
    safe_ggsave(file.path(output_dir, "qpcr_fit_v27.png"),              qpcr_plot,
                width = 8,  height = 6,  dpi = 150)
    safe_ggsave(file.path(output_dir, "detection_curve_v27.png"),       detection_plot,
                width = 8,  height = 6,  dpi = 150)
    safe_ggsave(file.path(output_dir, "edna_fit_v27.png"),              combined_edna,
                width = 10, height = 8,  dpi = 150)
    safe_ggsave(file.path(output_dir, "variance_decomposition_v27.png"), var_decomp_plot,
                width = 10, height = 6,  dpi = 150)
    if (!is.null(latent_field_plot))
      safe_ggsave(file.path(output_dir, "latent_fields_v27.png"),       latent_field_plot,
                  width = 12, height = 4,  dpi = 150)
    if (!is.null(loading_plot))
      safe_ggsave(file.path(output_dir, "loading_matrix_v27.png"),      loading_plot,
                  width = 8,  height = 6,  dpi = 150)
    if (!is.null(rho_factor_plot))
      safe_ggsave(file.path(output_dir, "rho_factor_posteriors_v27.png"), rho_factor_plot,
                  width = 10, height = 4,  dpi = 150)
    if (!is.null(combined_true_vs_pred))
      safe_ggsave(file.path(output_dir, "true_vs_predicted_v27.png"),   combined_true_vs_pred,
                  width = 12, height = 10, dpi = 150)
    utils::write.csv(top_rhat, file.path(output_dir, "top_rhat_v27.csv"), row.names = FALSE)
    cat("Saved outputs to:", output_dir, "\n")
  }

  list(
    fit              = fit,
    stan_data        = stan_data,
    formatted_data   = formatted_data,
    qpcr_params      = qpcr_params,
    factor_params    = factor_params,
    nugget_summary   = nugget_summary,
    conc_summary     = conc_summary,
    diagnostics = list(n_divergent = n_div, max_treedepth_hits = max_td,
                       max_rhat = max_rhat, n_high_rhat = n_high_rhat, top_rhat = top_rhat),
    heatmap_plots           = heatmap_plots,
    combined_heatmaps       = combined_heatmaps,
    qpcr_plot               = qpcr_plot,
    detection_plot          = detection_plot,
    edna_plots              = edna_plots,
    combined_edna           = combined_edna,
    var_decomp_plot         = var_decomp_plot,
    true_vs_pred_plots      = true_vs_pred_plots,
    combined_true_vs_pred   = combined_true_vs_pred,
    latent_field_plot       = latent_field_plot,
    loading_plot            = loading_plot,
    rho_factor_plot         = rho_factor_plot
  )
}


# =============================================================================
# SECTION 8: RUN WORKFLOW
# =============================================================================

#' Full simulate -> format -> Stan data -> fit -> plot pipeline for v27
#'
#' Combines all steps into a single function call.  Pass true_values through
#' automatically from the simulation for validation plots.
#'
#' @param stan_file       Path to quant_metabar_qpcr_v27.stan.
#' @param n_locations,n_species,reference_species,n_qpcr_per_location  Simulation dims.
#' @param K_f            Number of shared latent GP fields.
#' @param chains,iter,warmup,seed  Sampling settings.
#' @param m_x,m_y        HSGP basis function counts.
#' @param max_rho        Maximum rho (km) for HSGP domain.
#' @param rho_min        Hard lower bound on rho_factor (km).
#' @param rho_factor_prior_mean,rho_factor_prior_sd  Prior on rho_factor.
#' @param gp_alpha_factor Fixed GP amplitude per factor.
#' @param loading_prior_sd Regularizing prior on loading elements.
#' @param rho_factor_km  True isotropic length scale per factor for simulation.
#' @param gp_alpha_factor_true,loading_true,sigma_nugget_true  True spatial params.
#' @param log_lambda_mean,mu_species_true  Mean log-copies and per-species means.
#' @param sigma_qpcr_true  True qPCR measurement SD.
#' @param qpcr_detection_intercept_true,qpcr_detection_slope_true  True detection params.
#' @param total_reads_per_rep  Sequencing depth per PCR replicate.
#' @param beta0_true,beta1_true,gamma0_true,gamma1_true  True dispersion params.
#' @param sigma_qpcr_prior_sd,qpcr_detection_intercept_prior_sd,
#'   qpcr_detection_slope_prior_mean,qpcr_detection_slope_prior_sd  qPCR priors.
#' @param sigma_nugget_prior_mean,sigma_nugget_prior_sd,sigma_nugget_sigma_rate  Nugget priors.
#' @param log_conc_prior_sd,beta0_prior_sd,gamma0_prior_mean,gamma0_prior_sd,
#'   log_gamma1_prior_mean,log_gamma1_prior_sd  Remaining priors.
#' @param adapt_delta,max_treedepth  HMC tuning.
#' @param save_outputs,output_dir  Output persistence.
#' @return Named list with sim_data, formatted_data, stan_data, fit_successful, and
#'   all elements returned by fit_and_plot_v27.
run_workflow_v27 <- function(
    stan_file           = "code/spatial_qPCR_mb/quant_metabar_qpcr_v27.stan",
    n_locations         = 25,
    n_species           = 4,
    K_f                 = 3,
    reference_species   = 1,
    n_qpcr_per_location = 3,
    chains              = 3,
    iter                = 2000,
    warmup              = 1000,
    seed                = 42,
    m_x = 10, m_y = 8,
    max_rho               = 300,
    rho_min               = 20.0,
    rho_factor_prior_mean = 100.0,
    rho_factor_prior_sd   = 200.0,
    gp_alpha_factor       = 1.0,
    loading_prior_sd      = 1.0,
    rho_factor_km         = 50,
    gp_alpha_factor_true  = NULL,
    loading_true          = NULL,
    sigma_nugget_true     = 0.3,
    log_lambda_mean       = log(5),
    mu_species_true       = NULL,
    sigma_qpcr_true       = 0.3,
    qpcr_detection_intercept_true = -2.0,
    qpcr_detection_slope_true     = 1.5,
    total_reads_per_rep   = 10000,
    beta0_true            = 5.0,
    beta1_true            = -2.0,
    gamma0_true           = 2.0,
    gamma1_true           = exp(1.0),
    sigma_qpcr_prior_sd   = 1.0,
    qpcr_detection_intercept_prior_sd = 5.0,
    qpcr_detection_slope_prior_mean   = 2.0,
    qpcr_detection_slope_prior_sd     = 2.0,
    sigma_nugget_prior_mean  = log(0.1),
    sigma_nugget_prior_sd    = 1e-5,
    sigma_nugget_sigma_rate  = 1.0,
    log_conc_prior_sd        = 5.0,
    beta0_prior_sd           = 1.0,
    gamma0_prior_mean        = 2.0,
    gamma0_prior_sd          = 1.0,
    log_gamma1_prior_mean    = 1.0,
    log_gamma1_prior_sd      = 0.75,
    adapt_delta              = 0.99,
    max_treedepth            = 14,
    save_outputs             = TRUE,
    output_dir               = "code/spatial_qPCR_mb/outputs_v27"
) {
  cat("\n=== STEP 1: Simulating data (v27 factor-analytic) ===\n")
  sim_data <- simulate_qpcr_metabar_data_v27(
    n_locations                   = n_locations,
    n_species                     = n_species,
    K_f                           = K_f,
    n_qpcr_per_location           = n_qpcr_per_location,
    reference_species             = reference_species,
    rho_factor_km                 = rho_factor_km,
    gp_alpha_factor_true          = gp_alpha_factor_true,
    loading_true                  = loading_true,
    sigma_nugget_true             = sigma_nugget_true,
    log_lambda_mean               = log_lambda_mean,
    mu_species_true               = mu_species_true,
    sigma_qpcr_true               = sigma_qpcr_true,
    qpcr_detection_intercept_true = qpcr_detection_intercept_true,
    qpcr_detection_slope_true     = qpcr_detection_slope_true,
    total_reads_per_rep           = total_reads_per_rep,
    beta0_true                    = beta0_true,
    beta1_true                    = beta1_true,
    gamma0_true                   = gamma0_true,
    gamma1_true                   = gamma1_true,
    seed                          = seed
  )

  cat("\n=== STEP 2: Formatting data (v27) ===\n")
  formatted_data <- format_qpcr_metabar_data_v27(
    edna_data         = sim_data$edna_data,
    qpcr_data         = sim_data$qpcr_data,
    alpha_known       = sim_data$alpha_known,
    reference_species = reference_species,
    species_names     = sim_data$species_names,
    pcr_id_col        = "pcr_replicate_id"
  )

  cat("\n=== STEP 2.5: Building Stan data (v27) ===\n")
  stan_data <- make_stan_data_v27(
    formatted_data,
    m_x = m_x, m_y = m_y,
    max_rho = max_rho,
    K_f = K_f,
    gp_alpha_factor = gp_alpha_factor,
    rho_min = rho_min,
    rho_factor_prior_mean = rho_factor_prior_mean,
    rho_factor_prior_sd   = rho_factor_prior_sd,
    loading_prior_sd      = loading_prior_sd,
    log_conc_prior_sd     = log_conc_prior_sd,
    sigma_qpcr_prior_sd   = sigma_qpcr_prior_sd,
    qpcr_detection_intercept_prior_sd = qpcr_detection_intercept_prior_sd,
    qpcr_detection_slope_prior_mean   = qpcr_detection_slope_prior_mean,
    qpcr_detection_slope_prior_sd     = qpcr_detection_slope_prior_sd,
    sigma_nugget_prior_mean = sigma_nugget_prior_mean,
    sigma_nugget_prior_sd   = sigma_nugget_prior_sd,
    sigma_nugget_sigma_rate = sigma_nugget_sigma_rate,
    beta0_prior_sd          = beta0_prior_sd,
    gamma0_prior_mean       = gamma0_prior_mean,
    gamma0_prior_sd         = gamma0_prior_sd,
    log_gamma1_prior_mean   = log_gamma1_prior_mean,
    log_gamma1_prior_sd     = log_gamma1_prior_sd
  )

  # Align true concentrations with formatted locations (handles filtered locations)
  true_values_adj <- sim_data$true_values
  if (!is.null(true_values_adj$true_concentration)) {
    n_sim_locs <- nrow(true_values_adj$true_concentration)
    if (n_sim_locs != formatted_data$N_locations) {
      cat(sprintf("Matching true concentrations: %d sim locs -> %d formatted locs\n",
                  n_sim_locs, formatted_data$N_locations))
      true_values_adj$true_concentration <- match_true_to_formatted(
        true_values_adj$true_concentration, sim_data$locations, formatted_data$location_list)
    }
  }

  cat("\n=== STEP 3: Fit + Plots (v27) ===\n")
  fit_result <- tryCatch({
    fit_and_plot_v27(
      formatted_data           = formatted_data,
      stan_file                = stan_file,
      stan_data                = stan_data,
      chains                   = chains,
      iter                     = iter,
      warmup                   = warmup,
      seed                     = seed,
      beta0_prior_sd           = beta0_prior_sd,
      gamma0_prior_mean        = gamma0_prior_mean,
      gamma0_prior_sd          = gamma0_prior_sd,
      log_gamma1_prior_mean    = log_gamma1_prior_mean,
      log_gamma1_prior_sd      = log_gamma1_prior_sd,
      adapt_delta              = adapt_delta,
      max_treedepth            = max_treedepth,
      true_values              = true_values_adj,
      save_outputs             = save_outputs,
      output_dir               = output_dir
    )
  }, error = function(e) {
    cat("\nERROR:", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(fit_result))
    return(list(sim_data = sim_data, formatted_data = formatted_data,
                stan_data = stan_data, fit_successful = FALSE))

  c(fit_result, list(sim_data = sim_data, formatted_data = formatted_data, fit_successful = TRUE))
}


# =============================================================================
# SECTION 9: EXAMPLE USAGE (uncomment to run)
# =============================================================================

# --- Quick convergence test (K_f=3, N_species=6) ---
# result <- run_workflow_v27(
#   n_locations = 100, n_species = 6, K_f = 3,
#   chains = 3, iter = 1000, warmup = 500, seed = 42,
#   m_x = 4, m_y = 5, rho_min = 20,
#   rho_factor_km = 80, log_lambda_mean = log(100),
#   n_qpcr_per_location = 3,
#   output_dir = "code/spatial_qPCR_mb/outputs_v27"
# )
# print(result$fit, pars = c("rho_factor", "loading_diag", "loading_offdiag_raw"))
# result$factor_params
# result$latent_field_plot
