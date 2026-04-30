# =============================================================================
# Driver: RL2501 metabarcoding + ddPCR (Engraulis mordax anchor)
# =============================================================================
# Adds a binomial-on-droplets ddPCR likelihood to the metabar-only model.
# ddPCR anchors Engraulis concentration at locations where it was measured.
#
# Inputs (paths relative to project root):
#   Mod2_data/RL2501_3cetacean_3fish_combined_otu.csv   — species x PCR-rep reads
#   Mod2_data/RL2501_merged_sample_data.csv             — metadata, one row per PCR rep
#   Mod2_data/4_20_ddpcr_mosaic_wells_F065_2.csv        — ddPCR wells (Engraulis FAM + Sardinops HEX)
#   format_data_metabar_only.R                          — formatter for metabar-only model
#   quant_metabar_ddpcr.stan                            — new ddPCR-extended model
# =============================================================================

library(dplyr)
library(tidyr)

# ---- paths ------------------------------------------------------------------
otu_path   <- "Mod2_data/RL2501_3cetacean_3fish_combined_otu.csv"
meta_path  <- "Mod2_data/RL2501_merged_sample_data.csv"
ddpcr_path <- "Mod2_data/RL2501_ddpcr_engraulis_cleaned.csv"
funcs_path <- "format_data_metabar_only.R"
stan_path  <- "quant_metabar_ddpcr.stan"

# ---- 1. Load data -----------------------------------------------------------
otu  <- read.csv(otu_path,  row.names = 1, check.names = FALSE,
                 stringsAsFactors = FALSE)
meta <- read.csv(meta_path, stringsAsFactors = FALSE)
names(meta)[1] <- "pcr_replicate_id"

# ---- 2. Wide -> long for OTU ------------------------------------------------
otu_long <- otu %>%
  tibble::rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "pcr_replicate_id", values_to = "Nreads") %>%
  mutate(Nreads = as.integer(Nreads))

# Parse FilterID and FilterBioRep from PCR rep ID
# F006_2_1 -> Filter=F006, BioRep=2, TechRep=1, FilterBioRep=F006_2
otu_long <- otu_long %>%
  mutate(
    FilterID     = sub("_.*$", "", pcr_replicate_id),
    FilterBioRep = sub("_[^_]+$", "", pcr_replicate_id)
  )

# ---- 3. Build location_id by 1-km spatial merge of transect start points ----
# location unit = FilterBioRep
filt_meta <- meta %>%
  filter(pcr_replicate_id %in% colnames(otu)) %>%
  mutate(FilterBioRep = sub("_[^_]+$", "", pcr_replicate_id)) %>%
  distinct(FilterBioRep, Latitude.start, Longitude.start,
           Latitude.end,   Longitude.end) %>%
  filter(!is.na(Latitude.start), !is.na(Longitude.start))

stopifnot(nrow(filt_meta) == length(unique(filt_meta$FilterBioRep)))

hav_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  p1 <- lat1 * pi/180; p2 <- lat2 * pi/180
  dlat <- (lat2 - lat1) * pi/180
  dlon <- (lon2 - lon1) * pi/180
  a <- sin(dlat/2)^2 + cos(p1) * cos(p2) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

n <- nrow(filt_meta)
dmat <- matrix(0, n, n)
for (i in seq_len(n)) {
  dmat[i, ] <- hav_km(filt_meta$Latitude.start[i], filt_meta$Longitude.start[i],
                      filt_meta$Latitude.start,    filt_meta$Longitude.start)
}

adj <- dmat < 1.0
cluster <- rep(NA_integer_, n)
nxt <- 1L
for (i in seq_len(n)) {
  if (!is.na(cluster[i])) next
  q <- i; cluster[i] <- nxt
  while (length(q) > 0) {
    j <- q[1]; q <- q[-1]
    nbrs <- which(adj[j, ] & is.na(cluster))
    cluster[nbrs] <- nxt
    q <- c(q, nbrs)
  }
  nxt <- nxt + 1L
}
filt_meta$location_id <- sprintf("loc_%03d", cluster)
cat(sprintf("Spatial merge (<1 km starts): %d FilterBioReps -> %d locations\n",
            n, length(unique(filt_meta$location_id))))

# ---- 4. Join location_id back onto long-format OTU --------------------------
edna_data <- otu_long %>%
  left_join(filt_meta %>% select(FilterBioRep, location_id), by = "FilterBioRep") %>%
  filter(!is.na(location_id))

# ---- 5. Format eDNA data for the model --------------------------------------
source(funcs_path)

species_names <- c(
  "Engraulis mordax",            # reference -> index 1
  "Merluccius productus",
  "Symbolophorus californiensis",
  "Delphinus bairdii",
  "Megaptera novaeangliae",
  "Orcinus orca"
)
alpha_known       <- rep(1.0, length(species_names))
reference_species <- 1L

formatted <- format_data_metabar_only(
  edna_data         = edna_data,
  alpha_known       = alpha_known,
  reference_species = reference_species,
  species_names     = species_names,
  pcr_id_col        = "pcr_replicate_id",
  location_id_col   = "location_id",
  species_col       = "species",
  reads_col         = "Nreads"
)

# ---- 6. Load pre-cleaned ddPCR data -----------------------------------------
# The cleaned CSV has already been filtered to RL2501 + Engraulis only,
# with the F065_1 typo corrected and BioRep suffixes added. See cleaning script.
ddpcr_rl <- read.csv(ddpcr_path, stringsAsFactors = FALSE, check.names = FALSE)
names(ddpcr_rl) <- gsub("\\(copies/[^)]*\\)", "_copies_per_uL",
                        gsub("[^[:alnum:]_]+", "_", names(ddpcr_rl)))

# Standardize key column names
ddpcr_rl <- ddpcr_rl %>%
  rename_with(~ "AcceptedDroplets", matches("^Accepted_?Droplets$", ignore.case = TRUE)) %>%
  rename_with(~ "Positives",        matches("^Positives$",          ignore.case = TRUE)) %>%
  rename_with(~ "Replicate",        matches("^replicate$",          ignore.case = TRUE))

cat(sprintf("\nddPCR replicates loaded: %d (across %d FilterBioReps)\n",
            nrow(ddpcr_rl), length(unique(ddpcr_rl$FilterBioRep))))

# Sanity check: every ddPCR FilterBioRep should match a metabarcoding location
unmatched <- setdiff(unique(ddpcr_rl$FilterBioRep), filt_meta$FilterBioRep)
if (length(unmatched) > 0) {
  warning("ddPCR FilterBioReps not in metabarcoding metadata: ",
          paste(unmatched, collapse = ", "))
}

# Map FilterBioRep -> location_idx (the integer index used by the Stan model)
loc_idx_map        <- formatted$location_list %>% select(location_id, location_idx)
filtbiorep_to_locid <- filt_meta %>% select(FilterBioRep, location_id)

ddpcr_rl <- ddpcr_rl %>%
  left_join(filtbiorep_to_locid, by = "FilterBioRep") %>%
  left_join(loc_idx_map,         by = "location_id")

stopifnot(all(!is.na(ddpcr_rl$location_idx)))

cat("\nddPCR rows used in fit:\n")
print(ddpcr_rl %>% select(FilterBioRep, Replicate, Positives, AcceptedDroplets, location_idx))

# ---- 7. Build Stan data list ------------------------------------------------
# Start with the metabar-only stan_data, then add ddPCR fields
stan_data <- make_stan_data_metabar_only(formatted, N_pcr = 35)

# Augment with ddPCR fields for the new model
N_ddpcr_reps         <- nrow(ddpcr_rl)
ddpcr_rep_location   <- as.integer(ddpcr_rl$location_idx)
ddpcr_positives      <- as.integer(ddpcr_rl$Positives)
ddpcr_total_droplets <- as.integer(ddpcr_rl$AcceptedDroplets)

has_ddpcr <- rep(0L, formatted$N_locations)
has_ddpcr[unique(ddpcr_rep_location)] <- 1L

stan_data$N_ddpcr_reps         <- N_ddpcr_reps
stan_data$ddpcr_rep_location   <- ddpcr_rep_location
stan_data$ddpcr_positives      <- ddpcr_positives
stan_data$ddpcr_total_droplets <- ddpcr_total_droplets
stan_data$droplet_volume_uL    <- 0.000795          # Bio-Rad QX standard
stan_data$has_ddpcr            <- has_ddpcr

# ---- 7b. Tighten the log_concentration prior --------------------------------
# Default prior was sd=5.0, which allows log_conc from -8 to +12 (essentially flat
# over biologically reasonable values). With weak likelihood signal at non-ddPCR
# locations, this lets chains wander apart -> Rhat blows up. Tightening to sd=2.0
# is still very weakly informative but prevents the random walk to extreme values.
stan_data$log_conc_prior_sd <- 2.0
cat(sprintf("log_conc_prior: Normal(%.2f, %.2f)\n",
            stan_data$log_conc_prior_mean, stan_data$log_conc_prior_sd))

cat(sprintf("\nStan data summary: %d locations, %d species, %d eDNA obs, %d ddPCR reps\n",
            stan_data$N_locations, stan_data$N_species,
            stan_data$N_edna_obs, stan_data$N_ddpcr_reps))

# ---- 8. Initial values ------------------------------------------------------
# Data-driven init for log_concentration[loc, sp]:
#
#   (A) For Engraulis at the 6 ddPCR-anchored locations: use the binomial MLE
#       from droplet counts (most informative).
#   (B) For all other (loc, sp) cells: estimate log-concentration from the
#       observed read share at that location. Specifically, if Engraulis at the
#       location has read share p_ref and concentration anchor c_ref (~prior mean),
#       then for species sp with read share p_sp, set
#         log_conc[loc, sp] ~ log(c_ref) + log(p_sp / p_ref)
#       Cells where the species had zero reads at every replicate get a low floor.
#
# All cells are clamped >= 0.5 (Stan needs log_conc > 0 for log1m_exp).
init_fun <- function() {
  N_loc <- formatted$N_locations
  N_sp  <- formatted$N_species
  ref   <- reference_species

  # ---- (B) Build per-location read-share matrix from edna_obs_data ----
  # Sum reads per (location, species), then normalize to per-location proportions
  reads_loc <- matrix(0, N_loc, N_sp)
  for (obs in seq_len(formatted$N_edna_obs)) {
    loc <- formatted$edna_obs_location[obs]
    reads_loc[loc, ] <- reads_loc[loc, ] + formatted$edna_obs_data[obs, ]
  }
  prop_loc <- reads_loc / pmax(rowSums(reads_loc), 1)   # per-location proportions

  # Anchor concentration: prior mean (in log space)
  anchor_log_conc <- stan_data$log_conc_prior_mean        # ~2.0 -> ~7 copies/uL

  # log_conc[loc, sp] = anchor + log(p_sp / p_ref) at each location
  # If p_ref is zero at a location, fall back to anchor for ref, floor for others
  log_conc <- matrix(0.5, N_loc, N_sp)                    # default floor
  for (loc in seq_len(N_loc)) {
    p_ref <- prop_loc[loc, ref]
    if (p_ref > 0) {
      log_conc[loc, ref] <- anchor_log_conc
      for (sp in seq_len(N_sp)) {
        if (sp != ref && prop_loc[loc, sp] > 0) {
          log_conc[loc, sp] <- anchor_log_conc + log(prop_loc[loc, sp] / p_ref)
        }
      }
    } else {
      log_conc[loc, ref] <- anchor_log_conc
    }
  }

  # ---- (A) Override Engraulis at ddPCR-anchored locations with binomial MLE ----
  ddpcr_loc_summary <- aggregate(
    cbind(Positives, AcceptedDroplets) ~ location_idx,
    data = ddpcr_rl, FUN = sum
  )
  for (i in seq_len(nrow(ddpcr_loc_summary))) {
    loc  <- ddpcr_loc_summary$location_idx[i]
    pos  <- max(ddpcr_loc_summary$Positives[i], 1)        # avoid log(0)
    drps <- ddpcr_loc_summary$AcceptedDroplets[i]
    p_hat   <- pos / drps
    lam_hat <- -log(1 - p_hat) / 0.000795
    log_conc[loc, ref] <- log(lam_hat)
  }

  # ---- Clamp all cells to be safely positive (Stan needs log_conc > 0) ----
  log_conc <- pmax(log_conc, 0.5)

  # ---- Add small per-chain jitter so the 3 chains explore different regions ----
  log_conc <- log_conc + matrix(rnorm(N_loc * N_sp, 0, 0.2), N_loc, N_sp)
  log_conc <- pmax(log_conc, 0.5)                          # re-clamp post-jitter

  list(
    log_concentration = log_conc,
    beta0      = 0,
    beta1      = -5,
    gamma0     = 1,
    log_gamma1 = 0
  )
}

# ---- 9. Fit -----------------------------------------------------------------
cat("\n=== Compiling Stan model ===\n")
stan_mod <- rstan::stan_model(stan_path)

cat("\n=== Sampling ===\n")
fit <- rstan::sampling(
  stan_mod,
  data    = stan_data,
  chains  = 3,
  iter    = 2000,
  warmup  = 1000,
  seed    = 42,
  init    = init_fun,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

# ---- 10. Save and diagnose --------------------------------------------------
out_dir <- "results_Mod2_ddpcr_metabar"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(fit, file.path(out_dir, "fit_ddpcr_metabar.rds"))

cat("\n=== Diagnostics ===\n")
# Use rstan:: explicitly to avoid method dispatch issues, and limit pars to
# real-valued ones (skip integer counts like n_edna_locations that confuse summary).
diag_pars   <- c("log_concentration", "concentration",
                 "beta0", "beta1", "gamma0", "gamma1", "log_gamma1",
                 "psi", "species_offset", "mean_phi",
                 "amplification_bias_ratio", "lp__")
fit_summary <- rstan::summary(fit, pars = diag_pars)$summary
cat(sprintf("Max Rhat:  %.3f\n", max(fit_summary[, "Rhat"], na.rm = TRUE)))
cat(sprintf("Min n_eff: %.0f\n", min(fit_summary[, "n_eff"], na.rm = TRUE)))
sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
n_div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
cat(sprintf("Divergences: %d\n", n_div))

cat("\nFit saved to:", file.path(out_dir, "fit_ddpcr_metabar.rds"), "\n")

# ---- 11. ddPCR-specific posterior check -------------------------------------
# Compare observed positive droplets to posterior predicted at each ddPCR rep
posterior_pred <- rstan::extract(fit, pars = "ddpcr_pred_positives")$ddpcr_pred_positives
ddpcr_check <- data.frame(
  FilterBioRep      = ddpcr_rl$FilterBioRep,
  Replicate         = ddpcr_rl$Replicate,
  observed          = ddpcr_positives,
  pred_mean         = colMeans(posterior_pred),
  pred_q025         = apply(posterior_pred, 2, quantile, 0.025),
  pred_q975         = apply(posterior_pred, 2, quantile, 0.975),
  AcceptedDroplets  = ddpcr_total_droplets
)
cat("\nddPCR observed vs posterior predicted:\n")
print(ddpcr_check)
write.csv(ddpcr_check, file.path(out_dir, "ddpcr_posterior_check.csv"), row.names = FALSE)

# ---- 12. Plots --------------------------------------------------------------
cat("\n=== Generating plots ===\n")

# Use the metabar-only plotting helpers
bb_plot     <- plot_bb_params_metabar(fit)
edna_plots  <- plot_edna_fit_metabar(fit, formatted)
offset_plot <- plot_species_offsets_metabar(fit, formatted)

edna_clean <- Filter(Negate(is.null), edna_plots)
combined_edna <- if (length(edna_clean) > 0) {
  patchwork::wrap_plots(edna_clean, ncol = ceiling(sqrt(length(edna_clean)))) +
    patchwork::plot_annotation(title = "eDNA Model Fit (Proportions)")
} else NULL

# ---- True-vs-predicted concentration plot for Engraulis ---------------------
# Build a partial "true_concentration" matrix from ddPCR data:
# Engraulis at 6 ddPCR-anchored locations gets the binomial MLE concentration.
# Other cells stay NA (will be plotted with no point but with the posterior CI).
true_conc_partial <- matrix(NA_real_,
                            nrow = formatted$N_locations,
                            ncol = formatted$N_species)
ddpcr_loc_summary <- aggregate(
  cbind(Positives, AcceptedDroplets) ~ location_idx,
  data = ddpcr_rl, FUN = sum
)
for (i in seq_len(nrow(ddpcr_loc_summary))) {
  loc  <- ddpcr_loc_summary$location_idx[i]
  pos  <- ddpcr_loc_summary$Positives[i]
  drps <- ddpcr_loc_summary$AcceptedDroplets[i]
  # Binomial MLE; for pos==0, use upper bound from "see <= 0 in N trials" (1/N)
  if (pos > 0) {
    p_hat   <- pos / drps
    lam_hat <- -log(1 - p_hat) / 0.000795
    true_conc_partial[loc, reference_species] <- lam_hat
  }
  # pos==0 cells stay NA; we don't have a point estimate, only an upper bound
}

# Custom plot for Engraulis only — show posterior with ddPCR-anchored points where available
engraulis_idx <- reference_species
posterior     <- rstan::extract(fit, pars = "concentration")$concentration
sp_conc       <- posterior[, , engraulis_idx]
engraulis_df  <- data.frame(
  location_idx = seq_len(ncol(sp_conc)),
  pred_mean    = colMeans(sp_conc),
  pred_q025    = apply(sp_conc, 2, quantile, 0.025),
  pred_q975    = apply(sp_conc, 2, quantile, 0.975),
  true_ddpcr   = true_conc_partial[, engraulis_idx]
)

library(ggplot2)
ddpcr_truth_plot <- ggplot(engraulis_df, aes(x = location_idx)) +
  geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975), color = "steelblue", alpha = 0.4) +
  geom_point(aes(y = pred_mean), color = "steelblue", size = 1.5) +
  geom_point(aes(y = true_ddpcr), color = "darkorange", size = 3, shape = 18,
             na.rm = TRUE) +
  scale_y_log10() +
  labs(title    = "Engraulis mordax: posterior concentration vs ddPCR ground truth",
       subtitle = "Blue = posterior mean ± 95% CI from joint model; orange diamond = ddPCR MLE",
       x = "Location index", y = "Concentration (copies/uL, log scale)") +
  theme_bw()

# ---- Save plots -------------------------------------------------------------
if (!is.null(bb_plot))       ggplot2::ggsave(file.path(out_dir, "bb_params.png"),         bb_plot,         width = 12, height = 6,  dpi = 150)
if (!is.null(combined_edna)) ggplot2::ggsave(file.path(out_dir, "edna_fit.png"),         combined_edna,   width = 10, height = 8,  dpi = 150)
if (!is.null(offset_plot))   ggplot2::ggsave(file.path(out_dir, "species_offsets.png"),   offset_plot,     width = 10, height = 8,  dpi = 150)
ggplot2::ggsave(file.path(out_dir, "engraulis_vs_ddpcr.png"), ddpcr_truth_plot,            width = 12, height = 5,  dpi = 150)

# ddPCR observed vs predicted positive droplets (per replicate)
ddpcr_plot <- ggplot(ddpcr_check, aes(x = observed, y = pred_mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = pred_q025, ymax = pred_q975), width = 0, alpha = 0.4, color = "darkorange3") +
  geom_point(size = 3, color = "darkorange") +
  geom_text(aes(label = FilterBioRep), size = 3, vjust = -0.7, alpha = 0.7) +
  labs(title    = "ddPCR: Observed vs Posterior Predicted Positive Droplets",
       subtitle = "Engraulis mordax — error bars = 95% predictive interval",
       x = "Observed positives", y = "Posterior predicted positives") +
  theme_bw()
ggplot2::ggsave(file.path(out_dir, "ddpcr_obs_vs_pred.png"), ddpcr_plot,
                width = 8, height = 6, dpi = 150)

# ---- Top-Rhat table ---------------------------------------------------------
# Save the 30 worst-converged parameters as a CSV for diagnostic review
worst_pars <- fit_summary[order(fit_summary[, "Rhat"], decreasing = TRUE), ][1:30, ]
worst_df <- data.frame(
  parameter = rownames(worst_pars),
  mean      = worst_pars[, "mean"],
  sd        = worst_pars[, "sd"],
  n_eff     = worst_pars[, "n_eff"],
  Rhat      = worst_pars[, "Rhat"]
)
write.csv(worst_df, file.path(out_dir, "top_rhat.csv"), row.names = FALSE)
cat("\nTop 5 worst-converged parameters:\n")
print(head(worst_df, 5))

cat("\nDone. Outputs in", out_dir, "\n")
