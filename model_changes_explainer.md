# From `quant_metabar_only.stan` to `quant_metabar_ddpcr.stan`

This document describes what changed between the original metabarcoding-only Stan model (`quant_metabar_only.stan`) and the version that ran successfully (`quant_metabar_ddpcr.stan`), which adds a ddPCR likelihood term to anchor *Engraulis mordax* concentrations at the 6 RL2501 locations where ddPCR data exists.

The two files share the same eDNA likelihood structure — the differences are concentrated in the qPCR/ddPCR observation block and the parameters/data fields it depends on.

## Summary of changes

| Change | What was in `quant_metabar_only` | What's in `quant_metabar_ddpcr` |
|---|---|---|
| qPCR/ddPCR machinery | All commented out (placeholder) | Active ddPCR likelihood |
| Reference-species observation model | None (eDNA-only) | Binomial on positive droplets |
| Detection hurdle | qPCR-style (commented) | None (ddPCR has no LOD) |
| `lambda_K` ZTP mean | Declared but never assigned | Declared and properly assigned |
| Header docstring | Bare one-liner | Documents the ddPCR extension |

## 1. The active ddPCR observation model

This is the heart of the change — replacing the dormant qPCR hurdle with a working ddPCR likelihood.

### What was there (commented out in `quant_metabar_only`)

```stan
// for (rep in 1:N_qpcr_reps) {
//   int loc = qpcr_rep_location[rep];
//   real log_conc_ref = log_concentration[loc, reference_species];
//   real logit_p = qpcr_detection_intercept + qpcr_detection_slope * log_conc_ref;
//   target += bernoulli_logit_lpmf(qpcr_rep_detected[rep] | logit_p);
//   if (qpcr_rep_detected[rep] == 1)
//     target += normal_lpdf(qpcr_rep_log_conc[rep] | log_conc_ref, sigma_qpcr);
// }
```

A two-part **hurdle**: a Bernoulli detection submodel (`logit(p_detect) = intercept + slope × log_conc`) plus a Normal observation model on log-concentration when detected. This is the right structure for qPCR because qPCR has a real limit-of-detection: below some threshold the assay fails to amplify and you get no signal at all.

### What's there now (active in `quant_metabar_ddpcr`)

```stan
for (rep in 1:N_ddpcr_reps) {
  int loc = ddpcr_rep_location[rep];
  real lambda_ref = lambda[loc, reference_species];
  real p_positive = 1 - exp(-lambda_ref * droplet_volume_uL);
  ddpcr_positives[rep] ~ binomial(ddpcr_total_droplets[rep], p_positive);
}
```

A single binomial likelihood on positive droplets per well. The probability that any given droplet is positive comes from the Poisson zero formula: a droplet of volume `v_drop` (~0.795 nL on QX200) sampled from a sample with concentration `λ` copies/μL has probability `1 - exp(-λ × v_drop)` of containing at least one target molecule.

### Why this is the right change for ddPCR

ddPCR is fundamentally different from qPCR. Where qPCR gives you a continuous Cq value (or non-detection), ddPCR partitions the sample into ~16,000 droplets and counts how many are positive. There's no LOD-driven censoring: even one positive droplet gives you a (wide-uncertainty) quantitative concentration estimate, and zero positives is informative quantitative data — a Poisson zero that puts an upper bound on concentration. The binomial likelihood handles all of these regimes coherently with no special cases. The hurdle structure would have either crashed on `log(0)` for zero-positive wells or required mis-coding them as "non-detections," which is the wrong framing.

## 2. Parameters added and removed

### Added to the `data` block

```stan
int<lower=0> N_ddpcr_reps;
array[N_ddpcr_reps] int<lower=1, upper=N_locations> ddpcr_rep_location;
array[N_ddpcr_reps] int<lower=0>                     ddpcr_positives;
array[N_ddpcr_reps] int<lower=1>                     ddpcr_total_droplets;
real<lower=0>                                         droplet_volume_uL;
array[N_locations] int<lower=0, upper=1> has_ddpcr;
```

Each ddPCR replicate contributes `(positives, total_droplets, location)`. The droplet volume is fixed data because it's a property of the ddPCR system, not something to estimate.

### Removed from the data block

The original (commented) qPCR data fields and their priors:

- `qpcr_rep_detected` (binary 0/1)
- `qpcr_rep_log_conc` (continuous, sentinel 0 when not detected)
- `sigma_qpcr_prior_sd`, `qpcr_detection_intercept_prior_sd`, `qpcr_detection_slope_prior_mean`, `qpcr_detection_slope_prior_sd`

### Removed from the `parameters` block

```stan
// real log_sigma_qpcr;
// real qpcr_detection_intercept;
// real qpcr_detection_slope;
```

These three free parameters described qPCR's detection-vs-quantification regime. ddPCR doesn't have that regime, so they're gone. **No new free parameters are introduced** — the binomial likelihood ties directly to the existing `log_concentration[loc, reference_species]`. ddPCR is a purely observational anchor; it doesn't add modeling complexity.

### Counter changes in `transformed data`

Old: `n_qpcr_per_loc`, `n_detected_per_loc` (commented).
New: `n_ddpcr_per_loc` (number of ddPCR replicates per location).

## 3. The `lambda_K` fix (a real bug from the original model)

In the original `quant_metabar_only.stan`, the matrix `lambda_K` was declared but never populated:

```stan
matrix[N_locations, N_species] lambda_K;          // declared
matrix[N_locations, N_species] log_lambda_K;
// ... loop computes log_lambda_K but never lambda_K
```

`lambda_K` is the **zero-truncated Poisson (ZTP) mean** — the conditional mean of read counts given that the species is present (`λ / (1 − exp(−λ))`). It's distinct from `lambda` (the unconditional Poisson mean, which includes the zero mass). The original model used `log_lambda_K` correctly in the likelihood but left `lambda_K` itself as `NaN`, which caused `summary(fit)` to fail with "undefined values" warnings and prevented diagnostic computation.

The new model fixes this:

```stan
for (sp in 1:N_species) {
  lambda[loc, sp]       = exp(log_concentration[loc, sp]);
  log_lambda_K[loc, sp] = log_concentration[loc, sp] - log1m_exp(-log_concentration[loc, sp]);
  lambda_K[loc, sp]     = exp(log_lambda_K[loc, sp]);     // <-- new line
}
```

This is a strict bug fix — the model behavior is unchanged for sampling (the likelihood didn't depend on `lambda_K`), but post-hoc diagnostic and plotting code now works correctly.

## 4. Posterior predictive and log-likelihood for ddPCR

The `generated quantities` block also got the ddPCR analogues of what the qPCR machinery would have produced:

### Replaced (qPCR posterior predictive, was commented)

```stan
// array[N_qpcr_reps] int qpcr_pred_detected;
// vector[N_qpcr_reps] qpcr_pred_log_conc;
```

### Added (ddPCR posterior predictive)

```stan
array[N_ddpcr_reps] int ddpcr_pred_positives;
for (rep in 1:N_ddpcr_reps) {
  int loc = ddpcr_rep_location[rep];
  real lambda_ref = lambda[loc, reference_species];
  real p_positive = 1 - exp(-lambda_ref * droplet_volume_uL);
  ddpcr_pred_positives[rep] = binomial_rng(ddpcr_total_droplets[rep], p_positive);
}
```

This is what powers `ddpcr_obs_vs_pred.png` in the output.

### Replaced (qPCR log-likelihood, was commented)

```stan
// vector[N_qpcr_reps] log_lik_qpcr;
```

### Added (ddPCR log-likelihood)

```stan
vector[N_ddpcr_reps] log_lik_ddpcr;
for (rep in 1:N_ddpcr_reps) {
  int loc = ddpcr_rep_location[rep];
  real lambda_ref = lambda[loc, reference_species];
  real p_positive = 1 - exp(-lambda_ref * droplet_volume_uL);
  log_lik_ddpcr[rep] = binomial_lpmf(ddpcr_positives[rep]
                                      | ddpcr_total_droplets[rep], p_positive);
}
```

Plus a location-level aggregate `log_lik_ddpcr_loc` for LOO-CV at the location level.

The `expected_qpcr_loc` quantity was renamed to `expected_log_conc_ref` (same definition: `log_concentration[loc, reference_species]`) since "expected qPCR" no longer makes semantic sense.

## 5. What stayed the same

The *eDNA* observation model is byte-for-byte identical between the two files:

- Beta-Binomial dispersion structure (`beta0`, `beta1`, `gamma0`, `log_gamma1`)
- Poisson zero-inflation
- ZTP-corrected proportions via softmax over `log_lambda_K`
- All eDNA priors
- `concentration`, `total_concentration`, `true_proportions`, `mu_metabar_loc`, `species_offset`, `psi`, `amplification_bias_ratio` derived quantities
- eDNA posterior predictive (`edna_pred`, `edna_pred_loc`) and log-likelihood (`log_lik_edna`, `log_lik_edna_loc`)
- The 6 species, 48 locations, `reference_species = 1` (Engraulis mordax) configuration

This means **the metabarcoding analysis itself didn't change** — the new model adds an additional source of information (ddPCR droplet counts) that pulls `log_concentration[loc, ref]` toward what the droplet counts imply at the 6 anchored locations. Through the eDNA likelihood, that propagates to influence other species' concentrations at those same locations (because eDNA constrains relative concentrations).

## 6. What this means in practice

**At the 6 ddPCR-anchored locations** (F022_1, F031_1, F047_2, F065_2, F074_1, F095_1):

- Engraulis concentration is now constrained by both metabarcoding read shares and a Poisson process on droplet counts. Posterior CI should be much tighter than in the metabar-only fit.
- Other species' concentrations at these same locations gain indirect constraint via the metabarcoding likelihood: if Engraulis is anchored at concentration X, and the read share of (say) hake is some fraction of Engraulis at this location, then hake's concentration is also indirectly constrained.

**At the 42 locations without ddPCR data**:

- The model behaves like the original metabar-only model. Posterior reflects whatever metabarcoding alone can identify (which is mostly relative composition; absolute scale stays weakly informed by the `log_conc_prior_mean`).

**Globally**:

- Species-level offsets (`species_offset[sp]`) and the global mean (`psi`) are now identified more sharply, because the ddPCR anchor pins down absolute scale at 6 sites — which then informs the priors that the rest of the species/locations live within.

## 7. Note on the driver script and convergence fixes

In addition to the Stan model changes, the working run also depended on two changes to the R driver script (`run_metabar_ddpcr_RL2501.R`) that aren't part of the Stan model itself but were necessary for convergence:

1. **Tighter prior**: `log_conc_prior_sd` was reduced from 5.0 to 2.0. The default of 5.0 allowed log-concentration values from −8 to +12 (~`exp(-8)` to ~`exp(12)` copies/μL), which is so wide that the prior is nearly flat. Tightening to 2.0 keeps the prior weakly informative without allowing chains to wander into biologically nonsensical regions.

2. **Data-driven initialization**: instead of random `runif(1, 3)` starting values for all 288 `log_concentration` entries, each chain now starts from values derived from the observed read shares at each location, with the 6 ddPCR-anchored Engraulis cells initialized from the binomial MLE of the droplet counts. Small per-chain jitter (`Normal(0, 0.2)`) keeps chains slightly different so Rhat is still meaningful.

Together these two changes brought the fit from `Max Rhat = 36.8`, 512 divergences (essentially unconverged) to a successful run.
