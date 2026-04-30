# =============================================================================
# Driver: Format RL2501 OTU data and run quant_metabar_only.stan
# =============================================================================
# Inputs (edit paths as needed):
#   - RL2501_3cetacean_3fish_combined_otu.csv  : species (rows) x PCR replicates (cols)
#       column names parsed as F<filter>_<biorep>_<techrep>, e.g. F006_2_1
#   - RL2501_merged_sample_data.csv            : metadata, one row per PCR replicate,
#       first column = PCR replicate ID (matches OTU column name)
#   - format_data_metabar_only.R               : function library (provided)
#   - quant_metabar_only.stan                  : Stan model file
# =============================================================================

library(dplyr)
library(tidyr)

# ---- paths ------------------------------------------------------------------
otu_path   <- "Mod2_data/RL2501_3cetacean_3fish_combined_otu_EmLast.csv"
meta_path  <- "Mod2_data/RL2501_merged_sample_data.csv"
funcs_path <- "format_data_metabar_only.R"
stan_path  <- "quant_metabar_only.stan"   # change if your model lives elsewhere

# ---- 1. Load data -----------------------------------------------------------
otu  <- read.csv(otu_path,  row.names = 1, check.names = FALSE,
                 stringsAsFactors = FALSE)
meta <- read.csv(meta_path, stringsAsFactors = FALSE)
names(meta)[1] <- "pcr_replicate_id"   # first col is the PCR-rep ID

# ---- 2. Wide -> long --------------------------------------------------------
# rows are species, cols are PCR replicates -> one row per (rep, species)
otu_long <- otu %>%
  tibble::rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "pcr_replicate_id", values_to = "Nreads") %>%
  mutate(Nreads = as.integer(Nreads))

# Parse FilterID from PCR rep ID (e.g. F006_2_1 -> F006)
otu_long <- otu_long %>%
  mutate(FilterID = sub("_.*$", "", pcr_replicate_id))

# ---- 3. Build location_id by 1-km spatial merge of transect start points ---
# One row per Filter for clustering
filt_meta <- meta %>%
  filter(pcr_replicate_id %in% colnames(otu)) %>%
  distinct(FilterID, Latitude.start, Longitude.start,
           Latitude.end,   Longitude.end) %>%
  filter(!is.na(Latitude.start), !is.na(Longitude.start))

stopifnot(nrow(filt_meta) == length(unique(filt_meta$FilterID)))

# Haversine distance matrix (km)
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

# Single-linkage clustering: any two filters whose starts are <1 km share a location
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
cat(sprintf("Spatial merge (<1 km starts): %d filters -> %d locations\n",
            n, length(unique(filt_meta$location_id))))

# ---- 4. Join location_id back onto long-format OTU --------------------------
edna_data <- otu_long %>%
  left_join(filt_meta %>% select(FilterID, location_id), by = "FilterID") %>%
  filter(!is.na(location_id))     # drop any rep without coords

# ---- 5. Format for the model ------------------------------------------------
source(funcs_path)

species_names <- c(
  "Engraulis mordax",            # reference (most abundant) -> index 1
  "Merluccius productus",
  "Symbolophorus californiensis",
  "Delphinus bairdii",
  "Megaptera novaeangliae",
  "Orcinus orca"
)
alpha_known       <- rep(1.0, length(species_names))   # placeholder
reference_species <- 1L                                # Engraulis mordax

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

# ---- 6. Fit and plot --------------------------------------------------------
# Build the Stan data list explicitly so we can set N_pcr (PCR cycle count)
stan_data <- make_stan_data_metabar_only(formatted, N_pcr = 35)

results <- fit_and_plot_metabar_only(
  formatted_data = formatted,
  stan_file      = stan_path,
  stan_data      = stan_data,
  chains         = 3,
  iter           = 2000,
  warmup         = 1000,
  seed           = 42,
  adapt_delta    = 0.95,
  max_treedepth  = 12,
  save_outputs   = TRUE,
  output_dir     = "outputs_RL2501_metabar"
)

cat("\nDone. Outputs in outputs_RL2501_metabar/\n")
