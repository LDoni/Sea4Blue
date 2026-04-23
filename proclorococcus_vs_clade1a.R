library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)

# ------------------------------
# FUNZIONE: analisi distance decay per un taxon
# ------------------------------
distance_decay_taxon <- function(ps_obj, genus_name = NULL, clade_name = NULL) {
  # Estrai tax_table come character
  tax_df <- as.data.frame(as.matrix(tax_table(ps_obj)))
  
  if (!is.null(genus_name)) {
    # Filtra ASV con Genus == genus_name
    taxa_to_keep <- rownames(tax_df)[tax_df$Genus == genus_name]
    label <- genus_name
  } else if (!is.null(clade_name)) {
    # Filtra ASV con Genus == clade_name
    taxa_to_keep <- rownames(tax_df)[tax_df$Genus == clade_name]
    label <- clade_name
  } else {
    stop("Devi fornire genus_name o clade_name")
  }
  
  # Controllo se ci sono ASV
  if (length(taxa_to_keep) == 0) stop("Nessun ASV trovato per ", label)
  
  ps_sub <- prune_taxa(taxa_to_keep, ps_obj)
  
  cat("Number of ASVs for", label, ":", ntaxa(ps_sub), "\n")
  
  # Distance matrices
  bray_sub <- phyloseq::distance(ps_sub, method = "bray")
  unifrac_sub <- phyloseq::distance(ps_sub, method = "unifrac", weighted = FALSE)
  
  # Geographic distances
  coords <- data.frame(
    lon = sample_data(ps_sub)$Longitude.y,
    lat = sample_data(ps_sub)$Latitude.y
  )
  geo_dist <- as.dist(geosphere::distm(coords, fun = geosphere::distHaversine)/1000)
  
  # Environmental distances
  num_vars <- c("thetao", "so", "chl", "o2", "no3", "po4", "si")
  env_data <- sample_data(ps_sub)[, num_vars]
  env_dist <- dist(scale(env_data))
  
  # Dataframe per plotting
  dd_df <- tibble(
    geographic_km = as.vector(geo_dist),
    bray = as.vector(bray_sub),
    unifrac = as.vector(unifrac_sub),
    env_dist = as.vector(env_dist)
  ) %>% drop_na()
  
  # Linear models
  lm_geo_bray <- lm(bray ~ geographic_km, data = dd_df)
  lm_geo_unif <- lm(unifrac ~ geographic_km, data = dd_df)
  lm_env_bray <- lm(bray ~ env_dist, data = dd_df)
  lm_env_unif <- lm(unifrac ~ env_dist, data = dd_df)
  
  # Mantel tests
  mantel_geo_bray <- vegan::mantel(bray_sub, geo_dist, permutations = 9999)
  mantel_geo_unif <- vegan::mantel(unifrac_sub, geo_dist, permutations = 9999)
  mantel_env_bray <- vegan::mantel(bray_sub, env_dist, permutations = 9999)
  mantel_env_unif <- vegan::mantel(unifrac_sub, env_dist, permutations = 9999)
  
  # Plots
  p_bray_geo <- ggplot(dd_df, aes(geographic_km, bray)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(title = paste("Bray–Curtis vs Geo –", label),
         x = "Geographic distance (km)", y = "Bray–Curtis")
  
  p_unif_geo <- ggplot(dd_df, aes(geographic_km, unifrac)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(title = paste("UniFrac vs Geo –", label),
         x = "Geographic distance (km)", y = "UniFrac")
  
  list(
    label = label,
    n_ASVs = ntaxa(ps_sub),
    plots = list(p_bray_geo = p_bray_geo, p_unif_geo = p_unif_geo),
    lms = list(lm_geo_bray = lm_geo_bray, lm_geo_unif = lm_geo_unif,
               lm_env_bray = lm_env_bray, lm_env_unif = lm_env_unif),
    mantel = list(geo_bray = mantel_geo_bray, geo_unif = mantel_geo_unif,
                  env_bray = mantel_env_bray, env_unif = mantel_env_unif)
  )
}

# ------------------------------
# ANALISI TAXA
# ------------------------------

# 1️⃣ Prochlorococcus MIT9313
res_proch <- distance_decay_taxon(phyloseq_obj_css, genus_name = "Prochlorococcus MIT9313")

# 2️⃣ SAR11 Clade Ia
res_sar11 <- distance_decay_taxon(phyloseq_obj_css, genus_name = NULL, clade_name = "Clade Ia")

# ------------------------------
# VISUALIZZAZIONE
# ------------------------------

# Esempio: plot Bray–Curtis distance decay geografica
res_proch$plots$bray_geo
res_sar11$plots$bray_geo

# Esempio: plot UniFrac vs ambiente
res_proch$plots$unif_env
res_sar11$plots$unif_env

# ------------------------------
# OPZIONALE: stampare i valori dei Mantel
# ------------------------------
res_proch$mantel
res_sar11$mantel

# ------------------------------
# OPZIONALE: stampare R² dei linear models
summary(res_proch$lms$lm_geo_bray)
summary(res_sar11$lms$lm_geo_bray)
