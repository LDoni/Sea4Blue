source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

run_distance_decay <- function(cfg = sea4blue_config()) {
  load_packages(c("phyloseq", "dplyr", "tibble", "tidyr", "ggplot2", "geosphere", "vegan", "patchwork"))
  objs <- load_phyloseq_atlantic(cfg)
  ps <- objs$phyloseq_obj_css

  sample_data(ps)$Latitude.y <- as.numeric(sample_data(ps)$Latitude.y)
  sample_data(ps)$Longitude.y <- as.numeric(sample_data(ps)$Longitude.y)
  sample_data(ps)$thetao <- as.numeric(sample_data(ps)$thetao)

  bray_dist <- phyloseq::distance(ps, method = "bray")
  unifrac_dist <- phyloseq::distance(ps, method = "unifrac", weighted = FALSE)
  coords <- data.frame(lon = sample_data(ps)$Longitude.y, lat = sample_data(ps)$Latitude.y)
  geo_dist <- as.dist(geosphere::distm(coords, fun = geosphere::distHaversine) / 1000)
  temp_dist <- dist(sample_data(ps)$thetao)

  dd_df <- tibble::tibble(
    geographic_km = as.vector(geo_dist),
    delta_temp = as.vector(temp_dist),
    bray = as.vector(bray_dist),
    unifrac = as.vector(unifrac_dist)
  ) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(bray_sim = 1 - bray, unifrac_sim = 1 - unifrac)

  p_bray_geo <- ggplot(dd_df, aes(geographic_km, bray_sim)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(x = "Geographic distance (km)", y = "Bray-Curtis similarity") +
    theme_bw()

  p_uni_geo <- ggplot(dd_df, aes(geographic_km, unifrac_sim)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(x = "Geographic distance (km)", y = "UniFrac similarity") +
    theme_bw()

  p_bray_temp <- ggplot(dd_df, aes(delta_temp, bray_sim)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(x = "Temperature difference (°C)", y = "Bray-Curtis similarity") +
    theme_bw()

  p_uni_temp <- ggplot(dd_df, aes(delta_temp, unifrac_sim)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(x = "Temperature difference (°C)", y = "UniFrac similarity") +
    theme_bw()

  combined_plot <- (p_bray_geo + p_uni_geo) / (p_bray_temp + p_uni_temp)

  stats <- list(
    bray_geo_lm = summary(lm(bray ~ geographic_km, data = dd_df)),
    unifrac_geo_lm = summary(lm(unifrac ~ geographic_km, data = dd_df)),
    bray_geo_mantel = vegan::mantel(bray_dist, geo_dist, permutations = 9999),
    unifrac_geo_mantel = vegan::mantel(unifrac_dist, geo_dist, permutations = 9999),
    bray_partial = vegan::mantel.partial(bray_dist, geo_dist, temp_dist, permutations = 9999),
    unifrac_partial = vegan::mantel.partial(unifrac_dist, geo_dist, temp_dist, permutations = 9999)
  )

  list(data = dd_df, plots = list(bray_geo = p_bray_geo, unifrac_geo = p_uni_geo, bray_temp = p_bray_temp, unifrac_temp = p_uni_temp, combined = combined_plot), stats = stats)
}
