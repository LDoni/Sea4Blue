source('/home/userbio/Sea4Blue/scripts/utils/helpers.R')

run_mechanisms_distance <- function(cfg = sea4blue_config()) {
  load_packages(c('dplyr', 'ggplot2', 'reshape2', 'magrittr', 'geosphere', 'phyloseq'))
  objs <- load_phyloseq_atlantic(cfg)
  meta <- objs$sample_data
  meta$Sample_ID <- rownames(meta)

  bNTI <- as.matrix(read.csv(cfg$paths$btnti_csv, row.names = 1, check.names = FALSE))
  RC_BC <- as.matrix(read.csv(cfg$paths$rc_bray_csv, row.names = 1, check.names = FALSE))

  bNTI_mod <- magrittr::set_colnames(bNTI, rownames(bNTI))
  bNTI_mod[lower.tri(bNTI_mod)] <- NA
  diag(bNTI_mod) <- NA
  RC_BC_mod <- RC_BC
  rownames(RC_BC_mod) <- rownames(bNTI_mod)
  colnames(RC_BC_mod) <- colnames(bNTI_mod)
  RC_BC_mod[lower.tri(RC_BC_mod)] <- NA
  diag(RC_BC_mod) <- NA

  pair_df <- reshape2::melt(bNTI_mod) %>%
    cbind(reshape2::melt(RC_BC_mod)) %>%
    magrittr::set_colnames(c('From_Sample', 'To_Sample', 'bNTI', 'x1', 'x2', 'RC_BC')) %>%
    dplyr::select(From_Sample, To_Sample, bNTI, RC_BC) %>%
    dplyr::filter(!is.na(bNTI)) %>%
    dplyr::mutate(Mechanism = classify_mechanism(bNTI, RC_BC)) %>%
    dplyr::left_join(meta %>% dplyr::select(Sample_ID, Longitude.y, Latitude.y, Zona) %>% dplyr::rename(From_Lon = Longitude.y, From_Lat = Latitude.y, Zone_1 = Zona), by = c('From_Sample' = 'Sample_ID')) %>%
    dplyr::left_join(meta %>% dplyr::select(Sample_ID, Longitude.y, Latitude.y, Zona) %>% dplyr::rename(To_Lon = Longitude.y, To_Lat = Latitude.y, Zone_2 = Zona), by = c('To_Sample' = 'Sample_ID')) %>%
    dplyr::mutate(
      geographic_km = geosphere::distHaversine(cbind(From_Lon, From_Lat), cbind(To_Lon, To_Lat)) / 1000,
      distance_bin = cut(geographic_km, breaks = c(0, 250, 500, 1000, 2000, 3000, 5000, Inf), include.lowest = TRUE, right = FALSE),
      zone_pair = vapply(seq_len(dplyr::n()), function(i) paste(sort(c(as.character(Zone_1[i]), as.character(Zone_2[i]))), collapse = '-'), character(1))
    )

  summary_tbl <- pair_df %>%
    dplyr::filter(!is.na(distance_bin)) %>%
    dplyr::group_by(distance_bin, Mechanism) %>%
    dplyr::summarise(N = dplyr::n(), .groups = 'drop_last') %>%
    dplyr::mutate(Prop = N / sum(N)) %>%
    dplyr::ungroup()

  box_plot <- ggplot(pair_df, aes(Mechanism, geographic_km, fill = Mechanism)) +
    geom_boxplot(outlier.alpha = 0.2) +
    coord_flip() +
    labs(x = NULL, y = 'Geographic distance (km)') +
    theme_bw() +
    theme(legend.position = 'none')

  stacked_plot <- ggplot(summary_tbl, aes(distance_bin, Prop * 100, fill = Mechanism)) +
    geom_col(color = 'black', linewidth = 0.2) +
    labs(x = 'Geographic distance bin (km)', y = 'Mechanism proportion (%)') +
    theme_bw()

  zone_plot <- pair_df %>%
    dplyr::group_by(zone_pair, Mechanism) %>%
    dplyr::summarise(N = dplyr::n(), .groups = 'drop_last') %>%
    dplyr::mutate(Prop = N / sum(N)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(zone_pair, Prop * 100, fill = Mechanism)) +
    geom_col(color = 'black', linewidth = 0.2) +
    labs(x = 'Zone pair', y = 'Mechanism proportion (%)') +
    theme_bw()

  corr <- suppressWarnings(cor.test(pair_df$geographic_km, abs(pair_df$bNTI), method = 'spearman', exact = FALSE))
  corr_tbl <- tibble::tibble(metric = 'geo_km_vs_abs_bNTI', rho = unname(corr$estimate), p_value = corr$p.value)

  list(detail = pair_df, summary = summary_tbl, correlation = corr_tbl, box_plot = box_plot, stacked_plot = stacked_plot, zone_plot = zone_plot)
}
