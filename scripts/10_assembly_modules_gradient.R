source('/home/userbio/Sea4Blue/scripts/utils/helpers.R')
source('/home/userbio/Sea4Blue/scripts/05_stegen_mechanisms.R')
source('/home/userbio/Sea4Blue/scripts/07_module_environment_preferences.R')

ordered_pair <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  a <- ifelse(x <= y, x, y)
  b <- ifelse(x <= y, y, x)
  list(a = a, b = b)
}

pair_from_matrix <- function(mat, value_name) {
  idx <- which(lower.tri(mat), arr.ind = TRUE)
  out <- data.frame(
    Sample_1 = rownames(mat)[idx[, 1]],
    Sample_2 = colnames(mat)[idx[, 2]],
    stringsAsFactors = FALSE
  )
  pair <- ordered_pair(out$Sample_1, out$Sample_2)
  out[[value_name]] <- mat[idx]
  out$Pair_A <- pair$a
  out$Pair_B <- pair$b
  tibble::as_tibble(out)
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(estimate = NA_real_, p_value = NA_real_))
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = 'spearman', exact = FALSE))
  c(estimate = unname(test$estimate), p_value = test$p.value)
}

run_assembly_modules_gradient <- function(cfg = sea4blue_config(), top_n = 4) {
  load_packages(c('dplyr', 'tibble', 'tidyr', 'ggplot2', 'phyloseq', 'vegan', 'patchwork', 'forcats'))
  zone_colors <- sea4blue_zone_colors()
  objs <- load_phyloseq_atlantic(cfg)
  meta <- objs$sample_data
  meta$Sample_ID <- rownames(meta)
  meta$thetao <- as.numeric(meta$thetao)
  meta$Longitude.x <- as.numeric(meta$Longitude.x)
  meta$Zona <- factor(as.character(meta$Zona), levels = names(zone_colors))

  steg <- get_mechanism_proportions(cfg)
  mod <- run_module_environment_preferences(cfg)
  main_modules <- mod$diversity |>
    dplyr::arrange(dplyr::desc(mean_relative_abundance), module_sort_key(Module)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Module) |>
    as.character()

  module_mat <- t(mod$module_abundance[main_modules, , drop = FALSE])
  module_dist <- as.matrix(vegan::vegdist(module_mat, method = 'bray'))
  module_pair <- pair_from_matrix(module_dist, 'module_bray') |>
    dplyr::select(Pair_A, Pair_B, module_bray)

  step_pair <- ordered_pair(steg$pairwise$Sample_ID, steg$pairwise$To_Sample)
  pair_df <- steg$pairwise |>
    dplyr::rename(Sample_1 = Sample_ID, Sample_2 = To_Sample) |>
    dplyr::mutate(Pair_A = step_pair$a, Pair_B = step_pair$b) |>
    dplyr::left_join(module_pair, by = c('Pair_A', 'Pair_B')) |>
    dplyr::left_join(
      meta |>
        dplyr::select(Sample_ID, Zona, thetao, Longitude.x) |>
        dplyr::rename(Zone_1 = Zona, Temp_1 = thetao, Lon_1 = Longitude.x),
      by = c('Sample_1' = 'Sample_ID')
    ) |>
    dplyr::left_join(
      meta |>
        dplyr::select(Sample_ID, Zona, thetao, Longitude.x) |>
        dplyr::rename(Zone_2 = Zona, Temp_2 = thetao, Lon_2 = Longitude.x),
      by = c('Sample_2' = 'Sample_ID')
    ) |>
    dplyr::mutate(
      module_similarity = 1 - module_bray,
      temp_diff = abs(Temp_1 - Temp_2),
      lon_diff = abs(Lon_1 - Lon_2),
      zone_pair = vapply(
        seq_len(dplyr::n()),
        function(i) paste(sort(c(as.character(Zone_1[i]), as.character(Zone_2[i]))), collapse = '-'),
        character(1)
      ),
      temp_bin = cut(temp_diff, breaks = c(0, 1, 2, 4, 6, 8, Inf), include.lowest = TRUE, right = FALSE),
      selection_state = dplyr::case_when(
        Mechanism %in% c('Homogeneous Selection', 'Heterogeneous Selection') ~ 'Selection',
        Mechanism %in% c('Dispersal Limitation', 'Homogenising Dispersal') ~ 'Dispersal',
        TRUE ~ 'Drift'
      )
    )

  module_mech_summary <- pair_df |>
    dplyr::group_by(Mechanism) |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_module_similarity = mean(module_similarity, na.rm = TRUE),
      median_module_similarity = median(module_similarity, na.rm = TRUE),
      mean_temp_diff = mean(temp_diff, na.rm = TRUE),
      .groups = 'drop'
    ) |>
    dplyr::arrange(dplyr::desc(mean_module_similarity))

  zone_summary <- pair_df |>
    dplyr::filter(!is.na(zone_pair)) |>
    dplyr::group_by(zone_pair, Mechanism) |>
    dplyr::summarise(n = dplyr::n(), .groups = 'drop_last') |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()

  temp_summary <- pair_df |>
    dplyr::filter(!is.na(temp_bin)) |>
    dplyr::group_by(temp_bin, Mechanism) |>
    dplyr::summarise(n = dplyr::n(), .groups = 'drop_last') |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::ungroup()

  cor_tests <- tibble::tibble(
    metric = c('module_similarity_vs_bNTI', 'module_similarity_vs_RC_BC', 'module_similarity_vs_temp_diff', 'module_similarity_vs_lon_diff'),
    estimate = c(
      safe_spearman(pair_df$module_similarity, pair_df$bNTI)[['estimate']],
      safe_spearman(pair_df$module_similarity, pair_df$RC_BC)[['estimate']],
      safe_spearman(pair_df$module_similarity, pair_df$temp_diff)[['estimate']],
      safe_spearman(pair_df$module_similarity, pair_df$lon_diff)[['estimate']]
    ),
    p_value = c(
      safe_spearman(pair_df$module_similarity, pair_df$bNTI)[['p_value']],
      safe_spearman(pair_df$module_similarity, pair_df$RC_BC)[['p_value']],
      safe_spearman(pair_df$module_similarity, pair_df$temp_diff)[['p_value']],
      safe_spearman(pair_df$module_similarity, pair_df$lon_diff)[['p_value']]
    )
  )

  mech_similarity_plot <- ggplot(pair_df, aes(forcats::fct_reorder(Mechanism, module_similarity, .fun = median, .na_rm = TRUE), module_similarity, fill = Mechanism)) +
    geom_boxplot(outlier.alpha = 0.2) +
    coord_flip() +
    labs(x = NULL, y = 'Module composition similarity') +
    theme_bw() +
    theme(legend.position = 'none')

  mech_temp_plot <- ggplot(temp_summary, aes(temp_bin, prop * 100, fill = Mechanism)) +
    geom_col(color = 'black', linewidth = 0.2) +
    labs(x = 'Temperature difference bin (deg C)', y = 'Mechanism proportion (%)') +
    theme_bw()

  zone_plot <- ggplot(zone_summary, aes(zone_pair, prop * 100, fill = Mechanism)) +
    geom_col(color = 'black', linewidth = 0.2) +
    labs(x = 'Zone pair', y = 'Mechanism proportion (%)') +
    theme_bw()

  station_order <- meta |>
    dplyr::arrange(thetao, Longitude.x) |>
    dplyr::pull(Sample_ID)

  module_long <- as.data.frame(module_mat) |>
    tibble::rownames_to_column('Sample_ID') |>
    tidyr::pivot_longer(cols = dplyr::all_of(main_modules), names_to = 'Module', values_to = 'Abundance') |>
    dplyr::left_join(meta |>
      dplyr::select(Sample_ID, Zona, thetao, Longitude.x), by = 'Sample_ID') |>
    dplyr::mutate(Sample_ID = factor(Sample_ID, levels = station_order), Module = factor(Module, levels = main_modules))

  gradient_plot <- ggplot(module_long, aes(Sample_ID, Abundance, fill = Module)) +
    geom_col(width = 0.9) +
    facet_grid(~ Zona, scales = 'free_x', space = 'free_x') +
    labs(x = 'Stations ordered by temperature and longitude', y = 'Main-module relative abundance') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  combo_plot <- (mech_similarity_plot | mech_temp_plot) / (zone_plot | gradient_plot)

  list(
    pairwise = pair_df,
    module_mechanism_summary = module_mech_summary,
    zone_summary = zone_summary,
    temp_summary = temp_summary,
    correlations = cor_tests,
    main_modules = main_modules,
    plots = list(
      mechanism_similarity = mech_similarity_plot,
      mechanism_temperature = mech_temp_plot,
      mechanism_zone_pair = zone_plot,
      module_gradient = gradient_plot,
      integrated = combo_plot
    )
  )
}
