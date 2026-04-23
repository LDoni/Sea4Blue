source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

weighted_mean_module <- function(module_abund, meta_df, var) {
  vals <- meta_df[[var]]
  tibble::tibble(
    Module = rownames(module_abund),
    Variable = var,
    WeightedMean = apply(module_abund, 1, function(w) weighted.mean(vals, w = w, na.rm = TRUE))
  )
}

module_env_correlation <- function(module_abund_t, meta_df, env_vars) {
  out <- lapply(colnames(module_abund_t), function(mod) {
    lapply(env_vars, function(v) {
      test <- suppressWarnings(cor.test(module_abund_t[[mod]], meta_df[[v]], method = 'spearman', exact = FALSE))
      tibble::tibble(Module = mod, Variable = v, rho = unname(test$estimate), p_value = test$p.value)
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()
  out |> dplyr::mutate(fdr = p.adjust(p_value, method = 'fdr'))
}

run_module_environment_preferences <- function(cfg = sea4blue_config()) {
  load_packages(c('dplyr', 'tibble', 'ggplot2', 'tidyr', 'pheatmap', 'RColorBrewer', 'phyloseq', 'readr'))
  objs <- load_phyloseq_atlantic(cfg)
  ps <- objs$phyloseq_obj_css
  sample_df <- data.frame(sample_data(ps))
  zone_colors <- sea4blue_zone_colors()

  cluster_tbl <- readRDS(cfg$paths$sparcc_cluster_rds) |>
    dplyr::transmute(Module = as.character(Cluster), OTU_ID) |>
    dplyr::filter(OTU_ID %in% taxa_names(ps))

  otu_mat <- as(phyloseq::otu_table(ps), 'matrix')
  if (!phyloseq::taxa_are_rows(ps)) otu_mat <- t(otu_mat)
  otu_mat <- otu_mat[cluster_tbl$OTU_ID, , drop = FALSE]
  module_ids <- cluster_tbl$Module
  module_abund <- rowsum(otu_mat, group = module_ids)
  module_abund <- sweep(module_abund, 2, colSums(module_abund), FUN = '/')
  module_abund <- module_abund[order(module_sort_key(rownames(module_abund))), , drop = FALSE]

  env_vars <- c('thetao', 'so', 'chl', 'o2', 'no3', 'po4', 'si', 'Longitude.x', 'Latitude.x')
  sample_df <- sample_df |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(env_vars, 'Longitude.x', 'Latitude.x', 'thetao')), as.numeric)) |>
    dplyr::mutate(Zona = factor(as.character(Zona), levels = names(zone_colors)))

  signature_long <- lapply(env_vars, function(v) weighted_mean_module(module_abund, sample_df, v)) |>
    dplyr::bind_rows()

  signature_wide <- signature_long |>
    tidyr::pivot_wider(names_from = Variable, values_from = WeightedMean) |>
    dplyr::arrange(module_sort_key(Module))

  heatmap_mat <- signature_wide |>
    tibble::column_to_rownames('Module') |>
    as.matrix()
  heatmap_scaled <- scale(heatmap_mat)

  warm_rank <- signature_wide |>
    dplyr::arrange(dplyr::desc(thetao), module_sort_key(Module)) |>
    dplyr::mutate(Warm_Rank = dplyr::row_number())

  diversity_tbl <- cluster_tbl |>
    dplyr::count(Module, name = 'n_asv') |>
    dplyr::left_join(
      tibble::tibble(Module = rownames(module_abund), mean_relative_abundance = rowMeans(module_abund, na.rm = TRUE)),
      by = 'Module'
    ) |>
    dplyr::arrange(module_sort_key(Module))

  module_abund_t <- as.data.frame(t(module_abund))
  module_abund_t$Sample_ID <- rownames(module_abund_t)
  module_sample <- sample_df |>
    dplyr::mutate(Sample_ID = rownames(sample_df)) |>
    dplyr::left_join(module_abund_t, by = 'Sample_ID')

  module_cols <- rownames(module_abund)
  corr_tbl <- module_env_correlation(module_sample[, module_cols, drop = FALSE], module_sample, env_vars)
  corr_mat <- corr_tbl |>
    dplyr::select(Module, Variable, rho) |>
    tidyr::pivot_wider(names_from = Variable, values_from = rho) |>
    dplyr::arrange(module_sort_key(Module)) |>
    tibble::column_to_rownames('Module') |>
    as.matrix()

  corr_heatmap <- pheatmap::pheatmap(
    t(corr_mat),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdBu')))(100),
    breaks = seq(-1, 1, length.out = 101),
    main = 'Module-environment Spearman correlations',
    silent = TRUE
  )

  module_long <- module_sample |>
    dplyr::select(Sample_ID, Zona, Longitude.x, thetao, dplyr::all_of(module_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(module_cols), names_to = 'Module', values_to = 'Abundance') |>
    dplyr::mutate(Module = factor(Module, levels = rownames(module_abund)))

  temperature_plot <- ggplot(module_long, aes(thetao, Abundance, color = Zona)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = 'lm', se = TRUE, linewidth = 0.6) +
    scale_color_manual(values = zone_colors, drop = FALSE) +
    facet_wrap(~ Module, scales = 'free_y') +
    labs(x = 'Temperature (thetao)', y = 'Module relative abundance') +
    theme_bw()

  zone_plot <- ggplot(module_long, aes(Zona, Abundance, fill = Zona)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 1.4, alpha = 0.75) +
    scale_fill_manual(values = zone_colors, drop = FALSE) +
    facet_wrap(~ Module, scales = 'free_y') +
    labs(x = 'Zone', y = 'Module relative abundance') +
    theme_bw()

  zone_tests <- lapply(module_cols, function(mod) {
    fit <- kruskal.test(module_sample[[mod]] ~ module_sample$Zona)
    tibble::tibble(Module = mod, statistic = unname(fit$statistic), p_value = fit$p.value)
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(fdr = p.adjust(p_value, method = 'fdr')) |>
    dplyr::arrange(fdr, module_sort_key(Module))

  heatmap_plot <- pheatmap::pheatmap(
    t(heatmap_scaled),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdBu')))(100),
    main = 'Environmental preferences of prokaryotic modules',
    silent = TRUE
  )

  list(
    cluster_table = cluster_tbl,
    module_abundance = module_abund,
    module_sample = module_sample,
    signature = signature_wide,
    warm_rank = warm_rank,
    diversity = diversity_tbl,
    heatmap = heatmap_plot,
    heatmap_matrix = t(heatmap_scaled),
    correlations = corr_tbl,
    correlation_matrix = corr_mat,
    correlation_heatmap = corr_heatmap,
    temperature_plot = temperature_plot,
    zone_plot = zone_plot,
    zone_tests = zone_tests
  )
}
