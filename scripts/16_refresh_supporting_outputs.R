source('/home/userbio/Sea4Blue/scripts/03_alpha_beta_diversity.R')
source('/home/userbio/Sea4Blue/scripts/09_main_module_analyses.R')

save_plot_a4 <- function(plot_obj, path, landscape = FALSE) {
  width <- if (landscape) 11.69 else 8.27
  height <- if (landscape) 8.27 else 11.69
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, units = 'in', device = grDevices::cairo_pdf)
}

run_refresh_supporting_outputs <- function(cfg = sea4blue_config()) {
  load_packages(c('dplyr', 'ggplot2', 'phyloseq', 'scales'))
  ensure_result_dirs(cfg)

  alpha_beta <- run_alpha_beta_diversity(cfg)
  alpha_df <- alpha_beta$alpha_plot$data |>
    dplyr::transmute(
      Sample_ID = as.character(samples),
      Zona = as.character(Zona),
      Longitude.x = as.numeric(Longitude.x),
      Observed = as.numeric(value)
    ) |>
    dplyr::arrange(Longitude.x)

  lm_fit <- stats::lm(Observed ~ Longitude.x, data = alpha_df)
  spearman_fit <- stats::cor.test(alpha_df$Observed, alpha_df$Longitude.x, method = 'spearman', exact = FALSE)
  kruskal_fit <- stats::kruskal.test(Observed ~ Zona, data = alpha_df)

  alpha_tests <- dplyr::bind_rows(
    tibble::tibble(
      Analysis = 'Alpha diversity vs longitude (linear model)',
      Statistic = unname(summary(lm_fit)$coefficients['Longitude.x', 't value']),
      Estimate = unname(summary(lm_fit)$coefficients['Longitude.x', 'Estimate']),
      P_value = unname(summary(lm_fit)$coefficients['Longitude.x', 'Pr(>|t|)']),
      Extra = paste0('R2=', round(summary(lm_fit)$r.squared, 3))
    ),
    tibble::tibble(
      Analysis = 'Alpha diversity vs longitude (Spearman)',
      Statistic = unname(spearman_fit$statistic),
      Estimate = unname(spearman_fit$estimate),
      P_value = spearman_fit$p.value,
      Extra = paste0('n=', nrow(alpha_df))
    ),
    tibble::tibble(
      Analysis = 'Alpha diversity across zones (Kruskal-Wallis)',
      Statistic = unname(kruskal_fit$statistic),
      Estimate = NA_real_,
      P_value = kruskal_fit$p.value,
      Extra = paste0('df=', unname(kruskal_fit$parameter))
    )
  )

  utils::write.csv(alpha_df, file.path(cfg$tables_dir, 'alpha_diversity_observed.csv'), row.names = FALSE)
  utils::write.csv(alpha_tests, file.path(cfg$tables_dir, 'alpha_diversity_tests.csv'), row.names = FALSE)

  main_mod <- run_main_module_analyses(cfg)
  thermal_tbl <- main_mod$thermal_niche |>
    dplyr::arrange(thermal_optimum) |>
    dplyr::mutate(
      Module = factor(Module, levels = Module),
      mean_pct = mean_relative_abundance * 100,
      rho_label = paste0('rho = ', sprintf('%.2f', thetao_rho), ifelse(thetao_fdr < 0.05, ' *', ''))
    )

  thermal_plot_pretty <- ggplot2::ggplot(thermal_tbl, ggplot2::aes(y = Module, x = thermal_optimum)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = q10, xend = q90, yend = Module),
      linewidth = 1.2,
      color = '#bdbdbd',
      lineend = 'round'
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = thermal_optimum - thermal_breadth_sd, xmax = thermal_optimum + thermal_breadth_sd, color = thermal_optimum),
      linewidth = 4.6,
      lineend = 'round'
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = mean_pct, fill = thermal_optimum),
      shape = 21,
      color = 'black',
      stroke = 0.35,
      alpha = 0.98
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Module),
      fontface = 'bold',
      size = 4.1,
      color = 'black'
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = q90 + 0.35, label = rho_label),
      hjust = 0,
      size = 3.5,
      color = '#444444'
    ) +
    ggplot2::scale_fill_gradientn(colors = c('#2b8cbe', '#7bccc4', '#fdd49e', '#ef6548')) +
    ggplot2::scale_color_gradientn(colors = c('#2b8cbe', '#7bccc4', '#fdd49e', '#ef6548'), guide = 'none') +
    ggplot2::scale_size_continuous(range = c(8, 20), breaks = c(10, 25, 45), labels = function(x) paste0(x, '%')) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.14))) +
    ggplot2::labs(
      x = expression('Thermal optimum and breadth ('*degree*C*')'),
      y = 'Dominant module',
      size = 'Mean relative abundance'
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(face = 'bold'),
      axis.title.x = ggplot2::element_text(face = 'bold'),
      axis.text.y = ggplot2::element_text(face = 'bold', color = 'black', size = 12),
      legend.position = 'right'
    )

  save_plot_a4(thermal_plot_pretty, file.path(cfg$figures_dir, 'main_module_thermal_niche.pdf'), landscape = FALSE)
  ggplot2::ggsave(
    filename = file.path(cfg$project_root, 'paper', 'draft', 'figures', 'main_module_thermal_niche.png'),
    plot = thermal_plot_pretty,
    width = 8.27,
    height = 11.69,
    units = 'in',
    dpi = 300,
    bg = 'white'
  )

  invisible(list(alpha = alpha_df, alpha_tests = alpha_tests, thermal = thermal_tbl))
}

if (sys.nframe() == 0) {
  out <- run_refresh_supporting_outputs()
  print(out$alpha_tests)
}
