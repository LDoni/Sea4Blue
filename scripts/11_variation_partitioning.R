source('/home/userbio/Sea4Blue/scripts/utils/helpers.R')
source('/home/userbio/Sea4Blue/scripts/07_module_environment_preferences.R')

run_variation_partitioning <- function(cfg = sea4blue_config(), top_n = 4) {
  load_packages(c('dplyr', 'tibble', 'phyloseq', 'vegan', 'ggplot2'))
  objs <- load_phyloseq_atlantic(cfg)
  ps <- objs$phyloseq_obj_css
  meta <- objs$sample_data
  meta$Sample_ID <- rownames(meta)
  meta <- meta |>
    dplyr::mutate(
      thetao = as.numeric(thetao),
      o2 = as.numeric(o2),
      si = as.numeric(si),
      Longitude.x = as.numeric(Longitude.x),
      Latitude.x = as.numeric(Latitude.x)
    )

  keep <- stats::complete.cases(meta[, c('thetao', 'o2', 'si', 'Longitude.x', 'Latitude.x')])
  meta <- meta[keep, , drop = FALSE]
  ps <- phyloseq::prune_samples(meta$Sample_ID, ps)

  otu <- as(phyloseq::otu_table(ps), 'matrix')
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  otu <- otu[, meta$Sample_ID, drop = FALSE]
  comm <- t(otu)
  comm_hel <- vegan::decostand(comm, method = 'hellinger')

  env_df <- meta |>
    dplyr::select(thetao, o2, si)
  space_df <- meta |>
    dplyr::select(Longitude.x, Latitude.x)

  vp_comm <- vegan::varpart(comm_hel, env_df, space_df)

  comm_rda_env <- vegan::rda(comm_hel ~ thetao + o2 + si, data = meta)
  comm_rda_space <- vegan::rda(comm_hel ~ Longitude.x + Latitude.x, data = meta)
  comm_rda_env_space <- vegan::rda(comm_hel ~ thetao + o2 + si + Condition(Longitude.x + Latitude.x), data = meta)
  comm_rda_space_env <- vegan::rda(comm_hel ~ Longitude.x + Latitude.x + Condition(thetao + o2 + si), data = meta)

  comm_tests <- list(
    env = vegan::anova.cca(comm_rda_env, permutations = 999),
    space = vegan::anova.cca(comm_rda_space, permutations = 999),
    pure_env = vegan::anova.cca(comm_rda_env_space, permutations = 999),
    pure_space = vegan::anova.cca(comm_rda_space_env, permutations = 999)
  )

  mod <- run_module_environment_preferences(cfg)
  main_modules <- mod$diversity |>
    dplyr::arrange(dplyr::desc(mean_relative_abundance), module_sort_key(Module)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Module) |>
    as.character()
  module_mat <- t(mod$module_abundance[main_modules, meta$Sample_ID, drop = FALSE])
  module_hel <- vegan::decostand(module_mat, method = 'hellinger')

  vp_mod <- vegan::varpart(module_hel, env_df, space_df)

  mod_rda_env <- vegan::rda(module_hel ~ thetao + o2 + si, data = meta)
  mod_rda_space <- vegan::rda(module_hel ~ Longitude.x + Latitude.x, data = meta)
  mod_rda_env_space <- vegan::rda(module_hel ~ thetao + o2 + si + Condition(Longitude.x + Latitude.x), data = meta)
  mod_rda_space_env <- vegan::rda(module_hel ~ Longitude.x + Latitude.x + Condition(thetao + o2 + si), data = meta)

  mod_tests <- list(
    env = vegan::anova.cca(mod_rda_env, permutations = 999),
    space = vegan::anova.cca(mod_rda_space, permutations = 999),
    pure_env = vegan::anova.cca(mod_rda_env_space, permutations = 999),
    pure_space = vegan::anova.cca(mod_rda_space_env, permutations = 999)
  )

  frac_tbl <- dplyr::bind_rows(
    tibble::tibble(
      Response = 'Community',
      Fraction = c('Environment+Space', 'Environment|Space', 'Space|Environment'),
      Adj_R2 = c(
        vegan::RsquareAdj(comm_rda_env)$adj.r.squared + vegan::RsquareAdj(comm_rda_space_env)$adj.r.squared,
        vegan::RsquareAdj(comm_rda_env_space)$adj.r.squared,
        vegan::RsquareAdj(comm_rda_space_env)$adj.r.squared
      )
    ),
    tibble::tibble(
      Response = 'Main modules',
      Fraction = c('Environment+Space', 'Environment|Space', 'Space|Environment'),
      Adj_R2 = c(
        vegan::RsquareAdj(mod_rda_env)$adj.r.squared + vegan::RsquareAdj(mod_rda_space_env)$adj.r.squared,
        vegan::RsquareAdj(mod_rda_env_space)$adj.r.squared,
        vegan::RsquareAdj(mod_rda_space_env)$adj.r.squared
      )
    )
  )

  frac_plot <- ggplot(frac_tbl, aes(Fraction, Adj_R2, fill = Response)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    coord_flip() +
    labs(x = NULL, y = 'Adjusted R2') +
    theme_bw()

  list(
    community_varpart = vp_comm,
    module_varpart = vp_mod,
    community_tests = comm_tests,
    module_tests = mod_tests,
    main_modules = main_modules,
    fraction_table = frac_tbl,
    fraction_plot = frac_plot
  )
}
