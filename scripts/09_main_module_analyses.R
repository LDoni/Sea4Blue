source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")
source("/home/userbio/Sea4Blue/scripts/07_module_environment_preferences.R")

weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (!length(x)) {
    return(stats::setNames(rep(NA_real_, length(probs)), paste0("q", probs * 100)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  stats::setNames(
    vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1)),
    paste0("q", probs * 100)
  )
}

module_taxonomy_table <- function(cluster_tbl, tax_df, rank = "Phylum", top_n = 6) {
  cluster_tbl |>
    dplyr::left_join(
      tibble::rownames_to_column(tax_df, "OTU_ID"),
      by = "OTU_ID"
    ) |>
    dplyr::mutate(
      RankValue = .data[[rank]],
      RankValue = dplyr::if_else(is.na(RankValue) | RankValue == "", "Unclassified", RankValue)
    ) |>
    dplyr::count(Module, RankValue, name = "n_asv") |>
    dplyr::group_by(Module) |>
    dplyr::mutate(prop_asv = n_asv / sum(n_asv)) |>
    dplyr::arrange(Module, dplyr::desc(prop_asv), RankValue) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()
}

run_main_module_analyses <- function(cfg = sea4blue_config(), top_n = 4) {
  load_packages(c("dplyr", "tibble", "tidyr", "ggplot2", "vegan", "picante", "phyloseq", "ggrepel", "ape", "scales"))
  objs <- load_phyloseq_atlantic(cfg)
  mod <- run_module_environment_preferences(cfg)
  zone_colors <- sea4blue_zone_colors()

  ps <- objs$phyloseq_obj_css
  sample_df <- data.frame(phyloseq::sample_data(ps))
  sample_df$Sample_ID <- rownames(sample_df)
  sample_df <- sample_df |>
    dplyr::mutate(
      dplyr::across(dplyr::all_of(c("thetao", "so", "chl", "o2", "no3", "po4", "si", "Longitude.x", "Latitude.x")), as.numeric),
      Zona = factor(as.character(Zona), levels = names(zone_colors))
    )

  module_abund <- mod$module_abundance
  main_modules <- mod$diversity |>
    dplyr::arrange(dplyr::desc(mean_relative_abundance), module_sort_key(Module)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(Module)
  main_modules <- as.character(main_modules)

  thermal_tbl <- lapply(main_modules, function(mod_id) {
    w <- module_abund[mod_id, ]
    opt <- stats::weighted.mean(sample_df$thetao, w = w, na.rm = TRUE)
    breadth <- sqrt(stats::weighted.mean((sample_df$thetao - opt)^2, w = w, na.rm = TRUE))
    qs <- weighted_quantile(sample_df$thetao, w = w, probs = c(0.1, 0.5, 0.9))
    thetao_row <- mod$correlations |>
      dplyr::filter(Module == mod_id, Variable == "thetao") |>
      dplyr::slice(1)
    thetao_rho <- if (nrow(thetao_row) == 1) thetao_row$rho[[1]] else NA_real_
    thetao_fdr <- if (nrow(thetao_row) == 1) thetao_row$fdr[[1]] else NA_real_
    tibble::tibble(
      Module = mod_id,
      mean_relative_abundance = mean(module_abund[mod_id, ], na.rm = TRUE),
      thermal_optimum = opt,
      thermal_breadth_sd = breadth,
      q10 = qs[[1]],
      q50 = qs[[2]],
      q90 = qs[[3]],
      thetao_rho = thetao_rho,
      thetao_fdr = thetao_fdr
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(module_sort_key(Module))

  thermal_tbl$Module <- factor(thermal_tbl$Module, levels = rev(main_modules))
  thermal_plot <- ggplot(thermal_tbl, aes(y = Module)) +
    geom_segment(
      aes(x = q10, xend = q90, yend = Module),
      linewidth = 1.1,
      color = "grey55"
    ) +
    geom_linerange(
      aes(xmin = thermal_optimum - thermal_breadth_sd, xmax = thermal_optimum + thermal_breadth_sd),
      linewidth = 3,
      color = "#1f4e79"
    ) +
    geom_point(aes(x = thermal_optimum, size = mean_relative_abundance), color = "#d94801") +
    labs(
      x = "Temperature niche (thetao)",
      y = "Main module",
      size = "Mean relative abundance"
    ) +
    theme_bw()

  cluster_tbl <- mod$cluster_table |>
    dplyr::filter(Module %in% main_modules)
  tax_df <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)

  phylum_tbl <- module_taxonomy_table(cluster_tbl, tax_df, rank = "Phylum", top_n = 8) |>
    dplyr::mutate(Module = factor(Module, levels = main_modules))
  family_tbl <- module_taxonomy_table(cluster_tbl, tax_df, rank = "Family", top_n = 10) |>
    dplyr::mutate(Module = factor(Module, levels = main_modules))

  taxonomy_plot <- ggplot(phylum_tbl, aes(Module, prop_asv, fill = RankValue)) +
    geom_col(color = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Main module", y = "ASV composition", fill = "Phylum") +
    theme_bw()

  module_taxa_matrix <- matrix(
    0L,
    nrow = length(main_modules),
    ncol = phyloseq::ntaxa(ps),
    dimnames = list(main_modules, phyloseq::taxa_names(ps))
  )
  for (mod_id in main_modules) {
    module_taxa_matrix[mod_id, cluster_tbl$OTU_ID[cluster_tbl$Module == mod_id]] <- 1L
  }

  phy_tree <- phyloseq::phy_tree(ps)
  tree_dist <- ape::cophenetic.phylo(phy_tree)
  tree_dist <- tree_dist[colnames(module_taxa_matrix), colnames(module_taxa_matrix)]

  ses_mpd <- picante::ses.mpd(
    samp = module_taxa_matrix,
    dis = tree_dist,
    null.model = "taxa.labels",
    abundance.weighted = FALSE,
    runs = 999
  )
  ses_mntd <- picante::ses.mntd(
    samp = module_taxa_matrix,
    dis = tree_dist,
    null.model = "taxa.labels",
    abundance.weighted = FALSE,
    runs = 999
  )

  phylo_tbl <- tibble::tibble(
    Module = rownames(ses_mpd),
    n_taxa = rowSums(module_taxa_matrix),
    MPD_obs = ses_mpd$mpd.obs,
    MPD_z = ses_mpd$mpd.obs.z,
    MPD_p = ses_mpd$mpd.obs.p,
    MNTD_obs = ses_mntd$mntd.obs,
    MNTD_z = ses_mntd$mntd.obs.z,
    MNTD_p = ses_mntd$mntd.obs.p
  ) |>
    dplyr::mutate(
      NRI = -MPD_z,
      NTI = -MNTD_z,
      MPD_fdr = p.adjust(MPD_p, method = "fdr"),
      MNTD_fdr = p.adjust(MNTD_p, method = "fdr"),
      Structure = dplyr::case_when(
        NRI > 1.96 | NTI > 1.96 ~ "Clustered",
        NRI < -1.96 | NTI < -1.96 ~ "Overdispersed",
        TRUE ~ "Neutral"
      )
    ) |>
    dplyr::arrange(module_sort_key(Module))

  phylo_long <- phylo_tbl |>
    dplyr::select(Module, NRI, NTI) |>
    tidyr::pivot_longer(cols = c(NRI, NTI), names_to = "Metric", values_to = "Value") |>
    dplyr::mutate(Module = factor(Module, levels = main_modules))

  phylo_plot <- ggplot(phylo_long, aes(Module, Value, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey45") +
    scale_fill_manual(values = c(NRI = "#1f78b4", NTI = "#e31a1c")) +
    labs(x = "Main module", y = "Standardized effect size", fill = NULL) +
    theme_bw()

  module_comm <- t(module_abund[main_modules, , drop = FALSE])
  module_comm_hel <- vegan::decostand(module_comm, method = "hellinger")
  env_vars <- c("thetao", "o2", "no3", "po4", "si", "so", "chl", "Longitude.x", "Latitude.x")
  env_df <- sample_df |>
    dplyr::select(dplyr::all_of(env_vars), Zona) |>
    stats::na.omit()
  module_comm_hel <- module_comm_hel[rownames(env_df), , drop = FALSE]

  dbrda_mod <- vegan::capscale(
    module_comm_hel ~ thetao + o2 + no3 + po4 + si + so + chl + Longitude.x + Latitude.x,
    data = env_df,
    distance = "euclidean"
  )

  dbrda_terms <- vegan::anova.cca(dbrda_mod, by = "terms", permutations = 999)
  dbrda_axes <- vegan::anova.cca(dbrda_mod, by = "axis", permutations = 999)
  site_scores <- as.data.frame(vegan::scores(dbrda_mod, display = "sites", choices = 1:2))
  site_scores$Zona <- env_df$Zona
  site_scores$Sample_ID <- rownames(site_scores)
  bp_scores <- as.data.frame(vegan::scores(dbrda_mod, display = "bp", choices = 1:2))
  bp_scores$Variable <- rownames(bp_scores)
  sp_scores <- as.data.frame(vegan::scores(dbrda_mod, display = "species", choices = 1:2))
  sp_scores$Module <- rownames(sp_scores)
  colnames(site_scores)[1:2] <- c("CAP1", "CAP2")
  colnames(bp_scores)[1:2] <- c("CAP1", "CAP2")
  colnames(sp_scores)[1:2] <- c("CAP1", "CAP2")

  dbrda_plot <- ggplot(site_scores, aes(CAP1, CAP2, color = Zona)) +
    geom_point(size = 2.8, alpha = 0.9) +
    stat_ellipse(aes(group = Zona), linewidth = 0.6, alpha = 0.5, show.legend = FALSE) +
    geom_segment(
      data = bp_scores,
      aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.18, "cm")),
      color = "grey30"
    ) +
    ggrepel::geom_text_repel(
      data = bp_scores,
      aes(CAP1, CAP2, label = Variable),
      inherit.aes = FALSE,
      size = 3,
      color = "grey20"
    ) +
    ggrepel::geom_text_repel(
      data = sp_scores,
      aes(CAP1, CAP2, label = Module),
      inherit.aes = FALSE,
      size = 3.6,
      fontface = "bold",
      color = "black"
    ) +
    scale_color_manual(values = zone_colors, drop = FALSE) +
    labs(x = "dbRDA1", y = "dbRDA2", color = "Zone") +
    theme_bw()

  list(
    main_modules = main_modules,
    thermal_niche = thermal_tbl,
    thermal_plot = thermal_plot,
    phylum_composition = phylum_tbl,
    family_composition = family_tbl,
    taxonomy_plot = taxonomy_plot,
    phylogenetic_structure = phylo_tbl,
    phylo_plot = phylo_plot,
    dbrda = dbrda_mod,
    dbrda_terms = dbrda_terms,
    dbrda_axes = dbrda_axes,
    dbrda_plot = dbrda_plot
  )
}
