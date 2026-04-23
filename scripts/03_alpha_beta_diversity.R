source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

calc_vif_numeric <- function(df) {
  num_df <- df[, vapply(df, is.numeric, logical(1)), drop = FALSE]
  out <- lapply(names(num_df), function(v) {
    others <- setdiff(names(num_df), v)
    if (length(others) == 0) return(data.frame(variable = v, vif = NA_real_))
    mod <- lm(reformulate(others, response = v), data = num_df)
    r2 <- summary(mod)$r.squared
    data.frame(variable = v, vif = 1 / (1 - r2))
  })
  dplyr::bind_rows(out)
}

run_alpha_beta_diversity <- function(cfg = sea4blue_config()) {
  load_packages(c("phyloseq", "dplyr", "ggplot2", "ggpubr", "ggrepel", "ggforce", "ggside", "vegan", "tibble"))
  objs <- load_phyloseq_atlantic(cfg)
  phyloseq_obj_css <- objs$phyloseq_obj_css

  physeq_pa <- transform_sample_counts(phyloseq_obj_css, function(x) ifelse(x > 0, 1, 0))
  richness_plot_obj <- plot_richness(physeq_pa, x = "Longitude.x", color = "Zona", measures = c("Observed"))
  richness_plot_obj$data$Zona <- factor(as.character(richness_plot_obj$data$Zona), levels = c("ANW", "ANC", "ANE"))

  alpha_plot <- ggplot(richness_plot_obj$data, aes(Longitude.x, value)) +
    theme_bw() +
    geom_point(aes(colour = Zona)) +
    scale_color_manual(values = sea4blue_zone_colors(), drop = FALSE) +
    ggside::geom_ysideboxplot(aes(x = Zona, y = value, colour = Zona), orientation = "x") +
    theme(ggside.panel.scale.y = .4) +
    ggside::scale_ysidex_discrete() +
    geom_smooth(color = 'black', linewidth = 1.1, se = TRUE) +
    ylab("Richness")

  bray <- phyloseq::distance(phyloseq_obj_css, method = "bray")
  uniF <- phyloseq::distance(phyloseq_obj_css, method = "wunifrac")
  sampledf <- data.frame(sample_data(phyloseq_obj_css))

  bray_ord <- ordinate(physeq = phyloseq_obj_css, method = "PCoA", distance = "bray")
  wunifrac_ord <- ordinate(physeq = phyloseq_obj_css, method = "PCoA", distance = "wunifrac")

  bray_plot <- plot_ordination(phyloseq_obj_css, bray_ord, color = "Zona", axes = c(1, 2)) +
    scale_color_manual(values = sea4blue_zone_colors(), drop = FALSE) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = Giorno.di.navigazione.), max.overlaps = Inf, show.legend = FALSE) +
    ggforce::geom_mark_ellipse(aes(color = Zona), show.legend = FALSE) +
    theme_bw()

  wunifrac_plot <- plot_ordination(phyloseq_obj_css, wunifrac_ord, color = "Zona", axes = c(1, 2)) +
    scale_color_manual(values = sea4blue_zone_colors(), drop = FALSE) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = Giorno.di.navigazione.), max.overlaps = Inf, show.legend = FALSE) +
    ggforce::geom_mark_ellipse(aes(color = Zona), show.legend = FALSE) +
    theme_bw()

  env_vars <- c("Zona", "Giorno.di.navigazione.", "thetao", "so", "chl", "o2", "no3", "po4", "si", "Latitude.x", "Longitude.x")
  env_df <- sampledf[, env_vars]
  env_df <- dplyr::mutate(env_df, dplyr::across(where(is.character), as.factor))
  env_df <- dplyr::mutate(env_df, dplyr::across(c(thetao, so, chl, o2, no3, po4, si, Latitude.x, Longitude.x), as.numeric))

  bray_betadisper <- vegan::betadisper(bray, group = env_df$Zona)
  uni_betadisper <- vegan::betadisper(uniF, group = env_df$Zona)
  vif_tbl <- calc_vif_numeric(env_df)
  env_numeric <- dplyr::select(env_df, thetao, so, chl, o2, no3, po4, si, Latitude.x, Longitude.x)
  env_corr <- cor(env_numeric, use = "pairwise.complete.obs")
  db_rda_bray <- vegan::capscale(bray ~ thetao + so + chl + o2 + no3 + po4 + si + Longitude.x + Latitude.x + Condition(Zona), data = env_df)
  db_rda_uni <- vegan::capscale(uniF ~ thetao + so + chl + o2 + no3 + po4 + si + Longitude.x + Latitude.x + Condition(Zona), data = env_df)

  permanova <- list(
    bray_zone = vegan::adonis2(bray ~ Zona, data = env_df),
    bray_navigation = vegan::adonis2(bray ~ Zona * Giorno.di.navigazione., data = env_df),
    bray_full = vegan::adonis2(bray ~ Zona + Giorno.di.navigazione. + thetao + so + chl + o2 + no3 + po4 + si + Latitude.x + Longitude.x, data = env_df, by = "margin"),
    wunifrac_zone = vegan::adonis2(uniF ~ Zona, data = env_df),
    wunifrac_full = vegan::adonis2(uniF ~ Zona + Giorno.di.navigazione. + thetao + so + chl + o2 + no3 + po4 + si + Latitude.x + Longitude.x, data = env_df, by = "margin")
  )

  diagnostics <- list(
    bray_betadisper = bray_betadisper,
    bray_betadisper_anova = anova(bray_betadisper),
    uni_betadisper = uni_betadisper,
    uni_betadisper_anova = anova(uni_betadisper),
    vif = vif_tbl,
    env_corr = env_corr,
    dbrda_bray = db_rda_bray,
    dbrda_bray_anova = anova(db_rda_bray),
    dbrda_uni = db_rda_uni,
    dbrda_uni_anova = anova(db_rda_uni)
  )

  list(alpha_plot = alpha_plot, bray_plot = bray_plot, wunifrac_plot = wunifrac_plot, permanova = permanova, diagnostics = diagnostics, sample_data = env_df)
}
