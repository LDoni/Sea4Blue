source('/home/userbio/Sea4Blue/scripts/03_alpha_beta_diversity.R')
source('/home/userbio/Sea4Blue/scripts/04_distance_decay.R')
source('/home/userbio/Sea4Blue/scripts/05_stegen_mechanisms.R')
source('/home/userbio/Sea4Blue/scripts/06_mechanisms_temperature.R')
source('/home/userbio/Sea4Blue/scripts/07_module_environment_preferences.R')
source('/home/userbio/Sea4Blue/scripts/08_ocean_connectivity_associations.R')
source('/home/userbio/Sea4Blue/scripts/09_main_module_analyses.R')
source('/home/userbio/Sea4Blue/scripts/10_assembly_modules_gradient.R')
source('/home/userbio/Sea4Blue/scripts/11_variation_partitioning.R')

cfg <- sea4blue_config()
ensure_result_dirs(cfg)

save_plot_a4 <- function(plot_obj, path, landscape = TRUE) {
  width <- if (landscape) 11.69 else 8.27
  height <- if (landscape) 8.27 else 11.69
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, units = 'in', device = grDevices::cairo_pdf)
}

save_pheatmap_a4 <- function(pheatmap_obj, path, landscape = TRUE) {
  width <- if (landscape) 11.69 else 8.27
  height <- if (landscape) 8.27 else 11.69
  grDevices::cairo_pdf(filename = path, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(pheatmap_obj$gtable)
  grDevices::dev.off()
}

save_table <- function(x, path, row.names = FALSE) {
  utils::write.csv(x, path, row.names = row.names)
}

write_capture <- function(object, path) {
  utils::capture.output(print(object), file = path)
}

results <- list()

results$alpha_beta <- run_alpha_beta_diversity(cfg)
save_plot_a4(results$alpha_beta$alpha_plot, file.path(cfg$figures_dir, 'alpha_diversity_longitude.pdf'))
save_plot_a4(results$alpha_beta$bray_plot, file.path(cfg$figures_dir, 'beta_bray_ordination.pdf'))
save_plot_a4(results$alpha_beta$wunifrac_plot, file.path(cfg$figures_dir, 'beta_wunifrac_ordination.pdf'))
write_capture(results$alpha_beta$permanova$bray_zone, file.path(cfg$tables_dir, 'permanova_bray_zone.txt'))
write_capture(results$alpha_beta$permanova$bray_navigation, file.path(cfg$tables_dir, 'permanova_bray_navigation.txt'))
write_capture(results$alpha_beta$permanova$bray_full, file.path(cfg$tables_dir, 'permanova_bray_full.txt'))
write_capture(results$alpha_beta$permanova$wunifrac_zone, file.path(cfg$tables_dir, 'permanova_wunifrac_zone.txt'))
write_capture(results$alpha_beta$permanova$wunifrac_full, file.path(cfg$tables_dir, 'permanova_wunifrac_full.txt'))
write_capture(results$alpha_beta$diagnostics$bray_betadisper_anova, file.path(cfg$tables_dir, 'betadisper_bray_anova.txt'))
write_capture(results$alpha_beta$diagnostics$uni_betadisper_anova, file.path(cfg$tables_dir, 'betadisper_wunifrac_anova.txt'))
save_table(results$alpha_beta$diagnostics$vif, file.path(cfg$tables_dir, 'environment_vif.csv'))
save_table(as.data.frame(results$alpha_beta$diagnostics$env_corr), file.path(cfg$tables_dir, 'environment_correlation_matrix.csv'), row.names = TRUE)
write_capture(results$alpha_beta$diagnostics$dbrda_bray_anova, file.path(cfg$tables_dir, 'dbrda_bray_anova.txt'))
write_capture(results$alpha_beta$diagnostics$dbrda_uni_anova, file.path(cfg$tables_dir, 'dbrda_unifrac_anova.txt'))

results$distance_decay <- run_distance_decay(cfg)
save_plot_a4(results$distance_decay$plots$combined, file.path(cfg$figures_dir, 'distance_decay_combined.pdf'))
save_table(results$distance_decay$data, file.path(cfg$tables_dir, 'distance_decay_pairs.csv'))
write_capture(results$distance_decay$stats$bray_geo_lm, file.path(cfg$tables_dir, 'distance_decay_bray_geo_lm.txt'))
write_capture(results$distance_decay$stats$unifrac_geo_lm, file.path(cfg$tables_dir, 'distance_decay_unifrac_geo_lm.txt'))
write_capture(results$distance_decay$stats$bray_geo_mantel, file.path(cfg$tables_dir, 'distance_decay_bray_geo_mantel.txt'))
write_capture(results$distance_decay$stats$unifrac_geo_mantel, file.path(cfg$tables_dir, 'distance_decay_unifrac_geo_mantel.txt'))
write_capture(results$distance_decay$stats$bray_partial, file.path(cfg$tables_dir, 'distance_decay_bray_partial.txt'))
write_capture(results$distance_decay$stats$unifrac_partial, file.path(cfg$tables_dir, 'distance_decay_unifrac_partial.txt'))

results$stegen <- get_mechanism_proportions(cfg)
save_plot_a4(results$stegen$plot, file.path(cfg$figures_dir, 'stegen_mechanism_proportions.pdf'))
save_table(results$stegen$summary, file.path(cfg$tables_dir, 'stegen_mechanism_summary.csv'))
save_table(results$stegen$pairwise, file.path(cfg$tables_dir, 'stegen_pairwise.csv'))

results$mech_temp <- run_mechanisms_temperature(cfg)
save_plot_a4(results$mech_temp$plot, file.path(cfg$figures_dir, 'stegen_vs_temperature.pdf'))
save_table(results$mech_temp$summary, file.path(cfg$tables_dir, 'stegen_temperature_summary.csv'))

results$modules <- run_module_environment_preferences(cfg)
save_table(results$modules$signature, file.path(cfg$tables_dir, 'module_environment_signature.csv'))
save_table(results$modules$warm_rank, file.path(cfg$tables_dir, 'module_warm_rank.csv'))
save_table(results$modules$diversity, file.path(cfg$tables_dir, 'module_diversity_summary.csv'))
save_table(results$modules$correlations, file.path(cfg$tables_dir, 'module_environment_correlations.csv'))
save_table(results$modules$zone_tests, file.path(cfg$tables_dir, 'module_zone_tests.csv'))
save_plot_a4(results$modules$temperature_plot, file.path(cfg$figures_dir, 'module_vs_temperature.pdf'))
save_plot_a4(results$modules$zone_plot, file.path(cfg$figures_dir, 'module_vs_zone.pdf'))
save_pheatmap_a4(results$modules$heatmap, file.path(cfg$figures_dir, 'module_environment_heatmap.pdf'))
save_pheatmap_a4(results$modules$correlation_heatmap, file.path(cfg$figures_dir, 'module_environment_correlogram.pdf'))

results$main_modules <- run_main_module_analyses(cfg)
save_table(data.frame(Module = results$main_modules$main_modules), file.path(cfg$tables_dir, 'main_modules_selected.csv'))
save_table(results$main_modules$thermal_niche, file.path(cfg$tables_dir, 'main_module_thermal_niche.csv'))
save_table(results$main_modules$phylum_composition, file.path(cfg$tables_dir, 'main_module_phylum_composition.csv'))
save_table(results$main_modules$family_composition, file.path(cfg$tables_dir, 'main_module_family_composition.csv'))
save_table(results$main_modules$phylogenetic_structure, file.path(cfg$tables_dir, 'main_module_phylogenetic_structure.csv'))
write_capture(results$main_modules$dbrda_terms, file.path(cfg$tables_dir, 'main_module_dbrda_terms.txt'))
write_capture(results$main_modules$dbrda_axes, file.path(cfg$tables_dir, 'main_module_dbrda_axes.txt'))
save_plot_a4(results$main_modules$thermal_plot, file.path(cfg$figures_dir, 'main_module_thermal_niche.pdf'))
save_plot_a4(results$main_modules$taxonomy_plot, file.path(cfg$figures_dir, 'main_module_taxonomy.pdf'))
save_plot_a4(results$main_modules$phylo_plot, file.path(cfg$figures_dir, 'main_module_phylogenetic_structure.pdf'))
save_plot_a4(results$main_modules$dbrda_plot, file.path(cfg$figures_dir, 'main_module_dbrda.pdf'))

results$assembly_modules <- run_assembly_modules_gradient(cfg)
save_table(results$assembly_modules$pairwise, file.path(cfg$tables_dir, 'assembly_module_pairwise.csv'))
save_table(results$assembly_modules$module_mechanism_summary, file.path(cfg$tables_dir, 'assembly_module_summary.csv'))
save_table(results$assembly_modules$zone_summary, file.path(cfg$tables_dir, 'assembly_zone_summary.csv'))
save_table(results$assembly_modules$temp_summary, file.path(cfg$tables_dir, 'assembly_temperature_summary.csv'))
save_table(results$assembly_modules$correlations, file.path(cfg$tables_dir, 'assembly_module_correlations.csv'))
save_plot_a4(results$assembly_modules$plots$mechanism_similarity, file.path(cfg$figures_dir, 'assembly_module_similarity.pdf'))
save_plot_a4(results$assembly_modules$plots$mechanism_temperature, file.path(cfg$figures_dir, 'assembly_mechanism_temperature_bins.pdf'))
save_plot_a4(results$assembly_modules$plots$mechanism_zone_pair, file.path(cfg$figures_dir, 'assembly_mechanism_zone_pairs.pdf'))
save_plot_a4(results$assembly_modules$plots$module_gradient, file.path(cfg$figures_dir, 'assembly_module_gradient.pdf'))
save_plot_a4(results$assembly_modules$plots$integrated, file.path(cfg$figures_dir, 'assembly_module_gradient_integrated.pdf'))

results$varpart <- run_variation_partitioning(cfg)
save_table(results$varpart$fraction_table, file.path(cfg$tables_dir, 'variation_partitioning_summary.csv'))
write_capture(results$varpart$community_varpart, file.path(cfg$tables_dir, 'variation_partitioning_community.txt'))
write_capture(results$varpart$module_varpart, file.path(cfg$tables_dir, 'variation_partitioning_modules.txt'))
write_capture(results$varpart$community_tests$env, file.path(cfg$tables_dir, 'variation_partitioning_community_env.txt'))
write_capture(results$varpart$community_tests$pure_env, file.path(cfg$tables_dir, 'variation_partitioning_community_pure_env.txt'))
write_capture(results$varpart$community_tests$pure_space, file.path(cfg$tables_dir, 'variation_partitioning_community_pure_space.txt'))
write_capture(results$varpart$module_tests$env, file.path(cfg$tables_dir, 'variation_partitioning_modules_env.txt'))
write_capture(results$varpart$module_tests$pure_env, file.path(cfg$tables_dir, 'variation_partitioning_modules_pure_env.txt'))
write_capture(results$varpart$module_tests$pure_space, file.path(cfg$tables_dir, 'variation_partitioning_modules_pure_space.txt'))
save_plot_a4(results$varpart$fraction_plot, file.path(cfg$figures_dir, 'variation_partitioning_fraction_plot.pdf'))

results$ocean <- run_ocean_connectivity_associations(cfg)
save_plot_a4(results$ocean$plot, file.path(cfg$figures_dir, 'ocean_connectivity_vs_bray.pdf'))
save_table(results$ocean$pairwise, file.path(cfg$tables_dir, 'ocean_connectivity_pairs.csv'))
write_capture(results$ocean$stats$mantel_bray_ocean, file.path(cfg$tables_dir, 'ocean_mantel_bray.txt'))
write_capture(results$ocean$stats$mantel_unifrac_ocean, file.path(cfg$tables_dir, 'ocean_mantel_unifrac.txt'))
write_capture(results$ocean$stats$partial_bray_ocean_geo, file.path(cfg$tables_dir, 'ocean_partial_bray_geo.txt'))
write_capture(results$ocean$stats$partial_unifrac_ocean_geo, file.path(cfg$tables_dir, 'ocean_partial_unifrac_geo.txt'))
write_capture(results$ocean$stats$partial_bray_ocean_temp, file.path(cfg$tables_dir, 'ocean_partial_bray_temp.txt'))
write_capture(results$ocean$stats$partial_unifrac_ocean_temp, file.path(cfg$tables_dir, 'ocean_partial_unifrac_temp.txt'))
save_table(as.data.frame(results$ocean$connectivity), file.path(cfg$tables_dir, 'ocean_connectivity_symmetric.csv'), row.names = TRUE)

saveRDS(results, file.path(cfg$project_root, 'results', 'analysis_bundle.rds'))
cat('Final analyses completed
')
